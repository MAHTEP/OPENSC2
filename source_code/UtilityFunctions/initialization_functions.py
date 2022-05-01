import warnings
import numpy as np


def conductor_spatial_discretization(simulation, conductor):

    """
    ##############################################################################
    # CREATE CREATE THE MESH COORDINATES
    # gridO function to create the grid for the solver
    # INPUTS:
    # NELEMS = number of elements(input parameter)
    # XBREFI = starting coordinate of a refined zone, if any (input parameter)
    # XEREFI = end coordinate of a refined zone, if any (input parameter)
    # NELREF = number of elements in the refined zone, if any (input parameter)
    # grid_input["ITYMSH"] = mesh type: =0 fixed and uniform, =1 fixed refined,
    # =3 adapted with # initial refinement; =-1 read from file
    # -  XLENGTH = conductor length (input parameter, property of the conductor)
    # OBSOLETE IT is in conductor instance
    # -  NNODES = number of points in the computational grid, computed parameter # = NELEMS + 1 OBSOLETE IT is computed
    # Icond = conductor index OBSOLETE IT is computed
    # OUTPUT:
    # XCOORD = output vector containing the coordinates of the NNODES grid points

    #   def grid0(self, NELEMS, sizemin, sizemax, XBREFI, XEREFI, NELREF, grid_input["ITYMSH"]):
    #     numcond = len(conductor)
    #     NNODES = NELEMS + 1
    ##############################################################################
    # ITYMSH == 1 and ITYMSH == 3 by D.Placido PoliTo 06/2020.
    ##############################################################################
    """

    XLENGTH = conductor.dict_input["XLENGTH"]
    MAXNOD = conductor.dict_input["MAXNOD"]
    ITYMSH = conductor.dict_discretization["Grid_input"]["ITYMSH"]
    NELEMS = conductor.dict_discretization["Grid_input"]["NELEMS"]
    XBREFI = conductor.dict_discretization["Grid_input"]["XBREFI"]
    XEREFI = conductor.dict_discretization["Grid_input"]["XEREFI"]

    nnodes = NELEMS + 1
    # conductor spatial discretization initialization
    xcoord = np.zeros(nnodes)

    if ITYMSH == -1:
        # User defined mesh
        xcoord, nnodes = conductor.load_user_defined_quantity(
            simulation, "EXTERNAL_GRID", f"x_{conductor.name} [m]"
        )
        # Evaluate the number of elements from the number of nodes
        conductor.dict_discretization["Grid_input"]["NELEMS"] = nnodes - 1

    # COMPUTE THE COORDINATES IN THE FIRST TURN
    elif ITYMSH == 0 or ITYMSH == 2 or abs(XEREFI - XBREFI) <= 1e-3:
        # Consider the case of adaptive, not-initially-refined grid
        # UNIFORM SPACING
        xcoord = np.linspace(0.0, XLENGTH, nnodes)
    # !*LOCALLY REFINED MESH. COMPUTED ON A SINGLE TURN BASIS
    elif ITYMSH == 1 or ITYMSH == 3:

        NELREF = conductor.dict_discretization["Grid_input"]["NELREF"]
        SIZMIN = conductor.dict_discretization["Grid_input"]["SIZMIN"]
        SIZMAX = conductor.dict_discretization["Grid_input"]["SIZMAX"]
        DXINCRE = conductor.dict_discretization["Grid_input"]["DXINCRE"]

        # total number of elements to be used for coarse region of the mesh
        NELCOARS = NELEMS - NELREF
        # number of elements to be used in coarse region left to refined mesh zone
        NELLEFT = round((XBREFI - 0.0) / (XLENGTH - (XEREFI - XBREFI)) * NELCOARS)
        NELRIGHT = NELCOARS - NELLEFT

        NOD_ref = NELREF + 1

        # Refined zone discretization
        dx_ref = (XEREFI - XBREFI) / NELREF  # m refined zone discretization pitch
        if (dx_ref >= SIZMIN) and (dx_ref <= SIZMAX):
            # refined mesh
            xcoord[NELLEFT : NELLEFT + NELREF + 1] = np.linspace(
                XBREFI, XEREFI, NOD_ref
            )
        elif dx_ref < SIZMIN:
            raise ValueError("ERROR: GRID0 dx in refined zone < sizemin!!!\n")
        elif dx_ref > SIZMAX:
            raise ValueError("ERROR: GRID0 dx in refined zone > sizemax!!!\n")

        if NELLEFT > 0:
            # Discretization of coarse region left to refined zone
            dx_try = (XBREFI - 0.0) / NELLEFT
            dx1 = dx_ref  # dummy to not overwrite dx_ref
            ii = 0
            while (dx_try / dx1 > DXINCRE) and (ii <= NELLEFT):
                ii = ii + 1
                dx = dx1 * DXINCRE
                xcoord[NELLEFT - ii] = xcoord[NELLEFT + 1 - ii] - dx
                dx1 = dx
                dx_try = (xcoord[NELLEFT - ii] - 0.0) / (NELLEFT - ii)

            xcoord[0 : NELLEFT - ii + 1] = np.linspace(
                0.0, xcoord[NELLEFT - ii], NELLEFT - ii + 1
            )

        if NELRIGHT > 0:
            # Discretization of coarse region right to refined zone
            dx_try = (XLENGTH - XEREFI) / NELRIGHT
            dx1 = dx_ref  # dummy to not overwrite dx_ref
            ii = 0
            while (dx_try / dx1 > DXINCRE) and (ii <= NELRIGHT):
                ii = ii + 1
                dx = dx1 * DXINCRE
                xcoord[NELLEFT + NELREF + ii] = xcoord[NELLEFT + NELREF + ii - 1] + dx
                dx1 = dx
                dx_try = (XLENGTH - xcoord[NELLEFT + NELREF + ii]) / (NELRIGHT - ii)

            xcoord[NELLEFT + NELREF + ii : NELEMS + 1] = np.linspace(
                xcoord[NELLEFT + NELREF + ii], XLENGTH, NELRIGHT - ii + 1
            )

    if len(xcoord) > MAXNOD:
        warnings.warn(
            f"Number of discretization nodes larger than maximum allowed \
                  values:\n{len(xcoord)} > {MAXNOD}\n"
        )

    if xcoord[0] < 0.0 or xcoord[0] > 0.0:
        warnings.warn(
            f"From GRID0 XCOORD[0] ~= 0!\nXCOORD[0] = {xcoord[0]}, icond \
                  = {conductor.ICOND}"
        )

    if xcoord[-1] < XLENGTH or xcoord[-1] > XLENGTH:
        warnings.warn(
            f"WARNING> From GRID0:  XCOORD[-1] ~= XLENGTH!\nXCOORD[-1] = \
                  {xcoord[-1]}, icond = {conductor.ICOND}"
        )
        conductor.inHard = xcoord

    conductor.dict_discretization["N_nod"] = nnodes
    conductor.dict_discretization["xcoord"] = xcoord
    conductor.dict_discretization["Delta_x"] = xcoord[1:] - xcoord[:-1]
    conductor.dict_discretization["dx"] = conductor.dict_discretization["Delta_x"].max()
    # end function


# function to impose user defined grid
def ext_grid(path, NN):  # optimized and testetd: ok (cdp,06/2020)
    # read mesh from user defined file .dat made of two columns and NELEMS + 1 \
    # rows: first column is for node index, second row is for coordinate value
    xx = np.loadtxt(f"{path}/user_grid.dat")
    xx = xx[:, 1]  # get only x coordinates value
    if xx.size != NN:  # check on number of nodes
        raise ValueError(
            f"ERROR: number of rows in file user_grid.dat ({xx.size}) \
                     is not consistent with given input number of nodes \
                     ({NN})!!!.\n User should modify file user_grid.dat.\n"
        )
    return xx  # xcoord
