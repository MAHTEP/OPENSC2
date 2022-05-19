import logging
import numpy as np
import os
import pandas as pd
import warnings
from typing import Union

from fluid_component import FluidComponent
from jacket_component import JacketComponent

# from conductor import Conductor
from strand_component import StrandComponent
from strand_mixed_component import StrandMixedComponent
from strand_stabilizer_component import StrandStabilizerComponent
from strand_superconductor_component import StrandSuperconductorComponent
from cylindrical_helix import CylindricalHelix

logger_discretization = logging.getLogger("opensc2Logger.discretization")


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
    # =3 adapted with # initial refinement =-1 read from file
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

    XLENGTH = conductor.inputs["XLENGTH"]
    MAXNOD = conductor.inputs["MAXNOD"]
    ITYMSH = conductor.grid_input["ITYMSH"]
    NELEMS = conductor.grid_input["NELEMS"]
    XBREFI = conductor.grid_input["XBREFI"]
    XEREFI = conductor.grid_input["XEREFI"]

    nnodes = NELEMS + 1
    # conductor spatial discretization initialization
    xcoord = np.zeros(nnodes)

    if ITYMSH == -1:
        # User defined mesh
        xcoord, nnodes = conductor.load_user_defined_quantity(
            simulation, "EXTERNAL_GRID", f"x_{conductor.name} [m]"
        )
        # Evaluate the number of elements from the number of nodes
        conductor.grid_input["NELEMS"] = nnodes - 1

    # COMPUTE THE COORDINATES IN THE FIRST TURN
    elif ITYMSH == 0 or ITYMSH == 2 or abs(XEREFI - XBREFI) <= 1e-3:
        # Consider the case of adaptive, not-initially-refined grid
        # UNIFORM SPACING
        xcoord = np.linspace(0.0, XLENGTH, nnodes)
    # !*LOCALLY REFINED MESH. COMPUTED ON A SINGLE TURN BASIS
    elif ITYMSH == 1 or ITYMSH == 3:

        NELREF = conductor.grid_input["NELREF"]
        SIZMIN = conductor.grid_input["SIZMIN"]
        SIZMAX = conductor.grid_input["SIZMAX"]
        DXINCRE = conductor.grid_input["DXINCRE"]

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

    conductor.gird_features["N_nod"] = nnodes
    conductor.gird_features["xcoord"] = xcoord
    conductor.gird_features["Delta_x"] = xcoord[1:] - xcoord[:-1]
    conductor.gird_features["dx"] = conductor.gird_features["Delta_x"].max()
    # end function


def uniform_spatial_discretization(conductor: object, _=None) -> np.ndarray:
    """Evaluates straight uniform spatial discretization in z direction.

    Args:
        conductor (Conductor): conductor object, has all the information to evaluate the unifrom mesh.
        _ (_type_): not used input argument.

    Returns:
        np.ndarray: array with uniform spatial discretization along z direction of length conductor.gird_features["N_nod"] for straight conductor components.
    """
    return np.linspace(
        0.0, conductor.inputs["XLENGTH"], conductor.gird_features["N_nod"]
    )


def uniform_angular_discretization(
    conductor: object,
    comp: Union[
        StrandMixedComponent, StrandStabilizerComponent, StrandSuperconductorComponent
    ],
) -> np.ndarray:
    """Function that evaluates uniform angular discretization, used for helicoidal geometry.

    Args:
        conductor (Conductor): conductor object, has all the information to evaluate the unifrom mesh.
        comp (Union[StrandMixedComponent, StrandStabilizerComponent, StrandSuperconductorComponent]): generic object of for wich the uniform angular discretization should be evaluated.

    Returns:
        np.ndarray: array with uniform angular discretization of length conductor.gird_features["N_nod"] for helicoidal components.
    """
    return np.linspace(
        0.0,
        comp.cyl_helix.windings_number * 2 * np.pi,
        conductor.gird_features["N_nod"],
    )


def fixed_refined_spatial_discretization(
    conductor: object, zcoord: np.ndarray
) -> np.ndarray:
    """Function that evaluate fixed refined spatial discretization along z direction.

    Args:
        conductor (Conductor): conductor object, has all the information to evaluate the fixed refined grid.
        zcoord (np.ndarray): array of length conductor.gird_features["N_nod"] initialized to zeros.

    Raises:
        ValueError: if dz in refined region is lower than minimum dz.
        ValueError: if dz in refined region is larger than maximum dz.

    Returns:
        np.ndarray: array of length conductor.gird_features["N_nod"] with fixed refined spatial discretization for straight conductor components.
    """
    n_elem = dict()

    # total number of elements to be used for coarse region of the mesh
    n_elem["coarse"] = conductor.grid_input["NELEMS"] - conductor.grid_input["NELREF"]
    # number of elements to be used in coarse region left to refined mesh zone
    n_elem["left"] = round(
        (conductor.grid_input["XBREFI"] - 0.0)
        / (
            conductor.inputs["XLENGTH"]
            - (conductor.grid_input["XEREFI"] - conductor.grid_input["XBREFI"])
        )
        * n_elem["coarse"]
    )
    n_elem["right"] = n_elem["coarse"] - n_elem["left"]

    NOD_ref = conductor.grid_input["NELREF"] + 1

    # Refined zone discretization
    dx_ref = (
        conductor.grid_input["XEREFI"] - conductor.grid_input["XBREFI"]
    ) / conductor.grid_input[
        "NELREF"
    ]  # m refined zone discretization pitch
    if (dx_ref >= conductor.grid_input["SIZMIN"]) and (
        dx_ref <= conductor.grid_input["SIZMAX"]
    ):
        # refined mesh
        zcoord[
            n_elem["left"] : n_elem["left"] + conductor.grid_input["NELREF"] + 1
        ] = np.linspace(
            conductor.grid_input["XBREFI"], conductor.grid_input["XEREFI"], NOD_ref
        )
    elif dx_ref < conductor.grid_input["SIZMIN"]:
        raise ValueError(
            f"ERROR: {dx_ref=} m in refined zone < {conductor.grid_input['SIZMIN']=} m!!!\n"
        )
    elif dx_ref > conductor.grid_input["SIZMAX"]:
        raise ValueError(
            f"ERROR: {dx_ref=} m in refined zone > {conductor.grid_input['SIZMAX']} m!!!\n"
        )

    if n_elem["left"] > 0:
        # Discretization of coarse region left to refined zone
        dx_try = (conductor.grid_input["XBREFI"] - 0.0) / n_elem["left"]
        dx1 = dx_ref  # dummy to not overwrite dx_ref
        ii = 0
        while (dx_try / dx1 > conductor.grid_input["DXINCRE"]) and (
            ii <= n_elem["left"]
        ):
            ii = ii + 1
            dx = dx1 * conductor.grid_input["DXINCRE"]
            zcoord[n_elem["left"] - ii] = zcoord[n_elem["left"] + 1 - ii] - dx
            dx1 = dx
            dx_try = (zcoord[n_elem["left"] - ii] - 0.0) / (n_elem["left"] - ii)

        zcoord[0 : n_elem["left"] - ii + 1] = np.linspace(
            0.0, zcoord[n_elem["left"] - ii], n_elem["left"] - ii + 1
        )

    if n_elem["right"] > 0:
        # Discretization of coarse region right to refined zone
        dx_try = (
            conductor.inputs["XLENGTH"] - conductor.grid_input["XEREFI"]
        ) / n_elem["right"]
        dx1 = dx_ref  # dummy to not overwrite dx_ref
        ii = 0
        while (dx_try / dx1 > conductor.grid_input["DXINCRE"]) and (
            ii <= n_elem["right"]
        ):
            ii = ii + 1
            dx = dx1 * conductor.grid_input["DXINCRE"]
            zcoord[n_elem["left"] + conductor.grid_input["NELREF"] + ii] = (
                zcoord[n_elem["left"] + conductor.grid_input["NELREF"] + ii - 1] + dx
            )
            dx1 = dx
            dx_try = (
                conductor.inputs["XLENGTH"]
                - zcoord[n_elem["left"] + conductor.grid_input["NELREF"] + ii]
            ) / (n_elem["right"] - ii)

        zcoord[
            n_elem["left"]
            + conductor.grid_input["NELREF"]
            + ii : conductor.grid_input["NELEMS"]
            + 1
        ] = np.linspace(
            zcoord[n_elem["left"] + conductor.grid_input["NELREF"] + ii],
            conductor.inputs["XLENGTH"],
            n_elem["right"] - ii + 1,
        )
    return zcoord


def fixed_refined_angular_discretization(
    conductor: object,
    comp: Union[
        StrandMixedComponent, StrandStabilizerComponent, StrandSuperconductorComponent
    ],
    tau: np.ndarray,
) -> np.ndarray:
    """Function that evaluate fixed refined angular discretization, used for hlicoidal geometry.

    Args:
        conductor (Conductor): conductor object, has all the information to evaluate the fixed refined grid.
        comp (Union[ StrandMixedComponent, StrandStabilizerComponent, StrandSuperconductorComponent ]): generic object for wich the fixed refined angular discretization should be evaluated.
        tau (np.ndarray): array of length conductor.gird_features["N_nod"] initialized to zeros.

    Raises:
        ValueError: if dtau in refined region is lower than minimum dtau.
        ValueError: if dtau in refined region is larger than maximum dtau.

    Returns:
        np.ndarray: array of length conductor.gird_features["N_nod"] with fixed refined angular discretization for helicoidal conductor components.
    """
    n_elem = dict()
    n_elem["coarse"] = conductor.grid_input["NELEMS"] - conductor.grid_input["NELREF"]
    n_elem["left"] = round(
        (conductor.grid_input["XBREFI"] - 0.0)
        / (
            conductor.inputs["XLENGTH"]
            - (conductor.grid_input["XEREFI"] - conductor.grid_input["XBREFI"])
        )
        * n_elem["coarse"]
    )
    n_elem["right"] = n_elem["coarse"] - n_elem["left"]

    n_winding = dict()
    # number of windings left to the refined region
    n_winding["left"] = conductor.grid_input["XBREFI"] / (
        2 * np.pi * comp.cyl_helix.reduced_pitch
    )
    # number of windings right to the refined region
    n_winding["right"] = (
        conductor.inputs["XLENGTH"] - conductor.grid_input["EBREFI"]
    ) / (2 * np.pi * comp.cyl_helix.reduced_pitch)
    # number of windings in the refined region
    n_winding["ref"] = (
        conductor.grid_input["EBREFI"] - conductor.grid_input["XBREFI"]
    ) / (2 * np.pi * comp.cyl_helix.reduced_pitch)

    assert (
        comp.cyl_helix.winding_number
        - (n_winding["left"] + n_winding["ref"] + n_winding["right"])
        <= 1e-15
    )

    dtau_ref = 2 * np.pi * n_winding["ref"] / conductor.grid_input["NELREF"]

    if (dtau_ref >= conductor.grid_input["SIZMIN"]) and (
        dtau_ref <= conductor.grid_input["SIZMAX"]
    ):

        tau_beg = 2 * np.pi * n_winding["left"] + dtau_ref
        tau_end = 2 * np.pi * (n_winding["left"] + n_winding["ref"])
        tau[
            n_elem["left"] + 1 : (n_elem["left"] + 1) + (conductor.grid_input["NELREF"])
        ] = np.linspace(tau_beg, tau_end, conductor.grid_input["NELREF"] + 1)
    elif dtau_ref < conductor.grid_input["SIZMIN"]:
        raise ValueError(
            f"ERROR: {dtau_ref=} m in refined zone < {conductor.grid_input['SIZMIN']=} m!!!\n"
        )
    elif dtau_ref > conductor.grid_input["SIZMAX"]:
        raise ValueError(
            f"ERROR: {dtau_ref=} m in refined zone > {conductor.grid_input['SIZMAX']} m!!!\n"
        )

    if n_elem["left"] > 0:

        d_tau_try = 2 * np.pi * n_winding["left"] / n_elem["left"]
        d_tau1 = dtau_ref
        ii = 0
        while (d_tau_try / d_tau1 > conductor.grid_input["DXINCRE"]) and (
            ii <= n_elem["left"]
        ):
            d_tau = d_tau1 * conductor.grid_input["DXINCRE"]
            tau[n_elem["left"] - ii] = tau[n_elem["left"] + 1 - ii] - d_tau
            d_tau1 = d_tau
            d_tau_try = tau[n_elem["left"] - ii] / (n_elem["left"] - ii - 1)
            ii = ii + 1

        tau_beg = 0.0
        tau_end = tau(n_elem["left"] + 1 - ii)
        tau[1 : n_elem["left"] + 1 - ii] = np.linspace(
            tau_beg, tau_end, n_elem["left"] + 1 - ii
        )

    if n_elem["right"] > 0:

        d_tau_try = 2 * np.pi * n_winding["right"] / n_elem["right"]
        d_tau1 = dtau_ref
        ii = 0
        while (d_tau_try / d_tau1 > conductor.grid_input["DXINCRE"]) and (
            ii <= n_elem["right"]
        ):
            d_tau = d_tau1 * conductor.grid_input["DXINCRE"]
            tau[(n_elem["left"] + 1) + (conductor.grid_input["NELREF"]) + ii] = (
                tau[(n_elem["left"] + 1) + (conductor.grid_input["NELREF"]) + ii - 1]
                + d_tau
            )
            d_tau1 = d_tau
            d_tau_try = (
                2 * np.pi * comp.cyl_helix.winding_number
                - tau[(n_elem["left"] + 1) + (conductor.grid_input["NELREF"]) + ii]
            ) / (n_elem["left"] - ii - 1)
            ii = ii + 1

        tau_beg = tau[(n_elem["left"] + 1) + (conductor.grid_input["NELREF"]) + ii - 1]
        tau_end = 2 * np.pi * comp.cyl_helix.winding_number
        tau[
            (n_elem["left"] + 1) + (conductor.grid_input["NELREF"]) + ii - 1 : -1
        ] = np.linspace(tau_beg, tau_end, n_elem["right"] + 1 - ii + 1)

    return tau


def user_defined_grid(conductor: object):
    """Fuction that loads user defined spatial discretization of the generic conductor component.

    Args:
        conductor (Conductor): conductor object, has all the information to evaluate load and assign user defined grid.

    Notes: does not allow interpolation in space and or in time, i.e. only fixed spatial discretization is available.
    """
    # Build file path.
    file_path = os.path.join(conductor.BASE_PATH, conductor.file_input["EXTERNAL_GRID"])
    # Load all sheets of user defined grid auxiliary input file.
    coord_dfs = pd.read_excel(file_path, sheet_name=None)
    # Check user defined grid features.
    (
        conductor.grid_features["N_nod"],
        conductor.grid_features["zcoord"],
    ) = check_user_defined_grid(coord_dfs, file_path)
    # Evaluate the number of elements from the checked number of nodes.
    conductor.grid_input["NELEMS"] = conductor.grid_features["N_nod"] - 1
    # Assign spatial discretization coordinates to each conductor component.
    for comp in conductor.inventory["all_component"].collection:
        assign_user_defined_spatial_discretization(comp, coord_dfs[comp.ID])
        check_grid_features(conductor.inputs["XLENGHT"], comp)


def check_grid_features(
    zlenght: float,
    comp: Union[
        FluidComponent,
        JacketComponent,
        StrandMixedComponent,
        StrandStabilizerComponent,
        StrandSuperconductorComponent,
    ],
):
    """Function that cheks initial anf final coordinates of the discretization to be consistent with the input values.

    Args:
        zlenght (float): straight length of the cable
        comp (Union[ FluidComponent, JacketComponent, StrandMixedComponent, StrandStabilizerComponent, StrandSuperconductorComponent, ]): generic object for which the spatial discretization is evaluated.

    Raises:
        ValueError: if first coordinate is smaller 0.0 (below a tolerance)
        ValueError: if first coordinate is lager than 0.0 (above a tolerance)
        ValueError: if last coordinate is smaller or larger than zlenght (within a tollerance).
    """

    tol = 1e-6
    if (comp.coordinate["z"][0] - 0.0) > -tol:
        message = f"{comp.ID = }: comp.coordinate['z'][0] must be 0.0; {comp.coordinate['z'][0]=} (m) < {0.0} (m)."
        logger_discretization.error(message)
        raise ValueError(message)
    if (comp.coordinate["z"][0] - 0.0) > tol:
        message = f"{comp.ID = }: comp.coordinate['z'][0] must be 0.0; {comp.coordinate['z'][0]=} (m) > {0.0} (m)."
        logger_discretization.error(message)
        raise ValueError(message)

    if abs(comp.coordinate["z"][-1] - zlenght) > tol:
        message = f"{comp.ID = }: comp.coordinate['z'][-1] does not equal zlenght!\n{comp.coordinate['z'][-1]=} m; {zlenght =} m"
        logger_discretization.warning(message)
        raise ValueError(message)


def check_user_defined_grid(
    dfs: dict, conductor: object, file_path: str
) -> "tuple[int, np.ndarray]":
    """Function that checks the consistency of the user defined spatial discretization for all user defined components.

    Args:
        dfs (dict): dictionary of dataframes, each dataframes has the coordinate of the corresponding coductor component.
        conductor (Conductor): conductor object.
        file_path (str): path of the input file with user defined spatial discretization.

    Raises:
        ValueError: if the number of sheets in used defined auxiliary input file differs from the total number of conductor components defined.
        ValueError: if less than three spatial coordinates are provided (x, y, z).
        ValueError: if the number of nodes is not the same in at least one sheet of the file. The reference number of nodes is inferred from the first sheet of the file.
        ValueError: if the z component of the spatial discretization is not exactly the same for FluidComponent and JacketComponent objects.

    Returns:
        tuple[int, np.ndarray]: total number of nodes of the user defined spatial discretization; z component of the conuctor spatial discretization: it is the z component of the spatial discretization of the first deifined FluidComponent object.
    """
    if len(dfs) != conductor.inventory["all_component"].number:
        raise ValueError(
            f"The number of sheets in file {file_path} must be equal to the number of defined conductor components.\n{len(dfs)} != {conductor.inventory['all_component'].number}.\n"
        )

    comp_ref = conductor.inventory["FluidComponent"].collection[0]
    n_node_ref = check_max_node_number(
        dfs[comp_ref.ID].shape[0], conductor, file_path, comp_ref.ID
    )
    z_ref = dfs[comp_ref.ID]["z [m]"].to_numpy()

    for comp in conductor.inventory["all_component"].collection[1:]:
        if dfs[comp.ID].shape[1] < 3:
            raise ValueError(
                f"User must provide three coordinates in sheeT {comp.ID} of input file {file_path}.\n"
            )
        if dfs[comp.ID].shape[0] != n_node_ref:
            raise ValueError(
                f"Inconsistent number of user defined nodes. The number of nodes in sheet {comp.ID} of file {file_path} must be equal to the one defined in sheet {comp_ref.ID} of the same file."
            )

        if not np.array_equal(z_ref, dfs[comp.ID]["z [m]"].to_numpy()):
            raise ValueError(
                f"User must provide the same z component of the coordinate for FluidComponent, JacketComponent and StrandComponent objects. Please check column z [m] in sheet {comp.ID} of file {file_path}."
            )

    return n_node_ref, z_ref


def check_max_node_number(
    n_nod: int, conductor: object, file_path: str, *args: str
) -> int:
    """Function that checks if the number of nodes of the spatial discretization is lower than the maximum allowed value.

    Args:
        n_nod (int): user defined number of nodes. Different oringin according to flag ITYMSH; see file_path for further detail.
        conductor (Conductor): conductor object with useful informations to make the check.
        file_path (str): path of the file from which the total number of nodes is evaluated.\n
        * If flag ITYMSH > 0 n_nod is evaluated from the number of elements assigned to conductor in input file conductor_grid;
        * If flag ITYMSH < 0 n_nod is set equal to the number of rows of the first sheet of user defined auxiliary input file used for the spatial discretization.

    Raises:
        ValueError: raises error message if the number of nodes excheedes the maximum allowed number of nodes.

    Returns:
        (int): number of nodes used for the spatial discretization.
    """
    if n_nod > conductor.grid_input["MAXNOD"]:
        if conductor.grid_input["ITYMSH"] > 0:
            message = f"The number of nodes should not exceed the maximum value. {n_nod = } > {conductor.grid_input['MAXNOD'] = }.\nPlease check {conductor.ID} in file {file_path}.\n"
        else:
            message = f"The number of nodes should not exceed the maximum value. {n_nod = } > {conductor.grid_input['MAXNOD'] = }.\nPlease check sheet {args[0]} in file {file_path}.\n"
        raise ValueError(message)
    return n_nod


def assign_user_defined_spatial_discretization(
    comp: Union[
        FluidComponent,
        JacketComponent,
        StrandMixedComponent,
        StrandStabilizerComponent,
        StrandSuperconductorComponent,
    ],
    df: pd.DataFrame,
):
    """Function that assigns the user defined coordinates in sheet comp.ID of auxiliary input file with the spatial discretization to the corresponding conductor component.

    Args:
        comp (Union[ FluidComponent, JacketComponent, StrandMixedComponent, StrandStabilizerComponent, StrandSuperconductorComponent, ]): generic object for which user defines spatial discretization.
        dfs (pd.DataFrame): dataframe with the coordinate of the corresponding coductor component.
    """
    comp.coordinate["x"] = df["x [m]"].to_numpy()
    comp.coordinate["y"] = df["y [m]"].to_numpy()
    comp.coordinate["z"] = df["z [m]"].to_numpy()


def build_coordinates_of_barycenter(cond:object, comp: Union[
        FluidComponent,
        JacketComponent,
        StrandMixedComponent,
        StrandStabilizerComponent,
        StrandSuperconductorComponent,
    ]):
    """Function that builds the coordinates of the barycenter of conductor components.

    Args:
        cond (object): conductor object, has all the information to evaluate load and assign user defined grid.
        comp (Union[ FluidComponent, JacketComponent, StrandMixedComponent, StrandStabilizerComponent, StrandSuperconductorComponent, ]): generic object for which the coordinates of the barycenter should be evaluated.

    Raises:
        ValueError: if costetha is not equal to 1 for FluidComponent and JacketComponent.
    """

    if comp.inputs["COSTETA"] == 1:
        comp.coordinate["x"] = comp.inputs["X_barycenter"] * np.ones(
            cond.grid_features["N_nod"]
        )
        comp.coordinate["y"] = comp.inputs["Y_barycenter"] * np.ones(
            cond.grid_features["N_nod"]
        )
        if (
            cond.grid_input["ITYMSH"] == 0
            or cond.grid_input["ITYMSH"] == 2
            or abs(cond.grid_input["XEREFI"] - cond.grid_input["XBREFI"]) <= 1e-3
        ):
            comp.coordinate["z"] = uniform_spatial_discretization(cond)
        elif cond.grid_input["ITYMSH"] == 1 or cond.grid_input["ITYMSH"] == 3:
            comp.coordinate["z"] = fixed_refined_spatial_discretization(
                cond, zcoord=np.zeros(cond.grid_features["N_nod"])
            )
    else:
        if issubclass(comp, StrandComponent):
            comp.cyl_helix = CylindricalHelix(
                comp.inputs["X_barycenter"],
                comp.inputs["Y_barycenter"],
                cond.inputs["XLENGHT"],
                comp.inputs["COSTHETA"],
            )
            if (
                cond.grid_input["ITYMSH"] == 0
                or cond.grid_input["ITYMSH"] == 2
                or abs(cond.grid_input["XEREFI"] - cond.grid_input["XBREFI"]) <= 1e-3
            ):
                tau = uniform_angular_discretization(cond, comp)
            elif cond.grid_input["ITYMSH"] == 1 or cond.grid_input["ITYMSH"] == 3:
                tau = fixed_refined_angular_discretization(cond, comp)
            (
                comp.coordinate["x"],
                comp.coordinate["y"],
                comp.coordinate["z"],
            ) = comp.cyl_helix.helix_parametrization(tau)
        else:
            raise ValueError(
                r"$Cos(\theta)$"
                + f"for {comp.__class__.__name__} must be 1.0; current value {comp.inputs['COSTETA'] = }\n"
            )
    
    check_grid_features(cond.inputs["XLENGHT"],comp)
