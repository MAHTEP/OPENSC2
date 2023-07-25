# Functions invoked to solve conductors transient (cdp, 07/2020)

import numpy as np
import os
import re
from scipy.linalg import solve_banded
from utility_functions.auxiliary_functions import (
    get_from_xlsx,
)
from typing import Union

from conductor import Conductor
from utility_functions.step_matrix_construction import (
    matrix_initialization,
    array_initialization,
    build_amat,
    build_transport_coefficients,
    build_kmat_fluid,
    build_smat_fluid,
    build_smat_fluid_interface,
    build_smat_fluid_solid_interface,
    build_mmat_solid,
    build_kmat_solid,
    build_smat_solid_interface,
    build_smat_env_solid_interface,
    build_svec,
    build_svec_env_jacket_interface,
    build_elmmat,
    build_elamat,
    build_elkmat,
    build_elsmat,
    build_elslod,
    assemble_matrix,
    assemble_syslod,
    eval_system_matrix,
    build_known_therm_vector,
)

def get_time_step(conductor, transient_input, num_step):

    """
    ##############################################################################
      SUBROUTINE GETSTP(TIME,TEND  ,STPMIN,STPMAX, &
      PRVSTP,OPTSTP,ICOND,NCOND) # cza added NCOND (August 14, 2017)
    ##############################################################################
    # Translation from Fortran to Python: Placido D. PoliTo, 07/08/2020
    ##############################################################################
    """

    TINY = 1.0e-10
    FACTUP = 1.2
    FACTLO = 0.5

    if num_step == 1:
        # at the first step time_step is equal to STPMIN for all the conductors \
        # (cdo, 08/2020)
        conductor.time_step = transient_input["STPMIN"]
    else:
        # STORE THE PREVIOUS VALUE OF THE OPTIMAL TIME STEP
        #      PRVSTP=OPTSTP
        PRVSTP = conductor.time_step
        if transient_input["IADAPTIME"] == 0:
            conductor.time_step = transient_input["STPMIN"]
            # C * LIMIT THE TIME STEP IF PRINT-OUT OR STORAGE IS REQUIRED
            conductor.time_step = min(
                conductor.time_step, transient_input["TEND"] - conductor.cond_time[-1]
            )  # crb (March 9, 2011)
            return

        # Ad hoc to emulate time adaptivity for simulation whith feeder CS3U2.
        # Do not use STPMIN; it is tuned on the AC loss time evolution.
        if transient_input["IADAPTIME"] == 3:
            # Array of times.
            times = np.array([0.0, transient_input["TIME_SHIFT"], transient_input["TIME_SHIFT"]+10.0-1.0, transient_input["TIME_SHIFT"] + 69.0, transient_input["TIME_SHIFT"] + 69.8, transient_input["TIME_SHIFT"] + 69.9, transient_input["TIME_SHIFT"]+ 70.0, transient_input["TIME_SHIFT"] + 70.1, transient_input["TIME_SHIFT"] + 70.15, transient_input["TIME_SHIFT"] + 70.5, transient_input["TIME_SHIFT"] + 71.0,transient_input["TIME_SHIFT"] + 73.0])
            # Adapt the time step according to the simulation time.
            if conductor.cond_time[-1] >= times[0] and conductor.cond_time[-1] <= times[1]:
                # Phase in which there is only radiative heat that stabilizes 
                # at time_shift; for the next 9 s there is no magnetic fiel or 
                # current or ac loss.
                conductor.time_step = 1.0 # s
            elif conductor.cond_time[-1] >= times[1] and conductor.cond_time[-1] < times[2]:
                # Some AC loss appears
                conductor.time_step = 1e-1 # s
            elif conductor.cond_time[-1] >= times[2] and conductor.cond_time[-1] < times[3]:
                # Some AC loss appears
                conductor.time_step = 1e-2 # s
            elif conductor.cond_time[-1] >= times[3] and conductor.cond_time[-1] < times[4]:
                conductor.time_step = 1e-3 # s
            elif conductor.cond_time[-1] >= times[4] and conductor.cond_time[-1] < times[5]:
                conductor.time_step = 1e-4 # s
            elif conductor.cond_time[-1] >= times[5] and conductor.cond_time[-1] < times[6]:
                # High and localized AC loss peack
                conductor.time_step = 1e-5 # s
            elif conductor.cond_time[-1] >= times[6] and conductor.cond_time[-1] < times[7]:
                conductor.time_step = 1e-4 # s
            elif conductor.cond_time[-1] >= times[7] and conductor.cond_time[-1] < times[8]:
                conductor.time_step = 1e-3 # s
            elif conductor.cond_time[-1] >= times[8] and conductor.cond_time[-1] < times[9]:
                conductor.time_step = 5e-3 # s
            elif conductor.cond_time[-1] >= times[9] and conductor.cond_time[-1] < times[10]:
                conductor.time_step = 1e-2 # s
            else:
                conductor.time_step = 5e-2 # s
                # C * LIMIT THE TIME STEP IF PRINT-OUT OR STORAGE IS REQUIRED
                conductor.time_step = min(
                    conductor.time_step, transient_input["TEND"] - conductor.cond_time[-1])
            return
        
        # crb Differentiate the indexes depending on ischannel (December 16, 2015)
        t_step_comp = np.zeros(conductor.dict_N_equation["NODOFS"])
        for ii in range(conductor.inventory["FluidComponent"].number):
            # FluidComponent objects (cdp, 08/2020)
            # C * THE FOLLOWING STATEMENTS WOULD CONTROL THE ACCURACY OF MOMENTUM...
            if abs(transient_input["IADAPTIME"]) == 1:  # crb (Jan 20, 2011)
                # (cdp, 08/2020)
                t_step_comp[ii] = conductor.EIGTIM / (conductor.EQTEIG[ii] + TINY)
                t_step_comp[
                    ii + conductor.inventory["FluidComponent"].number
                ] = conductor.EIGTIM / (
                    conductor.EQTEIG[
                        ii + conductor.inventory["FluidComponent"].number
                    ]
                    + TINY
                )
            elif transient_input["IADAPTIME"] == 2:
                # C * ... BUT ARE SUBSTITUTED BY THESE
                t_step_comp[ii] = 1.0e10
                t_step_comp[
                    ii + conductor.inventory["FluidComponent"].number
                ] = 1.0e10
            # endif (iadaptime)
            t_step_comp[
                ii + 2 * conductor.inventory["FluidComponent"].number
            ] = conductor.EIGTIM / (
                conductor.EQTEIG[
                    ii + 2 * conductor.inventory["FluidComponent"].number
                ]
                + TINY
            )
            # TSTPVH2 = 1.0e+10
            # TSTPPH2 = 1.0e+10
            # TSTPTH2 = 1.0e+10
            # crb (June 26, 2015) # cod (July 23, 2015)
        for ii in range(conductor.inventory["SolidComponent"].number):
            # SolidComponent objects (cdp, 08/2020)
            t_step_comp[
                ii + conductor.dict_N_equation["FluidComponent"]
            ] = conductor.EIGTIM / (
                conductor.EQTEIG[ii + conductor.dict_N_equation["FluidComponent"]]
                + TINY
            )

        # C * STORED THE OPTIMAL TIME STEP (FROM ACCURACY POINT OF VIEW)
        OPTSTP = min(t_step_comp)  # cod (July 23, 2015)

        # CL* CONTROL THE TIME STEPPING ALSO BY THE VOLTAGE
        # cl* add following lines
        # cl* August 17, 2000 - start *************************************
        # è associata almodulo elettrico,ignorare per ora
        ##if FFLAG_IOP > 0:
        ##	TSTEP_VOLT = EFIELD_INCR(ICOND)/conductor.time_step
        ##	TSTEP_VOLT = DBLE(EIGTIM(ICOND)/TSTEP_VOLT)
        ##	OPTSTP=min(TSTPVB,TSTPVH1,TSTPVH2,TSTPPH1,TSTPPH2,TSTPPB,&# cod (July 23, 2015)
        ##	 &             TSTPTH1,TSTPTH2,TSTPTB,TSTPCO,TSTPJK,TSTEP_VOLT)

        # cl* August 17, 2000 - end ***************************************

        # C * CHANGE THE STEP SMOOTHLY
        if conductor.time_step < 0.5 * OPTSTP:
            conductor.time_step = conductor.time_step * FACTUP
        elif conductor.time_step > 1.0 * OPTSTP:
            conductor.time_step = conductor.time_step * FACTLO
        # endif

        # C * LIMIT THE TIME STEP IN THE WINDOW ALLOWED BY THE USER
        conductor.time_step = max(conductor.time_step, transient_input["STPMIN"])
        # caf June 26, 2015 moved these 2 lines to the end of the subroutine
        # C * LIMIT THE TIME STEP IF PRINT-OUT OR STORAGE IS REQUIRED
        conductor.time_step = min(
            conductor.time_step, transient_input["TEND"] - conductor.cond_time[-1]
        )
        # C --------------
        print(f"Selected conductor time step is: {conductor.time_step}\n")
        # C --------------
        # caf end *********************************************** June 26, 2015

def step(conductor, environment, qsource, num_step):

    """
    ##############################################################################
        SUBROUTINE STEP(XCOORD,TMPTCO ,TMPTJK ,PRSSREH1,PRSSREH2,DENSTYH1,
       &                DENSTYH2,TMPTHEH1,TMPTHEH2,PRSSREB,DENSTYB,
       &                TMPTHEB ,VELCTH1 ,VELCTH2 ,VELCTB,BFIELD ,TCSHRE ,
       &                JHFLXC , JHFLXJ ,EXFLXC ,EXFLXJ, REYNOH1 ,REYNOH2,
       &                REYNOB ,FRICTH1 ,FRICTH2,FRICTB, HTCH1, HTCH2   ,
       &                HTCB,HTHB1,HTHB2, TIME     ,GAMMABH1,
       &				  GAMMABH2, NCHKCON, hcjw  ,heqw1,heqw2,htjb,htjh1,htjh2, # crb (January 21, 2018)
       &                icond  ,ncond,
       &				  T ,P ,R , MDOT  ,PIN    ,POUT  ,TIN     ,TOUT    ,EPSI,
       &                tb_prvstp,th1_prvstp,th2_prvstp,    # cab (March 02, 2018) Add all prvstp variables for correct fixed_point iterations
       &                pb_prvstp,ph1_prvstp,ph2_prvstp,
       &                vb_prvstp,vh1_prvstp,vh2_prvstp,
       &                tc_prvstp,tj_prvstp)
    ##############################################################################
    # cod 13/07/2015
    """

    path = os.path.join(conductor.BASE_PATH, conductor.file_input["EXTERNAL_FLOW"])
    TINY = 1.0e-5

    # CLUCA ADDNOD = MAXNOD*(ICOND-1)

    # Collection of valid dictionary keys
    basic_mat_names = ("MMAT","AMAT","KMAT","SMAT")
    element_mat_names = ("ELMMAT","ELAMAT","ELKMAT","ELSMAT")
    final_mat_names = ("MASMAT","FLXMAT","DIFMAT","SORMAT")

    # Final matrices initialization, collected in dictionary final_mat
    final_mat = matrix_initialization(
        conductor.dict_band["Full"],
        conductor.dict_N_equation["Total"],
        final_mat_names,
    )

    # Stiffness matrix initialization. Not included in dictionary final_mat to 
    # simplify the code below, make it explicit and clear to read and maintain.
    SYSMAT = np.zeros(
        (conductor.dict_band["Full"],conductor.dict_N_equation["Total"])
    )
    # SYSVAR = np.zeros(conductor.dict_N_equation["Total"])
    ASCALING = np.zeros(conductor.dict_N_equation["Total"])
    UPWEQT = np.zeros(conductor.dict_N_equation["NODOFS"])
    # Known terms vector initilaization
    Known = np.zeros_like(ASCALING)
    
    if conductor.inputs["METHOD"] == "BE" or conductor.inputs["METHOD"] == "CN":
        # Backward Euler or Crank-Nicolson (cdp, 10/2020)
        if conductor.cond_num_step > 1:
            # Copy the load vector at the previous time step in the second column to \
            # correctly apply the theta method (cdp, 10/2020)
            conductor.dict_Step["SYSLOD"][:, 1] = conductor.dict_Step["SYSLOD"][
                :, 0
            ].copy()
            conductor.dict_Step["SYSLOD"][:, 0] = 0.0

    # qsource initialization to zeros (cdp, 07/2020)
    # questa inizializzazione è provvisoria, da capire cosa succede quando ci \
    # sono più conduttori tra di loro a contatto o non in contatto.
    # if numObj == 1:
    # 	qsource = dict_qsource["single_conductor"]

    UPWEQT[:conductor.dict_N_equation["FluidComponent"]] = 1.0
    # if conductor.inputs["METHOD"] == "CN":
    # 		UPWEQT[2*conductor.inventory["FluidComponent"].number:conductor.dict_N_equation["FluidComponent"]] = 0.0
    # else it is 1.0 by initialization; as far as SolidComponent are \
    # concerned UPWEQT is always 0 by initialization. (cdp, 08/2020)

    # LUMPED/CONSISTENT MASS PARAMETER
    # if LMPMAS:
    # 	ALFA = 1.0/6.0
    # else:
    # 	ALFA = 0.0

    # practical case, I comment the above more general if-else, when needed \
    # remember that LMPMAS came from mandm.for subroutine SETNUM (cdp, 07/2020)
    ALFA = 0.0

    # COMPUTE AND ASSEMBLE THE ELEMENT NON-LINEAR MATRICES AND LOADS

    # Build transport coefficients K', K'' and K'''.
    conductor = build_transport_coefficients(conductor)

    # cl* * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    # cl* add the turn-to-turn coupling
    # qturn1 = interturn(nod1,zcoord,TMPTJK,nnodes(icond),icond)
    # qturn2 = interturn(nod2,zcoord,TMPTJK,nnodes(icond),icond)
    # cl* * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    # ** MATRICES CONSTRUCTION **

    # riscrivere in forma array smart una volta risolti tutti i dubbi, se possibile (cdp, 07/2020)
    for elem_index in range(conductor.grid_input["NELEMS"]):
        
        # Basic matrices initialization to zeros at each Gauss point, collected 
        # in dictionary basic_mat.
        basic_mat = matrix_initialization(
            conductor.dict_N_equation["NODOFS"],
            conductor.dict_N_equation["NODOFS"],
            basic_mat_names,
        )
        # Basic source term vector initialization to zeros at each Gauss point. 
        # Not included in dictionary base_mat to simplify the code, make it 
        # explicit and easy to read and maintain.
        SVEC = array_initialization(
            conductor.dict_N_equation["NODOFS"],
            conductor.cond_num_step,
            col=2,
        )

        # Element matrices initialization to zeros at each Gauss point, 
        # collected in dictionary element_mat.
        element_mat = matrix_initialization(
            conductor.dict_N_equation["NODOFS2"],
            conductor.dict_N_equation["NODOFS2"],
            element_mat_names,
        )
        # Element source term vector initialization to zeros at each Gauss 
        # point. Not included in dictionary element_mat to simplify the code, 
        # make it explicit and easy to read and maintain.
        ELSLOD = array_initialization(
            conductor.dict_N_equation["NODOFS2"],
            conductor.cond_num_step,
        )
        
        # ** FORM THE M, A, K, S MATRICES AND S VECTOR AT THE GAUSS POINT, 
        # FLUID COMPONENTS EQUATIONS **

        # FORM THE M MATRIX AT THE GAUSS POINT (MASS AND CAPACITY)
        # FluidComponent equation: array smart
        basic_mat["MMAT"][
            :conductor.dict_N_equation["FluidComponent"],
            :conductor.dict_N_equation["FluidComponent"],
        ] = np.eye(conductor.dict_N_equation["FluidComponent"])
        # END M MATRIX: fluid components equations
        
        for fluid_comp_j in conductor.inventory["FluidComponent"].collection:

            # FORM THE A MATRIX AT THE GAUSS POINT (FLUX JACOBIAN)
            basic_mat["AMAT"] = build_amat(
                basic_mat["AMAT"],
                fluid_comp_j,
                elem_index,
                conductor.equation_index[fluid_comp_j.identifier]
            )

            # FORM THE K MATRIX AT THE GAUSS POINT (INCLUDING UPWIND)
            basic_mat["KMAT"] = build_kmat_fluid(
                basic_mat["KMAT"],
                UPWEQT,
                fluid_comp_j,
                conductor,
                elem_index,
            )

            # FORM THE S MATRIX AT THE GAUSS POINT (SOURCE JACOBIAN)
            basic_mat["SMAT"] = build_smat_fluid(
                basic_mat["SMAT"],
                fluid_comp_j,
                elem_index,
                conductor.equation_index[fluid_comp_j.identifier]
            )

        # FORM THE S MATRIX AT THE GAUSS POINT (SOURCE JACOBIAN)
        # Therms associated to fluid-fluid interfaces.
        basic_mat["SMAT"] = build_smat_fluid_interface(
            basic_mat["SMAT"],
            conductor,
            elem_index
        )
        # Therms associated to fluid-solid interfaces.
        basic_mat["SMAT"] = build_smat_fluid_solid_interface(
            basic_mat["SMAT"],
            conductor,
            elem_index,
        )
        # END S MATRIX: fluid components equations

        # * FORM THE M, A, K, S MATRICES AND S VECTOR AT THE GAUSS POINT, SOLID
        # COMPONENTS EQUATIONS *
        for s_comp_idx, s_comp in enumerate(
            conductor.inventory["SolidComponent"].collection
        ):
            # FORM THE M MATRIX AT THE GAUSS POINT (MASS AND CAPACITY)
            # SolidComponent equation.
            basic_mat["MMAT"] = build_mmat_solid(
                basic_mat["MMAT"],
                s_comp,
                elem_index,
                conductor.equation_index[s_comp.identifier]
            )
            # END M MATRIX: SolidComponent equation.

            # FORM THE A MATRIX AT THE GAUSS POINT (FLUX JACOBIAN)
            # No elements here.
            # END A MATRIX: SolidComponent equation.

            # FORM THE K MATRIX AT THE GAUSS POINT (INCLUDING UPWIND)
            basic_mat["KMAT"] = build_kmat_solid(
                basic_mat["KMAT"],
                s_comp,
                elem_index,
                conductor.equation_index[s_comp.identifier]
            )
            # END K MATRIX: SolidComponent equation.

            # FORM THE S VECTOR AT THE NODAL POINTS (SOURCE)
            SVEC = build_svec(
                SVEC,
                s_comp,
                elem_index,
                conductor.equation_index[s_comp.identifier],
                num_step=conductor.cond_num_step,
                qsource=qsource,
                comp_idx=s_comp_idx,
            )

        # FORM THE S MATRIX AT THE GAUSS POINT (SOURCE JACOBIAN)
        basic_mat["SMAT"] = build_smat_solid_interface(
            basic_mat["SMAT"],
            conductor,
            elem_index,
        )

        for interface in conductor.interface.env_solid:
            # Convective heating with the external environment (implicit 
            # treatment).
            basic_mat["SMAT"] = build_smat_env_solid_interface(
                basic_mat["SMAT"],
                conductor,
                interface,
                elem_index,
            )
            # END S MATRIX: solid components equation.

            SVEC = build_svec_env_jacket_interface(
                SVEC,
                conductor,
                interface,
                elem_index,
            )
            # END S VECTOR: solid components equation.

        # COMPUTE THE MASS AND CAPACITY MATRIX
        # array smart
        element_mat["ELMMAT"] = build_elmmat(
            element_mat["ELMMAT"],
            basic_mat["MMAT"],
            conductor,
            elem_index,
            ALFA,
        )

        # COMPUTE THE CONVECTION MATRIX
        # array smart
        element_mat["ELAMAT"] = build_elamat(
            element_mat["ELAMAT"],
            basic_mat["AMAT"],
            conductor,
        )

        # COMPUTE THE DIFFUSION MATRIX
        # array smart
        element_mat["ELKMAT"] = build_elkmat(
            element_mat["ELKMAT"],
            basic_mat["KMAT"],
            conductor,
            elem_index,
        )

        # COMPUTE THE SOURCE MATRIX
        # array smart
        element_mat["ELSMAT"] = build_elsmat(
            element_mat["ELSMAT"],
            basic_mat["SMAT"],
            conductor,
            elem_index,
        )

        # COMPUTE THE SOURCE VECTOR (ANALYTIC INTEGRATION)
        # array smart
        ELSLOD = build_elslod(
            ELSLOD,
            SVEC,
            conductor,
            elem_index,
        )
        
        # ASSEMBLE THE MATRICES AND THE LOAD VECTOR
        
        jump = conductor.dict_N_equation["NODOFS"] * elem_index
        
        # array smart
        final_mat = assemble_matrix(
            final_mat,
            element_mat,
            conductor,
            jump,
        )

        conductor.dict_Step["SYSLOD"][
            jump:jump + conductor.dict_band["Half"],:
        ] = assemble_syslod(ELSLOD,conductor,jump)

    # end for elem_index
    # ** END MATRICES CONSTRUCTION **

    # ** COMPUTE SYSTEM MATRIX **
    SYSMAT = eval_system_matrix(
        SYSMAT,
        final_mat,
        conductor,
    )

    # ADD THE LOAD CONTRIBUTION FROM PREVIOUS STEP
    # array smart
    Known = build_known_therm_vector(
        Known,
        final_mat,
        conductor,
    )

    # Call function to save ndarrays before the application of BC; this can be 
    # useful to create and make tests according to TDD approach.
    #   ** Decomment the following lines to save ndarrays **
    # save_ndarray(
    #     conductor,
    #     (
    #         final_mat["MASMAT"],final_mat["FLXMAT"],final_mat["DIFMAT"],
    #         final_mat["SORMAT"],SYSMAT,conductor.dict_Step["SYSVAR"],
    #         conductor.dict_Step["SYSLOD"]
    #     )
    # )

    # IMPOSE BOUNDARY CONDITIONS AT INLET/OUTLET
    for f_comp in conductor.inventory["FluidComponent"].collection:

        intial = abs(f_comp.coolant.operations["INTIAL"])
        # Apply boundary conditions according to the absloute value of flag
        # INTIAL.
        Known,SYSMAT = f_comp.apply_th_bc[intial](
            (Known,SYSMAT),
            conductor,
            path,
        )

    # DIAGONAL ROW SCALING

    # SELECT THE MAX FOR EACH ROW
    ASCALING = abs(SYSMAT).max(0)
    ind_ASCALING = np.nonzero(ASCALING == 0.0)
    # ind_ASCALING = np.nonzero(ASCALING <= 1e-6)
    if ind_ASCALING[0].shape == 0:
        raise ValueError(
            f"""ERROR: ASCALING[ind_ASCALING[0]] = 
           {ASCALING[ind_ASCALING[0]]}!\n"""
        )

    # SCALE THE SYSTEM MATRIX
    SYSMAT = SYSMAT / ASCALING

    # SCALE THE LOAD VECTOR
    Known = Known / ASCALING

    old_temperature_gauss = {
        obj.identifier: obj.coolant.dict_Gauss_pt["temperature"]
        for obj in conductor.inventory["FluidComponent"].collection
    }
    old_temperature_gauss.update(
        {
            obj.identifier: obj.dict_Gauss_pt["temperature"]
            for obj in conductor.inventory["SolidComponent"].collection
        }
    )

    SYSMAT = gredub(conductor, SYSMAT)
    # Compute the solution at current time stepand overwrite key SYSVAR of \
    # dict_Step
    conductor.dict_Step["SYSVAR"][:, 0] = gbacsb(conductor, SYSMAT, Known)

    # SYSVAR = solve_banded((15, 15), SYSMAT, Known)

    # COMPUTE THE NORM OF THE SOLUTION AND OF THE SOLUTION CHANGE (START)
    # array smart optimization

    CHG = np.zeros(conductor.dict_N_equation["Total"])
    EIG = np.zeros(conductor.dict_N_equation["Total"])

    # Evaluate the norm of the solution.
    conductor.dict_norm["Solution"] = eval_sub_array_norm(Known,conductor)

    # COMPUTE THE NORM OF THE SOLUTION CHANGE, THE EIGENVALUES AND RECOVER THE \
    # VARIABLES FROM THE SYSTEM SOLUTION (START)

    # Those are arrays
    # Solution change
    CHG = Known - conductor.dict_Step["SYSVAR"][:, 0]
    # Eigenvalues (sort of??)
    EIG = abs(CHG / conductor.time_step) / (abs(Known) + TINY)
    # Evaluate the norm of the solution change.
    conductor.dict_norm["Change"] = eval_sub_array_norm(CHG,conductor)
    # Evaluate the eigenvalues of the solution.
    conductor.EQTEIG = eval_eigenvalues(EIG,conductor)
    # Reorganize thermal hydraulic solution
    reorganize_th_solution(
        conductor,
        old_temperature_gauss,
    )
    
    # COMPUTE THE NORM OF THE SOLUTION CHANGE, THE EIGENVALUES AND RECOVER THE \
    # VARIABLES FROM THE SYSTEM SOLUTION (END)


def gredub(conductor, A):

    """
    ##############################################################################
    #    SUBROUTINE GREDUB(A     ,dict_N_equation["Total"],conductor.dict_band["Main_diag"] ,IERR  )
    ##############################################################################
    #
    # REDUCTION OF A NON-SINGULAR SYSTEM OF EQUATIONS. COEFFICIENTS IN A
    # STORED ONLY WITHIN THE BANDWIDTH conductor.dict_band["Main_diag"] AS FOLLOWS:
    # A(I-conductor.dict_band["Main_diag"],I),A(I-conductor.dict_band["Main_diag"]+1,I)...A(I,I)...A(I+conductor.dict_band["Main_diag"]-1,I),A(I+conductor.dict_band["Main_diag"]+1,I)
    #
    #   IERR =  0 IF NO ERROR IS DETECTED
    #   IERR =  1 THE VALUE OF THE HALF-BANDWIDTH GIVEN IS NOT CONSISTENT
    #             WITH THE MATRIX DIMENSION
    #   IERR <  0 MATRIX SINGULAR AT LINE -IERR
    #
    ##############################################################################
    #
    #  THIS VERSION USES THE DOUBLE PRECISION
    #
    ##############################################################################
    #
    # Translation form Frortran to Python: D.Placido PoliTo, 23/08/2020
    # Array smart optimization: D.Placido PoliTo, 27/08/2020
    #
    ##############################################################################
    """

    TINY = 1.0e-20
    IERR = 0
    if (
        conductor.dict_band["Main_diag"] < 0
        or conductor.dict_band["Main_diag"] > conductor.dict_N_equation["Total"]
    ):
        IERR = 1
        if conductor.dict_band["Main_diag"] < 0:
            raise ValueError(
                f"""ERROR! The value of the half-band width given is not 
      consistent with the matrix dimension:\n
      {conductor.dict_band["Main_diag"]} < 0\nReturned error: IERR = {IERR}"""
            )
        elif conductor.dict_band["Main_diag"] > conductor.dict_N_equation["Total"]:
            raise ValueError(
                f"""ERROR! The value of the half-band width given is not 
      consistent with the matrix dimension:\n
      {conductor.dict_band["Main_diag"]} > {conductor.dict_N_equation["Total"]}\nReturned error: IERR = {IERR}
      \nEnd of the program\n"""
            )

    K1 = conductor.dict_band["Main_diag"]  # 15
    K2 = conductor.dict_band["Main_diag"] + 1  # 16
    K21 = 2 * conductor.dict_band["Main_diag"]  # 30
    JJ = K21
    II = conductor.dict_N_equation["Total"] - conductor.dict_band["Main_diag"]
    for I in range(II, conductor.dict_N_equation["Total"]):
        A[JJ : K21 + 1, I] = np.zeros(K21 - JJ + 1)
        JJ = JJ - 1

    # Evaluate J1 and JK only once since they have always the same value \
    # (cdp, 08/2020)
    # remember that the stop value is not included (cdp, 08/2020)
    J1 = np.arange(start=1, stop=K2, step=1, dtype=int)
    JK = np.arange(start=conductor.dict_band["Main_diag"], stop=K21, step=1, dtype=int)

    for I in range(1, conductor.dict_N_equation["Total"]):
        # remember that the stop value is not included (cdp, 08/2020)
        II = np.arange(
            start=I - conductor.dict_band["Main_diag"], stop=I, step=1, dtype=int
        )
        # thake only positive values (cdp, 08/2020)
        II = II[II >= 0]
        # this is an array (cdp, 08/2020)
        ind1 = np.nonzero(abs(A[K1, II]) <= TINY)[0]
        # check if matrix is singular (cdp, 08/2020)
        if len(ind1) > 0:
            IERR = -II[ind1[0]]
            raise ValueError(
                f"""ERROR! Matrix is singular at line {II[ind1[0]]}: 
      {A[K1, II[ind1[0]]]} < {TINY}\nReturned error: IERR = {IERR}\n
      End of the program\n"""
            )

        J_min = conductor.dict_band["Main_diag"] - len(II)
        Q = np.zeros(II.shape)

        # da ragionare
        for ii in range(len(II)):
            Q[ii] = A[J_min + ii, I] / A[K1, II[ii]]
            A[J1[J_min + ii] : JK[J_min + ii] + 1, I] = (
                A[J1[J_min + ii] : JK[J_min + ii] + 1, I]
                - A[K2 : K21 + 1, II[ii]] * Q[ii]
            )
        # end for ii
        A[J_min : conductor.dict_band["Main_diag"], I] = Q
    # end for I
    return A
    # end of the function GREDUB


def gbacsb(conductor, A, B):

    """
    ##############################################################################
        SUBROUTINE GBACSB (A, B, X IERR)
    ##############################################################################
    #
    # BACK SUBSTITUION AND SOLUTION OF THE SYSTEM A X = B. THE MATRIX
    # A HAS BEEN REDUCED BY GREDUB. X AND B CAN BE COINCIDENT
    #
    #   IERR =  0 IF NO ERROR IS DETECTED
    #   IERR =  1 THE VALUE OF THE HALF-BANDWIDTH GIVEN IS NOT CONSISTENT
    #             WITH THE MATRIX DIMENSION
    #   IERR <  0 MATRIX SINGULAR AT LINE -IERR
    #
    ##############################################################################
    #
    #  THIS VERSION USES THE DOUBLE PRECISION
    #
    ##############################################################################
    #
    # Translation form Frortran to Python: D.Placido PoliTo, 23/08/2020
    # Array smart optimization: D.Placido PoliTo, 27/08/2020
    #
    ##############################################################################
    """

    TINY = 1.0e-20
    IERR = 0
    if (
        conductor.dict_band["Main_diag"] < 0
        or conductor.dict_band["Main_diag"] > conductor.dict_N_equation["Total"]
    ):
        IERR = 1
        if conductor.dict_band["Main_diag"] < 0:
            raise ValueError(
                f"""ERROR! The value of the half-band width given is not 
      consistent with the matrix dimension:\n
      {conductor.dict_band["Main_diag"]} < 0\nReturned error: IERR = {IERR}\n
      End of the program\n"""
            )
        elif conductor.dict_band["Main_diag"] > conductor.dict_N_equation["Total"]:
            raise ValueError(
                f"""ERROR! The value of the half-band width given is not 
      consistent with the matrix dimension:\n
      {conductor.dict_band["Main_diag"]} > {conductor.dict_N_equation["Total"]}\nReturned error: IERR = {IERR}
      \nEnd of the program\n"""
            )

    K1 = conductor.dict_band["Main_diag"]
    K2 = conductor.dict_band["Main_diag"] + 1
    for I in range(1, conductor.dict_N_equation["Total"]):
        # remember that the stop value is not included (cdp, 08/2020)
        II = np.arange(
            start=I - conductor.dict_band["Main_diag"], stop=I, step=1, dtype=int
        )
        # thake only positive values (cdp, 08/2020)
        II = II[II >= 0]
        J_min = conductor.dict_band["Main_diag"] - len(II)
        for ii in range(len(II)):
            B[I] = B[I] - B[II[ii]] * A[J_min + ii, I]

    if abs(A[K1, conductor.dict_N_equation["Total"] - 1]) <= TINY:
        IERR = -conductor.dict_N_equation["Total"] - 1
        raise ValueError(
            f"""ERROR! Matrix is singular at line 
    {conductor.dict_N_equation["Total"] - 1}: {A[K1, conductor.dict_N_equation["Total"] - 1]} < {TINY}\n
    Returned error: IERR = {IERR}\n
    End of the program\n"""
        )

    B[-1] = B[-1] / A[K1, -1]
    # Index array (cdp, 08/2020)
    II = np.arange(
        start=conductor.dict_N_equation["Total"] - 2, stop=-1, step=-1, dtype=int
    )

    # this is an array (cdp, 08/2020)
    ind = np.nonzero(abs(A[K1, II]) <= TINY)[0]
    if len(ind) > 0:
        IERR = -II[ind]
        raise ValueError(
            f"""ERROR! Matrix is singular at line {II[ind]}: 
    {A[K1, II[ind]]} < {TINY}\nReturned error: IERR = {IERR}\n
    End of the program\n"""
        )

    # Index array (cdp, 08/2020)
    JJ = II + conductor.dict_band["Main_diag"]
    JJ[JJ >= conductor.dict_N_equation["Total"]] = (
        conductor.dict_N_equation["Total"] - 1
    )
    # Index array (cdp, 08/2020)
    II1 = II + 1
    L_max = JJ - II1 + 1

    Q = B[II]

    for ii in range(len(II)):
        for jj in range(L_max[ii]):
            Q[ii] = Q[ii] - A[K2 + jj, II[ii]] * B[II1[ii] + jj]
        B[II[ii]] = Q[ii] / A[K1, II[ii]]

    X = np.zeros(conductor.dict_N_equation["Total"])
    X = B

    return X
    # end of the function GBACSB


def natural_sort(comp_a, comp_b):
    # Use the regexes to sort naturally (human like) the IDs of the components to be able to deal with all the interfaces in a general way.
    match_a = re.search(
        r"(?P<Fluid_component>CHAN)?(?P<Stack>STACK)?(?P<Mixed_sc_stab>STR_MIX)??(?P<StrandStabilizerComponent>STR_STAB)?(?P<JacketComponent>Z_JACKET)?_(\d+)",
        comp_a.identifier,
    )
    # r'((CHAN)?(STACK)?(STR_MIX)?(STR_STAB)?(Z_JACKET)?)_(\d)+'
    match_b = re.search(
        r"(?P<Fluid_component>CHAN)?(?P<Stack>STACK)?(?P<Mixed_sc_stab>STR_MIX)?(?P<StrandStabilizerComponent>STR_STAB)?(?P<JacketComponent>Z_JACKET)?_(\d+)",
        comp_b.identifier,
    )

    if match_a.group(comp_a.KIND) < match_b.group(comp_b.KIND):
        return f"{comp_a.identifier}_{comp_b.identifier}"
    elif match_a.group(comp_a.KIND) > match_b.group(comp_b.KIND):
        return f"{comp_b.identifier}_{comp_a.identifier}"
    else:
        # Equal string part, sort by number
        if int(match_a.group(6)) < int(match_b.group(6)):
            return f"{comp_a.identifier}_{comp_b.identifier}"
        else:
            return f"{comp_b.identifier}_{comp_a.identifier}"
        # end if
    # end if

def eval_sub_array_norm(
    array:np.ndarray,
    conductor:Conductor,
    )->np.ndarray:
    """Function that evaluates the euclidean norm of as many sub arrays as the number of unknowns of the thermal hydraulic problem stored insde input argument array. Being jj the j-th unknown (i.e. CHAN_1 temperature), the sub array is given by sub_arr = array[jj::ndf] if ndf is the number of unknowns (number of degrees of freedom). The euclidean norm is applied to this sub array. The final outcome is an array of eucliean norms with ndf elements.
    This function is used both to evaluate the norm of the solution and the norm of the solution change.

    Args:
        array (np.ndarray): array containing ndf sub arrays (each being the current thermal hydraulic solution or its change wrt the previous solution).
        conductor (Conductor): object with all the information of the conductor.

    Returns:
        np.ndarray: array of eucliean norms with ndf elements.
    """

    # Alias
    ndf = conductor.dict_N_equation["NODOFS"]
    sub_array_norm = np.zeros(ndf)
    # Collection of NamedTuple with fluid equation index (velocity, pressure 
    # and temperaure equations) and of integer for solid equation index.
    eq_idx = conductor.equation_index
    # Evaluate the square of the array.
    array **= 2.0

    # Evaluate the sub arrays euclidean norm.
    # Loop on FluidComponent.
    for f_comp in conductor.inventory["FluidComponent"].collection:
        # velocity
        sub_array_norm[eq_idx[f_comp.identifier].velocity] = np.sum(
            array[eq_idx[f_comp.identifier].velocity::ndf]
        )
        # pressure
        sub_array_norm[eq_idx[f_comp.identifier].pressure] = np.sum(
            array[eq_idx[f_comp.identifier].pressure::ndf]
        )
        # temperature
        sub_array_norm[eq_idx[f_comp.identifier].temperature] = np.sum(
            array[eq_idx[f_comp.identifier].temperature::ndf]
        )
    # Loop on SolidComponent.
    for s_comp in conductor.inventory["SolidComponent"].collection:
        # temperature
        sub_array_norm[eq_idx[s_comp.identifier]] = np.sum(
            array[eq_idx[s_comp.identifier]::ndf]
        )
    
    return np.sqrt(sub_array_norm)

def eval_eigenvalues(
    array:np.ndarray,
    conductor:Conductor,
    )->np.ndarray:
    """
    Function that evaluate an approximation of the eigenvalues of the solution of as many sub arrays as the number of unknowns of the thermal hydraulic problem stored insde input argument array. Being jj the j-th unknown (i.e. CHAN_1 temperature), the sub array is given by sub_arr = array[jj::ndf] if ndf is the number of unknowns (number of degrees of freedom). The eigenvalue is the maximum value of this sub array. The final outcome is an array of eigenvalues with ndf elements.

    Args:
        array (np.ndarray): array containing ndf sub arrays (each being an approximation of the eigenvalues of the thermal hydraulic solution).
        conductor (Conductor): object with all the information of the conductor.
    
    Returns:
        np.ndarray: array of eigenvalues with ndf elements.
    """

    # Alias
    ndf = conductor.dict_N_equation["NODOFS"]
    sub_array = np.zeros(ndf)
    # Collection of NamedTuple with fluid equation index (velocity, pressure 
    # and temperaure equations) and of integer for solid equation index.
    eq_idx = conductor.equation_index
    # COMPUTE THE EIGENVALUES
    for f_comp in conductor.inventory["FluidComponent"].collection:
        # velocity
        sub_array[eq_idx[f_comp.identifier].velocity] = max(
            array[eq_idx[f_comp.identifier].velocity::ndf]
        )
        # pressure
        sub_array[eq_idx[f_comp.identifier].pressure] = max(
            array[eq_idx[f_comp.identifier].velocity::ndf]
        )
        # temperature
        sub_array[eq_idx[f_comp.identifier].temperature] = max(
            array[eq_idx[f_comp.identifier].temperature::ndf]
        )
    # Loop on SolidComponent.
    for s_comp in conductor.inventory["SolidComponent"].collection:
        # temperature
        sub_array[eq_idx[s_comp.identifier]] = max(
            array[eq_idx[s_comp.identifier]::ndf]
        )
    
    return sub_array


def eval_array_by_fn(
    array:np.ndarray,
    conductor:Conductor,
    fn: Union[np.sum,np.max],
    )->np.ndarray:
    """Function that evaluates the euclidean norm of as many sub arrays as the number of unknowns of the thermal hydraulic problem stored insde input argument array -if fn is np.sum- or an approximation of the eigenvalues of the solution of as many sub arrays as the number of unknowns of the thermal hydraulic problem stored inside input argument array -if fn is np.sum-.

    Being jj the j-th unknown (i.e. CHAN_1 temperature), the sub array is given by sub_arr = array[jj::ndf] if ndf is the number of unknowns (number of degrees of freedom).

    If fn is np.sum, the euclidean norm is applied to this sub array. The final outcome is an array of eucliean norms with ndf elements. This can be used both to evaluate the norm of the solution and the norm of the solution change.

    If fn is np.max, the eigenvalue is the maximum value of this sub array. The final outcome is an array of eigenvalues with ndf elements.

    Args:
        array (np.ndarray): array containing ndf sub arrays (each being the current thermal hydraulic solution or its change wrt the previous solution -if fn is np.sum- or an approximation of the eigenvalues of the thermal hydraulic solution -if fn is np.max-).
        conductor (Conductor): object with all the information of the conductor.
        fn (Union[np.sum,np.max]): aggregation function (namely np.sum if euclidean norm should be evaluated, np.max if eigenvalues should be evaluated).

    Returns:
        np.ndarray: array of eucliean norms with ndf elements if fn is np.sum; array of eigenvalues with ndf elements if fn is np.max.
    """
    # Alias
    ndf = conductor.dict_N_equation["NODOFS"]
    sub_array = np.zeros(ndf)
    # Collection of NamedTuple with fluid equation index (velocity, pressure 
    # and temperaure equations) and of integer for solid equation index.
    eq_idx = conductor.equation_index

    # Exponent 2 is to evaluate the euclidean norm, exponent 1 is to evaluate 
    # the eigenvalue (the values does not change is powered to 1).
    pow_exp = {
        np.sum: 2.0,
        np.max: 1.0,
    }

    # Evaluate the sub arrays euclidean norm.
    # Loop on FluidComponent.
    for f_comp in conductor.inventory["FluidComponent"].collection:
        
        # Power velocity to 2 to evaluate the euclidean norm, to 1 to evaluate 
        # the eigenvalues.
        vv = array[eq_idx[f_comp.identifier].velocity::ndf] ** pow_exp[fn]
        # Power pressure to 2 to evaluate the euclidean norm, to 1 to evaluate 
        # the eigenvalues.
        pp = array[eq_idx[f_comp.identifier].pressure::ndf] ** pow_exp[fn]
        # Power temperature to 2 to evaluate the euclidean norm, to 1 to 
        # evaluate the eigenvalues.
        tt = array[eq_idx[f_comp.identifier].temperature::ndf] ** pow_exp[fn]
    
        # velocity
        sub_array[eq_idx[f_comp.identifier].velocity] = fn(vv)
        # pressure
        sub_array[eq_idx[f_comp.identifier].pressure] = fn(pp)
        # temperature
        sub_array[eq_idx[f_comp.identifier].temperature] = fn(tt)
    # Loop on SolidComponent.
    for s_comp in conductor.inventory["SolidComponent"].collection:
        # Power temperature to 2 to evaluate the euclidean norm, to 1 to 
        # evaluate the eigenvalues.
        tt = array[eq_idx[s_comp.identifier]::ndf] ** pow_exp[fn]
        # temperature
        sub_array[eq_idx[s_comp.identifier]] = fn(tt)
    
    if fn is np.sum:
        # Return the euclidean norm.
        return np.sqrt(sub_array)
    elif fn is np.max:
        # Return the eigenvalues.
        return sub_array

def reorganize_th_solution(
    conductor:Conductor,
    old_temperature:dict,
    ):
    """
    Function that reorganizes the thermal hydraulic solution into ndf arrays if ndf is the number of unknowns (number of degrees of freedom) according to the following rationale:
        for each FluidComponent object
            * velocity
            * pressure
            * temperature
        for each SolidComponent object
            * temperature
    
    Sub arrays (i.e. CHAN_1 temperature spatial distribution) are given by sub_arr = array[jj::ndf] if jj is the index of the j-th conductor component object (i.e. CHAN_1).

    Attribute f_comp.coolant.dict_node_pt (that stores fluid properties in nodal points) and s_comp.dict_node_pt (that stores solid properties) are updated inplace.

    Args:
        conductor (Conductor): object with all the information of the conductor.
        old_temperature (dict): collection of arrays of the temperature distribution at the previous time step for each conductor component.
    """

    # Alias
    ndf = conductor.dict_N_equation["NODOFS"]
    sysvar = conductor.dict_Step["SYSVAR"]
    # Collection of NamedTuple with fluid equation index (velocity, pressure 
    # and temperaure equations) and of integer for solid equation index.
    eq_idx = conductor.equation_index
    # Reorganize thermal hydraulic solution.
    for f_comp in conductor.inventory["FluidComponent"].collection:
        # velocity
        f_comp.coolant.dict_node_pt["velocity"] = sysvar[
            eq_idx[f_comp.identifier].velocity::ndf,0
        ].copy()
        # pressure
        f_comp.coolant.dict_node_pt["pressure"] = sysvar[
            eq_idx[f_comp.identifier].pressure::ndf,0
        ].copy()
        # temperature
        f_comp.coolant.dict_node_pt["temperature"] = sysvar[
            eq_idx[f_comp.identifier].temperature::ndf,0
        ].copy()
        # Get temperature change in Gauss points.
        f_comp.coolant.dict_Gauss_pt["temperature_change"] = (
            f_comp.coolant.dict_node_pt["temperature"][:-1]
            + f_comp.coolant.dict_node_pt["temperature"][1:]
        ) / 2. - old_temperature[f_comp.identifier]
    # Loop on SolidComponent.
    for s_comp in conductor.inventory["SolidComponent"].collection:
        # temperature
        s_comp.dict_node_pt["temperature"] = sysvar[
            eq_idx[s_comp.identifier]::ndf,0
        ].copy()
        # Get temperature change in Gauss points.
        s_comp.dict_Gauss_pt["temperature_change"] = (
            s_comp.dict_node_pt["temperature"][:-1]
            + s_comp.dict_node_pt["temperature"][1:]
        ) / 2. - old_temperature[s_comp.identifier]

def save_ndarray(conductor:Conductor,ndarrays:tuple):
    """Function that saves selected ndarrays to perform test according to the TDD approach. The path to the directory where to save the ndarrays is hardcoded in the function. This is not optimum but it is not trivial to do it in a different way.

    Args:
        conductor (Conductor): object with all the information of the conductor.
        ndarrays (tuple): collection of ndarrays in the following order: MASMAT, FLXMAT, DIFMAT, SORMAT, SYSMAT, SYSVAR, SYSLOD
    """
    
    if conductor.cond_num_step == 1 or np.isclose(conductor.Space_save[conductor.i_save],conductor.cond_time[-1]):

        path_matr = os.path.join("D:/refactoring/function_step", "matrices/before")
        os.makedirs(path_matr, exist_ok = True)

        if conductor.cond_num_step == 1:
            sfx = conductor.i_save - 1
        else:
            sfx = conductor.i_save
        # Collection of ndimensional array.
        nda_name = (
            "MASMAT","FLXMAT","DIFMAT","SORMAT","SYSMAT","SYSVAR","SYSLOD"
        )
        # Build path to save ndimensional array.
        paths = tuple(
            os.path.join(path_matr, f"{name}_{sfx}.tsv") for name in nda_name
        )
        # Loop to save ndimensional array at the given path
        for save_path,nda in zip(paths,ndarrays):
            with open(save_path, "w") as writer:
                np.savetxt(writer, nda, delimiter = "\t")
