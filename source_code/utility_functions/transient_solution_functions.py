# Functions invoked to solve conductors transient (cdp, 07/2020)

import numpy as np
import os
import re
from scipy.linalg import solve_banded
from utility_functions.auxiliary_functions import (
    get_from_xlsx,
    filter_component,
)
from collections import namedtuple
from typing import Union, NamedTuple

from fluid_component import FluidComponent
from conductor import Conductor

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

def matrix_initialization(row:int,col:int)->tuple:
    """Wrapper of function np.zeros that inizializes five identical rectangular matrices.

    Args:
        row (int): number of rows of the matrix.
        col (int): number of columns of the matrix.


    Returns:
        tuple: collection of initialized matrices.
    """
    
    matrix = np.zeros((row,col))
    return (
        matrix,
        matrix.copy(), # Use method copy to get hard copy of the matrix.
        matrix.copy(),
        matrix.copy(),
        matrix.copy(),
    )

def ndarray_initialization(dimension:int,num_step:int,col:int=0)->tuple:
    """Wrapper of function np.zeros that inizializes four identical square matrices and an additional variable that can be a column vector or a tuple with matrix of shape (dimension,col) according to the value of num_step (this is necessary to correctly apply the theta method).
    N.B. the application of the theta method should be completely rivisited in the whole code.

    Args:
        dimension (int): number of rows/columns of the matrices to be initialized and number of rows in the additional variable.
        num_step (int): time step number.
        col (int, optional): number of columns to be assigned to the additional variable. If col is 0, the array shape is (dimension,), else array shape is (dimension,col). Defaults to 0.

    Returns:
        tuple: collection of initialized ndarrays.
    """

    matrix = np.zeros((dimension,dimension))
    return (
        matrix,
        matrix.copy(), # Use method copy to get hard copy of the matrix.
        matrix.copy(),
        matrix.copy(),
        array_initialization(dimension, num_step,col)
    )

def array_initialization(dimension:int, num_step:int, col:int=0)-> Union[NamedTuple,np.ndarray]:
    """Wrapper of function np.zeros that initializes array of shape (dimension, col) according to the time step number.
    N.B. the application of the theta method should be completely rivisited in the whole code.

    Args:
        dimension (int): number of elements (rows) of the array to be initialized.
        num_step (int): time step number.
        col (int, optional): number of columns to be assigned to the array. If col is 0, the array shape is (dimension,), else array shape is (dimension,col). Defaults to 0.

    Returns:
        Union[NamedTuple,np.ndarray]: namedtuple with array if num_step is 1; np.ndarray in all other cases.
    """

    if num_step == 1:
        Array = namedtuple("Array",["previous","present"])
        # To correctly apply the theta method (to be rivisited in the whole 
        # code!).
        # previous is for the initialization (time step number is 0);
        # present is for the first time step after the initialization
        
        # Check on col to assign the correct shape to the array.
        if col: # used to define SVEC.
            return Array(
                previous=np.zeros((dimension, col)),
                present=np.zeros((dimension, col)),
            )
        else: # used to define ELSLOD.
            return Array(
                previous=np.zeros(dimension),
                present=np.zeros(dimension),
            )
    else:
        # Check on col to assign the correct shape to the array.
        if col: # used to define SVEC.
            return np.zeros((dimension, col))
        else: # used to define ELSLOD.
            return np.zeros(dimension)

def __build_fluid_eq_idx(conductor:Conductor)->dict:
    """Function that evaluates the index of the velocity, pressure and temperature equation of the FludiComponent objects, collecting them in a dictionary of a namedtuple.

    Args:
        conductor (Conductor): object with all the information of the conductor.

    Returns:
        dict: collection of namedtuple with the index of velocity, pressure and temperaure equation for FludiComponent objects.
    """
    
    # Constructor of the namedtuple to store the index of the equations for 
    # FludiComponent objects.
    Fluid_eq_idx = namedtuple(
        "Fluid_eq_idx",
        ["velocity","pressure","temperature"]
    )

    # Build dictionary of NamedTuple with the index of the equations for 
    # FludiComponent objects exploiting dictionary comprehension.
    return {
        fcomp.identifier:Fluid_eq_idx(
            # velocity equation index
            velocity=fcomp_idx,
            # pressure equation index
            pressure=fcomp_idx + conductor.inventory["FluidComponent"].number,
            # temperature equation index
            temperature=fcomp_idx + 2 * conductor.inventory["FluidComponent"].number
        )
        for fcomp_idx,fcomp in enumerate(
            conductor.inventory["FluidComponent"].collection
        )
    }


def __build_amat(
    matrix:np.ndarray,
    f_comp:FluidComponent,
    elem_idx:int,
    eq_idx:NamedTuple,
    )->np.ndarray:
    """Function that builds the A matrix (AMAT) at the Gauss point (flux Jacobian).

    Args:
        matrix (np.ndarray): initialized A matrix (np.zeros)
        f_comp (FluidComponent): fluid component object from which get all info to buld the coefficients.
        elem_idx (int): index of the i-th element of the spatial discretization.
        eq_idx (NamedTuple): collection of fluid equation index (velocity, pressure and temperaure equations).

    Returns:
        np.ndarray: matrix with updated elements.
    """
    
    # Build array to assign diagonal coefficients.
    diag_idx = np.array(eq_idx)
    
    # Set diagonal elements (exploit broadcasting).
    matrix[diag_idx,diag_idx] = f_comp.coolant.dict_Gauss_pt["velocity"][elem_idx]

    # Set off diagonal coefficients.
    # from velocity equation.
    matrix[eq_idx.velocity, eq_idx.pressure] = (
        1 / f_comp.coolant.dict_Gauss_pt["total_density"][elem_idx]
    )
    # from pressure equation.
    matrix[eq_idx.pressure, eq_idx.velocity] = (
        f_comp.coolant.dict_Gauss_pt["total_speed_of_sound"][elem_idx] ** 2
        * f_comp.coolant.dict_Gauss_pt["total_density"][elem_idx]
    )
    # from temperature equation.
    matrix[eq_idx.temperature, eq_idx.velocity] = (
        f_comp.coolant.dict_Gauss_pt["Gruneisen"][elem_idx]
        * f_comp.coolant.dict_Gauss_pt["temperature"][elem_idx]
    )

    return matrix

def build_transport_coefficients(conductor:Conductor)->Conductor:
    """Function that builds the transport coefficients K', K'' and K''' that appears in the source terms of the fluid components equations.

    Args:
        conductor (Conductor): object with all the information of the conductor.

    Returns:
        Conductor: updated version of the conductor object.
    """

    key_names = {"K1","K2","K3"}
    for key in key_names:
        # Dictionaties declaration.
        conductor.dict_Gauss_pt[key] = dict()

    # Loop in fluid-fluid interfaces.
    for interface in conductor.interface.fluid_fluid:
        # constuct recurrent coefficients of matrix S elements.
        for key in key_names:
            # Initialization.
            conductor.dict_Gauss_pt[key][interface.interf_name] = np.zeros_like(
                conductor.grid_features["zcoord_gauss"]
            )
        # Evaluate pressure difference bethween comp_1 and comp_2
        delta_p = np.abs(
            interface.comp_1.coolant.dict_Gauss_pt["pressure"]
            - interface.comp_2.coolant.dict_Gauss_pt["pressure"]
        )
        # Array smart
        delta_p[delta_p < conductor.Delta_p_min] = conductor.Delta_p_min

        # Find index such that # P_comp_2 < P_comp_1.
        ind_1 = np.nonzero(
            interface.comp_2.coolant.dict_Gauss_pt["pressure"]
            < interface.comp_1.coolant.dict_Gauss_pt["pressure"]
        )[0]
        # Find index such that P_comp_2 >= P_comp_1.
        ind_2 = np.nonzero(
            interface.comp_2.coolant.dict_Gauss_pt["pressure"]
            >= interface.comp_1.coolant.dict_Gauss_pt["pressure"]
        )[0]
        
        # Compute transport coefficients K', K'' and K'''
        conductor = eval_transport_coefficients(
            conductor,
            interface.interf_name,
            interface.comp_1,
            ind_1, # P_comp_2 < P_comp_1
            delta_p
        )
        conductor = eval_transport_coefficients(
            conductor,
            interface.interf_name,
            interface.comp_2,
            ind_2, # P_comp_2 >= P_comp_1
            delta_p
        )

    return conductor

def eval_transport_coefficients(conductor:Conductor,
    interf_name:str,
    comp:FluidComponent,
    index:np.ndarray,
    delta_p:np.ndarray
    )->Conductor:
    """Function that evaluates the transport coefficients K', K'' and K''' that appears in the source terms of the fluid components equations.

    Args:
        conductor (Conductor): object with all the information of the conductor.
        interf_name (str): name of the interface between fluid component objects.
        comp (FluidComponent): fluid component object of the interface with the dominant pressure (index of the gauss points where this is true are passed in inupt argument index.)
        index (np.ndarray): array with the index of the Gauss points where comp pressure is the dominant one (with respect to the pressure of the other fluid component in the interface)
        delta_p (np.ndarray): array with the pressure differece between the component of the interface).

    Returns:
        Conductor: conductor with updated values of K', K'' and K'''.
    """
    
    # K' evaluation [ms]:
    # K' = A_othogonal*sqrt(2*density/k_loc*abs(Delta_p))
    conductor.dict_Gauss_pt["K1"][interf_name][index] = (
        conductor.dict_interf_peri["ch_ch"]["Open"][interf_name]
        * np.sqrt(
            2.0
            * comp.coolant.dict_Gauss_pt["total_density"][index]
            / (conductor.k_loc * delta_p[index])
        )
    )
    # K'' evaluation [m^2]:
    # K'' = K'*lambda_v*velocity
    conductor.dict_Gauss_pt["K2"][interf_name][index] = (
        conductor.dict_Gauss_pt["K1"][interf_name][index]
        * comp.coolant.dict_Gauss_pt["velocity"][index]
        * conductor.lambda_v
    )
    # K''' evaluation [m^3/s]:
    # K''' = K'*(enthalpy + (velocity*lambda_v)^2/2)
    conductor.dict_Gauss_pt["K3"][interf_name][index] = (
        conductor.dict_Gauss_pt["K1"][interf_name][index]
        * (
            comp.coolant.dict_Gauss_pt["total_enthalpy"][index]
            + 
            0.5
            * (
                comp.coolant.dict_Gauss_pt["velocity"][index]
                * conductor.lambda_v
            )
            ** 2
        )
    )

    return conductor

def build_kmat_fluid(
    matrix:np.ndarray,
    upweqt:np.ndarray,
    velocity:float,
    delta_z:float,
    eq_idx:NamedTuple,
    )->np.ndarray:

    """Function that builds the K matrix (KMAT) at the Gauss point, UPWIND is included.
    UPWIND: differencing contribution a' la finite-difference, necessary to guarantee numarical stability that does not came from KMAT algebraic construction.

    Args:
        matrix (np.ndarray): initialized K matrix (np.zeros)
        upweqt (np.ndarray): array with the upwind numerical scheme.
        velocity (float): fluid velocity at present gauss point index.
        delta_z (float): length of the interval that includes the present gauss points.
        eq_idx (NamedTuple): collection of fluid equation index (velocity, pressure and temperaure equations).

    Returns:
        np.ndarray: matrix with updated elements.
    """

    velocity = np.abs(velocity)
    # Build array to assign diagonal coefficients.
    diag_idx = np.array(eq_idx)

    # For the fluid equation this matrix has only the diagonal elements in the 
    # velocity, pressure and temperature equations.
    # Diagonal therms definition: dz * upweqt * v / 2
    # Set diagonal elements (exploit broadcasting).
    matrix[diag_idx,diag_idx] = delta_z * upweqt[diag_idx] * velocity / 2.0

    return matrix

def build_smat_fluid(
    matrix:np.ndarray,
    f_comp:FluidComponent,
    elem_idx:int,
    eq_idx:NamedTuple,
    )->np.ndarray:

    """Function that builds the S matrix (SMAT) therms of the fluid at the Gauss point (SOURCE JACOBIAN).

    Args:
        matrix (np.ndarray): initialized S matrix (np.zeros)
        f_comp (FluidComponent): fluid component object from which get all info to build the coefficients.
        elem_idx (int): index of the i-th element of the spatial discretization.
        eq_idx (NamedTuple): collection of fluid equation index (velocity, pressure and temperaure equations).

    Returns:
        np.ndarray: matrix with updated elements.
    """

    # velocity equation: main diagonal elements construction
    # (j,j) [vel_j]
    matrix[eq_idx.velocity,eq_idx.velocity] = (
        2.0
        # dict_friction_factor[False]["total"]: total friction factor in Gauss 
        # points (see __init__ of class Channel for details).
        * f_comp.channel.dict_friction_factor[False]["total"][elem_idx]
        * np.abs(f_comp.coolant.dict_Gauss_pt["velocity"][elem_idx])
        / f_comp.channel.inputs["HYDIAMETER"]
    )
    
    # pressure equation: elements below main diagonal construction
    # (j+num_fluid_components,0:num_fluid_components) [Pres]
    matrix[eq_idx.pressure,eq_idx.velocity] = (
        - matrix[eq_idx.velocity,eq_idx.velocity]
        * f_comp.coolant.dict_Gauss_pt["Gruneisen"][elem_idx]
        * f_comp.coolant.dict_Gauss_pt["total_density"][elem_idx]
        * f_comp.coolant.dict_Gauss_pt["velocity"][elem_idx]
    )
    
    # temperature equation: elements below main diagonal construction
    # (j+2*num_fluid_components,0:num_fluid_components) [Temp]
    matrix[eq_idx.temperature,eq_idx.velocity] = (
        - matrix[eq_idx.velocity,eq_idx.velocity]
        / f_comp.coolant.dict_Gauss_pt["total_isochoric_specific_heat"][elem_idx]
        * f_comp.coolant.dict_Gauss_pt["velocity"][elem_idx]
    )

    return matrix

def build_smat_fluid_interface(
    matrix:np.ndarray,
    conductor:Conductor,
    elem_idx:int,
    eq_idx:dict,
    )->np.ndarray:

    """Function that builds the S matrix (SMAT) therms due to fluid component interfaces at the Gauss point (SOURCE JACOBIAN).

    Args:
        matrix (np.ndarray): S matrix after call to function buld_smat_fluid.
        conductor (Conductor): object with all the information of the conductor.
        elem_idx (int): index of the i-th element of the spatial discretization.
        eq_idx (dict): collection of NamedTyple with fluid equation index (velocity, pressure and temperaure equations).

    Returns:
        np.ndarray: matrix with updated elements.
    """

    # NOMENCLATURE
    # w: enthalpy
    # phi: Gruneisen
    # rho: density
    # c0: speed of sound
    # h: heat transfer coefficient (_o: open; _c:close)
    # P: contact perimeter (_o: open; _c:close)
    # A: cross section
    # v: velocity
    # T: temperature
    # c_v: isochoric specific heat

    for interface in conductor.interface.fluid_fluid:
        
        # VELOCITY EQUATION: above/below main diagonal elements construction:
        # (j,j+num_fluid_components) [Pres_j]

        # s_vj_pj = K1 * v - K2 / (A * rho)

        s_vj_pj = (
            conductor.dict_Gauss_pt["K1"][interface.interf_name][elem_idx]
            * interface.comp_1.coolant.dict_Gauss_pt["velocity"][elem_idx]
            - conductor.dict_Gauss_pt["K2"][interface.interf_name][elem_idx]
        ) / (
            interface.comp_1.channel.inputs["CROSSECTION"]
            * interface.comp_1.coolant.dict_Gauss_pt["total_density"][elem_idx]
        )

        matrix[
            eq_idx[interface.comp_1.identifier].velocity,
            eq_idx[interface.comp_1.identifier].pressure,
        ] = matrix[
            eq_idx[interface.comp_1.identifier].velocity,
            eq_idx[interface.comp_1.identifier].pressure,
        ] - s_vj_pj

        # (j,k + num_fluid_components:2*num_fluid_components) 
        # [Pres_k]
        matrix[
            eq_idx[interface.comp_1.identifier].velocity,
            eq_idx[interface.comp_2.identifier].pressure,
        ] = s_vj_pj

        # PRESSURE EQUATION: main diagonal elements construction:
        # (j+num_fluid_components,j+num_fluid_components) [Pres_j]

        # coef_grun_area = phi / A
        coef_grun_area = (
            interface.comp_1.coolant.dict_Gauss_pt["Gruneisen"][elem_idx]
            / interface.comp_1.channel.inputs["CROSSECTION"]
        )

        # s_pj_pj = phi/A * [K3 - vK2 - (w - v^2/2 - c0^2/phi)K1]
        #         = coef_grun_area * [K3 - vK2 - (w - v^2/2 - c0^2/phi)K1]
        s_pj_pj = (
            coef_grun_area
            * (
                conductor.dict_Gauss_pt["K3"][interface.interf_name][elem_idx]
                - interface.comp_1.coolant.dict_Gauss_pt["velocity"][elem_idx]
                * conductor.dict_Gauss_pt["K2"][interface.interf_name][elem_idx]
                - (
                    interface.comp_1.coolant.dict_Gauss_pt["total_enthalpy"][elem_idx]
                    - interface.comp_1.coolant.dict_Gauss_pt["velocity"][elem_idx]
                    ** 2
                    / 2.0
                    - interface.comp_1.coolant.dict_Gauss_pt["total_speed_of_sound"][elem_idx]
                    ** 2
                    / interface.comp_1.coolant.dict_Gauss_pt["Gruneisen"][elem_idx]
                )
                * conductor.dict_Gauss_pt["K1"][interface.interf_name][elem_idx]
            )
        )

        matrix[
            eq_idx[interface.comp_1.identifier].pressure,
            eq_idx[interface.comp_1.identifier].pressure,
        ] = matrix[
            eq_idx[interface.comp_1.identifier].pressure,
            eq_idx[interface.comp_1.identifier].pressure,
        ] + s_pj_pj

        # PRESSURE EQUATION: above/below main diagonal elements construction:
        # (j+num_fluid_components,\
        # k + num_fluid_components:2*num_fluid_components) [Pres_k]
        matrix[
            eq_idx[interface.comp_1.identifier].pressure,
            eq_idx[interface.comp_2.identifier].pressure,
        ] = - s_pj_pj

        # (j+num_fluid_components,j+2*num_fluid_components)
        # [Temp_j] I
        # coef_htc = P_o * h_o + P_c * h_c
        coef_htc = (
            conductor.dict_interf_peri["ch_ch"]["Open"][interface.interf_name]
            * conductor.dict_Gauss_pt["HTC"]["ch_ch"]["Open"][
                interface.interf_name
            ][elem_idx]
            + conductor.dict_interf_peri["ch_ch"]["Close"][
                interface.interf_name
            ]
            * conductor.dict_Gauss_pt["HTC"]["ch_ch"]["Close"][
                interface.interf_name
            ][elem_idx]
        )
        # s_pj_tc = phi/A * (P_o * h_o + P_c * h_c)
        #         = coef_frun_area * coef_htc
        s_pj_tj = coef_grun_area * coef_htc

        matrix[
            eq_idx[interface.comp_1.identifier].pressure,
            eq_idx[interface.comp_1.identifier].temperature,
        ] = matrix[
            eq_idx[interface.comp_1.identifier].pressure,
            eq_idx[interface.comp_1.identifier].temperature,
        ] + s_pj_tj

        # (j+num_fluid_components, 
        # k + 2*num_fluid_components:dict_N_equation
        # ["FluidComponent"]) [Temp_j]
        matrix[
            eq_idx[interface.comp_1.identifier].pressure,
            eq_idx[interface.comp_2.identifier].temperature,
        ] = - s_pj_tj

        # TEMPERATURE EQUATION: elements below main diagonal \
        # construction:
        # (j+2*num_fluid_components,j+num_fluid_components) [Pres_j]

        # coef_rho_cv_area = 1/(rho * c_v * A)
        coef_rho_cv_area = 1.0 / (
            interface.comp_1.coolant.dict_Gauss_pt["total_density"][elem_idx]
            * interface.comp_1.coolant.dict_Gauss_pt[
                "total_isochoric_specific_heat"
            ][elem_idx]
            * interface.comp_1.channel.inputs["CROSSECTION"]
        )
        # s_tj_pj = 1/(rho * c_v * A) * [K3 - vK2 - (w - v^2/2 - phi*c_v*T)K1]
        #         = coef_rho_cv_area * [K3 - vK2 - (w - v^2/2 - phi*c_v*T)K1]
        s_tj_pj = (
            coef_rho_cv_area
            * (
                conductor.dict_Gauss_pt["K3"][interface.interf_name][elem_idx]
                - interface.comp_1.coolant.dict_Gauss_pt["velocity"][elem_idx]
                * conductor.dict_Gauss_pt["K2"][interface.interf_name][elem_idx]
                - (
                    interface.comp_1.coolant.dict_Gauss_pt["total_enthalpy"][elem_idx]
                    - interface.comp_1.coolant.dict_Gauss_pt["velocity"][elem_idx]
                    ** 2
                    / 2.0
                    - interface.comp_1.coolant.dict_Gauss_pt["Gruneisen"][elem_idx]
                    * interface.comp_1.coolant.dict_Gauss_pt["total_isochoric_specific_heat"][elem_idx]
                    * interface.comp_1.coolant.dict_Gauss_pt["temperature"][elem_idx]
                )
                * conductor.dict_Gauss_pt["K1"][interface.interf_name][elem_idx]
            )
        )

        matrix[
            eq_idx[interface.comp_1.identifier].temperature,
            eq_idx[interface.comp_1.identifier].pressure,
        ] = matrix[
            eq_idx[interface.comp_1.identifier].temperature,
            eq_idx[interface.comp_1.identifier].pressure,
        ] + s_tj_pj

        # (j+2*num_fluid_components,\
        # k + num_fluid_components:2*num_fluid_components) [Pres_k]
        matrix[
            eq_idx[interface.comp_1.identifier].temperature,
            eq_idx[interface.comp_2.identifier].pressure,
        ] = - s_tj_pj

        # TEMPERATURE EQUATION: main diagonal element construction:
        # (j+2*num_fluid_components,j+2*num_fluid_components) 
        # [Temp_j] I

        # s_tj_tj = 1/(rho * c_v * A) * (P_o * h_o + P_c * h_c)
        #         = coef_rho_cv_area * coef_htc
        s_tj_tj = coef_rho_cv_area * coef_htc

        matrix[
            eq_idx[interface.comp_1.identifier].temperature,
            eq_idx[interface.comp_1.identifier].temperature,
        ] = matrix[
            eq_idx[interface.comp_1.identifier].temperature,
            eq_idx[interface.comp_1.identifier].temperature,
        ] + s_tj_tj

        # TEMPERATURE EQUATION: above/below main diagonal elements 
        # construction:
        # (j+2*num_fluid_components,k + 2*num_fluid_components) 
        # [Temp_k]
        matrix[
            eq_idx[interface.comp_1.identifier].temperature,
            eq_idx[interface.comp_2.identifier].temperature,
        ] = - s_tj_tj

def build_smat_fluid_solid_interface(
    matrix:np.ndarray,
    conductor:Conductor,
    elem_idx:int,
    eq_idx:dict,
    )->np.ndarray:

    """Function that builds the S matrix (SMAT) therms due to fluid-solid component interfaces at the Gauss point (SOURCE JACOBIAN).

    Args:
        matrix (np.ndarray): S matrix after call to function buld_smat_fluid_interface.
        conductor (Conductor): object with all the information of the conductor.
        elem_idx (int): index of the i-th element of the spatial discretization.
        eq_idx (dict): collection of NamedTuple with fluid equation index (velocity, pressure and temperaure equations) and of scalar for solid equation index.

    Returns:
        np.ndarray: matrix with updated elements.
    """
    
    # NOMENCLATURE
    # phi: Gruneisen
    # rho: density
    # h: heat transfer coefficient (_o: open; _c:close)
    # P: contact perimeter (_o: open; _c:close)
    # A: cross section
    # c_v: isochoric specific heat

    for interface in enumerate(conductor.interface.fluid_solid):
        # pressure equation: above main diagonal elements construction.
        # (j+num_fluid_components,j+2*num_fluid_components) [Temp_j] II + III

        # coef_grun_area = phi / A
        coef_grun_area = (
            interface.comp_1.coolant.dict_Gauss_pt["Gruneisen"][elem_idx]
            / interface.comp_1.channel.inputs["CROSSECTION"]
        )

        # coef_htc = P * h
        coef_htc = (
            conductor.dict_interf_peri["ch_sol"][interface.interf_name]
            * conductor.dict_Gauss_pt["HTC"]["ch_sol"][interface.interf_name][
                elem_idx
            ]
        )
        
        # s_pj_tj = phi / A * P * h
        #         = coef_grun_area * coef_htc
        s_pj_tj = coef_grun_area * coef_htc

        matrix[
            eq_idx[interface.comp_1.identifier].pressure,
            eq_idx[interface.comp_1.identifier].temperature
        ] = matrix[
            eq_idx[interface.comp_1.identifier].pressure,
            eq_idx[interface.comp_1.identifier].temperature,
        ] + s_pj_tj
        
        # (j+num_fluid_components,l + dict_N_equation["FluidComponent"]) [Temp_l]
        matrix[
            eq_idx[interface.comp_1.identifier].pressure,
            eq_idx[interface.comp_2.identifier],
        ] = - s_pj_tj

        # temperature equation: main diagonal element construction
        # (j+2*num_fluid_components,j+2*num_fluid_components) [Temp_j] II + III 
        
        # coef_rho_cv_area = 1/(rho * c_v * A)
        coef_rho_cv_area = 1.0 / (
            interface.comp_1.coolant.dict_Gauss_pt["total_density"][elem_idx]
            * interface.comp_1.coolant.dict_Gauss_pt[
                "total_isochoric_specific_heat"
            ][elem_idx]
            * interface.comp_1.channel.inputs["CROSSECTION"]
        )

        # s_tj_tj = 1/(rho * c_v * A) * P * h
        #         = coef_rho_cv_area * coef_htc
        s_tj_tj = coef_rho_cv_area * coef_htc
        
        matrix[
            eq_idx[interface.comp_1.identifier].temperature,
            eq_idx[interface.comp_1.identifier].temperature,
        ] = matrix[
            eq_idx[interface.comp_1.identifier].temperature,
            eq_idx[interface.comp_1.identifier].temperature,
        ] + s_tj_tj
        
        # temperature equation: above main diagonal elements construction
        # (j+2*num_fluid_components,l + dict_N_equation["FluidComponent"]) [Temp_l]
        matrix[
            eq_idx[interface.comp_1.identifier].temperature,
            eq_idx[interface.comp_2.identifier],
        ] = - s_tj_tj

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

    # Matrices initialization.
    MASMAT,FLXMAT,SORMAT,DIFMAT,SYSMAT = matrix_initialization(
        conductor.dict_band["Full"],
        conductor.dict_N_equation["Total"]
    )
    
    if conductor.inputs["METHOD"] == "BE" or conductor.inputs["METHOD"] == "CN":
        # Backward Euler or Crank-Nicolson (cdp, 10/2020)
        if conductor.cond_num_step > 1:
            # Copy the load vector at the previous time step in the second column to \
            # correctly apply the theta method (cdp, 10/2020)
            conductor.dict_Step["SYSLOD"][:, 1] = conductor.dict_Step["SYSLOD"][
                :, 0
            ].copy()
            conductor.dict_Step["SYSLOD"][:, 0] = 0.0

    # SYSVAR = np.zeros(conductor.dict_N_equation["Total"])
    ASCALING = np.zeros(conductor.dict_N_equation["Total"])
    UPWEQT = np.zeros(conductor.dict_N_equation["NODOFS"])
    # Known terms vector (cdp, 10/2020)
    Known = np.zeros_like(ASCALING)

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

    # ** MATRICES CONSTRUCTION (cdp, 07/2020) **

    # Build equations index for all FluidComponent objects
    fluid_eq_idx = __build_fluid_eq_idx(conductor)

    # riscrivere in forma array smart una volta risolti tutti i dubbi, se possibile (cdp, 07/2020)
    for ii in range(conductor.grid_input["NELEMS"]):

        # Auxiliary matrices initialization to zeros at each Gauss point.
        MMAT,AMAT,KMAT,SMAT,SVEC = ndarray_initialization(
            conductor.dict_N_equation["NODOFS"],
            conductor.cond_num_step,
            col=2
        )
        
        ELMMAT,ELAMAT,ELKMAT,ELSMAT,ELSLOD = ndarray_initialization(
            2 * conductor.dict_N_equation["NODOFS"],
            conductor.cond_num_step
        )
        
        jump = conductor.dict_N_equation["NODOFS"] * ii

        # (cdp, 07/2020)
        # ** FORM THE M, A, K, S MATRICES AND S VECTOR AT THE GAUSS POINT, FLUID \
        # COMPONENTS EQUATIONS **
        # FORM THE M MATRIX AT THE GAUSS POINT (MASS AND CAPACITY)
        # FluidComponent equation: array smart (cdp, 07/2020)
        MMAT[
            :conductor.dict_N_equation["FluidComponent"],
            :conductor.dict_N_equation["FluidComponent"],
        ] = np.eye(conductor.dict_N_equation["FluidComponent"])
        # END M MATRIX: fluid components equations (cdp, 07/2020)
        for jj, fluid_comp_j in enumerate(conductor.inventory["FluidComponent"].collection):

            # FORM THE A MATRIX AT THE GAUSS POINT (FLUX JACOBIAN)
            AMAT = __build_amat(
                AMAT,
                fluid_comp_j,
                ii,
                fluid_eq_idx[fluid_comp_j.identifier]
            )

            # FORM THE K MATRIX AT THE GAUSS POINT (INCLUDING UPWIND)
            KMAT = build_kmat_fluid(
                KMAT,
                UPWEQT,
                fluid_comp_j.coolant.dict_Gauss_pt["velocity"][ii],
                conductor.grid_features["delta_z"][ii],
                fluid_eq_idx[fluid_comp_j.identifier]
            )

            # FORM THE S MATRIX AT THE GAUSS POINT (SOURCE JACOBIAN)
            SMAT = build_smat_fluid(
                SMAT,
                fluid_comp_j,
                ii,
                fluid_eq_idx[fluid_comp_j.identifier]
            )

        # FORM THE S MATRIX AT THE GAUSS POINT (SOURCE JACOBIAN)
        # Therms associated to fluid-fluid interfaces.
        SMAT = build_smat_fluid_interface(
            SMAT,
            conductor,
            ii,
            fluid_eq_idx
        )

            for ll, s_comp in enumerate(conductor.inventory["SolidComponent"].collection):
                # chan_sol_topology is equivalent to \
                # conductor.dict_topology["ch_sol"][fluid_comp_r.identifier][s_comp.identifier] \
                # but it is shorter so I decide to use it here (cdp, 09/2020)
                chan_sol_topology = f"{fluid_comp_j.identifier}_{s_comp.identifier}"
                # Check for valid interface.
                if chan_sol_topology in conductor.dict_interf_peri["ch_sol"]:
                    # Perform calculation only if there is an interface, this \
                    # will reduce the computational time (cdp, 09/2020)
                    # pressure equation: above main diagonal elements
                    # construction (cdp, 07/2020)
                    # (j+num_fluid_components,j+2*num_fluid_components) [Temp_j] \
                    # II + III (cdp, 07/2020)
                    SMAT[
                        jj + conductor.inventory["FluidComponent"].number,
                        jj
                        + 2 * conductor.inventory["FluidComponent"].number,
                    ] = SMAT[
                        jj + conductor.inventory["FluidComponent"].number,
                        jj
                        + 2 * conductor.inventory["FluidComponent"].number,
                    ] + (
                        fluid_comp_j.coolant.dict_Gauss_pt["Gruneisen"][ii]
                        / fluid_comp_j.channel.inputs["CROSSECTION"]
                    ) * (
                        conductor.dict_interf_peri["ch_sol"][chan_sol_topology]
                        * conductor.dict_Gauss_pt["HTC"]["ch_sol"][chan_sol_topology][
                            ii
                        ]
                    )
                    # (j+num_fluid_components,l + dict_N_equation["FluidComponent"]) [Temp_l] (cdp, 07/2020)
                    SMAT[
                        jj + conductor.inventory["FluidComponent"].number,
                        ll + conductor.dict_N_equation["FluidComponent"],
                    ] = -(
                        fluid_comp_j.coolant.dict_Gauss_pt["Gruneisen"][ii]
                        / fluid_comp_j.channel.inputs["CROSSECTION"]
                    ) * (
                        conductor.dict_interf_peri["ch_sol"][chan_sol_topology]
                        * conductor.dict_Gauss_pt["HTC"]["ch_sol"][chan_sol_topology][
                            ii
                        ]
                    )
                    # temperature equation: main diagonal element construction \
                    # (cdp, 07/2020)
                    # (j+2*num_fluid_components,j+2*num_fluid_components) [Temp_j] \
                    # II + III (cdp, 07/2020)
                    SMAT[
                        jj
                        + 2 * conductor.inventory["FluidComponent"].number,
                        jj
                        + 2 * conductor.inventory["FluidComponent"].number,
                    ] = SMAT[
                        jj
                        + 2 * conductor.inventory["FluidComponent"].number,
                        jj
                        + 2 * conductor.inventory["FluidComponent"].number,
                    ] + 1.0 / (
                        fluid_comp_j.coolant.dict_Gauss_pt["total_density"][ii]
                        * fluid_comp_j.coolant.dict_Gauss_pt[
                            "total_isochoric_specific_heat"
                        ][ii]
                        * fluid_comp_j.channel.inputs["CROSSECTION"]
                    ) * (
                        conductor.dict_interf_peri["ch_sol"][chan_sol_topology]
                        * conductor.dict_Gauss_pt["HTC"]["ch_sol"][chan_sol_topology][
                            ii
                        ]
                    )
                    # temperature equation: above main diagonal elements
                    # construction (cdp, 07/2020)
                    # (j+2*num_fluid_components,l + dict_N_equation["FluidComponent"]) [Temp_l] (cdp, 07/2020)
                    SMAT[
                        jj
                        + 2 * conductor.inventory["FluidComponent"].number,
                        ll + conductor.dict_N_equation["FluidComponent"],
                    ] = (
                        -1.0
                        / (
                            fluid_comp_j.coolant.dict_Gauss_pt["total_density"][ii]
                            * fluid_comp_j.coolant.dict_Gauss_pt[
                                "total_isochoric_specific_heat"
                            ][ii]
                            * fluid_comp_j.channel.inputs["CROSSECTION"]
                        )
                        * (
                            conductor.dict_interf_peri["ch_sol"][chan_sol_topology]
                            * conductor.dict_Gauss_pt["HTC"]["ch_sol"][
                                chan_sol_topology
                            ][ii]
                        )
                    )
                # end if flag_chan_sol_contact (cdp, 09/2020)
            # end for ll (cdp, 07/2020)
            # END S MATRIX: fluid components equations (cdp, 07/2020)

            # FORM THE S VECTOR AT THE NODAL POINTS (SOURCE)
            # Set to zeros in initialization (cdp, 08/2020)
            # END S VECTOR: fluid components equations (cdp, 07/2020)
        # end for jj (cdp, 07/2020)

        # (cdp, 07/2020)
        # * FORM THE M, A, K, S MATRICES AND S VECTOR AT THE GAUSS POINT, SOLID \
        # COMPONENTS EQUATIONS *
        for ll, s_comp_l in enumerate(conductor.inventory["SolidComponent"].collection):
            neq = conductor.dict_N_equation["FluidComponent"] + ll
            # FORM THE M MATRIX AT THE GAUSS POINT (MASS AND CAPACITY)
            # SolidComponent equation (cdp, 07/2020)
            # A_{s_comp}*rho_{s_comp,homo}*cp_{s_comp,homo}/cos(theta); \
            # homo = homogenized (cdp, 07/2020)
            MMAT[neq, neq] = (
                s_comp_l.inputs["CROSSECTION"]
                * s_comp_l.dict_Gauss_pt["total_density"][ii]
                * s_comp_l.dict_Gauss_pt["total_isobaric_specific_heat"][ii]
                / s_comp_l.inputs["COSTETA"]
            )
            # END M MATRIX: solid components equation (cdp, 07/2020)

            # FORM THE A MATRIX AT THE GAUSS POINT (FLUX JACOBIAN)
            # No elements here (cdp, 07/2020)
            # END A MATRIX: solid components equation (cdp, 07/2020)

            # FORM THE K MATRIX AT THE GAUSS POINT (INCLUDING UPWIND)
            # A_{s_comp}*k_{s_comp,homo}; homo = homogenized (cdp, 07/2020)
            KMAT[neq, neq] = (
                s_comp_l.inputs["CROSSECTION"]
                * s_comp_l.dict_Gauss_pt["total_thermal_conductivity"][ii]
                / s_comp_l.inputs["COSTETA"]
            )
            # END K MATRIX: solid components equation (cdp, 07/2020)

            # FORM THE S MATRIX AT THE GAUSS POINT (SOURCE JACOBIAN)
            s_comp_filter = filter_component(
                conductor.inventory["SolidComponent"].collection,
                s_comp_l
            )
            for s_comp_m in s_comp_filter:
                # s_comp_topology is equivalent to 
                # conductor.dict_topology["sol_sol"][s_comp_m.identifier]
                # [s_comp_l.identifier] but it is shorter so I decide to 
                # use it here.
                s_comp_topology = natural_sort(s_comp_l, s_comp_m.obj)
                # Check for valid interface.
                if s_comp_topology in conductor.dict_interf_peri["sol_sol"]:
                    # Perform calculation only if there is an interface, 
                    # this will reduce the computational time.

                    # SOLID COMPONENTS CONDUCTION EQUATION: main diagonal 
                    # element construction:
                    # (l + dict_N_equation["FluidComponent"],l 
                    # + dict_N_equation["FluidComponent"]) [Temp_l] II 
                    # + III
                    SMAT[neq, neq] = (
                        SMAT[neq, neq]
                        + conductor.dict_interf_peri["sol_sol"][s_comp_topology]
                        * conductor.dict_Gauss_pt["HTC"]["sol_sol"]["cond"][
                            s_comp_topology
                        ][ii]
                    )
                    
                    # SOLID COMPONENTS CONDUCTION EQUATION: above/below 
                    # main diagonal elements construction:
                    # (l + dict_N_equation["FluidComponent"],m 
                    # + dict_N_equation["FluidComponent"]) [Temp_m]
                    SMAT[neq, s_comp_m.idx + conductor.dict_N_equation["FluidComponent"]] = (
                        -conductor.dict_interf_peri["sol_sol"][s_comp_topology]
                        * conductor.dict_Gauss_pt["HTC"]["sol_sol"]["cond"][
                            s_comp_topology
                        ][ii]
                    )
                # end if flag_sol_sol_contact (cdp, 09/2020)
            # end for s_comp_m (cdp, 07/2020)
            for jj, fluid_comp_j in enumerate(conductor.inventory["FluidComponent"].collection):
                chan_sol_topology = f"{fluid_comp_j.identifier}_{s_comp_l.identifier}"
                # Check for valid interface.
                if chan_sol_topology in conductor.dict_interf_peri["ch_sol"]:
                    # Perform calculation only if there is an interface, this \
                    # will reduce the computational time (cdp, 09/2020)
                    # SOLID COMPONENTS CONDUCTION EQUATION: main diagonal element \
                    # construction (cdp, 07/2020)
                    # (l + dict_N_equation["FluidComponent"],l + dict_N_equation["FluidComponent"]) [Temp_l] I (cdp, 07/2020)
                    SMAT[neq, neq] = (
                        SMAT[neq, neq]
                        + conductor.dict_interf_peri["ch_sol"][chan_sol_topology]
                        * conductor.dict_Gauss_pt["HTC"]["ch_sol"][chan_sol_topology][
                            ii
                        ]
                    )
                    # SOLID COMPONENTS CONDUCTION EQUATION: below main diagonal elements
                    # construction (cdp, 07/2020)
                    # (l + dict_N_equation["FluidComponent"],l + 2*num_fluid_components) [Temp_j] (cdp, 07/2020)
                    SMAT[
                        neq,
                        jj
                        + 2 * conductor.inventory["FluidComponent"].number,
                    ] = (
                        -conductor.dict_interf_peri["ch_sol"][chan_sol_topology]
                        * conductor.dict_Gauss_pt["HTC"]["ch_sol"][chan_sol_topology][
                            ii
                        ]
                    )
                # end if flag_chan_sol_contact (cdp, 09/2020)
            # end for jj (cdp, 07/2020)
            # Convective heating with the external environment (implicit treatment).
            if s_comp_l.name == conductor.inventory["JacketComponent"].name:
                if (
                    conductor.dict_df_coupling["contact_perimeter_flag"].at[
                        environment.KIND, s_comp_l.identifier
                    ]
                    == 1
                ):
                    if (
                        conductor.dict_df_coupling["HTC_choice"].at[
                            environment.KIND, s_comp_l.identifier
                        ]
                        == 2
                        and conductor.inputs["Is_rectangular"]
                    ):
                        # Rectangular duct.
                        SMAT[neq, neq] = (
                            SMAT[neq, neq]
                            + 2
                            * conductor.inputs["Height"]
                            * conductor.dict_Gauss_pt["HTC"]["env_sol"][
                                f"{environment.KIND}_{s_comp_l.identifier}"
                            ]["conv"]["side"][ii]
                            + conductor.inputs["width"]
                            * (
                                conductor.dict_Gauss_pt["HTC"]["env_sol"][
                                    f"{environment.KIND}_{s_comp_l.identifier}"
                                ]["conv"]["bottom"][ii]
                                + conductor.dict_Gauss_pt["HTC"]["env_sol"][
                                    f"{environment.KIND}_{s_comp_l.identifier}"
                                ]["conv"]["top"][ii]
                            )
                        )
                    else:
                        SMAT[neq, neq] = (
                            SMAT[neq, neq]
                            + conductor.dict_interf_peri["env_sol"][
                                f"{environment.KIND}_{s_comp_l.identifier}"
                            ]
                            * conductor.dict_Gauss_pt["HTC"]["env_sol"][
                                f"{environment.KIND}_{s_comp_l.identifier}"
                            ]["conv"][ii]
                        )
            # End s_comp_l.name

            # END S MATRIX: solid components equation (cdp, 07/2020)

            # FORM THE S VECTOR AT THE NODAL POINTS (SOURCE)
            # cl modify august 24 2019
            if s_comp_l.name != conductor.inventory["JacketComponent"].name:
                # StrandComponent objects (cdp, 08/2020)
                # This is independent from the solution method thanks to the \
                # escamotage of the dummy steady state corresponding to the \
                # initialization (cdp, 10/2020)
                if conductor.cond_num_step == 1:
                    # Current time step (cdp, 10/2020)
                    SVEC.present[neq, 0] = s_comp_l.dict_Gauss_pt["Q1"][ii, 0]
                    SVEC.present[neq, 1] = s_comp_l.dict_Gauss_pt["Q2"][ii, 0]
                    # Previous time step (cdp, 10/2020)
                    SVEC.previous[neq, 0] = s_comp_l.dict_Gauss_pt["Q1"][ii, 1]
                    SVEC.previous[neq, 1] = s_comp_l.dict_Gauss_pt["Q2"][ii, 1]
                else:
                    # Compute only at the current time step (cdp, 10/2020)
                    SVEC[neq, 0] = s_comp_l.dict_Gauss_pt["Q1"][ii, 0]
                    SVEC[neq, 1] = s_comp_l.dict_Gauss_pt["Q2"][ii, 0]
            else:
                # JacketComponents objects (cdp, 08/2020)
                # This is independent from the solution method thanks to the \
                # escamotage of the dummy steady state corresponding to the \
                # initialization (cdp, 10/2020)
                if conductor.cond_num_step == 1:
                    # Current time step (cdp, 10/2020)
                    SVEC.present[neq, 0] = (
                        s_comp_l.dict_Gauss_pt["Q1"][ii, 0]
                        - qsource[ii, ll - conductor.dict_N_equation["StrandComponent"] - 1]
                    )
                    SVEC.present[neq, 1] = (
                        s_comp_l.dict_Gauss_pt["Q2"][ii, 0]
                        - qsource[ii + 1, ll - conductor.dict_N_equation["StrandComponent"] - 1]
                    )
                    # Previous time step (cdp, 10/2020)
                    SVEC.previous[neq, 0] = (
                        s_comp_l.dict_Gauss_pt["Q1"][ii, 1]
                        - qsource[ii, ll - conductor.dict_N_equation["StrandComponent"] - 1]
                    )
                    SVEC.previous[neq, 1] = (
                        s_comp_l.dict_Gauss_pt["Q2"][ii, 1]
                        - qsource[ii + 1, ll - conductor.dict_N_equation["StrandComponent"] - 1]
                    )
                    if (
                        conductor.dict_df_coupling["contact_perimeter_flag"].at[
                            environment.KIND, s_comp_l.identifier
                        ]
                        == 1
                    ):
                        # Add the contribution of the external heating by convection to the known term vector.
                        if (
                            conductor.dict_df_coupling["HTC_choice"].at[
                                environment.KIND, s_comp_l.identifier
                            ]
                            == 2
                            and conductor.inputs["Is_rectangular"]
                        ):
                            # Rectangular duct.
                            coef = 2 * conductor.inputs[
                                "Height"
                            ] * conductor.dict_Gauss_pt["HTC"]["env_sol"][
                                f"{environment.KIND}_{s_comp_l.identifier}"
                            ][
                                "conv"
                            ][
                                "side"
                            ][
                                ii
                            ] + conductor.inputs[
                                "width"
                            ] * (
                                conductor.dict_Gauss_pt["HTC"]["env_sol"][
                                    f"{environment.KIND}_{s_comp_l.identifier}"
                                ]["conv"]["bottom"][ii]
                                + conductor.dict_Gauss_pt["HTC"]["env_sol"][
                                    f"{environment.KIND}_{s_comp_l.identifier}"
                                ]["conv"]["top"][ii]
                            )
                        else:
                            coef = (
                                conductor.dict_interf_peri["env_sol"][
                                    f"{environment.KIND}_{s_comp_l.identifier}"
                                ]
                                * conductor.dict_Gauss_pt["HTC"]["env_sol"][
                                    f"{environment.KIND}_{s_comp_l.identifier}"
                                ]["conv"][ii]
                            )
                        # End if conductor.inputs["Is_rectangular"]

                        SVEC.present[neq, 0] = (
                            SVEC.present[neq, 0]
                            + coef * environment.inputs["Temperature"]
                        )  # W/m
                        SVEC.present[neq, 1] = (
                            SVEC.present[neq, 1]
                            + coef * environment.inputs["Temperature"]
                        )  # W/m
                        SVEC.previous[neq, 0] = (
                            SVEC.previous[neq, 0]
                            + coef * environment.inputs["Temperature"]
                        )  # W/m
                        SVEC.previous[neq, 1] = (
                            SVEC.previous[neq, 1]
                            + coef * environment.inputs["Temperature"]
                        )  # W/m
                    # End if conductor.dict_df_coupling["contact_perimeter_flag"].at[environment.KIND, s_comp_l.identifier] == 1
                else:
                    # Compute only at the current time step (cdp, 10/2020)
                    SVEC[neq, 0] = (
                        s_comp_l.dict_Gauss_pt["Q1"][ii, 0]
                        - qsource[ii, ll - conductor.dict_N_equation["StrandComponent"] - 1]
                    )
                    SVEC[neq, 1] = (
                        s_comp_l.dict_Gauss_pt["Q2"][ii, 0]
                        - qsource[ii + 1, ll - conductor.dict_N_equation["StrandComponent"] - 1]
                    )
                    if (
                        conductor.dict_df_coupling["contact_perimeter_flag"].at[
                            environment.KIND, s_comp_l.identifier
                        ]
                        == 1
                    ):
                        # Add the contribution of the external heating by convection to the known term vector.
                        if (
                            conductor.dict_df_coupling["HTC_choice"].at[
                                environment.KIND, s_comp_l.identifier
                            ]
                            == 2
                            and conductor.inputs["Is_rectangular"]
                        ):
                            # Rectangular duct.
                            coef = 2 * conductor.inputs[
                                "Height"
                            ] * conductor.dict_Gauss_pt["HTC"]["env_sol"][
                                f"{environment.KIND}_{s_comp_l.identifier}"
                            ][
                                "conv"
                            ][
                                "side"
                            ][
                                ii
                            ] + conductor.inputs[
                                "width"
                            ] * (
                                conductor.dict_Gauss_pt["HTC"]["env_sol"][
                                    f"{environment.KIND}_{s_comp_l.identifier}"
                                ]["conv"]["bottom"][ii]
                                + conductor.dict_Gauss_pt["HTC"]["env_sol"][
                                    f"{environment.KIND}_{s_comp_l.identifier}"
                                ]["conv"]["top"][ii]
                            )
                        else:
                            coef = (
                                conductor.dict_interf_peri["env_sol"][
                                    f"{environment.KIND}_{s_comp_l.identifier}"
                                ]
                                * conductor.dict_Gauss_pt["HTC"]["env_sol"][
                                    f"{environment.KIND}_{s_comp_l.identifier}"
                                ]["conv"][ii]
                            )
                        # End if conductor.inputs["Is_rectangular"]

                        SVEC[neq, 0] = (
                            SVEC[neq, 0] + coef * environment.inputs["Temperature"]
                        )  # W/m
                        SVEC[neq, 1] = (
                            SVEC[neq, 1] + coef * environment.inputs["Temperature"]
                        )  # W/m
                    # End if conductor.dict_df_coupling["contact_perimeter_flag"].at[environment.KIND, s_comp_l.identifier] == 1
            # cl end august 24 2019
            # END S VECTOR: solid components equation (cdp, 07/2020)
        # end for ll (cdp, 07/2020)

        # COMPUTE THE MASS AND CAPACITY MATRIX
        # array smart (cdp, 07/2020)
        ELMMAT[
            :conductor.dict_N_equation["NODOFS"],
            :conductor.dict_N_equation["NODOFS"],
        ] = (conductor.grid_features["delta_z"][ii] * (1.0 / 3.0 + ALFA) * MMAT)
        ELMMAT[
            :conductor.dict_N_equation["NODOFS"],
            conductor.dict_N_equation["NODOFS"] : 2
            * conductor.dict_N_equation["NODOFS"],
        ] = (conductor.grid_features["delta_z"][ii] * (1.0 / 6.0 - ALFA) * MMAT)
        ELMMAT[
            conductor.dict_N_equation["NODOFS"] : 2
            * conductor.dict_N_equation["NODOFS"],
            :conductor.dict_N_equation["NODOFS"],
        ] = (conductor.grid_features["delta_z"][ii] * (1.0 / 6.0 - ALFA) * MMAT)
        ELMMAT[
            conductor.dict_N_equation["NODOFS"] : 2
            * conductor.dict_N_equation["NODOFS"],
            conductor.dict_N_equation["NODOFS"] : 2
            * conductor.dict_N_equation["NODOFS"],
        ] = (conductor.grid_features["delta_z"][ii] * (1.0 / 3.0 + ALFA) * MMAT)

        # COMPUTE THE CONVECTION MATRIX
        # array smart (cdp, 07/2020)
        ELAMAT[
            :conductor.dict_N_equation["NODOFS"],
            :conductor.dict_N_equation["NODOFS"],
        ] = (-AMAT / 2.0)
        ELAMAT[
            :conductor.dict_N_equation["NODOFS"],
            conductor.dict_N_equation["NODOFS"] : 2
            * conductor.dict_N_equation["NODOFS"],
        ] = (AMAT / 2.0)
        ELAMAT[
            conductor.dict_N_equation["NODOFS"] : 2
            * conductor.dict_N_equation["NODOFS"],
            :conductor.dict_N_equation["NODOFS"],
        ] = (-AMAT / 2.0)
        ELAMAT[
            conductor.dict_N_equation["NODOFS"] : 2
            * conductor.dict_N_equation["NODOFS"],
            conductor.dict_N_equation["NODOFS"] : 2
            * conductor.dict_N_equation["NODOFS"],
        ] = (AMAT / 2.0)

        # COMPUTE THE DIFFUSION MATRIX
        # array smart (cdp, 07/2020)
        ELKMAT[
            :conductor.dict_N_equation["NODOFS"],
            :conductor.dict_N_equation["NODOFS"],
        ] = (KMAT / conductor.grid_features["delta_z"][ii])
        ELKMAT[
            :conductor.dict_N_equation["NODOFS"],
            conductor.dict_N_equation["NODOFS"] : 2
            * conductor.dict_N_equation["NODOFS"],
        ] = (-KMAT / conductor.grid_features["delta_z"][ii])
        ELKMAT[
            conductor.dict_N_equation["NODOFS"] : 2
            * conductor.dict_N_equation["NODOFS"],
            :conductor.dict_N_equation["NODOFS"],
        ] = (-KMAT / conductor.grid_features["delta_z"][ii])
        ELKMAT[
            conductor.dict_N_equation["NODOFS"] : 2
            * conductor.dict_N_equation["NODOFS"],
            conductor.dict_N_equation["NODOFS"] : 2
            * conductor.dict_N_equation["NODOFS"],
        ] = (KMAT / conductor.grid_features["delta_z"][ii])

        # COMPUTE THE SOURCE MATRIX
        # array smart (cdp, 07/2020)
        ELSMAT[
            :conductor.dict_N_equation["NODOFS"],
            :conductor.dict_N_equation["NODOFS"],
        ] = (SMAT * conductor.grid_features["delta_z"][ii] / 3.0)
        ELSMAT[
            :conductor.dict_N_equation["NODOFS"],
            conductor.dict_N_equation["NODOFS"] : 2
            * conductor.dict_N_equation["NODOFS"],
        ] = (SMAT * conductor.grid_features["delta_z"][ii] / 6.0)
        ELSMAT[
            conductor.dict_N_equation["NODOFS"] : 2
            * conductor.dict_N_equation["NODOFS"],
            :conductor.dict_N_equation["NODOFS"],
        ] = (SMAT * conductor.grid_features["delta_z"][ii] / 6.0)
        ELSMAT[
            conductor.dict_N_equation["NODOFS"] : 2
            * conductor.dict_N_equation["NODOFS"],
            conductor.dict_N_equation["NODOFS"] : 2
            * conductor.dict_N_equation["NODOFS"],
        ] = (SMAT * conductor.grid_features["delta_z"][ii] / 3.0)

        # COMPUTE THE SOURCE VECTOR (ANALYTIC INTEGRATION)
        # array smart (cdp, 07/2020)
        # This is independent from the solution method thanks to the escamotage of \
        # the dummy steady state corresponding to the initialization (cdp, 10/2020)
        if conductor.cond_num_step == 1:
            # Current time step (cdp, 10/2020)
            ELSLOD.present[:conductor.dict_N_equation["NODOFS"]] = (
                conductor.grid_features["delta_z"][ii]
                / 6.0
                * (2.0 * SVEC.present[:, 0] + SVEC.present[:, 1])
            )
            ELSLOD.present[
                conductor.dict_N_equation["NODOFS"] : 2
                * conductor.dict_N_equation["NODOFS"]
            ] = (
                conductor.grid_features["delta_z"][ii]
                / 6.0
                * (SVEC.present[:, 0] + 2.0 * SVEC.present[:, 1])
            )
            # Previous time step (cdp, 10/2020)
            ELSLOD.previous[:conductor.dict_N_equation["NODOFS"]] = (
                conductor.grid_features["delta_z"][ii]
                / 6.0
                * (2.0 * SVEC.previous[:, 0] + SVEC.previous[:, 1])
            )
            ELSLOD.previous[
                conductor.dict_N_equation["NODOFS"] : 2
                * conductor.dict_N_equation["NODOFS"]
            ] = (
                conductor.grid_features["delta_z"][ii]
                / 6.0
                * (SVEC.previous[:, 0] + 2.0 * SVEC.previous[:, 1])
            )
        else:
            # Compute only at the current time step (cdp, 10/2020)
            ELSLOD[:conductor.dict_N_equation["NODOFS"]] = (
                conductor.grid_features["delta_z"][ii]
                / 6.0
                * (2.0 * SVEC[:, 0] + SVEC[:, 1])
            )
            ELSLOD[
                conductor.dict_N_equation["NODOFS"] : 2
                * conductor.dict_N_equation["NODOFS"]
            ] = (
                conductor.grid_features["delta_z"][ii]
                / 6.0
                * (SVEC[:, 0] + 2.0 * SVEC[:, 1])
            )
        # ASSEMBLE THE MATRICES AND THE LOAD VECTOR
        # array smart (cdp, 07/2020) check ok
        for iii in range(conductor.dict_band["Half"]):
            MASMAT[
                conductor.dict_band["Half"]
                - iii
                - 1 : conductor.dict_band["Full"]
                - iii,
                jump + iii,
            ] = (
                MASMAT[
                    conductor.dict_band["Half"]
                    - iii
                    - 1 : conductor.dict_band["Full"]
                    - iii,
                    jump + iii,
                ]
                + ELMMAT[iii, :]
            )
            FLXMAT[
                conductor.dict_band["Half"]
                - iii
                - 1 : conductor.dict_band["Full"]
                - iii,
                jump + iii,
            ] = (
                FLXMAT[
                    conductor.dict_band["Half"]
                    - iii
                    - 1 : conductor.dict_band["Full"]
                    - iii,
                    jump + iii,
                ]
                + ELAMAT[iii, :]
            )
            DIFMAT[
                conductor.dict_band["Half"]
                - iii
                - 1 : conductor.dict_band["Full"]
                - iii,
                jump + iii,
            ] = (
                DIFMAT[
                    conductor.dict_band["Half"]
                    - iii
                    - 1 : conductor.dict_band["Full"]
                    - iii,
                    jump + iii,
                ]
                + ELKMAT[iii, :]
            )
            SORMAT[
                conductor.dict_band["Half"]
                - iii
                - 1 : conductor.dict_band["Full"]
                - iii,
                jump + iii,
            ] = (
                SORMAT[
                    conductor.dict_band["Half"]
                    - iii
                    - 1 : conductor.dict_band["Full"]
                    - iii,
                    jump + iii,
                ]
                + ELSMAT[iii, :]
            )
        if (
            conductor.inputs["METHOD"] == "BE"
            or conductor.inputs["METHOD"] == "CN"
        ):
            # Backward Euler or Crank-Nicolson (cdp, 10, 2020)
            if conductor.cond_num_step == 1:
                # Construct key SYSLOD of dictionary dict_Step (cdp, 10/2020)
                # Current time step (cdp, 10/2020)
                conductor.dict_Step["SYSLOD"][
                    jump : jump + conductor.dict_band["Half"], 0
                ] = (
                    ELSLOD.present
                    + conductor.dict_Step["SYSLOD"][
                        jump : jump + conductor.dict_band["Half"], 0
                    ]
                )
                # Previous time step (cdp, 10/2020)
                conductor.dict_Step["SYSLOD"][
                    jump : jump + conductor.dict_band["Half"], 1
                ] = (
                    ELSLOD.previous
                    + conductor.dict_Step["SYSLOD"][
                        jump : jump + conductor.dict_band["Half"], 1
                    ]
                )
            else:
                # Update only the first column, that correspond to the current time \
                # step (cdp, 10/2020)
                conductor.dict_Step["SYSLOD"][
                    jump : jump + conductor.dict_band["Half"], 0
                ] = (
                    ELSLOD
                    + conductor.dict_Step["SYSLOD"][
                        jump : jump + conductor.dict_band["Half"], 0
                    ]
                )
        elif conductor.inputs["METHOD"] == "AM4":
            # Adams-Moulton order 4 (cdp, 10/2020)
            if conductor.cond_num_step == 1:
                # Construct key SYSLOD of dictionary dict_Step (cdp, 10/2020)
                # Current time step (cdp, 10/2020)
                conductor.dict_Step["SYSLOD"][
                    jump : jump + conductor.dict_band["Half"], 0
                ] = (
                    ELSLOD.present
                    + conductor.dict_Step["SYSLOD"][
                        jump : jump + conductor.dict_band["Half"], 0
                    ]
                )
                for cc in range(conductor.dict_Step["SYSLOD"].shape[1]):
                    # Dummy initial steady state (cdp, 10/2020)
                    conductor.dict_Step["SYSLOD"][
                        jump : jump + conductor.dict_band["Half"], cc
                    ] = (
                        ELSLOD.previous
                        + conductor.dict_Step["SYSLOD"][
                            jump : jump + conductor.dict_band["Half"], cc
                        ]
                    )
            else:
                # Shift the colums by one towards right and compute the new first \
                # column at the current time step (cdp, 10/2020)
                conductor.dict_Step["SYSLOD"][:, 1:4] = conductor.dict_Step["SYSLOD"][
                    :, 0:3
                ].copy()
                conductor.dict_Step["SYSLOD"][
                    jump : jump + conductor.dict_band["Half"], 0
                ] = (
                    ELSLOD
                    + conductor.dict_Step["SYSLOD"][
                        jump : jump + conductor.dict_band["Half"], 0
                    ]
                )
            # end if conductor.cond_num_step (cdp, 10/2020)
        # end conductor.inputs["METHOD"] (cdp, 10/2020)

    # end for ii (cdp, 07/2020)
    # ** END MATRICES CONSTRUCTION (cdp, 07/2020) **

    # SCRIPT TO SAVE MATRICES MASMAT, FLXMAT, DIFMAT, SORMAT, SYSVAR, SYSLOD.

    if conductor.cond_num_step == 1 or np.isclose(conductor.Space_save[conductor.i_save],conductor.cond_time[-1]):

        path_matr = os.path.join("D:/refactoring/function_step", "matrices/before")
        os.makedirs(path_matr, exist_ok = True)

        if conductor.cond_num_step == 1:
            sfx = conductor.i_save - 1
        else:
            sfx = conductor.i_save

        masmat_fname = os.path.join(
            path_matr, f"MASMAT_{sfx}.tsv")
        flxmat_fname = os.path.join(
            path_matr, f"FLXMAT_{sfx}.tsv")
        difmat_fname = os.path.join(
            path_matr, f"DIFMAT_{sfx}.tsv")
        sormat_fname = os.path.join(
            path_matr, f"SORMAT_{sfx}.tsv")
        sysvar_fname = os.path.join(
            path_matr, f"SYSVAR_{sfx}.tsv")
        syslod_fname = os.path.join(
            path_matr, f"SYSLOD_{sfx}.tsv")
        with open(masmat_fname, "w") as writer:
            np.savetxt(writer, MASMAT, delimiter = "\t")
        with open(flxmat_fname, "w") as writer:
            np.savetxt(writer, FLXMAT, delimiter = "\t")
        with open(difmat_fname, "w") as writer:
            np.savetxt(writer, DIFMAT, delimiter = "\t")
        with open(sormat_fname, "w") as writer:
            np.savetxt(writer, SORMAT, delimiter = "\t")
        with open(sysvar_fname, "w") as writer:
            np.savetxt(writer, conductor.dict_Step["SYSVAR"], delimiter = "\t")
        with open(syslod_fname, "w") as writer:
            np.savetxt(writer, conductor.dict_Step["SYSLOD"], delimiter = "\t")

    # ** COMPUTE SYSTEM MATRIX **
    if conductor.inputs["METHOD"] == "BE" or conductor.inputs["METHOD"] == "CN":
        # Backward Euler or Crank-Nicolson (cdp, 10, 2020)
        SYSMAT = MASMAT / conductor.time_step + conductor.theta_method * (
            FLXMAT + DIFMAT + SORMAT
        )
    elif conductor.inputs["METHOD"] == "AM4":
        # Adams-Moulton order 4 (cdp, 10, 2020)
        if conductor.cond_num_step == 1:
            # This is due to the dummy initial steady state (cdp, 10/2020)
            for cc in range(conductor.dict_Step["AM4_AA"].shape[2]):
                conductor.dict_Step["AM4_AA"][:, :, cc] = FLXMAT + DIFMAT + SORMAT
        else:
            # Shift the matrices by one towards right and compute the new first \
            # matrix at the current time step (cdp, 10/2020)
            conductor.dict_Step["AM4_AA"][:, :, 1:4] = conductor.dict_Step["AM4_AA"][
                :, :, 0:3
            ]
            conductor.dict_Step["AM4_AA"][:, :, 0] = FLXMAT + DIFMAT + SORMAT
        # end if conductor.cond_num_step (cdp, 10/2020)
        # compute SYSMAT
        SYSMAT = (
            MASMAT / conductor.time_step
            + 9 / 24 * conductor.dict_Step["AM4_AA"][:, :, 0]
        )
    # end conductor.inputs["METHOD"] (cdp, 10/2020)
    
    # lines of code to save SYSMAT and SYSLOD in .tsv files
    if conductor.cond_num_step == 1 or np.isclose(conductor.Space_save[conductor.i_save],conductor.cond_time[-1]):
        SYSMAT_f_name = os.path.join(
            path_matr, f"SYSMAT_{sfx}.tsv")
        with open(SYSMAT_f_name, "w") as writer:
            np.savetxt(writer, SYSMAT, delimiter = "\t")

    # ADD THE LOAD CONTRIBUTION FROM PREVIOUS STEP
    # array smart optimization (cdp, 07/2020)
    # check optimization (cdp, 09/2020)
    for I in range(conductor.dict_N_equation["Total"]):
        if I <= conductor.dict_band["Half"] - 1:
            # remember that arange stops before the stop value: \
            # last value = stop - step (cdp, 07/2020)
            J = np.arange(
                start=0, stop=conductor.dict_band["Half"] + I, step=1, dtype=int
            )
        elif I >= conductor.dict_N_equation["Total"] - (
            conductor.dict_band["Half"] - 1
        ):
            J = np.arange(
                start=I - (conductor.dict_band["Half"] - 1),
                stop=conductor.dict_N_equation["Total"],
                step=1,
                dtype=int,
            )
        else:
            J = np.arange(
                start=I - (conductor.dict_band["Half"] - 1),
                stop=I + conductor.dict_band["Half"],
                step=1,
                dtype=int,
            )
        JJ = J - I + conductor.dict_band["Half"] - 1
        if (
            conductor.inputs["METHOD"] == "BE"
            or conductor.inputs["METHOD"] == "CN"
        ):
            # Backward Euler or Crank-Nicolson (cdp, 10, 2020)
            # Matrix vector product contribution (cdp, 10, 2020)
            Known[I] = np.sum(
                (
                    MASMAT[JJ, I] / conductor.time_step
                    - (1.0 - conductor.theta_method)
                    * (FLXMAT[JJ, I] + DIFMAT[JJ, I] + SORMAT[JJ, I])
                )
                * conductor.dict_Step["SYSVAR"][J, 0]
            )
        elif conductor.inputs["METHOD"] == "AM4":
            # Adams-Moulton order 4 (cdp, 10, 2020)
            # Matrices vectors product contribution (cdp, 10, 2020)
            Known[I] = np.sum(
                (
                    MASMAT[JJ, I] / conductor.time_step
                    - 19 / 24 * conductor.dict_Step["AM4_AA"][JJ, I, 1]
                )
                * conductor.dict_Step["SYSVAR"][J, 0]
                + 5
                / 24
                * conductor.dict_Step["AM4_AA"][JJ, I, 2]
                * conductor.dict_Step["SYSVAR"][J, 1]
                - 1
                / 24
                * conductor.dict_Step["AM4_AA"][JJ, I, 3]
                * conductor.dict_Step["SYSVAR"][J, 2]
            )
    # end for I (cdp, 07/2020)
    if conductor.inputs["METHOD"] == "BE" or conductor.inputs["METHOD"] == "CN":
        # Backward Euler or Crank-Nicolson (cdp, 10, 2020)
        # External sources (SYSLOD) contribution (cdp, 10, 2020)
        Known = (
            Known
            + conductor.theta_method * conductor.dict_Step["SYSLOD"][:, 0]
            + (1.0 - conductor.theta_method) * conductor.dict_Step["SYSLOD"][:, 1]
        )
    elif conductor.inputs["METHOD"] == "AM4":
        # Adams-Moulton order 4 (cdp, 10, 2020)
        # External sources (SYSLOD) contribution (cdp, 10, 2020)
        Known = (
            Known
            + 9 / 24 * conductor.dict_Step["SYSLOD"][:, 0]
            + 19 / 24 * conductor.dict_Step["SYSLOD"][:, 1]
            - 5 / 24 * conductor.dict_Step["SYSLOD"][:, 2]
            + 1 / 24 * conductor.dict_Step["SYSLOD"][:, 3]
        )

    # lines of code to save SYSLOD in .tsv files
    if conductor.cond_num_step == 1 or np.isclose(conductor.Space_save[conductor.i_save],conductor.cond_time[-1]):
        SYSLOD_f_name = os.path.join(
            path_matr, f"SYSLOD_{sfx}.tsv")
        with open(SYSLOD_f_name, "w") as writer:
            np.savetxt(writer, conductor.dict_Step["SYSLOD"], delimiter = "\t")

    # IMPOSE BOUNDARY CONDITIONS AT INLET/OUTLET
    for jj, fluid_comp_j in enumerate(conductor.inventory["FluidComponent"].collection):
        INTIAL = fluid_comp_j.coolant.operations["INTIAL"]
        # index for inlet BC (cdp, 08/2020)
        Iiv_inl = dict(
            forward=jj,
            backward=conductor.dict_N_equation["Total"]
            - conductor.dict_N_equation["NODOFS"]
            + jj,
        )
        Iip_inl = dict(
            forward=jj + conductor.inventory["FluidComponent"].number,
            backward=conductor.dict_N_equation["Total"]
            - conductor.dict_N_equation["NODOFS"]
            + jj
            + conductor.inventory["FluidComponent"].number,
        )
        Iit_inl = dict(
            forward=jj + 2 * conductor.inventory["FluidComponent"].number,
            backward=conductor.dict_N_equation["Total"]
            - conductor.dict_N_equation["NODOFS"]
            + jj
            + 2 * conductor.inventory["FluidComponent"].number,
        )

        # index for outlet BC (cdp, 08/2020)
        Iiv_out = dict(
            forward=jj
            + conductor.dict_N_equation["Total"]
            - conductor.dict_N_equation["NODOFS"],
            backward=jj,
        )
        Iip_out = dict(
            forward=jj
            + conductor.dict_N_equation["Total"]
            - conductor.dict_N_equation["NODOFS"]
            + conductor.inventory["FluidComponent"].number,
            backward=jj + conductor.inventory["FluidComponent"].number,
        )
        Iit_out = dict(
            forward=jj
            + conductor.dict_N_equation["Total"]
            - conductor.dict_N_equation["NODOFS"]
            + 2 * conductor.inventory["FluidComponent"].number,
            backward=jj + 2 * conductor.inventory["FluidComponent"].number,
        )
        if abs(INTIAL) == 1:
            if INTIAL == 1:
                # inlet pressure (cdp, 07/2020)
                p_inl = fluid_comp_j.coolant.operations["PREINL"]
                # inlet temperature (cdp, 07/2020)
                T_inl = fluid_comp_j.coolant.operations["TEMINL"]
                # outlet pressure (cdp, 07/2020)
                p_out = fluid_comp_j.coolant.operations["PREOUT"]
                # outlet temperature: to be assigned if outlet velocity is negative \
                # (cdp, 08/2020)
                T_out = fluid_comp_j.coolant.operations["TEMOUT"]
            else:
                # get from file with interpolation in time (cdp,07/2020)
                [flow_par, flagSpecfield] = get_from_xlsx(
                    conductor, path, fluid_comp_j, "INTIAL", INTIAL
                )
                print(
                    f"""flagSpecfield == {flagSpecfield}: still to be decided if 
                it is useful and if yes still to be defined\n"""
                )
                # inlet pressure (cdp, 07/2020)
                p_inl = flow_par[2]
                # inlet temperature (cdp, 07/2020)
                T_inl = flow_par[0]
                # outlet pressure (cdp, 07/2020)
                p_out = flow_par[3]
                # outlet temperature: to be assigned if outlet velocity is negative \
                # (cdp, 08/2020)
                T_out = flow_par[1]
            # Assign BC: (cdp, 08/2020)
            # p_inl
            SYSMAT[
                :conductor.dict_band["Full"],
                Iip_inl[fluid_comp_j.channel.flow_dir[0]],
            ] = 0.0
            # main diagonal (cdp, 08/2020)
            SYSMAT[
                conductor.dict_band["Half"] - 1,
                Iip_inl[fluid_comp_j.channel.flow_dir[0]],
            ] = 1.0
            Known[Iip_inl[fluid_comp_j.channel.flow_dir[0]]] = p_inl
            # p_out
            SYSMAT[
                :conductor.dict_band["Full"],
                Iip_out[fluid_comp_j.channel.flow_dir[0]],
            ] = 0.0
            # main diagonal (cdp, 08/2020)
            SYSMAT[
                conductor.dict_band["Half"] - 1,
                Iip_out[fluid_comp_j.channel.flow_dir[0]],
            ] = 1.0
            Known[Iip_out[fluid_comp_j.channel.flow_dir[0]]] = p_out
            if fluid_comp_j.channel.flow_dir[0] == "forward":
                # T_inl
                if fluid_comp_j.coolant.dict_node_pt["velocity"][0] > 0:
                    SYSMAT[
                        :conductor.dict_band["Full"],
                        Iit_inl[fluid_comp_j.channel.flow_dir[0]],
                    ] = 0.0
                    # main diagonal (cdp, 08/2020)
                    SYSMAT[
                        conductor.dict_band["Half"] - 1,
                        Iit_inl[fluid_comp_j.channel.flow_dir[0]],
                    ] = 1.0
                    Known[Iit_inl[fluid_comp_j.channel.flow_dir[0]]] = T_inl
                # T_out
                if fluid_comp_j.coolant.dict_node_pt["velocity"][-1] < 0:
                    SYSMAT[
                        :conductor.dict_band["Full"],
                        Iit_out[fluid_comp_j.channel.flow_dir[0]],
                    ] = 0.0
                    # main diagonal (cdp, 08/2020)
                    SYSMAT[
                        conductor.dict_band["Half"] - 1,
                        Iit_out[fluid_comp_j.channel.flow_dir[0]],
                    ] = 1.0
                    Known[Iit_out[fluid_comp_j.channel.flow_dir[0]]] = T_out
            elif fluid_comp_j.channel.flow_dir[0] == "backward":
                # T_inl
                if fluid_comp_j.coolant.dict_node_pt["velocity"][-1] < 0:
                    SYSMAT[
                        :conductor.dict_band["Full"],
                        Iit_inl[fluid_comp_j.channel.flow_dir[0]],
                    ] = 0.0
                    # main diagonal (cdp, 08/2020)
                    SYSMAT[
                        conductor.dict_band["Half"] - 1,
                        Iit_inl[fluid_comp_j.channel.flow_dir[0]],
                    ] = 1.0
                    Known[Iit_inl[fluid_comp_j.channel.flow_dir[0]]] = T_inl
                # T_out
                if fluid_comp_j.coolant.dict_node_pt["velocity"][0] > 0:
                    SYSMAT[
                        :conductor.dict_band["Full"],
                        Iit_out[fluid_comp_j.channel.flow_dir[0]],
                    ] = 0.0
                    # main diagonal (cdp, 08/2020)
                    SYSMAT[
                        conductor.dict_band["Half"] - 1,
                        Iit_out[fluid_comp_j.channel.flow_dir[0]],
                    ] = 1.0
                    Known[Iit_out[fluid_comp_j.channel.flow_dir[0]]] = T_out
            # End if fluid_comp_j.channel.flow_dir[0] == "forward"
        elif abs(INTIAL) == 2:
            # INLET AND OUTLET RESERVOIRS, INLET CONDITIONS AND FLOW SPECIFIED
            if INTIAL == 2:
                # inlet mass flow rate (cdp, 07/2020)
                MDTIN = fluid_comp_j.coolant.operations["MDTIN"]
                # inlet pressure (cdp, 07/2020)
                # p_inl = fluid_comp_j.coolant.operations["PREINL"]
                # inlet temperature (cdp, 07/2020)
                T_inl = fluid_comp_j.coolant.operations["TEMINL"]
                # outlet pressure (cdp, 10/2020)
                p_out = fluid_comp_j.coolant.operations["PREOUT"]
                # outlet temperature: to be assigned if outlet velocity is negative \
                # (cdp, 08/2020)
                T_out = fluid_comp_j.coolant.operations["TEMOUT"]
            else:  # N.B va aggiustato per renderlo conforme caso positivo! \
                # (cdp, 10/2020)
                # all values from flow_dummy.xlsx: call get_from_xlsx (cdp, 07/2020)
                [flow_par, flagSpecfield] = get_from_xlsx(
                    conductor, path, fluid_comp_j, "INTIAL", INTIAL
                )
                print(
                    f"""flagSpecfield == {flagSpecfield}: still to be decided if it 
              useful and if yes still to be defined\n"""
                )
                MDTIN = flow_par[3]  # inlet mass flow rate (cdp, 07/2020)
                p_inl = flow_par[2]  # inlet pressure (cdp, 07/2020)
                T_inl = flow_par[0]  # inlet temperature (cdp, 07/2020)
                # outlet temperature: to be assigned if outlet velocity is negative \
                # (cdp, 08/2020)
                T_out = flow_par[1]
            # Assign BC: (cdp, 08/2020)
            # v_inl
            SYSMAT[
                :conductor.dict_band["Full"],
                Iiv_inl[fluid_comp_j.channel.flow_dir[0]],
            ] = 0.0
            # main diagonal (cdp, 08/2020)
            SYSMAT[
                conductor.dict_band["Half"] - 1,
                Iiv_inl[fluid_comp_j.channel.flow_dir[0]],
            ] = 1.0
            if fluid_comp_j.channel.flow_dir[0] == "forward":
                # Flow direction from x = 0 to x = L.
                Known[Iiv_inl[fluid_comp_j.channel.flow_dir[0]]] = (
                    MDTIN
                    / fluid_comp_j.coolant.dict_node_pt["total_density"][0]
                    / fluid_comp_j.channel.inputs["CROSSECTION"]
                )
            elif fluid_comp_j.channel.flow_dir[0] == "backward":
                # Flow direction from x = L to x = 0.
                Known[Iiv_inl[fluid_comp_j.channel.flow_dir[0]]] = (
                    MDTIN
                    / fluid_comp_j.coolant.dict_node_pt["total_density"][-1]
                    / fluid_comp_j.channel.inputs["CROSSECTION"]
                )
            ## p_inl
            # SYSMAT[0:conductor.dict_band["Full"], Iip_inl] = 0.0
            ## main diagonal (cdp, 08/2020)
            # SYSMAT[conductor.dict_band["Half"] - 1, Iip_inl] = 1.0

            # p_out
            SYSMAT[
                :conductor.dict_band["Full"],
                Iip_out[fluid_comp_j.channel.flow_dir[0]],
            ] = 0.0
            # main diagonal (cdp, 08/2020)
            SYSMAT[
                conductor.dict_band["Half"] - 1,
                Iip_out[fluid_comp_j.channel.flow_dir[0]],
            ] = 1.0
            Known[Iip_out[fluid_comp_j.channel.flow_dir[0]]] = p_out
            # Known[Iip_inl] = p_inl
            if fluid_comp_j.channel.flow_dir[0] == "forward":
                # T_inl
                if fluid_comp_j.coolant.dict_node_pt["velocity"][0] > 0:
                    SYSMAT[
                        :conductor.dict_band["Full"],
                        Iit_inl[fluid_comp_j.channel.flow_dir[0]],
                    ] = 0.0
                    # main diagonal (cdp, 08/2020)
                    SYSMAT[
                        conductor.dict_band["Half"] - 1,
                        Iit_inl[fluid_comp_j.channel.flow_dir[0]],
                    ] = 1.0
                    Known[Iit_inl[fluid_comp_j.channel.flow_dir[0]]] = T_inl
                # T_out (T_inl if MDTIN < 0)
                if fluid_comp_j.coolant.dict_node_pt["velocity"][-1] < 0:
                    SYSMAT[
                        :conductor.dict_band["Full"],
                        Iit_out[fluid_comp_j.channel.flow_dir[0]],
                    ] = 0.0
                    # main diagonal (cdp, 08/2020)
                    SYSMAT[
                        conductor.dict_band["Half"] - 1,
                        Iit_out[fluid_comp_j.channel.flow_dir[0]],
                    ] = 1.0
                    Known[Iit_out[fluid_comp_j.channel.flow_dir[0]]] = T_out
            elif fluid_comp_j.channel.flow_dir[0] == "backward":
                # T_inl
                if fluid_comp_j.coolant.dict_node_pt["velocity"][-1] < 0:
                    SYSMAT[
                        :conductor.dict_band["Full"],
                        Iit_inl[fluid_comp_j.channel.flow_dir[0]],
                    ] = 0.0
                    # main diagonal (cdp, 08/2020)
                    SYSMAT[
                        conductor.dict_band["Half"] - 1,
                        Iit_inl[fluid_comp_j.channel.flow_dir[0]],
                    ] = 1.0
                    Known[Iit_inl[fluid_comp_j.channel.flow_dir[0]]] = T_inl
                # T_out
                if fluid_comp_j.coolant.dict_node_pt["velocity"][0] > 0:
                    SYSMAT[
                        :conductor.dict_band["Full"],
                        Iit_out[fluid_comp_j.channel.flow_dir[0]],
                    ] = 0.0
                    # main diagonal (cdp, 08/2020)
                    SYSMAT[
                        conductor.dict_band["Half"] - 1,
                        Iit_out[fluid_comp_j.channel.flow_dir[0]],
                    ] = 1.0
                    Known[Iit_out[fluid_comp_j.channel.flow_dir[0]]] = T_out
            # End if fluid_comp_j.channel.flow_dir[0] == "forward"
        elif INTIAL == 3:
            # INLET RESERVOIR AND CLOSED OUTLET (SYMMETRY)
            # p_inl
            SYSMAT[
                :conductor.dict_band["Full"],
                Iip_inl[fluid_comp_j.channel.flow_dir[0]],
            ] = 0.0
            # main diagonal (cdp, 08/2020)
            SYSMAT[
                conductor.dict_band["Half"] - 1,
                Iip_inl[fluid_comp_j.channel.flow_dir[0]],
            ] = 1.0
            Known[
                Iip_inl[fluid_comp_j.channel.flow_dir[0]]
            ] = fluid_comp_j.coolant.operations["PREINL"]
            # T_inl
            if fluid_comp_j.coolant.dict_node_pt["velocity"][0] > 0:
                SYSMAT[
                    :conductor.dict_band["Full"],
                    Iit_inl[fluid_comp_j.channel.flow_dir[0]],
                ] = 0.0
                # main diagonal (cdp, 08/2020)
                SYSMAT[
                    conductor.dict_band["Half"] - 1,
                    Iit_inl[fluid_comp_j.channel.flow_dir[0]],
                ] = 1.0
                Known[
                    Iit_inl[fluid_comp_j.channel.flow_dir[0]]
                ] = fluid_comp_j.coolant.operations["TEMINL"]
            # v_out
            SYSMAT[
                :conductor.dict_band["Full"],
                Iiv_out[fluid_comp_j.channel.flow_dir[0]],
            ] = 0.0
            SYSMAT[
                conductor.dict_band["Half"] - 1,
                Iiv_out[fluid_comp_j.channel.flow_dir[0]],
            ] = 1.0
            Known[
                Iiv_out[fluid_comp_j.channel.flow_dir[0]]
            ] = 0.0  # closed outlet (cdp, 08/2020)
            # T_out
            if fluid_comp_j.coolant.dict_node_pt["velocity"][-1] < 0:
                SYSMAT[
                    :conductor.dict_band["Full"],
                    Iit_out[fluid_comp_j.channel.flow_dir[0]],
                ] = 0.0
                # main diagonal (cdp, 08/2020)
                SYSMAT[
                    conductor.dict_band["Half"] - 1,
                    Iit_out[fluid_comp_j.channel.flow_dir[0]],
                ] = 1.0
                Known[
                    Iit_out[fluid_comp_j.channel.flow_dir[0]]
                ] = fluid_comp_j.coolant.operations["TEMOUT"]
        elif INTIAL == 4:
            # v_inl
            SYSMAT[
                :conductor.dict_band["Full"],
                Iiv_inl[fluid_comp_j.channel.flow_dir[0]],
            ] = 0.0
            # main diagonal (cdp, 08/2020)
            SYSMAT[
                conductor.dict_band["Half"] - 1,
                Iiv_inl[fluid_comp_j.channel.flow_dir[0]],
            ] = 1.0
            Known[
                Iiv_inl[fluid_comp_j.channel.flow_dir[0]]
            ] = 0.0  # closed inlet (cdp, 08/2020)
            # v_out
            SYSMAT[
                :conductor.dict_band["Full"],
                Iiv_out[fluid_comp_j.channel.flow_dir[0]],
            ] = 0.0
            # main diagonal (cdp, 08/2020)
            SYSMAT[
                conductor.dict_band["Half"] - 1,
                Iiv_out[fluid_comp_j.channel.flow_dir[0]],
            ] = 1.0
            Known[
                Iiv_out[fluid_comp_j.channel.flow_dir[0]]
            ] = 0.0  # closed outlet (cdp, 08/2020)
        elif abs(INTIAL) == 5:
            if INTIAL == 5:
                # inlet mass flow rate (cdp, 07/2020)
                MDTIN = fluid_comp_j.coolant.operations["MDTIN"]
                # inlet temperature (cdp, 07/2020)
                T_inl = fluid_comp_j.coolant.operations["TEMINL"]
                # outlet pressure (cdp, 07/2020)
                p_out = fluid_comp_j.coolant.operations["PREOUT"]
                # outlet temperature: to be assigned if outlet velocity is negative \
                # (cdp, 08/2020)
                T_out = fluid_comp_j.coolant.operations["TEMOUT"]
            else:
                # all values from flow_dummy.xlsx: call get_from_xlsx (cdp, 07/2020)
                [flow_par, flagSpecfield] = get_from_xlsx(
                    conductor, path, fluid_comp_j, "INTIAL", INTIAL
                )
                print(
                    f"""flagSpecfield == {flagSpecfield}: still to be decided if it
              is useful and if yes still to be defined\n"""
                )
                MDTIN = flow_par[3]  # inlet mass flow rate (cdp, 07/2020)
                T_inl = flow_par[0]  # inlet temperature (cdp, 07/2020)
                p_out = flow_par[2]  # outlet pressure (cdp, 07/2020)
                # outlet temperature: to be assigned if outlet velocity is negative \
                # (cdp, 08/2020)
                T_out = flow_par[1]
            # Assign BC: (cdp, 08/2020)
            # v_inl
            SYSMAT[
                :conductor.dict_band["Full"],
                Iiv_inl[fluid_comp_j.channel.flow_dir[0]],
            ] = 0.0
            # main diagonal (cdp, 08/2020)
            SYSMAT[
                conductor.dict_band["Half"] - 1,
                Iiv_inl[fluid_comp_j.channel.flow_dir[0]],
            ] = 1.0
            if fluid_comp_j.channel.flow_dir[0] == "forward":
                # Flow direction from x = 0 to x = L.
                Known[Iiv_inl[fluid_comp_j.channel.flow_dir[0]]] = (
                    MDTIN
                    / fluid_comp_j.coolant.dict_node_pt["total_density"][0]
                    / fluid_comp_j.channel.inputs["CROSSECTION"]
                )
            elif fluid_comp_j.channel.flow_dir[0] == "backward":
                # Flow direction from x = L to x = 0.
                Known[Iiv_inl[fluid_comp_j.channel.flow_dir[0]]] = (
                    MDTIN
                    / fluid_comp_j.coolant.dict_node_pt["total_density"][-1]
                    / fluid_comp_j.channel.inputs["CROSSECTION"]
                )

            # p_out
            SYSMAT[
                :conductor.dict_band["Full"],
                Iip_out[fluid_comp_j.channel.flow_dir[0]],
            ] = 0.0
            # main diagonal (cdp, 08/2020)
            SYSMAT[
                conductor.dict_band["Half"] - 1,
                Iip_out[fluid_comp_j.channel.flow_dir[0]],
            ] = 1.0
            Known[Iip_out[fluid_comp_j.channel.flow_dir[0]]] = p_out
            # Known[Iip_inl] = p_inl
            if fluid_comp_j.channel.flow_dir[0] == "forward":
                # T_inl
                if fluid_comp_j.coolant.dict_node_pt["velocity"][0] > 0:
                    SYSMAT[
                        :conductor.dict_band["Full"],
                        Iit_inl[fluid_comp_j.channel.flow_dir[0]],
                    ] = 0.0
                    # main diagonal (cdp, 08/2020)
                    SYSMAT[
                        conductor.dict_band["Half"] - 1,
                        Iit_inl[fluid_comp_j.channel.flow_dir[0]],
                    ] = 1.0
                    Known[Iit_inl[fluid_comp_j.channel.flow_dir[0]]] = T_inl
                # T_out (T_inl if MDTIN < 0)
                if fluid_comp_j.coolant.dict_node_pt["velocity"][-1] < 0:
                    SYSMAT[
                        :conductor.dict_band["Full"],
                        Iit_out[fluid_comp_j.channel.flow_dir[0]],
                    ] = 0.0
                    # main diagonal (cdp, 08/2020)
                    SYSMAT[
                        conductor.dict_band["Half"] - 1,
                        Iit_out[fluid_comp_j.channel.flow_dir[0]],
                    ] = 1.0
                    Known[Iit_out[fluid_comp_j.channel.flow_dir[0]]] = T_out
            elif fluid_comp_j.channel.flow_dir[0] == "backward":
                # T_inl
                if fluid_comp_j.coolant.dict_node_pt["velocity"][-1] < 0:
                    SYSMAT[
                        :conductor.dict_band["Full"],
                        Iit_inl[fluid_comp_j.channel.flow_dir[0]],
                    ] = 0.0
                    # main diagonal (cdp, 08/2020)
                    SYSMAT[
                        conductor.dict_band["Half"] - 1,
                        Iit_inl[fluid_comp_j.channel.flow_dir[0]],
                    ] = 1.0
                    Known[Iit_inl[fluid_comp_j.channel.flow_dir[0]]] = T_inl
                # T_out
                if fluid_comp_j.coolant.dict_node_pt["velocity"][0] > 0:
                    SYSMAT[
                        :conductor.dict_band["Full"],
                        Iit_out[fluid_comp_j.channel.flow_dir[0]],
                    ] = 0.0
                    # main diagonal (cdp, 08/2020)
                    SYSMAT[
                        conductor.dict_band["Half"] - 1,
                        Iit_out[fluid_comp_j.channel.flow_dir[0]],
                    ] = 1.0
                    Known[Iit_out[fluid_comp_j.channel.flow_dir[0]]] = T_out
            # End if fluid_comp_j.channel.flow_dir[0] == "forward"

    # end for jj

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
    # dict_Step (cdp, 10/2020)
    conductor.dict_Step["SYSVAR"][:, 0] = gbacsb(conductor, SYSMAT, Known)

    # SYSVAR = solve_banded((15, 15), SYSMAT, Known)

    # COMPUTE THE NORM OF THE SOLUTION AND OF THE SOLUTION CHANGE (START)
    # array smart optimization (cdp, 08/2020)

    SOL = np.zeros(conductor.dict_N_equation["Total"])
    CHG = np.zeros(conductor.dict_N_equation["Total"])
    EIG = np.zeros(conductor.dict_N_equation["Total"])

    for jj in range(conductor.inventory["FluidComponent"].number):
        # velocity (cdp, 08/2020)
        SOL[
            jj : conductor.dict_N_equation["Total"] : conductor.dict_N_equation[
                "NODOFS"
            ]
        ] = Known[
            jj : conductor.dict_N_equation["Total"] : conductor.dict_N_equation[
                "NODOFS"
            ]
        ]
        conductor.dict_norm["Solution"][jj] = np.sqrt(
            np.sum(
                SOL[
                    jj : conductor.dict_N_equation["Total"] : conductor.dict_N_equation[
                        "NODOFS"
                    ]
                ]
                ** 2
            )
        )
        # pressure (cdp, 08/2020)
        SOL[
            jj
            + conductor.inventory["FluidComponent"].number : conductor.dict_N_equation["Total"] : conductor.dict_N_equation["NODOFS"]
        ] = Known[
            jj
            + conductor.inventory["FluidComponent"].number : conductor.dict_N_equation["Total"] : conductor.dict_N_equation["NODOFS"]
        ]
        conductor.dict_norm["Solution"][
            jj + conductor.inventory["FluidComponent"].number
        ] = np.sqrt(
            np.sum(
                SOL[
                    jj
                    + conductor.inventory["FluidComponent"].number : conductor.dict_N_equation["Total"] : conductor.dict_N_equation[
                        "NODOFS"
                    ]
                ]
                ** 2
            )
        )
        # temperature (cdp, 08/2020)
        SOL[
            jj
            + 2
            * conductor.inventory["FluidComponent"].number : conductor.dict_N_equation["Total"] : conductor.dict_N_equation["NODOFS"]
        ] = Known[
            jj
            + 2
            * conductor.inventory["FluidComponent"].number : conductor.dict_N_equation["Total"] : conductor.dict_N_equation["NODOFS"]
        ]
        conductor.dict_norm["Solution"][
            jj + 2 * conductor.inventory["FluidComponent"].number
        ] = np.sqrt(
            np.sum(
                SOL[
                    jj
                    + 2
                    * conductor.inventory["FluidComponent"].number : conductor.dict_N_equation["Total"] : conductor.dict_N_equation[
                        "NODOFS"
                    ]
                ]
                ** 2
            )
        )
    for ll in range(conductor.inventory["SolidComponent"].number):
        # temperature (cdp, 08/2020)
        SOL[
            ll
            + conductor.dict_N_equation["FluidComponent"] : conductor.dict_N_equation[
                "Total"
            ] : conductor.dict_N_equation["NODOFS"]
        ] = Known[
            ll
            + conductor.dict_N_equation["FluidComponent"] : conductor.dict_N_equation[
                "Total"
            ] : conductor.dict_N_equation["NODOFS"]
        ]
        conductor.dict_norm["Solution"][
            ll + conductor.dict_N_equation["FluidComponent"]
        ] = np.sqrt(
            np.sum(
                SOL[
                    ll
                    + conductor.dict_N_equation[
                        "FluidComponent"
                    ] : conductor.dict_N_equation["Total"] : conductor.dict_N_equation[
                        "NODOFS"
                    ]
                ]
                ** 2
            )
        )

    # COMPUTE THE NORM OF THE SOLUTION (END)

    # COMPUTE THE NORM OF THE SOLUTION CHANGE, THE EIGENVALUES AND RECOVER THE \
    # VARIABLES FROM THE SYSTEM SOLUTION (START)

    # Those are arrays (cdp, 08/2020)
    CHG = Known - conductor.dict_Step["SYSVAR"][:, 0]
    EIG = abs(CHG / conductor.time_step) / (abs(SOL) + TINY)

    for jj, fluid_comp in enumerate(conductor.inventory["FluidComponent"].collection):
        # velocity (cdp, 08/2020)
        conductor.dict_norm["Change"][jj] = np.sqrt(
            np.sum(
                CHG[
                    jj : conductor.dict_N_equation["Total"] : conductor.dict_N_equation[
                        "NODOFS"
                    ]
                ]
                ** 2
            )
        )
        conductor.EQTEIG[jj] = max(
            EIG[
                jj : conductor.dict_N_equation["Total"] : conductor.dict_N_equation[
                    "NODOFS"
                ]
            ]
        )
        fluid_comp.coolant.dict_node_pt["velocity"] = conductor.dict_Step["SYSVAR"][
            jj : conductor.dict_N_equation["Total"] : conductor.dict_N_equation[
                "NODOFS"
            ],
            0,
        ].copy()
        # pressure (cdp, 08/2020)
        conductor.dict_norm["Change"][
            jj + conductor.inventory["FluidComponent"].number
        ] = np.sqrt(
            np.sum(
                CHG[
                    jj
                    + conductor.inventory["FluidComponent"].number : conductor.dict_N_equation["Total"] : conductor.dict_N_equation[
                        "NODOFS"
                    ]
                ]
                ** 2
            )
        )
        conductor.EQTEIG[
            jj + conductor.inventory["FluidComponent"].number
        ] = max(
            EIG[
                jj
                + conductor.inventory["FluidComponent"].number : conductor.dict_N_equation["Total"] : conductor.dict_N_equation[
                    "NODOFS"
                ]
            ]
        )
        fluid_comp.coolant.dict_node_pt["pressure"] = conductor.dict_Step["SYSVAR"][
            jj
            + conductor.inventory["FluidComponent"].number : conductor.dict_N_equation["Total"] : conductor.dict_N_equation[
                "NODOFS"
            ],
            0,
        ].copy()
        # temperature (cdp, 08/2020)
        conductor.dict_norm["Change"][
            jj + 2 * conductor.inventory["FluidComponent"].number
        ] = np.sqrt(
            np.sum(
                CHG[
                    jj
                    + 2
                    * conductor.inventory["FluidComponent"].number: conductor.dict_N_equation["Total"] : conductor.dict_N_equation[
                        "NODOFS"
                    ]
                ]
                ** 2
            )
        )
        conductor.EQTEIG[
            jj + 2 * conductor.inventory["FluidComponent"].number
        ] = max(
            EIG[
                jj
                + 2
                * conductor.inventory["FluidComponent"].number : conductor.dict_N_equation["Total"] : conductor.dict_N_equation[
                    "NODOFS"
                ]
            ]
        )
        fluid_comp.coolant.dict_node_pt["temperature"] = conductor.dict_Step["SYSVAR"][
            jj
            + 2
            * conductor.inventory["FluidComponent"].number : conductor.dict_N_equation["Total"] : conductor.dict_N_equation[
                "NODOFS"
            ],
            0,
        ].copy()
        fluid_comp.coolant.dict_Gauss_pt["temperature_change"] = (
            fluid_comp.coolant.dict_node_pt["temperature"][:-1]
            + fluid_comp.coolant.dict_node_pt["temperature"][1:]
        ) / 2.0 - old_temperature_gauss[fluid_comp.identifier]
    for ll, comp in enumerate(conductor.inventory["SolidComponent"].collection):
        # temperature (cdp, 08/2020)
        conductor.dict_norm["Change"][
            ll + conductor.dict_N_equation["FluidComponent"]
        ] = np.sqrt(
            np.sum(
                CHG[
                    ll
                    + conductor.dict_N_equation[
                        "FluidComponent"
                    ] : conductor.dict_N_equation["Total"] : conductor.dict_N_equation[
                        "NODOFS"
                    ]
                ]
                ** 2
            )
        )
        conductor.EQTEIG[ll + conductor.dict_N_equation["FluidComponent"]] = max(
            EIG[
                ll
                + conductor.dict_N_equation[
                    "FluidComponent"
                ] : conductor.dict_N_equation["Total"] : conductor.dict_N_equation[
                    "NODOFS"
                ]
            ]
        )
        comp.dict_node_pt["temperature"] = conductor.dict_Step["SYSVAR"][
            ll
            + conductor.dict_N_equation["FluidComponent"] : conductor.dict_N_equation[
                "Total"
            ] : conductor.dict_N_equation["NODOFS"],
            0,
        ].copy()

        comp.dict_Gauss_pt["temperature_change"] = (
            comp.dict_node_pt["temperature"][:-1] + comp.dict_node_pt["temperature"][1:]
        ) / 2.0 - old_temperature_gauss[comp.identifier]

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
