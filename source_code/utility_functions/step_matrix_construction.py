import numpy as np

from collections import namedtuple
from typing import Union, NamedTuple

from fluid_component import FluidComponent
from solid_component import SolidComponent
from conductor import Conductor

def matrix_initialization(row:int,col:int)->tuple:
    """Wrapper of function np.zeros that inizializes five identical rectangular matrices.

    Args:
        row (int): number of rows of the matrix.
        col (int): number of columns of the matrix.

    Returns:
        tuple: collection of initialized matrices.
    """
    
    return tuple(np.zeros((row,col)) for _ in range(5))

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
    return (
        *tuple(np.zeros((dimension,dimension)) for _ in range(4)),
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
        Array = namedtuple("Array",("previous","present"))
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

def build_amat(
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
    density = f_comp.coolant.dict_Gauss_pt["total_density"][elem_idx]
    
    # Build array to assign diagonal coefficients.
    diag_idx = np.array(eq_idx)
    
    # Set diagonal elements (exploit broadcasting).
    matrix[diag_idx,diag_idx] = f_comp.coolant.dict_Gauss_pt["velocity"][elem_idx]

    # Set off diagonal coefficients.
    # from velocity equation.
    matrix[eq_idx.velocity, eq_idx.pressure] = 1. / density
    # from pressure equation.
    matrix[eq_idx.pressure, eq_idx.velocity] = (
        density
        * f_comp.coolant.dict_Gauss_pt["total_speed_of_sound"][elem_idx] ** 2.
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
    # The loop is done on all the type of interfaces between fluids since in 
    # the evaluation of K', K'' and K''' the relevant parameter is the open 
    # cross section per unit length (open contact perimeter): K', K'' and K''' 
    # are not 0 only if the open contact perimeter is not 0. In this way, a 
    # check on the interface kind between fluids is avoided for each interface 
    # and for each time step.
    for interface in conductor.interface.fluid_fluid:
        comp_1_pressure = interface.comp_1.coolant.dict_Gauss_pt["pressure"]
        comp_2_pressure = interface.comp_2.coolant.dict_Gauss_pt["pressure"]

        # constuct recurrent coefficients of matrix S elements.
        for key in key_names:
            # Initialization.
            conductor.dict_Gauss_pt[key][interface.interf_name] = np.zeros_like(
                conductor.grid_features["zcoord_gauss"]
            )
        # Evaluate pressure difference bethween comp_1 and comp_2
        delta_p = np.abs(comp_1_pressure - comp_2_pressure)
        # Array smart
        delta_p[delta_p < conductor.Delta_p_min] = conductor.Delta_p_min

        # Find index such that # P_comp_2 < P_comp_1.
        ind_1 = np.nonzero(comp_2_pressure < comp_1_pressure)[0]
        # Find index such that P_comp_2 >= P_comp_1.
        ind_2 = np.nonzero(comp_2_pressure >= comp_1_pressure)[0]
        
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
    velocity = comp.coolant.dict_Gauss_pt["velocity"][index]

    # K' evaluation [ms]:
    # K' = A_othogonal*sqrt(2*density/k_loc*abs(Delta_p))
    K1 = (
        conductor.dict_interf_peri["ch_ch"]["Open"][interf_name]
        * np.sqrt(
            2.
            * comp.coolant.dict_Gauss_pt["total_density"][index]
            / (conductor.k_loc * delta_p[index])
        )
    )
    
    # K'' evaluation [m^2]:
    # K'' = K'*lambda_v*velocity
    K2 = K1 * conductor.lambda_v * velocity

    # K''' evaluation [m^3/s]:
    # K''' = K'*(enthalpy + (velocity*lambda_v)^2/2)
    K3 = (
        K1
        * (
            comp.coolant.dict_Gauss_pt["total_enthalpy"][index]
            + .5 * (velocity * conductor.lambda_v) ** 2.
        )
    )

    # Assing evaluated K1, K2 and K3 to correspondig key in conductor attribute 
    # dict_Gauss_pt.
    for key, value in zip(("K1","K2","K3"),(K1,K2,K3)):
        conductor.dict_Gauss_pt[key][interf_name][index] = value

    return conductor

def build_kmat_fluid(
    matrix:np.ndarray,
    upweqt:np.ndarray,
    f_comp:FluidComponent,
    conductor:Conductor,
    elem_idx:int,
    )->np.ndarray:

    """Function that builds the K matrix (KMAT) at the Gauss point, UPWIND is included.
    UPWIND: differencing contribution a' la finite-difference, necessary to guarantee numarical stability that does not came from KMAT algebraic construction.

    Args:
        matrix (np.ndarray): initialized K matrix (np.zeros)
        upweqt (np.ndarray): array with the upwind numerical scheme.
        f_comp (FluidComponent): fluid component object from which get all info to buld the coefficients.
        conductor (Conductor): object with all the information of the conductor.
        elem_idx (int): index of the i-th element of the spatial discretization.

    Returns:
        np.ndarray: matrix with updated elements.
    """

    # Alias
    # Fluid velocity at present gauss point index.
    velocity = np.abs(f_comp.coolant.dict_Gauss_pt["velocity"][elem_idx])
    # Length of the interval that includes the present gauss points.
    delta_z = conductor.grid_features["delta_z"][elem_idx]
    # Collection of fluid equation index (velocity, pressure and temperaure 
    # equations).
    eq_idx = conductor.equation_index[f_comp.identifier]

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

    # Reference value for f_comp.coolant.dict_Gauss_pt["velocity"][elem_idx] 
    # (shallow copy).
    velocity = f_comp.coolant.dict_Gauss_pt["velocity"][elem_idx]
    # velocity equation: main diagonal elements construction
    # (j,j) [vel_j]
    matrix[eq_idx.velocity,eq_idx.velocity] = (
        2.0
        # dict_friction_factor[False]["total"]: total friction factor in Gauss 
        # points (see __init__ of class Channel for details).
        * f_comp.channel.dict_friction_factor[False]["total"][elem_idx]
        * np.abs(velocity) / f_comp.channel.inputs["HYDIAMETER"]
    )
    
    # pressure equation: elements below main diagonal construction
    # (j+num_fluid_components,0:num_fluid_components) [Pres]
    matrix[eq_idx.pressure,eq_idx.velocity] = (
        - matrix[eq_idx.velocity,eq_idx.velocity]
        * f_comp.coolant.dict_Gauss_pt["Gruneisen"][elem_idx]
        * f_comp.coolant.dict_Gauss_pt["total_density"][elem_idx]
        * velocity
    )
    
    # temperature equation: elements below main diagonal construction
    # (j+2*num_fluid_components,0:num_fluid_components) [Temp]
    matrix[eq_idx.temperature,eq_idx.velocity] = (
        - matrix[eq_idx.velocity,eq_idx.velocity]
        / f_comp.coolant.dict_Gauss_pt["total_isochoric_specific_heat"][elem_idx]
        * velocity
    )

    return matrix

def build_smat_fluid_interface(
    matrix:np.ndarray,
    conductor:Conductor,
    elem_idx:int,
    )->np.ndarray:

    """Function that builds the S matrix (SMAT) therms due to fluid component interfaces at the Gauss point (SOURCE JACOBIAN).

    Args:
        matrix (np.ndarray): S matrix after call to function buld_smat_fluid.
        conductor (Conductor): object with all the information of the conductor.
        elem_idx (int): index of the i-th element of the spatial discretization.

    Returns:
        np.ndarray: matrix with updated elements.
    """

    # NOMENCLATURE
    # h: heat transfer coefficient (_o: open; _c:close)
    # P: contact perimeter (_o: open; _c:close)
    
    # Alias
    # Collection of NamedTuple with fluid equation index (velocity, pressure 
    # and temperaure equations).
    eq_idx = conductor.equation_index

    for interface in conductor.interface.fluid_fluid:
        
        K1 = conductor.dict_Gauss_pt["K1"][interface.interf_name][elem_idx]
        K2 = conductor.dict_Gauss_pt["K2"][interface.interf_name][elem_idx]
        K3 = conductor.dict_Gauss_pt["K3"][interface.interf_name][elem_idx]

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

        # Fill rows of comp_1, columns involving comp_1 and comp_2.
        matrix = __smat_fluid_interface(
            matrix,
            interface.comp_1,
            interface.comp_2,
            elem_idx,
            eq_idx,
            K1=K1,
            K2=K2,
            K3=K3,
            coef_htc=coef_htc,
        )
        # Fill rows of comp_2, columns involving comp_2 and comp_1.
        matrix = __smat_fluid_interface(
            matrix,
            interface.comp_2,
            interface.comp_1,
            elem_idx,
            eq_idx,
            K1=K1,
            K2=K2,
            K3=K3,
            coef_htc=coef_htc,
        )

    return matrix

def __smat_fluid_interface(
    matrix:np.ndarray,
    comp_1:FluidComponent,
    comp_2:FluidComponent,
    elem_idx:int,
    eq_idx:dict,
    **kwargs
    )->np.ndarray:
    """Function that evaluates the S matrix (SMAT) terms due to fluid component interfaces at the Gauss point (SOURCE JACOBIAN) by rows (horizontally) corresponding to the equations of comp_1 for the columns of comp_1 including the contributions given by comp_2.

    Args:
        matrix (np.ndarray): S matrix after call to function buld_smat_fluid.
        comp_1 (FluidComponent): fluid component object from which get all info to build the coefficients.
        comp_2 (FluidComponent): fluid component object from which get all info to build the coefficients.
        elem_idx (int): index of the i-th element of the spatial discretization.
        eq_idx (dict): collection of NamedTuple with fluid equation index (velocity, pressure and temperaure equations).
    
    Kwargs:
        K1 (np.ndarray): transport coefficient K'.
        K2 (np.ndarray): transport coefficient K''.
        K3 (np.ndarray): transport coefficient K'''.
        coef_htc (float): heat transfer coefficient per unit of length evauated as the sum of the heat transfer coefficient of the open interface time the corresponding contact perimeter and the heat transfer coefficient of the close intefrace and the corresponding contact perimeter 
            P_o * h_o + P_c * h_c.

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

    # Alias
    comp_1_v = comp_1.coolant.dict_Gauss_pt["velocity"][elem_idx]
    comp_1_rho = comp_1.coolant.dict_Gauss_pt["total_density"][elem_idx]
    comp_1_A = comp_1.channel.inputs["CROSSECTION"]
    comp_1_phi = comp_1.coolant.dict_Gauss_pt["Gruneisen"][elem_idx]
    comp_1_enthalpy = comp_1.coolant.dict_Gauss_pt["total_enthalpy"][elem_idx]
    comp_1_cv = comp_1.coolant.dict_Gauss_pt["total_isochoric_specific_heat"][
        elem_idx
    ]
    K1 = kwargs["K1"]
    K2 = kwargs["K2"]
    K3 = kwargs["K3"]
    coef_htc = kwargs["coef_htc"]

    # VELOCITY EQUATION: above/below main diagonal elements construction:
    # (j,j+num_fluid_components) [Pres_j]

    # s_vj_pj = (K1 * v - K2) / (A * rho)

    s_vj_pj = (K1 * comp_1_v - K2) / (comp_1_A * comp_1_rho)

    matrix[
        eq_idx[comp_1.identifier].velocity,
        eq_idx[comp_1.identifier].pressure,
    ] -= s_vj_pj

    # (j,k + num_fluid_components:2*num_fluid_components) 
    # [Pres_k]
    matrix[
        eq_idx[comp_1.identifier].velocity,
        eq_idx[comp_2.identifier].pressure,
    ] = s_vj_pj

    # PRESSURE EQUATION: main diagonal elements construction:
    # (j+num_fluid_components,j+num_fluid_components) [Pres_j]

    # coef_grun_area = phi / A
    coef_grun_area = comp_1_phi / comp_1_A

    # s_pj_pj = phi/A * [K3 - vK2 - (w - v^2/2 - c0^2/phi)K1]
    #         = coef_grun_area * [K3 - vK2 - (w - v^2/2 - c0^2/phi)K1]
    s_pj_pj = (
        coef_grun_area
        * (K3 - comp_1_v * K2 - (comp_1_enthalpy - comp_1_v ** 2. / 2.
        - comp_1.coolant.dict_Gauss_pt["total_speed_of_sound"][elem_idx] ** 2. / comp_1_phi) * K1
        )
    )

    matrix[
        eq_idx[comp_1.identifier].pressure,
        eq_idx[comp_1.identifier].pressure,
    ] += s_pj_pj

    # PRESSURE EQUATION: above/below main diagonal elements construction:
    # (j+num_fluid_components,\
    # k + num_fluid_components:2*num_fluid_components) [Pres_k]
    matrix[
        eq_idx[comp_1.identifier].pressure,
        eq_idx[comp_2.identifier].pressure,
    ] = - s_pj_pj

    # (j+num_fluid_components,j+2*num_fluid_components)
    # [Temp_j] I
    # s_pj_tc = phi/A * (P_o * h_o + P_c * h_c)
    #         = coef_frun_area * coef_htc
    s_pj_tj = coef_grun_area * coef_htc

    matrix[
        eq_idx[comp_1.identifier].pressure,
        eq_idx[comp_1.identifier].temperature,
    ] += s_pj_tj

    # (j+num_fluid_components, 
    # k + 2*num_fluid_components:dict_N_equation
    # ["FluidComponent"]) [Temp_j]
    matrix[
        eq_idx[comp_1.identifier].pressure,
        eq_idx[comp_2.identifier].temperature,
    ] = - s_pj_tj

    # TEMPERATURE EQUATION: elements below main diagonal \
    # construction:
    # (j+2*num_fluid_components,j+num_fluid_components) [Pres_j]

    # coef_rho_cv_area = 1/(rho * c_v * A)
    coef_rho_cv_area = 1. / (comp_1_rho * comp_1_cv * comp_1_A)
    # s_tj_pj = 1/(rho * c_v * A) * [K3 - vK2 - (w - v^2/2 - phi*c_v*T)K1]
    #         = coef_rho_cv_area * [K3 - vK2 - (w - v^2/2 - phi*c_v*T)K1]
    s_tj_pj = (
        coef_rho_cv_area
        * (
            K3 - comp_1_v * K2
            - (comp_1_enthalpy - comp_1_v ** 2. / 2.
                - comp_1_phi * comp_1_cv
                * comp_1.coolant.dict_Gauss_pt["temperature"][elem_idx]
            )
            * K1
        )
    )

    matrix[
        eq_idx[comp_1.identifier].temperature,
        eq_idx[comp_1.identifier].pressure,
    ] += s_tj_pj

    # (j+2*num_fluid_components,\
    # k + num_fluid_components:2*num_fluid_components) [Pres_k]
    matrix[
        eq_idx[comp_1.identifier].temperature,
        eq_idx[comp_2.identifier].pressure,
    ] = - s_tj_pj

    # TEMPERATURE EQUATION: main diagonal element construction:
    # (j+2*num_fluid_components,j+2*num_fluid_components) 
    # [Temp_j] I

    # s_tj_tj = 1/(rho * c_v * A) * (P_o * h_o + P_c * h_c)
    #         = coef_rho_cv_area * coef_htc
    s_tj_tj = coef_rho_cv_area * coef_htc

    matrix[
        eq_idx[comp_1.identifier].temperature,
        eq_idx[comp_1.identifier].temperature,
    ] += s_tj_tj

    # TEMPERATURE EQUATION: above/below main diagonal elements 
    # construction:
    # (j+2*num_fluid_components,k + 2*num_fluid_components) 
    # [Temp_k]
    matrix[
        eq_idx[comp_1.identifier].temperature,
        eq_idx[comp_2.identifier].temperature,
    ] = - s_tj_tj

    return matrix

def build_smat_fluid_solid_interface(
    matrix:np.ndarray,
    conductor:Conductor,
    elem_idx:int,
    )->np.ndarray:

    """Function that builds the S matrix (SMAT) therms due to fluid-solid component interfaces at the Gauss point (SOURCE JACOBIAN).

    Args:
        matrix (np.ndarray): S matrix after call to function buld_smat_fluid_interface.
        conductor (Conductor): object with all the information of the conductor.
        elem_idx (int): index of the i-th element of the spatial discretization.
        eq_idx (dict): collection of NamedTuple with fluid equation index (velocity, pressure and temperaure equations) and of integer for solid equation index.

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

    # Alias
    # Collection of NamedTyple with fluid equation index (velocity, pressure 
    # and temperaure equations).
    eq_idx = conductor.equation_index

    for interface in conductor.interface.fluid_solid:

        comp_1_A = interface.comp_1.channel.inputs["CROSSECTION"]
        # pressure equation: above main diagonal elements construction.
        # (j+num_fluid_components,j+2*num_fluid_components) [Temp_j] II + III

        # coef_grun_area = phi / A
        coef_grun_area = (
            interface.comp_1.coolant.dict_Gauss_pt["Gruneisen"][elem_idx]
            / comp_1_A
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
        ] += s_pj_tj
        
        # (j+num_fluid_components,l + dict_N_equation["FluidComponent"]) [Temp_l]
        matrix[
            eq_idx[interface.comp_1.identifier].pressure,
            eq_idx[interface.comp_2.identifier],
        ] = - s_pj_tj

        # temperature equation: main diagonal element construction
        # (j+2*num_fluid_components,j+2*num_fluid_components) [Temp_j] II + III 
        
        # coef_rho_cv_area = 1/(rho * c_v * A)
        coef_rho_cv_area = 1. / (
            interface.comp_1.coolant.dict_Gauss_pt["total_density"][elem_idx]
            * interface.comp_1.coolant.dict_Gauss_pt[
                "total_isochoric_specific_heat"
            ][elem_idx]
            * comp_1_A
        )

        # s_tj_tj = 1/(rho * c_v * A) * P * h
        #         = coef_rho_cv_area * coef_htc
        s_tj_tj = coef_rho_cv_area * coef_htc
        
        matrix[
            eq_idx[interface.comp_1.identifier].temperature,
            eq_idx[interface.comp_1.identifier].temperature,
        ] += s_tj_tj
        
        # temperature equation: above main diagonal elements construction
        # (j+2*num_fluid_components,l + dict_N_equation["FluidComponent"]) [Temp_l]
        matrix[
            eq_idx[interface.comp_1.identifier].temperature,
            eq_idx[interface.comp_2.identifier],
        ] = - s_tj_tj

        # SOLID COMPONENTS CONDUCTION EQUATION: main diagonal element
        # construction.
        # (l + dict_N_equation["FluidComponent"],l + dict_N_equation["FluidComponent"]) [Temp_l] I
        matrix[
            eq_idx[interface.comp_2.identifier],
            eq_idx[interface.comp_2.identifier],
        ] += coef_htc
        
        # SOLID COMPONENTS CONDUCTION EQUATION: below main diagonal elements
        # construction.
        # (l + dict_N_equation["FluidComponent"],l + 2*num_fluid_components) [Temp_j]
        matrix[
            eq_idx[interface.comp_2.identifier],
            eq_idx[interface.comp_1.identifier].temperature,
        ] = -coef_htc

    return matrix

def build_mmat_solid(
    matrix:np.ndarray,
    s_comp:SolidComponent,
    elem_idx:int,
    eq_idx:int,
    )->np.ndarray:

    """Function that updates the M matrix (MMAT) at the Gauss point, for the SolidComponent equation.

    Args:
        matrix (np.ndarray): M matrix with the element from the fluid equations.
        s_comp (SolidComponent): solid component object from which get all info to build the coefficients.
        elem_idx (int): index of the i-th element of the spatial discretization.
        eq_idx (int): solid component equation index.

    Returns:
        np.ndarray: matrix with updated elements.
    """

    # FORM THE M MATRIX AT THE GAUSS POINT (MASS AND CAPACITY)
    # SolidComponent (homogenized) equation.
    # A * rho *cp / cos(theta)
    matrix[eq_idx, eq_idx] = (
        s_comp.inputs["CROSSECTION"]
        * s_comp.dict_Gauss_pt["total_density"][elem_idx]
        * s_comp.dict_Gauss_pt["total_isobaric_specific_heat"][elem_idx]
        / s_comp.inputs["COSTETA"]
    )

    return matrix

def build_kmat_solid(
    matrix:np.ndarray,
    s_comp:SolidComponent,
    elem_idx:int,
    eq_idx:int,
    )->np.ndarray:

    """Function that updates the K matrix (KMAT) at the Gauss point, for the SolidComponent equation.

    Args:
        matrix (np.ndarray): K matrix after call to build_kmat_fluid.
        s_comp (SolidComponent): solid component object from which get all info to build the coefficients.
        elem_idx (int): index of the i-th element of the spatial discretization.
        eq_idx (int): solid component equation index.

    Returns:
        np.ndarray: matrix with updated elements.
    """

    # FORM THE K MATRIX AT THE GAUSS POINT (INCLUDING UPWIND)
    # A_{s_comp}*k_{s_comp,homo}; homo = homogenized
    matrix[eq_idx,eq_idx] = (
        s_comp.inputs["CROSSECTION"]
        * s_comp.dict_Gauss_pt["total_thermal_conductivity"][elem_idx]
        / s_comp.inputs["COSTETA"]
    )

    return matrix

def build_smat_solid_interface(
    matrix:np.ndarray,
    conductor:Conductor,
    elem_idx:int,
    )->np.ndarray:

    """Function that builds the S matrix (SMAT) therms due to solid component interfaces at the Gauss point (SOURCE JACOBIAN).

    Args:
        matrix (np.ndarray): S matrix after call to function build_smat_fluid_solid_interface.
        conductor (Conductor): object with all the information of the conductor.
        elem_idx (int): index of the i-th element of the spatial discretization.

    Returns:
        np.ndarray: matrix with updated elements.
    """
    
    # NOMENCLATURE
    # P: contact perimeter
    # h_conv: convective heat transfer coefficient.

    # Alias
    # Collection of integer solid equation index.
    eq_idx = conductor.equation_index

    for interface in conductor.interface.solid_solid:

        # coef_htc = P * h_conv W / m / K
        coef_htc = (
            conductor.dict_interf_peri["sol_sol"][interface.interf_name]
            * conductor.dict_Gauss_pt["HTC"]["sol_sol"]["cond"][
                interface.interf_name
            ][elem_idx]
        )

        # SOLID COMPONENTS CONDUCTION EQUATION: main diagonal element 
        # construction:
        # (l + dict_N_equation["FluidComponent"],l 
        # + dict_N_equation["FluidComponent"]) [Temp_l] II + III
        matrix[
            eq_idx[interface.comp_1.identifier],
            eq_idx[interface.comp_1.identifier],
         ] += coef_htc
        
        # SOLID COMPONENTS CONDUCTION EQUATION: above/below main diagonal 
        # elements construction:
        # (l + dict_N_equation["FluidComponent"],m 
        # + dict_N_equation["FluidComponent"]) [Temp_m]
        matrix[
            eq_idx[interface.comp_1.identifier],
            eq_idx[interface.comp_2.identifier],
        ] = - coef_htc
    
    return matrix

def build_smat_env_solid_interface(
    matrix:np.ndarray,
    conductor:Conductor,
    interface:NamedTuple,
    elem_idx:int,
    )->np.ndarray:

    """Function that builds the S matrix (SMAT) therms due to environment and solid component interfaces at the Gauss point (SOURCE JACOBIAN).

    Args:
        matrix (np.ndarray): S matrix after call to function build_smat_solid_interface.
        conductor (Conductor): object with all the information of the conductor.
        interface (NamedTuple): collection of interface information like interface name and components that constitute the interface.
        elem_idx (int): index of the i-th element of the spatial discretization.

    Returns:
        np.ndarray: matrix with updated elements.
    """

    # Alias.
    h_conv = conductor.dict_Gauss_pt["HTC"]["env_sol"][
                interface.interf_name
            ]["conv"]
    # Collection of integer solid equation index.
    eq_idx = conductor.equation_index

    # Convective heating with the external environment (implicit treatment).
    if (
        conductor.dict_df_coupling["HTC_choice"].at[
            interface.comp_1.KIND,
            interface.comp_2.identifier,
        ]
        == 2
        and conductor.inputs["Is_rectangular"]
    ):
        # Rectangular duct.
        coef_htc = (
            + 2. * conductor.inputs["Height"] * h_conv["side"][elem_idx]
            + conductor.inputs["Width"]
            * (
                h_conv["bottom"][elem_idx] + h_conv["top"][elem_idx]
            )
        )
    else:
        coef_htc = (
            h_conv[elem_idx]
            * conductor.dict_interf_peri["env_sol"][
                interface.interf_name
            ]
        )
    # Update matrix coefficients.
    matrix[
            eq_idx[interface.comp_2.identifier],
            eq_idx[interface.comp_2.identifier],
        ] += coef_htc

    return matrix

def build_svec(
    array:np.ndarray,
    s_comp: SolidComponent,
    elem_idx:int,
    eq_idx:int,
    **kwargs,
    )->np.ndarray:
    """Function that builds the source vector (SVEC) elements at the Gauss point due to heat generation in strand and or jacket component objects and to thermal contact beween jacket compoments belonging to different conductors (qsource). For strand component objects the latter contribution is always zero.
    N.B. This function is a merge of the if statement if isinstance(scomp,StrandComponent) is true do not account for qsource else, account for qsource. Since, as mentioned qsourse = 0 for StrandComponent, the check chan be avoided. This should improve readability and maintainability of the function.

    Args:
        array (np.ndarray): initialized SVEC array.
        s_comp (SolidComponent): solid component object from which get all info to build the coefficients.
        elem_idx (int): index of the i-th element of the spatial discretization.
        eq_idx (int): solid equation index.

    Kwargs:
        num_step (int): present time step counter value.
        qsource (np.ndarray): matrix with heat due to thermal contact between jacket components of different conductors.
        comp_idx (int): component index, used to correctly assign the heat source term due to thermal contact between solid components of different conductors.

    Returns:
        np.ndarray: array with updated elements.
    """
    
    # Alias.
    qsource = kwargs["qsource"]
    comp_idx = kwargs["comp_idx"]
    Q1 = s_comp.dict_Gauss_pt["Q1"]
    Q2 = s_comp.dict_Gauss_pt["Q2"]

    # N.B. qsource has non zero values only in nodes and columns that represent 
    # the contact between jacket components of different conductors.

    # This is independent from the solution method thanks to the escamotage of 
    # the dummy steady state corresponding to the initialization.
    if kwargs["num_step"] == 1:
        # Present time step.
        array.present[eq_idx,0] = (
            Q1[elem_idx,0] - qsource[elem_idx,comp_idx]
        )
        array.present[eq_idx,1] = (
            Q2[elem_idx,0] - qsource[elem_idx + 1, comp_idx]
        )
        # Previous time step.
        array.previous[eq_idx,0] = (
            Q1[elem_idx,1] - qsource[elem_idx,comp_idx]
        )
        array.previous[eq_idx,1] = (
            Q2[elem_idx,1] - qsource[elem_idx + 1, comp_idx]
        )
    else:
        # Compute only at the current time step.
        array[eq_idx,0] = (
            Q1[elem_idx,0] - qsource[elem_idx, comp_idx]
        )
        array[eq_idx,1] = (
            Q2[elem_idx,0] - qsource[elem_idx + 1, comp_idx]
        )
    
    return array

def build_svec_env_jacket_interface(
    array:np.ndarray,
    conductor: Conductor,
    interface:NamedTuple,
    elem_idx:int,
    )->np.ndarray:
    """Function that builds the source vector (SVEC) terms at the Gauss point due to heat transfer by convection and/or radiation between environment and jacket component objects.

    Args:
        array (np.ndarray): SVEC array after call to function build_svec.
        interface (NamedTuple): collection of interface information like interface name and components that constitute the interface.
        elem_idx (int): index of the i-th element of the spatial discretization.

    Returns:
        np.ndarray: array with updated elements.
    """

    # Alias.
    h_conv = conductor.dict_Gauss_pt["HTC"]["env_sol"][
                interface.interf_name
            ]["conv"]
    height = conductor.inputs["Height"]
    width = conductor.inputs["Width"]
    env = interface.comp_1
    s_comp = interface.comp_2
    # Collection of integer solid equation index.
    eq_idx = conductor.equation_index
    
    # Add the contribution of the external heating by convection to the 
    # known term vector.
    if (
        conductor.dict_df_coupling["HTC_choice"].at[
            env.KIND, s_comp.identifier
        ]
        == 2
        and conductor.inputs["Is_rectangular"]
    ):
        # Rectangular duct.
        coef = 2. * height * h_conv["side"][elem_idx]
        + width* (h_conv["bottom"][elem_idx] + h_conv["top"][elem_idx])
    else:
        coef = (
            conductor.dict_interf_peri["env_sol"][
                interface.interf_name
            ]
            * h_conv[elem_idx]
        )
    
        # Linear heat flux from environment W/m
        env_heat = coef * env.inputs["Temperature"]

        if conductor.cond_num_step == 1:
            # Present time step.
            array.present[eq_idx[s_comp.identifier],0] += env_heat
            array.present[eq_idx[s_comp.identifier],1] += env_heat
            # Previous time step.
            array.previous[eq_idx[s_comp.identifier],0] += env_heat
            array.previous[eq_idx[s_comp.identifier],1] += env_heat
        else:
            # Present time step.
            array[eq_idx[s_comp.identifier],0] += env_heat
            array[eq_idx[s_comp.identifier],1] += env_heat

    return array

def build_elmmat(
    matrix:np.ndarray,
    mmat:np.ndarray,
    conductor:Conductor,
    elem_idx:int,
    alpha:float=0,
    )->np.ndarray:
    """Function that builds the mass and capacity matrix (ELMMAT) at the Gauss point exploiting slicing.

    Args:
        matrix (np.ndarray): Initialized ELM matrix
        mmat (np.ndarray): mass and capacity matrix MMAT after call to function build_mmat_solid.
        conductor (Conductor): object with all the information of the conductor.
        elem_idx (int): index of the i-th element of the spatial discretization.
        alpha (float, optional): Lumped/consistent mass parameter. Defaults to 0.

    Returns:
        np.ndarray: matrix with updated elements.
    """
    
    # Alias
    # Length of the present element of the spatial discretization.
    dz = conductor.grid_features["delta_z"][elem_idx]
    # Number of degrees of freedom, used to slice ELMMAT matrix.
    ndf = conductor.dict_N_equation["NODOFS"]
    # Twice the number of degrees of freedom, used to slice ELMMAT matrix.
    ndf2 = conductor.dict_N_equation["NODOFS2"]

    # Build diagonal block of the matrix.
    diag_block = dz * (1. / 3. + alpha) * mmat
    # Build off diagonal block of the matrix.
    off_diag_block = dz * (1. / 6. - alpha) * mmat
    
    # COMPUTE THE MASS AND CAPACITY MATRIX
    # array smart
    matrix[:ndf,:ndf] = diag_block
    matrix[:ndf,ndf:ndf2] = off_diag_block
    matrix[ndf:ndf2,:ndf] = off_diag_block
    matrix[ndf:ndf2,ndf:ndf2] = diag_block

    return matrix

def build_elamat(
    matrix:np.ndarray,
    amat:np.ndarray,
    conductor:Conductor,
    )->np.ndarray:
    """Function that builds the convection matrix (ELAMAT) at the Gauss point exploiting slicing.

    Args:
        matrix (np.ndarray): Initialized ELA matrix.
        amat (np.ndarray): matrix AMAT after call to function build_amat.
        conductor (Conductor): object with all the information of the conductor.

    Returns:
        np.ndarray: matrix with updated elements.
    """

    # Alias
    # Number of degrees of freedom, used to slice ELAMAT matrix.
    ndf = conductor.dict_N_equation["NODOFS"]
    # Twice the number of degrees of freedom, used to slice ELAMAT matrix.
    ndf2 = conductor.dict_N_equation["NODOFS2"]

    block = amat / 2.

    matrix[:ndf,:ndf] = -block
    matrix[:ndf,ndf:ndf2] = block
    matrix[ndf:ndf2,:ndf] = -block
    matrix[ndf:ndf2,ndf:ndf2] = block

    return matrix

def build_elkmat(
    matrix:np.ndarray,
    kmat:np.ndarray,
    conductor:Conductor,
    elem_idx:int,
    )->np.ndarray:
    """
    Function that builds the diffusion matrix (ELKMAT) at the Gauss point exploiting slicing.

    Args:
        matrix (np.ndarray): Initialized ELK matrix.
        kmat (np.ndarray): mass and capacity matrix KMAT after call to function build_kmat_solid.
        conductor (Conductor): object with all the information of the conductor.
        elem_idx (int): index of the i-th element of the spatial discretization.

    Returns:
        np.ndarray: matrix with updated elements.
    """

    # Alias
    # Length of the present element of the spatial discretization.
    dz = conductor.grid_features["delta_z"][elem_idx]
    # Number of degrees of freedom, used to slice ELKMAT matrix.
    ndf = conductor.dict_N_equation["NODOFS"]
    # Twice the number of degrees of freedom, used to slice ELKMAT matrix.
    ndf2 = conductor.dict_N_equation["NODOFS2"]

    block = kmat / dz
    
    # COMPUTE THE DIFFUSION MATRIX
    # array smart
    matrix[:ndf,:ndf] = block
    matrix[:ndf,ndf:ndf2] = - block
    matrix[ndf:ndf2,:ndf] = - block
    matrix[ndf:ndf2,ndf:ndf2] = block

    return matrix

def build_elsmat(
    matrix:np.ndarray,
    smat:np.ndarray,
    conductor:Conductor,
    elem_idx:int,
    )->np.ndarray:
    """
    Function that builds the source matrix (ELSMAT) at the Gauss point exploiting slicing.

    Args:
        matrix (np.ndarray): Initialized ELS matrix.
        smat (np.ndarray): source matrix SMAT after call to function build_smat_env_solid_interface.
        conductor (Conductor): object with all the information of the conductor.
        elem_idx (int): index of the i-th element of the spatial discretization.

    Returns:
        np.ndarray: matrix with updated elements.
    """

    # Alias
    # Length of the present element of the spatial discretization.
    dz = conductor.grid_features["delta_z"][elem_idx]
    # Number of degrees of freedom, used to slice ELSMAT matrix.
    ndf = conductor.dict_N_equation["NODOFS"]
    # Twice the number of degrees of freedom, used to slice ELSMAT matrix.
    ndf2 = conductor.dict_N_equation["NODOFS2"]

    diag_block = smat * dz / 3.
    off_diag_block = diag_block / 2.

    # COMPUTE THE SOURCE MATRIX
    # array smart
    matrix[:ndf,:ndf] = diag_block
    matrix[:ndf,ndf:ndf2] = off_diag_block
    matrix[ndf:ndf2,:ndf] = off_diag_block
    matrix[ndf:ndf2,ndf:ndf2] = diag_block

    return matrix

def build_elslod(array:np.ndarray,
    svec:np.ndarray,
    conductor:Conductor,
    elem_idx:int,
    )->np.ndarray:
    """
    Function that builds the source matrix (ELSLOD) at the Gauss point exploiting slicing.

    Args:
        array (np.ndarray): Initialized ELSLOD array.
        svec (np.ndarray): source array SVEC after call to function build_svec_env_jacket_interface.
        conductor (Conductor): object with all the information of the conductor.
        elem_idx (int): index of the i-th element of the spatial discretization.

    Returns:
        np.ndarray: array with updated elements.
    """

    # Alias
    # Length of the present element of the spatial discretization.
    dz = conductor.grid_features["delta_z"][elem_idx]
    # Number of degrees of freedom, used to slice ELSMAT matrix.
    ndf = conductor.dict_N_equation["NODOFS"]
    # Twice the number of degrees of freedom, used to slice ELSMAT matrix.
    ndf2 = conductor.dict_N_equation["NODOFS2"]

    dz_6 = dz / 6.

    # COMPUTE THE SOURCE VECTOR (ANALYTIC INTEGRATION)
    # array smart
    # This is independent from the solution method thanks to the escamotage of
    # the dummy steady state corresponding to the initialization.
    if conductor.cond_num_step == 1:
        # To correctly apply the theta method (to be rivisited in the whole 
        # code!).
        # Current time step
        array.present[:ndf] = 2. * svec.present[:,0] + svec.present[:,1]
        array.present[ndf:ndf2] = svec.present[:,0] + 2. * svec.present[:,1]
        # Not use simply array.present: raises AttributeError
        array.present[:] *= dz_6
        # Previous time step
        array.previous[:ndf] = 2. * svec.previous[:,0] + svec.previous[:,1]
        array.previous[ndf:ndf2] = svec.previous[:,0] + 2. * svec.previous[:,1]
        # Not use simply array.previous: raises AttributeError
        array.previous[:] *= dz_6
    else:
        # Compute only at the current time step
        array[:ndf] = 2. * svec[:,0] + svec[:,1]
        array[ndf:ndf2] = svec[:,0] + 2. * svec[:,1]
        array *= dz_6
    
    return array

def assemble_matrix(
    big_matrices:tuple,
    small_matrices:tuple,
    conductor:Conductor,
    jump_idx:int,
    )->tuple:
    """Function that assembles matrices MASMAT, FLXMAT, DIFMAT and SORMAT starting from the values of matrices ELMMAT, ELAMAT, ELKMAT and ELSMAT respectively. The matrices match is as follows:
        * ELMMAT builds MATMAT
        * ELAMAT builds FLXMAT
        * ELKMAT builds DIFMAT
        * ELSMAT builds SORMAT

    Args:
        big_matrices (tuple): collection of matrices MASMAT, FLXMAT, DIFMAT and SORMAT.
        small_matrices (tuple): collection of matrices ELMMAT, ELAMAT, ELKMAT and ELSMAT.
        conductor (Conductor): object with all the information of the conductor.
        jump_idx (int): index to jump over NODOFS * elem_idx position in the big matrices.

    Returns:
        tuple: collection of matrices MASMAT, FLXMAT, DIFMAT and SORMAT with updated elements.
    """
    
    # Alias
    half = conductor.dict_band["Half"]
    full = conductor.dict_band["Full"]

    # ASSEMBLE THE MATRICES AND THE LOAD VECTOR
    # array smart
    for hbw_idx in range(half):
        # hbw_idx: half band widht index.
        # Lower bound for row slicing.
        row_lb = half - hbw_idx - 1
        # Upper bound for row slicing.
        row_ub = full - hbw_idx
        # Column index.
        col = jump_idx + hbw_idx
        # Loop to update matrices masmat, flxmat, difmat and sormat stored in 
        # tuple big_matrices exploiting matrices elmmat, elamat, elkmat and 
        # elsmat respectively, stored in tuple small_matrices.
        for big_matrix, small_matrix in zip(big_matrices,small_matrices):
            # Sum elements from row_lb to row_ub of column col in big_matrix 
            # with all the elements in the row hbw_idx of small_matrix.
            big_matrix[row_lb:row_ub,col] +=small_matrix[hbw_idx,:]
    
    return big_matrices

def assemble_syslod(
    array:np.ndarray,
    conductor:Conductor,
    jump_idx:int,
)->np.ndarray:
    """Function that assembles the source term vector syslod exploiting the information stored in array ELSLOD.

    Args:
        array (np.ndarray): ELSLOD array after call to funciton build_elslod.
        conductor (Conductor): object with all the information of the conductor.
        jump_idx (int): index to jump over NODOFS * elem_idx position in syslod array, used to slice syslod.

    Returns:
        np.ndarray: array with updated elements (the sliced updated part f syslod).
    """

    # Alias
    method = conductor.inputs["METHOD"]
    num_step = conductor.cond_num_step
    half = conductor.dict_band["Half"]
    syslod = conductor.dict_Step["SYSLOD"][jump_idx:jump_idx + half,:]
    
    if method == "BE" or method == "CN":
        # Backward Euler or Crank-Nicolson
        if num_step == 1:
            # Construct key SYSLOD of dictionary dict_Step
            # Current time step
            syslod[:,0] += array.present
            # Previous time step
            syslod[:,1] += array.previous
        else:
            # Update only the first column, that correspond to the current time
            # step
            syslod[:,0] += array
    elif method == "AM4":
        # Adams-Moulton order 4
        # The implementation of higher order numerical schemes for time 
        # integration should be completely reviewed!
        if num_step == 1:
            # Construct key SYSLOD of dictionary dict_Step
            # Current time step
            syslod[:,0] += array.present
            # This loop should be checked carefully!
            for cc in range(syslod.shape[1]):
                # Dummy initial steady state
                syslod[:,cc] += array.previous
        else:
            # Shift the colums by one towards right and compute the new first 
            # column at the current time step
            syslod[:,1:] = syslod[:,:3].copy()
            syslod[:,0] += array
        # end if conductor.cond_num_step
    # end conductor.inputs["METHOD"]

    return syslod

def eval_system_matrix(
    matrix:np.ndarray,
    aux_matrices:tuple,
    conductor:Conductor,
    )->np.ndarray:
    """Function that evaluates the system matrix using the values of the matrix MASMAT, FLXMAT, DIFMAT and SORMAT, according to the selected method for time integration.

    Args:
        matrix (np.ndarray): initialized SYSMAT matrix
        aux_matrices (tuple): collection of matrix MASMAT, FLXMAT, DIFMAT and SORMAT after call to function assemble_matrix.
        conductor (Conductor): object with all the information of the conductor.

    Returns:
        np.ndarray: matrix with updated elements.
    """
    
    # Alias
    method = conductor.inputs["METHOD"]
    # Unpack auxiliary matrices (MASMAT,FLXMAT,DIFMAT,SORMAT)
    masmat,flxmat,difmat,sormat = aux_matrices
    # ** COMPUTE SYSTEM MATRIX **
    if method == "BE" or method == "CN":
        # Backward Euler or Crank-Nicolson
        matrix = (
            masmat / conductor.time_step
            + conductor.theta_method * (flxmat + difmat + sormat)
        )
        
    elif method == "AM4":
        # Adams-Moulton order 4
        # Alias
        am4_aa = conductor.dict_Step["AM4_AA"] # shallow copy
        if conductor.cond_num_step == 1:
            # This is due to the dummy initial steady state
            for cc in range(am4_aa.shape[2]):
                am4_aa[cc,:,:] = flxmat + difmat + sormat
        else:
            # Shift the matrices by one towards right and compute the new first
            # matrix at the current time step.
            am4_aa[1:4,:,:] = am4_aa[0:3,:,:]
            am4_aa[0,:,:] = flxmat + difmat + sormat
        # compute SYSMAT
        matrix = masmat / conductor.time_step + 9. / 24. * am4_aa[0,:,:]
    
    return matrix

def build_known_therm_vector(
    array:np.ndarray,
    aux_matrices:tuple,
    conductor:Conductor
)->np.ndarray:
    """Function that builds the known therm vector for the thermal hydraulic problem according to the selected method for time integration.

    Args:
        array (np.ndarray): initialized array Known.
        aaux_matrices (tuple): collection of matrix MASMAT, FLXMAT, DIFMAT and SORMAT after call to function assemble_matrix.
        conductor (Conductor): object with all the information of the conductor.

    Returns:
        np.ndarray: array with updated elements.
    """

    # Alias
    total = conductor.dict_N_equation["Total"]
    half = conductor.dict_band["Half"]
    half_1 = half - 1
    method = conductor.inputs["METHOD"]
    syslod = conductor.dict_Step["SYSLOD"] # shallow copy
    sysvar = conductor.dict_Step["SYSVAR"] # shallow copy

    if method == "AM4":
        # Alias
        am4_aa = conductor.dict_Step["AM4_AA"] # shallow copy
        am4_coef = np.array((9.,19.,5.,- 1.)) / 24.
    
    # Unpack auxiliary matrices (MASMAT,FLXMAT,DIFMAT,SORMAT)
    masmat,flxmat,difmat,sormat = aux_matrices
    
    # ADD THE LOAD CONTRIBUTION FROM PREVIOUS STEP
    # c_mat_idx: column index of the auxiliary matrices (MASMAT,FLXMAT,DIFMAT,
    # SORMAT); used also as row index of the known term vector.
    for c_mat_idx in range(total):
        if c_mat_idx <= half_1:
            # remember that arange stops before the stop value:
            # last value = stop - step
            # r_arr_idx: row index of the sysvar array
            r_arr_idx = np.arange(
                start=0,
                stop=half + c_mat_idx,
                step=1,
                dtype=int,
            )
        elif c_mat_idx >= total - half_1:
            r_arr_idx = np.arange(
                start=c_mat_idx - half_1,
                stop=total,
                step=1,
                dtype=int,
            )
        else:
            r_arr_idx = np.arange(
                start=c_mat_idx - half_1,
                stop=c_mat_idx + half,
                step=1,
                dtype=int,
            )
        # r_mat_idx: row index of the auxiliary matrices (MASMAT,FLXMAT,DIFMAT,
        # SORMAT)
        r_mat_idx = r_arr_idx - c_mat_idx + half_1
        if method == "BE" or method == "CN":
            # Backward Euler or Crank-Nicolson
            # Matrix vector product contribution
            array[c_mat_idx] = np.sum(
                (
                    masmat[r_mat_idx,c_mat_idx] / conductor.time_step
                    - (1.0 - conductor.theta_method)
                    * (
                        flxmat[r_mat_idx, c_mat_idx]
                        + difmat[r_mat_idx, c_mat_idx]
                        + sormat[r_mat_idx, c_mat_idx]
                    )
                )
                * sysvar[r_arr_idx,0]
            )
        elif method == "AM4":
            # Adams-Moulton order 4
            # Matrices vectors product contribution
            # np.sum(am4_coef[2:] * am4_aa[2:,r_mat_idx,c_mat_idx].T * sysvar[r_arr_idx,1:3],1) should be equivalent to 5. / 24.* am4_aa[2,r_mat_idx,c_mat_idx] * sysvar[r_arr_idx, 1] - 1. / 24. * am4_aa[3,r_mat_idx,c_mat_idx] * sysvar[r_arr_idx, 2]
            array[c_mat_idx] = np.sum(
                (
                    masmat[r_mat_idx, c_mat_idx] / conductor.time_step
                    - am4_coef[1] * am4_aa[1,r_mat_idx,c_mat_idx]
                ) * sysvar[r_arr_idx,0] # array of shape (r_arr_idx.shape[0],)
                + np.sum(
                    am4_coef[2:] * am4_aa[2:,r_mat_idx,c_mat_idx].T
                    * sysvar[r_arr_idx,1:3],1
                ) # array of shape (r_arr_idx.shape[0],)
            ) # array of shape (1,)

    if method == "BE" or method == "CN":
        # Backward Euler or Crank-Nicolson
        # External sources (SYSLOD) contribution
        array += (
            + conductor.theta_method * syslod[:,0]
            + (1.0 - conductor.theta_method) * syslod[:,1]
        )
    elif method == "AM4":
        # Adams-Moulton order 4

        # Chance coefficient sign to exploit sum (array smart).
        am4_coef[2:] = - am4_coef[2:]
        # External sources (SYSLOD) contribution
        array += np.sum(am4_coef * syslod,1)

    return array