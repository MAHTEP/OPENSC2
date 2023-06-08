import numpy as np

from collections import namedtuple
from typing import Union, NamedTuple

from utility_functions.auxiliary_functions import get_from_xlsx

from fluid_component import FluidComponent
from conductor import Conductor

def impose_pressure_drop(
    ndarrays:tuple,
    f_comp:FluidComponent,
    conductor:Conductor,
    path:str,
    )->tuple:
    
    """Function that imposes a pressure drop as boundary conditions for the thermal hydraulic problem. Imposed quantities are:
        * inlet pressure
        * inlet temperature (or outlet temperature if back flow)
        * outlet pressure
    The function suitably accounts for forward (left to right)/backward (right to left) flow as well as for back flow (part of the fluid moves towards the inlet due to a pressure build up).

    Args:
        ndarrays (tuple): collection of known term vector (Known) and system matrix (SYSMAT) np.ndarray on wich impose the boundary conditions.
        f_comp (FluidComponent): fluid component object from which get all info to apply the boundary conditions.
        conductor (Conductor): object with all the information of the conductor.
        path (str): path of the auxiliary input file with the values for the boundary conditions (if INTIAL = -1).

    Returns:
        tuple: collection of known term vector (Known) and system matrix (SYSMAT) np.ndarray with imposed boundary conditions.
    """
    
    known,sysmat = ndarrays
    # Get boundary condition values.
    if f_comp.coolant.operations["INTIAL"] == 1:
        # inlet pressure.
        p_inl = f_comp.coolant.operations["PREINL"]
        # inlet temperature.
        T_inl = f_comp.coolant.operations["TEMINL"]
        # outlet pressure.
        p_out = f_comp.coolant.operations["PREOUT"]
        # outlet temperature: to be assigned if outlet velocity is negative.
        T_out = f_comp.coolant.operations["TEMOUT"]
    else:
        # get from file with interpolation in time.
        [flow_par, flagSpecfield] = get_from_xlsx(
            conductor,
            path,
            f_comp,
            "INTIAL",
            f_comp.coolant.operations["INTIAL"]
        )
        print(
            f"""flagSpecfield == {flagSpecfield}: still to be decided if 
        it is useful and if yes still to be defined\n"""
        )
        # inlet pressure.
        p_inl = flow_par[2]
        # inlet temperature.
        T_inl = flow_par[0]
        # outlet pressure.
        p_out = flow_par[3]
        # outlet temperature: to be assigned if outlet velocity is negative.
        T_out = flow_par[1]
    
    # ALIAS
    _, inl_p_idx, inl_t_idx = f_comp.inl_idx # tuple unpack
    _, out_p_idx, out_t_idx = f_comp.out_idx # tuple unpack
    main_d_idx = conductor.dict_band["Main_diag"]
    flow_dir = f_comp.channel.flow_dir[0]
    velocity = f_comp.coolant.dict_node_pt["velocity"]
    
    # Assign BC
    if flow_dir == "forward":
        # p_inl
        sysmat[:,inl_p_idx.forward] = 0.0
        # main diagonal.
        sysmat[main_d_idx,inl_p_idx.forward] = 1.0
        known[inl_p_idx.forward] = p_inl
        # p_out
        sysmat[:,out_p_idx.forward] = 0.0
        # main diagonal.
        sysmat[main_d_idx,out_p_idx.forward] = 1.0
        known[out_p_idx.forward] = p_out
        # T_inl
        if velocity[0] > 0:
            sysmat[:,inl_t_idx.forward] = 0.0
            # main diagonal.
            sysmat[main_d_idx,inl_t_idx.forward] = 1.0
            known[inl_t_idx.forward] = T_inl
        # T_out
        if velocity[-1] < 0:
            sysmat[:,out_t_idx.forward] = 0.0
            # main diagonal.
            sysmat[main_d_idx,out_t_idx.forward] = 1.0
            known[out_t_idx.forward] = T_out
    elif flow_dir == "backward":
        # p_inl
        sysmat[:,inl_p_idx.backward] = 0.0
        # main diagonal.
        sysmat[main_d_idx,inl_p_idx.backward] = 1.0
        known[inl_p_idx.backward] = p_inl
        # p_out
        sysmat[:,out_p_idx.backward] = 0.0
        # main diagonal.
        sysmat[main_d_idx,out_p_idx.backward] = 1.0
        known[out_p_idx.backward] = p_out
        # T_inl
        if velocity[-1] < 0:
            sysmat[:,inl_t_idx.backward] = 0.0
            # main diagonal.
            sysmat[main_d_idx,inl_t_idx.backward] = 1.0
            known[inl_t_idx.backward] = T_inl
        # T_out
        if velocity[0] > 0:
            sysmat[:,out_t_idx.backward] = 0.0
            # main diagonal.
            sysmat[main_d_idx,out_t_idx.backward] = 1.0
            known[out_t_idx.backward] = T_out

    return known,sysmat

def impose_inl_p_out_v(
    ndarrays:tuple,
    f_comp:FluidComponent,
    conductor:Conductor,
    path:str,
    )->tuple:
    
    """Function that imposes a inlet pressure and outlet velocity as boundary conditions for the thermal hydraulic problem. Imposed quantities are:
        * inlet pressure
        * inlet temperature (or outlet temperature if backflow)
        * outlet velocity
    The function suitably accounts for forward (left to right)/backward (right to left) flow as well as for back flow (part of the fluid moves towards the inlet due to a pressure build up).

    Args:
        ndarrays (tuple): collection of known term vector (Known) and system matrix (SYSMAT) np.ndarray on wich impose the boundary conditions.
        f_comp (FluidComponent): fluid component object from which get all info to apply the boundary conditions.
        conductor (Conductor): object with all the information of the conductor.
        path (str): path of the auxiliary input file with the values for the boundary conditions (if INTIAL = -2).

    Returns:
        tuple: collection of known term vector (Known) and system matrix (SYSMAT) np.ndarray with imposed boundary conditions.
    """

    # INLET AND OUTLET RESERVOIRS, INLET CONDITIONS AND FLOW SPECIFIED
    known,sysmat = ndarrays
    # Get boundary condition values.
    if f_comp.coolant.operations["INTIAL"] == 2:
        # outlet mass flow rate.
        mfr_out = f_comp.coolant.operations["MDTOUT"]
        # inlet pressure.
        p_inl = f_comp.coolant.operations["PREINL"]
        # inlet temperature.
        T_inl = f_comp.coolant.operations["TEMINL"]
        # outlet temperature: to be assigned if outlet velocity is negative.
        T_out = f_comp.coolant.operations["TEMOUT"]
    else:  # N.B va aggiustato per renderlo conforme caso positivo!
        # all values from flow_dummy.xlsx: call get_from_xlsx.
        [flow_par, flagSpecfield] = get_from_xlsx(
            conductor,
            path,
            f_comp,
            "INTIAL",
            f_comp.coolant.operations["INTIAL"],
        )
        print(
            f"""flagSpecfield == {flagSpecfield}: still to be decided if it 
        useful and if yes still to be defined\n"""
        )
        mfr_out = flow_par[3]  # inlet mass flow rate.
        p_inl = flow_par[2]  # inlet pressure.
        T_inl = flow_par[0]  # inlet temperature.
        # outlet temperature: to be assigned if outlet velocity is negative \
        #.
        T_out = flow_par[1]

    # ALIAS
    _, inl_p_idx, inl_t_idx = f_comp.inl_idx # tuple unpack
    out_v_idx, _, out_t_idx = f_comp.out_idx # tuple unpack
    main_d_idx = conductor.dict_band["Main_diag"]
    flow_dir = f_comp.channel.flow_dir[0]
    velocity = f_comp.coolant.dict_node_pt["velocity"]
    density = f_comp.coolant.dict_node_pt["total_density"]
    cross_section = f_comp.channel.inputs["CROSSECTION"]
    
    # Assign BC
    if flow_dir == "forward":
        # Flow direction from x = 0 to x = L.
        # v_out
        sysmat[:,out_v_idx.forward] = 0.0
        # main diagonal.
        sysmat[main_d_idx,out_v_idx.forward] = 1.0
        known[out_v_idx.forward] = (mfr_out / density[-1] / cross_section)
        # p_inl
        sysmat[0:, inl_p_idx.forward] = 0.0
        # main diagonal.
        sysmat[main_d_idx, inl_p_idx.forward] = 1.0
        known[inl_p_idx.forward] = p_inl
        # T_inl
        if velocity[0] > 0:
            sysmat[:,inl_t_idx.forward] = 0.0
            # main diagonal.
            sysmat[main_d_idx,inl_t_idx.forward] = 1.0
            known[inl_t_idx.forward] = T_inl
        # T_out (T_inl if mfr_out < 0)
        if velocity[-1] < 0:
            sysmat[:,out_t_idx.forward] = 0.0
            # main diagonal.
            sysmat[main_d_idx,out_t_idx.forward] = 1.0
            known[out_t_idx.forward] = T_out
    elif flow_dir == "backward":
        # Flow direction from x = L to x = 0.
        # v_out
        sysmat[:,out_v_idx.backward] = 0.0
        # main diagonal.
        sysmat[main_d_idx,out_v_idx.backward] = 1.0
        known[out_v_idx.backward] = (mfr_out / density[0] / cross_section)
        # p_inl
        sysmat[0:, inl_p_idx.backward] = 0.0
        # main diagonal.
        sysmat[main_d_idx, inl_p_idx.backward] = 1.0
        known[inl_p_idx.backward] = p_inl
        # T_inl
        if velocity[-1] < 0:
            sysmat[:,inl_t_idx.backward] = 0.0
            # main diagonal.
            sysmat[main_d_idx,inl_t_idx.backward] = 1.0
            known[inl_t_idx.backward] = T_inl
        # T_out
        if velocity[0] > 0:
            sysmat[:,out_t_idx.backward] = 0.0
            # main diagonal.
            sysmat[main_d_idx,out_t_idx.backward] = 1.0
            known[out_t_idx.backward] = T_out
    
    return known,sysmat

def impose_inl_v_out_p(
    ndarrays:tuple,
    f_comp:FluidComponent,
    conductor:Conductor,
    path:str,
    )->tuple:
    
    """Function that imposes a inlet velocity and outlet pressure as boundary conditions for the thermal hydraulic problem. Imposed quantities are:
        * inlet velocity
        * inlet temperature (or outlet temperature if backflow)
        * outlet pressure
    The function suitably accounts for forward (left to right)/backward (right to left) flow as well as for back flow (part of the fluid moves towards the inlet due to a pressure build up).

    Args:
        ndarrays (tuple): collection of known term vector (Known) and system matrix (SYSMAT) np.ndarray on wich impose the boundary conditions.
        f_comp (FluidComponent): fluid component object from which get all info to apply the boundary conditions.
        conductor (Conductor): object with all the information of the conductor.
        path (str): path of the auxiliary input file with the values for the boundary conditions (if INTIAL = -3).

    Returns:
        tuple: collection of known term vector (Known) and system matrix (SYSMAT) np.ndarray with imposed boundary conditions.
    """

    known,sysmat = ndarrays
    # Get boundary condition values.
    if f_comp.coolant.operations["INTIAL"] == 3:
        # outlet mass flow rate.
        mfr_inl = f_comp.coolant.operations["MDTIN"]
        # outlet pressure.
        p_out = f_comp.coolant.operations["PREOUT"]
        # inlet temperature.
        T_inl = f_comp.coolant.operations["TEMINL"]
        # outlet temperature: to be assigned if outlet velocity is negative.
        T_out = f_comp.coolant.operations["TEMOUT"]
    else:  # N.B va aggiustato per renderlo conforme caso positivo!
        # all values from flow_dummy.xlsx: call get_from_xlsx.
        [flow_par, flagSpecfield] = get_from_xlsx(
            conductor,
            path,
            f_comp,
            "INTIAL",
            f_comp.coolant.operations["INTIAL"],
        )
        print(
            f"""flagSpecfield == {flagSpecfield}: still to be decided if it 
        useful and if yes still to be defined\n"""
        )
        mfr_inl = flow_par[3]  # inlet mass flow rate.
        p_out = flow_par[2]  # outlet pressure.
        T_inl = flow_par[0]  # inlet temperature.
        # outlet temperature: to be assigned if outlet velocity is negative.
        T_out = flow_par[1]

    # ALIAS
    inl_v_idx, _, inl_t_idx = f_comp.inl_idx # tuple unpack
    _, out_p_idx, out_t_idx = f_comp.out_idx # tuple unpack
    main_d_idx = conductor.dict_band["Main_diag"]
    flow_dir = f_comp.channel.flow_dir[0]
    velocity = f_comp.coolant.dict_node_pt["velocity"]
    density = f_comp.coolant.dict_node_pt["total_density"]
    cross_section = f_comp.channel.inputs["CROSSECTION"]
    
    # Assign BC
    if flow_dir == "forward":
        # Flow direction from x = 0 to x = L.
        # v_inl
        sysmat[:,inl_v_idx.forward] = 0.0
        # main diagonal.
        sysmat[main_d_idx,inl_v_idx.forward] = 1.0
        known[inl_v_idx.forward] = (mfr_inl / density[0] / cross_section)
        # p_out
        sysmat[0:, out_p_idx.forward] = 0.0
        # main diagonal.
        sysmat[main_d_idx, out_p_idx.forward] = 1.0
        known[out_p_idx.forward] = p_out
        # T_inl
        if velocity[0] > 0:
            sysmat[:,inl_t_idx.forward] = 0.0
            # main diagonal.
            sysmat[main_d_idx,inl_t_idx.forward] = 1.0
            known[inl_t_idx.forward] = T_inl
        # T_out (T_inl if mfr_inl < 0)
        if velocity[-1] < 0:
            sysmat[:,out_t_idx.forward] = 0.0
            # main diagonal.
            sysmat[main_d_idx,out_t_idx.forward] = 1.0
            known[out_t_idx.forward] = T_out
    elif flow_dir == "backward":
        # Flow direction from x = L to x = 0.
        # v_inl
        sysmat[:,inl_v_idx.backward] = 0.0
        # main diagonal.
        sysmat[main_d_idx,inl_v_idx.backward] = 1.0
        known[inl_v_idx.backward] = (mfr_inl / density[-1] / cross_section)
        # p_out
        sysmat[0:, out_p_idx.backward] = 0.0
        # main diagonal.
        sysmat[main_d_idx, out_p_idx.backward] = 1.0
        known[out_p_idx.backward] = p_out
        # T_inl
        if velocity[-1] < 0:
            sysmat[:,inl_t_idx.backward] = 0.0
            # main diagonal.
            sysmat[main_d_idx,inl_t_idx.backward] = 1.0
            known[inl_t_idx.backward] = T_inl
        # T_out
        if velocity[0] > 0:
            sysmat[:,out_t_idx.backward] = 0.0
            # main diagonal.
            sysmat[main_d_idx,out_t_idx.backward] = 1.0
            known[out_t_idx.backward] = T_out
    
    return known,sysmat