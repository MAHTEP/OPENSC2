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