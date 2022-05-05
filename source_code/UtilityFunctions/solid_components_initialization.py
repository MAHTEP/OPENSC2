# Placido Daniele 11/2020

import warnings
import numpy as np
import os


def solid_components_temperature_initialization(cond):
    """
    Function that initializes Solid Components temperature spatial distribution according to conductor topology and to the value of flag INTIAL (cdp, 12/2020)
    """
    T_min = np.zeros(cond.dict_obj_inventory["FluidComponents"]["Number"])
    # Loop on FluidComponents to get the minimum temperature among all the \
    # channels (cdp, 12/2020)
    for rr in range(cond.dict_obj_inventory["FluidComponents"]["Number"]):
        fluid_comp = cond.dict_obj_inventory["FluidComponents"]["Objects"][rr]
        T_min[rr] = fluid_comp.coolant.dict_node_pt["temperature"].min()
    # end for rr (cdp, 12/2020)
    # For each solid component evaluate temperature (cdp, 07/2020)
    # If needed read only the sub matrix describing channel - solid objects \
    # contact (cdp, 07/2020)
    # nested loop on channel - solid objects (cpd 07/2020)
    for cc in range(cond.dict_obj_inventory["SolidComponents"]["Number"]):
        s_comp = cond.dict_obj_inventory["SolidComponents"]["Objects"][cc]
        s_comp.dict_node_pt = dict()  # dictionary declaration (cdp, 07/2020)
        s_comp.dict_Gauss_pt = dict()  # dictionary declaration (cdp, 07/2020)
        # s_comp temperature initialization to 0 (cdp, 12/2020)
        s_comp.dict_node_pt["temperature"] = np.zeros(cond.dict_discretization["N_nod"])
        if s_comp.operations["INTIAL"] == 0:
            # Not user defined temperature initialization (cdp, 12/2020)
            # get contact perimeter flags reading the sub matrix channel - solid by \
            # column, array smart, in this way is possible to determine if s_comp is \
            # in thermal contact with channels or not (cdp, 12/2020)
            contact_flag = cond.dict_df_coupling["contact_perimeter"].iloc[
                0 : cond.dict_obj_inventory["FluidComponents"]["Number"],
                cc + cond.dict_obj_inventory["FluidComponents"]["Number"],
            ]
            if np.sum(contact_flag) > 0:
                # get weight reading the sub matrix channel - solid by column, array \
                # smart (cdp, 12/2020)
                weight = cond.dict_df_coupling["contact_perimeter"].iloc[
                    0 : cond.dict_obj_inventory["FluidComponents"]["Number"],
                    cc + cond.dict_obj_inventory["FluidComponents"]["Number"],
                ]
                # evaluate SolidComponents temperature as the weighted average on \
                # conctat_perimeter with channels (cpd 07/2020)
                for rr in range(cond.dict_obj_inventory["FluidComponents"]["Number"]):
                    fluid_comp = cond.dict_obj_inventory["FluidComponents"]["Objects"][
                        rr
                    ]
                    s_comp.dict_node_pt["temperature"] = s_comp.dict_node_pt[
                        "temperature"
                    ] + fluid_comp.coolant.dict_node_pt["temperature"] * weight[
                        rr
                    ] / np.sum(
                        weight
                    )
            else:
                # The s_comp object is not in thermal contact with channels: the \
                # temperature spatial distribution is initialized at the minimum \
                # temperature among the channels (cdp, 12/2020)
                s_comp.dict_node_pt["temperature"] = (
                    np.ones(cond.dict_discretization["N_nod"]) * T_min.min()
                )
            # end if np.sum(contact_flag) (cdp, 12/2020)
        elif abs(s_comp.operations["INTIAL"]) == 1:
            # User imposed initial temperature spatial distribution (cdp, 12/2020)
            if s_comp.operations["INTIAL"] == 1:
                # linear spatial temperature distribution (cdp, 12/2020)
                s_comp.dict_node_pt["temperature"] = np.interp(
                    cond.dict_discretization["xcoord"],
                    [0.0, cond.inputs["XLENGTH"]],
                    [s_comp.operations["TEMINL"], s_comp.operations["TEMOUT"]],
                )

            elif s_comp.operations["INTIAL"] == -1:
                print("still to do\n")
        else:
            # raise error due to not available INTIAL value (cdp, 12/2020)
            raise ValueError(
                f"""INTIAL value not available. Please check check INTIAL value in sheet {s_comp.NAME} of file {cond.file_input["OPERATION"]}.\n"""
            )
        # end if s_comp.operations["INTIAL"] (cdp, 12/2020)
    # end for cc (cdp, 12/2020)


# end function SolidComponents_T_initialization (cdp, 12/2020)
