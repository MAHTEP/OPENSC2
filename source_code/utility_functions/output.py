import numpy as np
import pandas as pd
import os

from collections import namedtuple

from fluid_component import FluidComponent


def save_properties(conductor, f_path):

    """Functions that save .tsv files with suitable file names, FluidComponent and SolidComponent initialization and final solution, together with spatial coordinate discretization. Channels saved variables are: temperature, pressure, density, viscosity, specific heat at constant pressure, thermal conductivity, velocity, Reynolds number and Prandtl number. StrandComponent saved variables are: temperature, density, specific heat at constant pressure, thermal conductivity, magnetic field, electrical resistivity, current sharing temperature; jackets saved variables are temperature, specific heat at constant pressure, thermal conductivity, magnetic field, electrical resistivity."""

    # Check if FluidComponent collection is not empty.
    if conductor.inventory["FluidComponent"].collection:
        # FludiComponent collection is not empty.
        list_prop_chan = list(
            conductor.inventory["FluidComponent"].collection[0].coolant.dict_node_pt.keys()
        )
        list_prop_chan.append("friction_factor")
        list_prop_chan.insert(0,"zcoord")
        list_prop_chan = tuple(list_prop_chan)
        list_units = (
            "(m)",
            "(Pa)",
            "(K)",
            "(kg/m^3)",
            "(m/s)",
            "(1/K)",
            "(1/Pa)",
            "(~)",
            "(Pa*s)",
            "(J/kg)",
            "(J/kg/K)",
            "(J/kg/K)",
            "(m/s)",
            "(W/m/K)",
            "(~)",
            "(~)",
            "(kg/s)",
            "(~)",
        )

        # Build a tuple of string with property name and property units 
        # exploiting generator expression; used to build header_chan.
        prop_unit_chan = tuple(f"{prop_name} {list_units[prop_idx]}\t" for prop_idx, prop_name in enumerate(list_prop_chan))
        # Build header_chan concatenating strings in prop_unit_chan; the 
        # trailing \t caracter is removed with rstrip.
        header_chan = "".join(prop for prop in prop_unit_chan).rstrip("\t")

    header_st = "zcoord (m)\ttemperature (K)\tB_field (T)\tT_cur_sharing (K)"
    header_stab = "zcoord (m)\ttemperature (K)\tB_field (T)"
    header_jk = "zcoord (m)\ttemperature (K)"
    for fluid_comp in conductor.inventory["FluidComponent"].collection:
        A_chan = np.zeros(
            (
                conductor.grid_features["N_nod"],
                len(fluid_comp.coolant.dict_node_pt) + 2,
            )
        )
        file_path = os.path.join(f_path, f"{fluid_comp.identifier}.tsv")
        A_chan[:, 0] = conductor.grid_features["zcoord"]
        for ii, prop_value in enumerate(fluid_comp.coolant.dict_node_pt.values(), 1):
            A_chan[:, ii] = prop_value
        # Save total friction factor
        A_chan[:, -1] = fluid_comp.channel.dict_friction_factor[True]["total"]
        with open(file_path, "w") as writer:
            np.savetxt(writer, A_chan, delimiter="\t", header=header_chan, comments="")
    for strand in conductor.inventory["StrandComponent"].collection:
        file_path = os.path.join(f_path, f"{strand.identifier}.tsv")
        if strand.name != conductor.inventory["StrandStabilizerComponent"].name:
            A_strand = np.zeros((conductor.grid_features["N_nod"], 4))
            A_strand[:, 3] = strand.dict_node_pt["T_cur_sharing"]
        else:
            A_strand = np.zeros((conductor.grid_features["N_nod"], 3))
        A_strand[:, 0] = conductor.grid_features["zcoord"]
        A_strand[:, 1] = strand.dict_node_pt["temperature"]
        A_strand[:, 2] = strand.dict_node_pt["B_field"]
        with open(file_path, "w") as writer:
            if strand.name != conductor.inventory["StrandStabilizerComponent"].name:
                np.savetxt(
                    writer, A_strand, delimiter="\t", header=header_st, comments=""
                )
            else:
                np.savetxt(
                    writer, A_strand, delimiter="\t", header=header_stab, comments=""
                )
    for jacket in conductor.inventory["JacketComponent"].collection:
        file_path = os.path.join(f_path, f"{jacket.identifier}.tsv")
        A_jacket = np.zeros((conductor.grid_features["N_nod"], 2))
        A_jacket[:, 0] = conductor.grid_features["zcoord"]
        A_jacket[:, 1] = jacket.dict_node_pt["temperature"]
        with open(file_path, "w") as writer:
            np.savetxt(writer, A_jacket, delimiter="\t", header=header_jk, comments="")


# end function Save_properties


def save_sd_nodal_comp(conductor, f_path:str, comp_type:str):
    """Funcition that saves spatial distributions of selected properties at user defined time step. This version of the function deals with spatial distribution in nodal points.

    Args:
        conductor (Conductor): object with all the information of the conductor.
        f_path (str): folder path where to store spatial distributions.
        comp_type (str): class name of the conductor component objects for which spatial distributions are be saved. Possible values:
            * FludiComponent
            * JacketComponent
            * StackComponent
            * StrandMixedComponent
            * StrantStabilizerComponent
    """

    if (
        comp_type == "FluidComponent" 
        or comp_type == "JacketComponent" 
        or comp_type == "StrandStabilizerComponent"):
        # Alias
        n_prop = conductor.relevant_prop_sd_num["node"][comp_type]
        props = conductor.relevant_prop_sd["node"][comp_type]
        headers = conductor.header_sd["node"][comp_type]

        for obj in conductor.inventory[comp_type].collection:
            file_path = os.path.join(
                f_path,
                f"{obj.identifier}_({conductor.cond_num_step})_sd.tsv"
            )
            AA = np.zeros((conductor.grid_features["N_nod"], n_prop))
            AA[:, 0] = conductor.store_sd_node["zcoord"]["t_save"]
            for prop_idx, prop_name in enumerate(props,1):
                AA[:, prop_idx] = obj.store_sd_node[prop_name]["t_save"]
            with open(file_path, "w") as writer:
                np.savetxt(
                    writer, AA, delimiter="\t", header=headers, comments=""
                )
    elif (
        comp_type == "StackComponent"
        or comp_type == "StrandMixedComponent"
        ):
        # Alias
        n_prop = conductor.relevant_prop_sd_num["node"][comp_type]
        props = conductor.relevant_prop_sd["node"][comp_type]
        headers = conductor.header_sd["node"][comp_type]

        for obj in conductor.inventory[comp_type].collection:
            obj_id = obj.identifier
            file_path = os.path.join(
                f_path,
                f"{obj_id}_({conductor.cond_num_step})_sd.tsv"
            )
            AA = np.zeros((conductor.grid_features["N_nod"], n_prop[obj_id]))
            AA[:, 0] = conductor.store_sd_node["zcoord"]["t_save"]
            for prop_idx, prop_name in enumerate(props[obj_id],1):
                AA[:, prop_idx] = obj.store_sd_node[prop_name]["t_save"]
            with open(file_path, "w") as writer:
                np.savetxt(
                    writer, AA, delimiter="\t", header=headers[obj_id], comments=""
                )

def save_sd_gauss_comp(conductor, f_path:str):
    """Funcition that saves spatial distributions of selected properties at user defined time step. This version of the function deals with spatial distribution in Gauss nodal points.

    Args:
        conductor (Conductor): object with all the information of the conductor.
        f_path (str): folder path where to store spatial distributions.
    """

    # Alias
    n_prop = conductor.relevant_prop_sd_num["gauss"]
    props = conductor.relevant_prop_sd["gauss"]["SolidComponent"]
    headers = conductor.header_sd["gauss"]

    for obj in conductor.inventory["SolidComponent"].collection:
        file_path = os.path.join(
            f_path, f"{obj.identifier}_({conductor.cond_num_step})_gauss_sd.tsv"
        )
        AA = np.zeros((conductor.grid_input["NELEMS"], n_prop))
        AA[:, 0] = conductor.store_sd_gauss["zcoord_gauss"]["t_save"]
        for prop_idx, prop_name in enumerate(props,1):
            AA[:, prop_idx] = obj.store_sd_gauss[prop_name]["t_save"]
        with open(file_path, "w") as writer:
            np.savetxt(writer, AA, delimiter="\t", header=headers, comments="")

def save_conductor_sd(conductor, f_path:str):
    """
    Function that saves spatial distributions of heat transfer coefficients between conductor components and heat exchanged between conductor components at user defined time step (on temporary files).
    List of saved quantities:
        * heat transfer coefficients of the open fraction between fluid components
        * heat transfer coefficients of the closed fraction between fluid components
        * heat transfer coefficients between fluid and solid components
        * conductive heat transfer coefficients between solid components
        * radiative heat transfer coefficients between solid components
        * convective heat transfer coefficients between environment and solid components
        * radiative heat transfer coefficients between environment and solid components
        * heat exchanged by radiation between jackets
        * heat exchanged by convection and/or radiation between outer surface of the conductor and the environment

    Args:
        conductor (Conductor): object with all the information of the conductor.
        f_path (str): folder path where to store spatial distributions.
    """

    # Loop on items of dict file_htc_sd_pref to avoid following error.
    # ValueError: The truth value of an array with more than one element is 
    # ambiguous. Use a.any() or a.all().
    # The error arised if the loop is on items of dict store_sd_node 
    # because keyword "zcoord" stores array and not a dictionary. Another 
    # possible solution is to loop on store_sd_node and check that key is not 
    # zcoord. The chosen soltution does not require this check.

    # Save heat transfer coefficients between conductor components.
    for key, val in conductor.file_htc_sd_pref.items():
        if conductor.store_sd_node[key]["t_save"]:
            file_name = f"{val}_({conductor.cond_num_step})_sd.tsv"
            file_path = os.path.join(f_path, file_name)
            # Build the dataframe from dictionary and save it as tsv file.
            pd.DataFrame.from_dict(
                conductor.store_sd_node[key]["t_save"],
                dtype=float,
            ).to_csv(file_path, sep="\t", index=False, header=True)
    
    # Save heat exchanged between conductor components.
    for key, val in conductor.file_heat_sd_pref.items():
        if conductor.store_sd_gauss[key]["t_save"]:
            file_name = f"{val}_({conductor.cond_num_step})_sd.tsv"
            file_path = os.path.join(f_path, file_name)
            # Build the dataframe from dictionary and save it as tsv file.
            pd.DataFrame.from_dict(
                conductor.store_sd_gauss[key]["t_save"],
                dtype=float,
            ).to_csv(file_path, sep="\t", index=False, header=True)

def save_simulation_space(conductor, f_path:str):

    """
    Function that save on files with suitable file names transient solution,
    spatial coordinate discretization, time and time step for each conductor.

    Args:
        conductor (Conductor): object with all the information of the conductor.
        f_path (str): folder path where to store spatial distributions.
    """

    # Save transient solution for each Conductor objects at each 
    # iteration_store iteration. For each conductor component, solution is
    # stored in a dedicated file. FluidComponent stored variables are velocity, 
    # pressure and temperature together with spatial discretization, while for 
    # SolidComponent only temperature and spatial discretization are stored. 
    # Hypothesis: for different input files parameters there are different 
    # simulation names

    # Alias
    conductor.num_step_save[conductor.i_save] = conductor.cond_num_step
    
    component_type = {
        "FluidComponent",
        "JacketComponent",
        "StackComponent",
        "StrandMixedComponent",
        "StrandStabilizerComponent",
    }

    for comp_type in component_type:
        save_sd_nodal_comp(conductor, f_path, comp_type)
    
    # Save linear power due to electric resistance along the SOs (available in
    # gauss nodal points).
    save_sd_gauss_comp(conductor,f_path)

    # Call function to save spatial distributions of heat transfer coefficients 
    # and heat exchanged between conductor components in temporary files.
    save_conductor_sd(conductor, f_path)

# end function save_simulation_space


def reorganize_spatial_distribution(cond, f_path, n_digit_time):
    """
    Function that reorganizes the files of the spatial distribution collecting in a single file for each property the spatial distribution at user defined times. In this way the file format is like the ones of the time evolution and this should simplify plots and furter data analysis. (cdp, 11/2020)
    """
    list_ch_key = [
        "velocity",
        "pressure",
        "temperature",
        "total_density",
        "friction_factor",
    ]
    # list_sol_key = ["temperature", "total_density", "total_isobaric_specific_heat", "total_thermal_conductivity", \
    # "EXTFLX", "JHTFLX"]

    Prop_list = namedtuple("Prop_list",["full","reduced"])

    sol_key = dict(
        sc=Prop_list(
            full=("temperature", "T_cur_sharing", "J_critical"),
            reduced=("temperature", "J_critical"),
        ),
        stab=("temperature",),
        jk=("temperature",),
    )
    list_sol_key_gauss = ("current_along", "delta_voltage_along", "P_along")
    # lists all the file .tsv in subfolder Spatial_distribution (cdp, 11/2020)
    # Round the time to save to n_digit_time digits only once
    time = np.around(cond.Space_save, n_digit_time)

    # declare dictionary to store the spatial diccretizations only once.
    dict_zcoord = dict()
    # Loop to save spatial coordinates.
    for ii,_ in enumerate(cond.Space_save):
        # Check if FluidComponent collection is not empty.
        if cond.inventory["FluidComponent"].collection:
            # FluidComponent collection is not empty.
            comp = cond.inventory["FluidComponent"].collection[0]
        else:
            # FluidComponent collection is empty: use first item in 
            # SolidComponent collection.
            comp = cond.inventory["SolidComponent"].collection[0]

        file_name = f"{comp.identifier}_({cond.num_step_save[ii]})_sd.tsv"
        file_load = os.path.join(f_path, file_name)
        # Load dataframe.
        df = pd.read_csv(file_load, delimiter="\t")
        # store the spatial discretizations at each required time step in file 
        # zcoord.tsv.
        dict_zcoord[f"time = {time[ii]} (s)"] = df["zcoord (m)"]
    # convert the dictionary to a DataFrame
    df_zcoord = pd.DataFrame(dict_zcoord)
    # build file name
    file_name = f"zcoord.tsv"
    path_save = os.path.join(f_path, file_name)
    # save the DataFrame as file zcoord.tsv
    df_zcoord.to_csv(path_save, sep="\t", index=False)

    # loop on FluidComponent (cdp, 11/2020)
    for fluid_comp in cond.inventory["FluidComponent"].collection:
        # create a list of files that have the fluid_comp.identifier and User in the name \
        # exploiting list compreension: these files are the ones that will be \
        # reorganized by this function (cdp, 11/2020)
        # list_ch_file = [ff for ff in list_file if (fluid_comp.identifier in ff and "User" in ff)]
        # declare the dictionary of data frame (cdp, 11/2020)
        dict_df = dict()
        dict_df_new = dict()
        for ii, _ in enumerate(cond.Space_save):
            file_name = f"{fluid_comp.identifier}_({cond.num_step_save[ii]})_sd.tsv"
            file_load = os.path.join(f_path, file_name)
            # Load file file_name as data frame as a value of dictionary \
            # corresponding to key file_name (cdp, 11/2020)
            dict_df[file_name] = pd.read_csv(
                filepath_or_buffer=file_load, delimiter="\t"
            )
            # Delete the old file format.
            os.remove(file_load)
            if ii == 0:
                # get columns names only the first time (cdp, 11/2020)
                header = list(dict_df[file_name].columns.values.tolist())
                for jj, prop in enumerate(list_ch_key):
                    # decompose the data frame in four dataframes (cdp, 11/2020)
                    dict_df_new[prop] = dict_df[file_name].filter(
                        items=[header[jj + 1]]
                    )
                    # rename data frames columns (cdp, 11/2020)
                    dict_df_new[prop].rename(
                        columns={header[jj + 1]: f"time = {time[ii]} (s)"}, inplace=True
                    )
                # end for jj (cdp, 11/2020)
            else:
                for jj, prop in enumerate(list_ch_key):
                    # construct the new data frames with concat method (cdp, 11/2020)
                    dict_df_new[prop] = pd.concat(
                        [
                            dict_df_new[prop],
                            dict_df[file_name].filter(items=[header[jj + 1]]),
                        ],
                        axis=1,
                    )
                    dict_df_new[prop].rename(
                        columns={header[jj + 1]: f"time = {time[ii]} (s)"}, inplace=True
                    )
                # end for jj (cdp, 11/2020)
            # end if ii (cdp, 11/2020)
        # end for ii (cdp, 11/2020)
        # for loop to save the new data frame (cdp, 11/2020)
        for prop in list_ch_key:
            # build file name (cdp, 11/2020)
            file_name = f"{fluid_comp.identifier}_{prop}_sd.tsv"
            # build path to save the file (cdp, 11/2020)
            path_save = os.path.join(f_path, file_name)
            # save the data frame, without the row index name (cdp, 11/2020)
            dict_df_new[prop].to_csv(path_save, sep="\t", index=False)
        # end for prop (cdp, 11/2020)
    # end for fluid_comp (cdp, 11/2020)
    # loop on SolidComponent (cdp, 11/2020)
    for s_comp in cond.inventory["SolidComponent"].collection:
        # declare the dictionary of data frame (cdp, 11/2020)
        dict_df = dict()
        dict_df_new = dict()
        for ii, _ in enumerate(cond.Space_save):
            file_name = f"{s_comp.identifier}_({cond.num_step_save[ii]})_sd.tsv"
            file_name_gauss = (
                f"{s_comp.identifier}_({cond.num_step_save[ii]})_gauss_sd.tsv"
            )
            file_load = os.path.join(f_path, file_name)
            file_load_gauss = os.path.join(f_path, file_name_gauss)
            # Load file file_name as data frame as a value of dictionary \
            # corresponding to key file_name (cdp, 11/2020)
            dict_df[file_name] = pd.read_csv(
                filepath_or_buffer=file_load, delimiter="\t"
            )
            dict_df[file_name_gauss] = pd.read_csv(
                filepath_or_buffer=file_load_gauss, delimiter="\t"
            )
            # Delete the old file format.
            os.remove(file_load)
            os.remove(file_load_gauss)
            if ii == 0:
                # get columns names only the first time (cdp, 11/2020)
                header = list(dict_df[file_name].columns.values.tolist())
                if s_comp.KIND == "Mixed_sc_stab" or s_comp.KIND == "Stack":
                    # Check if current sharing temperature is evaluated at each
                    # thermal time step.
                    if s_comp.operations["TCS_EVALUATION"]:
                        list_sol_key = sol_key["sc"].full
                    else:
                        list_sol_key = sol_key["sc"].reduced
                elif s_comp.KIND == "StrandStabilizerComponent":
                    list_sol_key = sol_key["stab"]
                else:  # Jacket
                    list_sol_key = sol_key["jk"]

                for jj, prop in enumerate(list_sol_key):
                    # decompose the data frame in several dataframes (cdp, 11/2020)
                    dict_df_new[prop] = dict_df[file_name].filter(
                        items=[header[jj + 1]]
                    )
                    # rename data frames columns (cdp, 11/2020)
                    dict_df_new[prop].rename(
                        columns={header[jj + 1]: f"time = {time[ii]} (s)"}, inplace=True
                    )
                header_gauss = list(dict_df[file_name_gauss].columns.values.tolist())
                for jj, prop in enumerate(list_sol_key_gauss):
                    # decompose the data frame in four dataframes (cdp, 11/2020)
                    dict_df_new[prop] = dict_df[file_name_gauss].filter(
                        items=[header_gauss[jj + 1]]
                    )
                    # rename data frames columns (cdp, 11/2020)
                    dict_df_new[prop].rename(
                        columns={header_gauss[jj + 1]: f"time = {time[ii]} (s)"},
                        inplace=True,
                    )
            else:
                for jj, prop in enumerate(list_sol_key):
                    # construct the new data frames with concat method (cdp, 11/2020)
                    dict_df_new[prop] = pd.concat(
                        [
                            dict_df_new[prop],
                            dict_df[file_name].filter(items=[header[jj + 1]]),
                        ],
                        axis=1,
                    )
                    dict_df_new[prop].rename(
                        columns={header[jj + 1]: f"time = {time[ii]} (s)"}, inplace=True
                    )
                for jj, prop in enumerate(list_sol_key_gauss):
                    # construct the new data frames with concat method (cdp, 11/2020)
                    dict_df_new[prop] = pd.concat(
                        [
                            dict_df_new[prop],
                            dict_df[file_name_gauss].filter(
                                items=[header_gauss[jj + 1]]
                            ),
                        ],
                        axis=1,
                    )
                    dict_df_new[prop].rename(
                        columns={header_gauss[jj + 1]: f"time = {time[ii]} (s)"},
                        inplace=True,
                    )
            # end if ii (cdp, 11/2020)
        # end for ii (cdp, 11/2020)
        # for loop to save the new data frame (cdp, 11/2020)
        for prop in list_sol_key:
            # build file name (cdp, 11/2020)
            file_name = f"{s_comp.identifier}_{prop}_sd.tsv"
            # build path to save the file (cdp, 11/2020)
            path_save = os.path.join(f_path, file_name)
            # save the data frame, without the row index name (cdp, 11/2020)
            dict_df_new[prop].to_csv(path_save, sep="\t", index=False)
        for prop in list_sol_key_gauss:
            # build file name (cdp, 11/2020)
            file_name = f"{s_comp.identifier}_{prop}_sd.tsv"
            # build path to save the file (cdp, 11/2020)
            path_save = os.path.join(f_path, file_name)
            # save the data frame, without the row index name (cdp, 11/2020)
            dict_df_new[prop].to_csv(path_save, sep="\t", index=False)
    # end for s_comp (cdp, 11/2020)

    # Manage files with heat exhanged between inner jackets by radiation.
    reorganize_heat_sd(cond, f_path, "Heat_rad_inner", "Heat_rad", n_digit_time)
    # Manage files with heat exhanged between outer conductor surface and environment by convection and/or radiation.
    reorganize_heat_sd(cond, f_path, "Heat_exch_env", "Heat_exch", n_digit_time)

    # Manage files with open heat transfer coefficients between fluid components.
    reorganize_heat_sd(cond, f_path, "HTC_ch_ch_o", "HTC_open", n_digit_time)
    # Manage files with close heat transfer coefficients between fluid components.
    reorganize_heat_sd(cond, f_path, "HTC_ch_ch_c", "HTC_close", n_digit_time)
    # Manage files with heat transfer coefficient between fluid and solid components.
    reorganize_heat_sd(cond, f_path, "HTC_ch_sol", "HTC", n_digit_time)
    # Manage files with conductive heat transfer coefficients between solid components.
    reorganize_heat_sd(cond, f_path, "HTC_sol_sol_cond", "HTC_cond", n_digit_time)
    # Manage files with radiative heat transfer coefficients between solid components.
    reorganize_heat_sd(cond, f_path, "HTC_sol_sol_rad", "HTC_rad", n_digit_time)

    # Manage files with convective heat transfer coefficients between 
    # environment and solid components.
    reorganize_heat_sd(cond, f_path, "HTC_env_sol_conv", "HTC_conv", n_digit_time)
    # Manage files with radiative heat transfer coefficients between 
    # environment and solid components.
    reorganize_heat_sd(cond, f_path, "HTC_env_sol_rad", "HTC_rad", n_digit_time)

# end function Reorganize_spatial_distribution (cdp, 11/2020)


def reorganize_heat_sd(cond, f_path, radix_old, radix_new, n_digit_time):
    """Function that reorganizes the files with the spatial distribution of the heat exchanged between inner jackets by radiation and between the outer surface of the conductor and the environment by convection and/or radiation.

    N.B. Questa funzione potrebbe essere adattata anche per riorganizzare i file delle distribuzioni spaziali dei componenti (deriva da questa con qualche semplificazione). Mi sembra troppo complicata: vedere se possibile semplificare.

    Args:
        cond ([type]): [description]
        f_path ([type]): [description]
        radix_old ([type]): [description]
        radix_new ([type]): [description]
    """
    old = dict()
    new = dict()
    cols = list()
    time = np.around(cond.Space_save, n_digit_time)
    for ii, _ in enumerate(cond.Space_save):
        file_name = f"{radix_old}_({cond.num_step_save[ii]})_sd.tsv"
        file_load = os.path.join(f_path, file_name)
        # Check if file exist and if True load it.
        if os.path.isfile(file_load):
            old[file_name] = pd.read_csv(file_load, delimiter="\t")
            # Delete the old file format.
            os.remove(file_load)
            if ii == 0:
                # get columns names only the first time.
                cols = old[file_name].columns.values.tolist()
                for col in cols:
                    # decompose the old dataframe in a sub set of dataframes.
                    new[col] = old[file_name].filter(items=[col])
                    # rename dataframes columns.
                    new[col].rename(
                        columns={col: f"time = {time[ii]} (s)"}, inplace=True
                    )
            else:
                for col in cols:
                    # construct the new dataframes with concat method (cdp, 11/2020)
                    new[col] = pd.concat(
                        [new[col], old[file_name].filter(items=[col])], axis="columns"
                    )
                    new[col].rename(
                        columns={col: f"time = {time[ii]} (s)"}, inplace=True
                    )
            # end if ii (cdp, 11/2020)
        # End os.path.isfile.
    # end for ii (cdp, 11/2020)
    # for loop to save the new data frame (cdp, 11/2020)
    for col in cols:
        # build file name (cdp, 11/2020)
        file_name = f"{radix_new}_{col}_sd.tsv"
        # build path to save the file (cdp, 11/2020)
        path_save = os.path.join(f_path, file_name)
        # save the data frame, without the row index name (cdp, 11/2020)
        new[col].to_csv(path_save, sep="\t", index=False)


# end function reorganize_heat_sd.


def save_time_evolution_init(simulation:object, conductor:object)-> tuple:
    """Function that performs the initialization steps needed to carry out the saving of time evolutions at the user defined spatial coordinates.
    The function:
        1. initializes the dictionaryies assiciated to each quantity of interest with empty lists, one for each spatial coordinates ad which saving the time evolutions. Those are updated inplace as they are attributes of conductor component instances.
        2. creates the files where the time evolutions will be stored, writing the header of the file only once.

    Args:
        simulation (object): object with all information on the simulation
        conductor (object): object with all information on the conductor

    Returns:
        tuple: collection of strings representing the coordinates at which saving the time evolution of each quantities of interest for each conductor component.
    """

    # Alias
    base_path = simulation.dict_path[
        f"Output_Time_evolution_{conductor.identifier}_dir"
    ]

    # Build keys of the dictionaries that will store the values of the time 
    # evolutions of each quantity of interest. Each key will store a list 
    # (array) with the values of a quantity evaluated at a given spatial 
    # coordinate. The first key of is collection is "time (s)" that will store 
    # the time at which the quantities are saved.
    key_zcoord = (
        "time (s)",
        *(
            f"zcoord = {conductor.Time_save[ii]} (m)"
            for ii in range(conductor.Time_save.size)
        ),
    )

    # Construct the header of the file with input and output quantities of 
    # interest only once.
    headers_inl_out = (
        "time (s)",
        "velocity_inl (m/s)",
        "pressure_inl (Pa)",
        "temperature_inl (K)",
        "total_density_inl (kg/m^3)",
        "mass_flow_rate_inl (kg/s)",
        "velocity_out (m/s)",
        "pressure_out (Pa)",
        "temperature_out (K)",
        "total_density_out (kg/m^3)",
        "mass_flow_rate_out (kg/s)",
    )

    # For each conductor component, initialize dictionaries and save heading 
    # for all the quantities of interest.

    for f_comp in conductor.inventory["FluidComponent"].collection:
        # Loop on velocity, pressure, temperature and total density.
        for prop_name in f_comp.coolant.time_evol.keys():
            # Inizialize the dictionary corresponding to key prop_name to a 
            # dictionary of empty lists for the first time.
            f_comp.coolant.time_evol[prop_name].update(
                {key: list() for key in key_zcoord}
            )
            # Save the headings only ones.
            pd.DataFrame(columns=key_zcoord).to_csv(
                os.path.join(
                    base_path,
                    f"{f_comp.identifier}_{prop_name}_te.tsv"
                ),
                sep="\t",
                index=False,
                header=True,
            )

        # Inizialize the dictionary corresponding to key friction_factor to a 
        # dictionary of empty lists for the first time.
        f_comp.channel.time_evol["friction_factor"].update(
            {key: list() for key in key_zcoord}
        )
        
        # Save the headings only ones.
        pd.DataFrame(columns=key_zcoord).to_csv(
            os.path.join(
                base_path,
                f"{f_comp.identifier}_friction_factor_te.tsv",
            ),
            sep="\t",
            index=False,
            header=True,
        )
        # Save the headings for inlet and outlet quantities only ones.
        pd.DataFrame(columns=headers_inl_out).to_csv(
            os.path.join(
                base_path,
                f"{f_comp.identifier}_inlet_outlet_te.tsv",
            ),
            sep="\t",
            index=False,
            header=True,
        )

    for s_comp in conductor.inventory["SolidComponent"].collection:
        # Loop on temperature, magnetic field and current sharing temperature
        # (if available).
        for prop_name in s_comp.time_evol.keys():
            # Inizialize the dictionary corresponding to key prop_name to a 
            # dictionary of empty lists for the first time.
            s_comp.time_evol[prop_name].update(
                {key: list() for key in key_zcoord}
            )
            # Save the headings only ones.
            pd.DataFrame(columns=key_zcoord).to_csv(
                os.path.join(
                    base_path,
                    f"{s_comp.identifier}_{prop_name}_te.tsv"
                ),
                sep="\t",
                index=False,
                header=True,
            )

        # Loop on current, voltage difference and linear joule power (if 
        # available)
        for prop_name in s_comp.time_evol_gauss.keys():
            # Inizialize the dictionary corresponding to key prop_name to a 
            #  dictionary of empty lists for the first time.
            s_comp.time_evol_gauss[prop_name].update(
                {key: list() for key in key_zcoord}
            )
            # Save the headings only ones.
            pd.DataFrame(columns=key_zcoord).to_csv(
                os.path.join(
                    base_path,
                    f"{s_comp.identifier}_{prop_name}_gauss_te.tsv",
                ),
                sep="\t",
                index=False,
                header=True,
            )

    return key_zcoord

def save_time_evolution(simulation:object, conductor:object):
    """Function that saves the time evolution of the quantities of interest at user defined spatial coordinates (sensor location). The value of the quantities of interest in these coordinates are obtained by interpolation on the mesh. The function updates a dictionary of list of values for each quantity of interest and for each conductor component. When the lenght of these list becomes equal to the CHUNCK_SIZE parameter, the content of these list is written in the corresponding files. This allows to reduce the number of writing process during the simulation.
    The dictionary of conductor componet are updated inplace; therefore this function does not have a return.

    Args:
        simulation (object): object with all information on the simulation
        conductor (object): object with all information on the conductor
    """

    # Alias
    z_sensor = conductor.Time_save
    zcoord = conductor.grid_features["zcoord"]
    zcoord_gauss = conductor.grid_features["zcoord_gauss"]
    base_path = simulation.dict_path[
        f"Output_Time_evolution_{conductor.identifier}_dir"
    ]
    tend = simulation.transient_input["TEND"]

    # Convert conductor.cond_time list to np.array
    time = np.array(conductor.cond_time[-1])

    prop_te = np.zeros(conductor.n_sensor_tot)
    prop_te[0] = time

    # FluidComponent objects
    for f_comp in conductor.inventory["FluidComponent"].collection:
        
        # Loop on velocity, pressure, temperature and total density.
        for prop_name in f_comp.coolant.time_evol.keys():
            
            f_comp.coolant.time_evol[prop_name] = update_time_evolution(
                f_comp.coolant.time_evol[prop_name],
                prop_te,
                z_sensor,
                zcoord,
                f_comp.coolant.dict_node_pt[prop_name]
            )

            # Write the content of the dictionary to file, if conditions are 
            # satisfied.
            f_comp.coolant.time_evol[prop_name] = save_te_on_file(
                conductor,
                f_comp.coolant.time_evol[prop_name],
                os.path.join(
                    base_path,
                    f"{f_comp.identifier}_{prop_name}_te.tsv",
                ),
                tend,
            )

        # Save friction factor time evolution.
        # Update the contend of the dictionary of lists with propertiy values 
        # at selected zcoord and current time.
        f_comp.channel.time_evol["friction_factor"] = update_time_evolution(
            f_comp.channel.time_evol["friction_factor"],
            prop_te,
            z_sensor,
            zcoord,
            f_comp.channel.dict_friction_factor[True]["total"]
        )

        # Write the content of the dictionary to file, if conditions are 
        # satisfied.
        f_comp.channel.time_evol["friction_factor"] = save_te_on_file(
            conductor,
            f_comp.channel.time_evol["friction_factor"],
            os.path.join(
                base_path,
                f"{f_comp.identifier}_friction_factor_te.tsv",
            ),
            tend,
        )

        f_comp.coolant.time_evol_io = update_time_evolution_io(
            f_comp,
            time,
        )
        
        f_comp.coolant.time_evol_io = save_te_on_file_io(
            conductor,
            f_comp.coolant.time_evol_io,
            os.path.join(
                base_path,
                f"{f_comp.identifier}_inlet_outlet_te.tsv",
            ),
            tend
        )

    # SolidComponent objects
    for s_comp in conductor.inventory["SolidComponent"].collection:

        for prop_name in s_comp.time_evol.keys():

            s_comp.time_evol[prop_name] = update_time_evolution(
                s_comp.time_evol[prop_name],
                prop_te,
                z_sensor,
                zcoord,
                s_comp.dict_node_pt[prop_name]
            )

            # Write the content of the dictionary to file, if conditions are 
            # satisfied.
            s_comp.time_evol[prop_name] = save_te_on_file(
                conductor,
                s_comp.time_evol[prop_name],
                os.path.join(
                    base_path,
                    f"{s_comp.identifier}_{prop_name}_te.tsv",
                ),
                tend,
            )

        for prop_name in s_comp.time_evol_gauss.keys():
            # Update the contend of the dictionary of lists with propertiy
            # values at selected zcoord and current time.
            if prop_name == "linear_power_el_resistance":

                s_comp.time_evol_gauss[prop_name] = update_time_evolution(
                    s_comp.time_evol_gauss[prop_name],
                    prop_te,
                    z_sensor,
                    zcoord_gauss,
                    s_comp.dict_Gauss_pt[prop_name][:, 0]
                )

            else:
                s_comp.time_evol_gauss[prop_name] = update_time_evolution(
                    s_comp.time_evol_gauss[prop_name],
                    prop_te,
                    z_sensor,
                    zcoord_gauss,
                    s_comp.dict_Gauss_pt[prop_name]
                )
            # Write the content of the dictionary to file, if conditions are
            # satisfied.
            s_comp.time_evol_gauss[prop_name] = save_te_on_file(
                conductor,
                s_comp.time_evol_gauss[prop_name],
                os.path.join(
                    base_path,
                    f"{s_comp.identifier}_{prop_name}_te.tsv",
                ),
                tend,
            )


    if np.isclose(time, tend):
        # TEND is reached: save the conductor time in file Time.tsv exploiting 
        # pandas series.
        pd.Series(conductor.cond_time, name="time (s)", dtype=float).to_csv(
            os.path.join(
                base_path,
                "Time.tsv",
            ),
            sep="\t",
            header=True,
            index=False,
        )

def update_time_evolution(
        t_evol:dict,
        prop_te:np.ndarray,
        zz:np.ndarray,
        zp:np.ndarray,
        prop:np.ndarray
    )-> dict:
    """Function that updates the dictionary that stores the time evolutions of the property at each user defined sensor coordinates. Time evolution of the properties in those coordiantes are evaluated by means of linear interpolation.

    Args:
        t_evol (dict): dictionary that collects the time evolution in all the user defined sensor coordinates of the property
        prop_te (np.ndarray): array that stores the interpolated time evolution of the properties; it is filled in this function.
        zz (np.ndarray): the coordinates at which evaluate the interpolated time evolution of the properties.
        zp (np.ndarray): the coordinates of the data points with which carry out the interpolation.
        prop (np.ndarray): the vaules of the properties used to carry out the interpolation.

    Returns:
        dict: dictionary with the updated time evolution of the properties in all the user defined sensor coordinates; the time at which this values are evaluated is also stored in the dictionary.
    """

    # Get the time evolution of the quantity of interest (prop) at sensor 
    # location (z_sensor) interpolatin on the mesh. Remember that in index 0 is 
    # stored the value of the time at which time evolution is saved.
    prop_te[1:] = np.interp(zz,zp,prop)
    
    # Update each list with the corresponding value of the time evolution 
    # stored in array prop_te.
    for ii, te_val in enumerate(t_evol.values()):
        te_val.append(prop_te[ii])

    return t_evol

def save_simulation_time(simulation, conductor):

    """
    Function to save time evolution of velocity, pressure, temperature, inlet
    and outlet mass flowrate of channels; temperature, magnetic field and
    current sharing temperature of strands and jackets temperature. (cdp, 08/2020)
    """

    # At each time step find the index corresponding to the maximum node \
    # lower or equal to the user defined coordinate. This is done at each \
    # time step since the spatial discretization may change and/or user may \
    # modify the spatial coordinates wrt which saving the variables (cdp, 08/2020)

    ind_zcoord = {
        f"zcoord = {conductor.Time_save[ii]} (m)": np.max(
            np.nonzero(
                conductor.grid_features["zcoord"]
                <= round(conductor.Time_save[ii], conductor.n_digit_z)
            )
        )
        for ii in range(conductor.Time_save.size)
    }
    ind_zcoord_gauss = {f"zcoord_g = {conductor.Time_save[0]} (m)": 0}
    ind_zcoord_gauss.update(
        {
            f"zcoord_g = {conductor.Time_save[ii]} (m)": np.max(
                np.nonzero(
                    conductor.grid_features["zcoord_gauss"]
                    <= round(conductor.Time_save[ii], conductor.n_digit_z)
                )
            )
            for ii in range(1, conductor.Time_save.size)
        }
    )
    # construct file header only once (cdp, 08/2020)
    if simulation.num_step == 0:
        headers = ["time (s)"]
        headers.extend([str(key) for key in ind_zcoord.keys()])
        headers_gauss = ["time (s)"]
        headers_gauss.extend([str(key) for key in ind_zcoord_gauss.keys()])
        headers_inl_out = [
            "time (s)",
            "velocity_inl (m/s)",
            "pressure_inl (Pa)",
            "temperature_inl (K)",
            "total_density_inl (kg/m^3)",
            "mass_flow_rate_inl (kg/s)",
            "velocity_out (m/s)",
            "pressure_out (Pa)",
            "temperature_out (K)",
            "total_density_out (kg/m^3)",
            "mass_flow_rate_out (kg/s)",
        ]
        for f_comp in conductor.inventory["FluidComponent"].collection:
            # Loop on velocity, pressure, temperature and total density.
            for key, value in f_comp.coolant.time_evol.items():
                # Inizialize dictionary corresponding to key to a dictionary of empty lists for the first time.
                f_comp.coolant.time_evol[key] = initialize_dictionaty_te(
                    value, ind_zcoord
                )
                # Save the headings only ones.
                pd.DataFrame(columns=headers).to_csv(
                    os.path.join(
                        simulation.dict_path[
                            f"Output_Time_evolution_{conductor.identifier}_dir"
                        ],
                        f"{f_comp.identifier}_{key}_te.tsv",
                    ),
                    sep="\t",
                    index=False,
                    header=True,
                )
            # End for key.
            # Inizialize dictionary corresponding to key to a dictionary of empty lists for the first time.
            f_comp.channel.time_evol["friction_factor"] = initialize_dictionaty_te(
                f_comp.channel.time_evol["friction_factor"], ind_zcoord
            )
            # Save the headings only ones.
            pd.DataFrame(columns=headers).to_csv(
                os.path.join(
                    simulation.dict_path[
                        f"Output_Time_evolution_{conductor.identifier}_dir"
                    ],
                    f"{f_comp.identifier}_friction_factor_te.tsv",
                ),
                sep="\t",
                index=False,
                header=True,
            )
            # Save the headings only ones.
            pd.DataFrame(columns=headers_inl_out).to_csv(
                os.path.join(
                    simulation.dict_path[
                        f"Output_Time_evolution_{conductor.identifier}_dir"
                    ],
                    f"{f_comp.identifier}_inlet_outlet_te.tsv",
                ),
                sep="\t",
                index=False,
                header=True,
            )
        # End for f_comp.
        for s_comp in conductor.inventory["SolidComponent"].collection:
            # Loop on velocity, pressure, temperature and total density.
            for key, value in s_comp.time_evol.items():
                # Inizialize dictionary corresponding to key to a dictionary of empty lists for the first time.
                s_comp.time_evol[key] = initialize_dictionaty_te(value, ind_zcoord)
                # Save the headings only ones.
                pd.DataFrame(columns=headers).to_csv(
                    os.path.join(
                        simulation.dict_path[
                            f"Output_Time_evolution_{conductor.identifier}_dir"
                        ],
                        f"{s_comp.identifier}_{key}_te.tsv",
                    ),
                    sep="\t",
                    index=False,
                    header=True,
                )
            # End for key.
            for key, value in s_comp.time_evol_gauss.items():
                # Inizialize dictionary corresponding to key to a dictionary of empty lists for the first time.
                s_comp.time_evol_gauss[key] = initialize_dictionaty_te(
                    value, ind_zcoord_gauss
                )
                # Save the headings only ones.
                pd.DataFrame(columns=headers_gauss).to_csv(
                    os.path.join(
                        simulation.dict_path[
                            f"Output_Time_evolution_{conductor.identifier}_dir"
                        ],
                        f"{s_comp.identifier}_{key}_te.tsv",
                    ),
                    sep="\t",
                    index=False,
                    header=True,
                )
            # End for key.
        # End for s_comp.
    # End if simulation.num_step (cdp, 10/2020)

    # convert conductor.cond_time list to np.array (cdp, 10/2020)
    time = np.array(conductor.cond_time[-1])

    # FluidComponent objects (cdp, 08/2020)
    for fluid_comp in conductor.inventory["FluidComponent"].collection:
        # Loop on velocity, pressure, temperature and total density.
        for key, value in fluid_comp.coolant.time_evol.items():
            # Update the contend of the dictionary of lists with propertiy values at selected zcoord and current time.
            fluid_comp.coolant.time_evol[key] = update_values(
                value, fluid_comp.coolant.dict_node_pt[key], time, ind_zcoord
            )
            # Write the content of the dictionary to file, if conditions are satisfied.
            fluid_comp.coolant.time_evol[key] = save_te_on_file(
                conductor,
                fluid_comp.coolant.time_evol[key],
                os.path.join(
                    simulation.dict_path[
                        f"Output_Time_evolution_{conductor.identifier}_dir"
                    ],
                    f"{fluid_comp.identifier}_{key}_te.tsv",
                ),
                simulation.transient_input["TEND"],
                ind_zcoord,
            )
        # End for key.

        # Save friction factor time evolution.
        # Update the contend of the dictionary of lists with propertiy values at selected zcoord and current time.
        fluid_comp.channel.time_evol["friction_factor"] = update_values(
            fluid_comp.channel.time_evol["friction_factor"],
            fluid_comp.channel.dict_friction_factor[True]["total"],
            time,
            ind_zcoord,
        )
        # Write the content of the dictionary to file, if conditions are satisfied.
        fluid_comp.channel.time_evol["friction_factor"] = save_te_on_file(
            conductor,
            fluid_comp.channel.time_evol["friction_factor"],
            os.path.join(
                simulation.dict_path[
                    f"Output_Time_evolution_{conductor.identifier}_dir"
                ],
                f"{fluid_comp.identifier}_friction_factor_te.tsv",
            ),
            simulation.transient_input["TEND"],
            ind_zcoord,
        )

        if fluid_comp.channel.flow_dir[0] == "forward":
            index_inl = 0
            index_out = -1
        elif fluid_comp.channel.flow_dir[0] == "backward":
            index_inl = -1
            index_out = 0

        # Inlet and outlet quantities (cdp, 08/2020)
        file_name_io = os.path.join(
            simulation.dict_path[f"Output_Time_evolution_{conductor.identifier}_dir"],
            f"{fluid_comp.identifier}_inlet_outlet_te.tsv",
        )
        fluid_comp.coolant.time_evol_io["time (s)"].append(time)
        # Append inlet properties to list; use dict.update to avoid error (do not understood why with fluid_comp.coolant.time_evol_io.update does not work).
        dict.update(
            {
                key: value.append(
                    fluid_comp.coolant.dict_node_pt[key.split("_inl")[0]][index_inl]
                )
                for key, value in fluid_comp.coolant.time_evol_io.items()
                if "inl" in key
            }
        )
        # Append outlet properties to list.
        dict.update(
            {
                key: value.append(
                    fluid_comp.coolant.dict_node_pt[key.split("_out")[0]][index_out]
                )
                for key, value in fluid_comp.coolant.time_evol_io.items()
                if "out" in key
            }
        )
        # Write the content of the dictionary to file, if conditions are satisfied.
        if len(fluid_comp.coolant.time_evol_io["time (s)"]) == conductor.CHUNCK_SIZE:
            pd.DataFrame(
                fluid_comp.coolant.time_evol_io,
                columns=list(fluid_comp.coolant.time_evol_io.keys()),
                dtype=float,
            ).to_csv(
                file_name_io,
                sep="\t",
                mode="a",
                chunksize=conductor.CHUNCK_SIZE,
                index=False,
                header=False,
            )
            # Initialize empty dictionary.
            fluid_comp.coolant.time_evol_io.update(
                {key: list() for key in fluid_comp.coolant.time_evol_io.keys()}
            )
        elif (
            abs(conductor.cond_time[-1] - simulation.transient_input["TEND"])
            / simulation.transient_input["TEND"]
            <= 1e-6
        ):
            pd.DataFrame(
                fluid_comp.coolant.time_evol_io,
                columns=list(fluid_comp.coolant.time_evol_io.keys()),
                dtype=float,
            ).to_csv(
                file_name_io,
                sep="\t",
                mode="a",
                chunksize=conductor.CHUNCK_SIZE,
                index=False,
                header=False,
            )
        # End if len().
    # End for fluid_comp.

    # SolidComponent objects (cdp, 08/2020)
    for s_comp in conductor.inventory["SolidComponent"].collection:
        for key, value in s_comp.time_evol.items():
            # Update the contend of the dictionary of lists with propertiy values at selected zcoord and current time.
            s_comp.time_evol[key] = update_values(
                value, s_comp.dict_node_pt[key], time, ind_zcoord
            )
            # Write the content of the dictionary to file, if conditions are satisfied.
            s_comp.time_evol[key] = save_te_on_file(
                conductor,
                s_comp.time_evol[key],
                os.path.join(
                    simulation.dict_path[
                        f"Output_Time_evolution_{conductor.identifier}_dir"
                    ],
                    f"{s_comp.identifier}_{key}_te.tsv",
                ),
                simulation.transient_input["TEND"],
                ind_zcoord,
            )
        # End for key.
        for key, value in s_comp.time_evol_gauss.items():
            # Update the contend of the dictionary of lists with propertiy
            # values at selected zcoord and current time.
            if key == "linear_power_el_resistance":
                s_comp.time_evol_gauss[key] = update_values(
                    value, s_comp.dict_Gauss_pt[key][:, 0], time, ind_zcoord_gauss
                )
            else:
                s_comp.time_evol_gauss[key] = update_values(
                    value, s_comp.dict_Gauss_pt[key], time, ind_zcoord_gauss
                )
            # Write the content of the dictionary to file, if conditions are
            # satisfied.
            s_comp.time_evol_gauss[key] = save_te_on_file(
                conductor,
                s_comp.time_evol_gauss[key],
                os.path.join(
                    simulation.dict_path[
                        f"Output_Time_evolution_{conductor.identifier}_dir"
                    ],
                    f"{s_comp.identifier}_{key}_te.tsv",
                ),
                simulation.transient_input["TEND"],
                ind_zcoord_gauss,
            )
        # End for key.
    # End for s_comp.

    if (
        abs(conductor.cond_time[-1] - simulation.transient_input["TEND"])
        / simulation.transient_input["TEND"]
        <= 1e-6
    ):
        # TEND is reached: save the conductor time in file Time.tsv exploiting pandas series
        pd.Series(conductor.cond_time, name="time (s)", dtype=float).to_csv(
            os.path.join(
                simulation.dict_path[
                    f"Output_Time_evolution_{conductor.identifier}_dir"
                ],
                "Time.tsv",
            ),
            sep="\t",
            header=True,
            index=False,
        )
    # End if abs.


# end function Save_simulation_time (cdp, 08/2020)


def initialize_dictionaty_te(val, ind_zcoord):

    val = {"time (s)": list()}
    val.update({key: list() for key in ind_zcoord.keys()})
    return val


# End function initialize_dictionaty_te.


def update_values(val, prop, time, ind_zcoord):

    val["time (s)"].append(time)
    # Use dict.update to avoid error (do not understood why with val.update does not work).
    dict.update(
        {
            key: value.append(prop[ind_zcoord[key]])
            for key, value in val.items()
            if "zcoord" in key
        }
    )
    return val


# End function update_values.


def save_te_on_file(
        conductor:object,
        t_evol:dict,
        file_name:str,
        tend:float
    )->dict:
    """Function that saves the time evolution of selectet variables at given saptial coordinates. If arrays have reached a lenght equal to CHUNK_SIZE, their values are stored in appropriate files and the conten of the array is cleared (new empty list are created) to perform a new storage cycle. This allows to reduce the number of times the code reads and writes files.

    Args:
        conductor (object): object with all information on the conductor.
        t_evol (dict): dictionary that collects the time evolution in all the user defined sensor coordinates of the property.
        file_name (str): name of the file where to save time evolutions at user defined sensor coordinates.
        tend (float): end time of the simulation.

    Returns:
        dict: re-initialized dictionary t_evol if list lenght is equal to the CHUNCK_SIZE; else the not modified dictionary t_evol.
    """

    # Alias
    key_zcoord = conductor.key_zcoord
    chunck_size = conductor.CHUNCK_SIZE
    time = conductor.cond_time[-1]

    if len(t_evol["time (s)"]) == chunck_size:
        pd.DataFrame(t_evol, columns=list(t_evol.keys()), dtype=float).to_csv(
            file_name,
            sep="\t",
            mode="a",
            chunksize=chunck_size,
            index=False,
            header=False,
        )
        # Initiazlie t_evol with empty list to start a new saving cycle.
        t_evol = {z_keys: list() for z_keys in key_zcoord}

    elif np.isclose(time, tend):
        pd.DataFrame(t_evol, columns=list(t_evol.keys()), dtype=float).to_csv(
            file_name,
            sep="\t",
            mode="a",
            chunksize=chunck_size,
            index=False,
            header=False,
        )

    return t_evol

def update_time_evolution_io(obj:FluidComponent, time:float)->dict:
    """Function that updates the dictionary time_evol_io that stores the time evolutions of the property at the inlet and at the outlet. Since these coordinates are fixed and correspont to indez 0 and -1, there is no need to carry out interpolation to evaluate the values. This function can be used only with instances of class FluidComponent.

    Args:
        obj (FluidComponent): instance of class FluidComponent
        time (float): time at which time evolution should be saved

    Raises:
        TypeError: if obj is not an instance of class FluidComponent
        KeyError: if not valid keys are assigned to dictionary time_evol_io.

    Returns:
        dict: updated dictionary time_evol_io with time evolution of quantities of interest at the inlet and at the outlet.
    """

    if not isinstance(obj,FluidComponent):
        raise TypeError("obj should be an instance of class FluidComponent")

    if obj.channel.flow_dir[0] == "forward":
        index_inl = 0
        index_out = -1
    elif obj.channel.flow_dir[0] == "backward":
        index_inl = -1
        index_out = 0

    # Inlet and outlet quantities

    for key, val in obj.coolant.time_evol_io.items():
        if "inl" in key:
            val.append(
                obj.coolant.dict_node_pt[key.split("_inl")[0]][index_inl]
            )
        elif "out" in key:
            val.append(
                obj.coolant.dict_node_pt[key.split("_out")[0]][index_out]
            )
        elif "time" in key:
            val.append(time)
        else:
            raise KeyError(f"Not valid key {key} in dictionary time_evol_io.")

    return obj.coolant.time_evol_io


def save_geometry_discretization(collection: list, file_path: str):
    """Function used to save the coordinates of the barycenter of each conductor component in file with .tsv extension.
    Cartesian reference frame is used.

    Args:
        collection (list): list with all the conductor component objects
        file_path (str): path where to save the file with the geometry discretization.
    """

    [
        pd.DataFrame(comp.coordinate).to_csv(
            os.path.join(file_path, f"{comp.identifier}_barycenter.tsv"),
            sep="\t",
            index=False,
            header=True,
        )
        for comp in collection
    ]
