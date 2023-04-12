import numpy as np
import pandas as pd
import os


def save_properties(conductor, f_path):

    """Functions that save .tsv files with suitable file names, FluidComponent and SolidComponent initialization and final solution, together with spatial coordinate discretization. Channels saved variables are: temperature, pressure, density, viscosity, specific heat at constant pressure, thermal conductivity, velocity, Reynolds number and Prandtl number. StrandComponent saved variables are: temperature, density, specific heat at constant pressure, thermal conductivity, magnetic field, electrical resistivity, current sharing temperature; jackets saved variables are temperature, specific heat at constant pressure, thermal conductivity, magnetic field, electrical resistivity."""

    # Check if FluidComponent collection is not empty.
    if conductor.inventory["FluidComponent"].collection:
        # FludiComponent collection is not empty.
        list_prop_chan = list(
            conductor.inventory["FluidComponent"].collection[0].coolant.dict_node_pt.keys()
        )
        list_prop_chan.append("friction_factor")
        list_units = [
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
        ]
        header_chan = "zcoord (m)"
        for jj in range(len(list_prop_chan)):
            header_chan = f"{header_chan}\t{list_prop_chan[jj]} {list_units[jj]}"
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


def save_simulation_space(conductor, f_path, n_digit_time):

    """
    Function that save on files with suitable file names transient solution,
    spatial coordinate discretization, time and time step for each conductor
    (cdp, 08/2020)
    """

    # Save transient solution for each Conductor objects at each \
    # iteration_store iteration. For each conductor component, solution is \
    # stored in a dedicated file. FluidComponent stored variables are velocity, \
    # pressure and temperature together with spatial discretization, while for \
    # SolidComponent only temperature and spatial discretization are stored. \
    # (cdp, 08/2020)
    # Hypothesis: for different input files parameters there are different \
    # simulation names (cdp, 08/2020)

    # get the time at which the saving is made
    time = round(conductor.Space_save[conductor.i_save], n_digit_time)
    conductor.num_step_save[conductor.i_save] = conductor.cond_num_step
    # end if len(tt[0]) (cdp, 12/2020)
    prop_chan = [
        "zcoord",
        "velocity",
        "pressure",
        "temperature",
        "total_density",
        "friction_factor",
    ]
    header_chan = "zcoord (m)\tvelocity (m/s)\tpressure (Pa)\ttemperature (K)\ttotal_density (kg/m^3)\tfriction_factor (~)"
    for fluid_comp in conductor.inventory["FluidComponent"].collection:
        file_path = os.path.join(
            f_path, f"{fluid_comp.identifier}_({conductor.cond_num_step})_sd.tsv"
        )
        A_chan = np.zeros((conductor.grid_features["N_nod"], len(prop_chan)))
        for ii in range(len(prop_chan)):
            if prop_chan[ii] == "zcoord":
                A_chan[:, ii] = conductor.grid_features[prop_chan[ii]]
            elif prop_chan[ii] != "friction_factor":
                A_chan[:, ii] = fluid_comp.coolant.dict_node_pt[prop_chan[ii]]
            else:
                # Save friction factor
                A_chan[:, ii] = fluid_comp.channel.dict_friction_factor[True]["total"]

            # end if prop_chan[ii] (cdp, 01/2021)
        # end for ii (cdp, 01/2021)
        with open(file_path, "w") as writer:
            np.savetxt(writer, A_chan, delimiter="\t", header=header_chan, comments="")
    # end for fluid_comp (cdp, 10/2020)
    # headers_s_comp = "zcoord (m)\ttemperature (K)\tdensity (kg/m^3)\tspec_heat_p (J/kg/K)\tther_cond (W/m/K)\tEXFTLX (W/m)\tJHTFLX (W/m^2)"
    # prop_s_comp = ["zcoord", "temperature", "total_density", "total_isobaric_specific_heat", \
    # "total_thermal_conductivity", "EXTFLX", "JHTFLX"]
    headers_full = "zcoord (m)\ttemperature (K)\tcurrent_sharing_temperature (K)"
    headers_reduced = "zcoord (m)\ttemperature (K)"
    prop_full = ["zcoord", "temperature", "T_cur_sharing"]
    prop_reduced = ["zcoord", "temperature"]
    for strand in conductor.inventory["StrandComponent"].collection:
        file_path = os.path.join(
            f_path, f"{strand.identifier}_({conductor.cond_num_step})_sd.tsv"
        )
        if strand.KIND != "StrandStabilizerComponent":
            # Check if current sharing temperature is evaluated at each
            # thermal time step.
            if strand.operations["TCS_EVALUATION"]:
                headers_strand = headers_full
                prop_strand = prop_full
            else:
                headers_strand = headers_reduced
                prop_strand = prop_reduced
        else:
            headers_strand = headers_reduced
            prop_strand = prop_reduced

        A_strand = np.zeros((conductor.grid_features["N_nod"], len(prop_strand)))
        for ii in range(len(prop_strand)):
            if prop_strand[ii] == "zcoord":
                A_strand[:, ii] = conductor.grid_features[prop_strand[ii]]
            else:
                A_strand[:, ii] = strand.dict_node_pt[prop_strand[ii]]
            # end if prop_strand[ii]
        # end for ii
        with open(file_path, "w") as writer:
            np.savetxt(
                writer, A_strand, delimiter="\t", header=headers_strand, comments=""
            )

    headers_jk = "zcoord (m)\ttemperature (K)"
    prop_jk = ["zcoord", "temperature"]
    # Loop to save jacket properties spatial distribution.
    for jk in conductor.inventory["JacketComponent"].collection:
        file_path = os.path.join(
            f_path, f"{jk.identifier}_({conductor.cond_num_step})_sd.tsv"
        )
        A_jk = np.zeros((conductor.grid_features["N_nod"], len(prop_jk)))
        for ii in range(len(prop_jk)):
            if prop_jk[ii] == "zcoord":
                A_jk[:, ii] = conductor.grid_features[prop_jk[ii]]
            else:
                A_jk[:, ii] = jk.dict_node_pt[prop_jk[ii]]
            # end if prop_s_comp[ii] (cdp, 01/2021)
        # end for ii (cdp, 01/2021)
        with open(file_path, "w") as writer:
            np.savetxt(writer, A_jk, delimiter="\t", header=headers_jk, comments="")
    # end for s_comp (cdp, 10/2020)
    # Save linear power due to electric resistance along the SOs (available in
    # gauss nodal points).
    headers_s_comp = (
        "zcoord_gauss (m)\tcurrent_along (A)\tvoltage_drop_along (V)\tP_along (W/m)"
    )
    prop_s_comp = [
        "zcoord_gauss",
        "current_along",
        "voltage_drop_along",
        "linear_power_el_resistance",
    ]
    for s_comp in conductor.inventory["SolidComponent"].collection:
        file_path = os.path.join(
            f_path, f"{s_comp.identifier}_({conductor.cond_num_step})_gauss_sd.tsv"
        )
        A_s_comp = np.zeros((conductor.grid_input["NELEMS"], len(prop_s_comp)))
        for ii in range(len(prop_s_comp)):
            if prop_s_comp[ii] == "zcoord_gauss":
                A_s_comp[:, ii] = conductor.grid_features[prop_s_comp[ii]]
            else:
                if prop_s_comp[ii] == "linear_power_el_resistance":
                    A_s_comp[:, ii] = s_comp.dict_Gauss_pt[prop_s_comp[ii]][:, 0]
                else:
                    A_s_comp[:, ii] = s_comp.dict_Gauss_pt[prop_s_comp[ii]]
            # end if prop_s_comp[ii] (cdp, 01/2021)
        # end for ii (cdp, 01/2021)
        with open(file_path, "w") as writer:
            np.savetxt(
                writer, A_s_comp, delimiter="\t", header=headers_s_comp, comments=""
            )

    # Check if dictionary conductor.heat_rad_jk is not empty in order to save the content in a file.
    if bool(conductor.heat_rad_jk):
        # Build path to save temporary file with the spatial distribution of the heat exchanged by radiation between jackets at each required time step.
        f_path_heat_rad = os.path.join(
            f_path, f"Heat_rad_inner_({conductor.cond_num_step})_sd.tsv"
        )
        # Build the dataframe from dictionary and save it as tsv file.
        pd.DataFrame.from_dict(conductor.heat_rad_jk, dtype=float).to_csv(
            f_path_heat_rad, sep="\t", index=False, header=True
        )

    if bool(conductor.heat_exchange_jk_env):
        # Build path to save temporary file with the spatial distribution of the heat exchanged by convection and/or radiation between outer surface of the conductor and the environment at each required time step.
        f_path_ex_jk_env = os.path.join(
            f_path, f"Heat_exch_env_({conductor.cond_num_step})_sd.tsv"
        )
        # Build the dataframe from dictionary and save it as tsv file.
        pd.DataFrame.from_dict(conductor.heat_exchange_jk_env, dtype=float).to_csv(
            f_path_ex_jk_env, sep="\t", index=False, header=True
        )

    if bool(conductor.dict_node_pt["HTC"]["ch_ch"]["Open"]):
        # Path to save temporary file with the open heat transfer coefficients between fluid components.
        f_path_htc_ch_ch_o = os.path.join(
            f_path, f"HTC_ch_ch_o_({conductor.cond_num_step})_sd.tsv"
        )
        # Build the dataframe from dictionary and save it as tsv file.
        pd.DataFrame.from_dict(
            conductor.dict_node_pt["HTC"]["ch_ch"]["Open"],
            dtype=float,
        ).to_csv(f_path_htc_ch_ch_o, sep="\t", index=False, header=True)

    if conductor.dict_node_pt["HTC"]["ch_ch"]["Close"]:
        # Path to save temporary file with the close heat transfer coefficients between fluid components.
        f_path_htc_ch_ch_c = os.path.join(
            f_path, f"HTC_ch_ch_c_({conductor.cond_num_step})_sd.tsv"
        )
        # Build the dataframe from dictionary and save it as tsv file.
        pd.DataFrame.from_dict(
            conductor.dict_node_pt["HTC"]["ch_ch"]["Close"],
            dtype=float,
        ).to_csv(f_path_htc_ch_ch_c, sep="\t", index=False, header=True)

    if conductor.dict_node_pt["HTC"]["ch_sol"]:
        # Path to save temporary file with the heat transfer coefficients between fluid and solid components.
        f_path_htc_ch_sol = os.path.join(
            f_path, f"HTC_ch_sol_({conductor.cond_num_step})_sd.tsv"
        )
        # Build the dataframe from dictionary and save it as tsv file.
        pd.DataFrame.from_dict(
            conductor.dict_node_pt["HTC"]["ch_sol"],
            dtype=float,
        ).to_csv(f_path_htc_ch_sol, sep="\t", index=False, header=True)

    if conductor.dict_node_pt["HTC"]["sol_sol"]["cond"]:
        # Path to save temporary file with the conductive heat transfer coefficients between solid components.
        f_path_htc_sol_sol_cond = os.path.join(
            f_path, f"HTC_sol_sol_cond_({conductor.cond_num_step})_sd.tsv"
        )
        # Build the dataframe from dictionary and save it as tsv file.
        pd.DataFrame.from_dict(
            conductor.dict_node_pt["HTC"]["sol_sol"]["cond"],
            dtype=float,
        ).to_csv(f_path_htc_sol_sol_cond, sep="\t", index=False, header=True)
    if conductor.dict_node_pt["HTC"]["sol_sol"]["rad"]:
        # Path to save temporary file with the radiative heat transfer coefficients between solid components.
        f_path_htc_sol_sol_rad = os.path.join(
            f_path, f"HTC_sol_sol_rad_({conductor.cond_num_step})_sd.tsv"
        )
        # Build the dataframe from dictionary and save it as tsv file.
        pd.DataFrame.from_dict(
            conductor.dict_node_pt["HTC"]["sol_sol"]["rad"],
            dtype=float,
        ).to_csv(f_path_htc_sol_sol_rad, sep="\t", index=False, header=True)

    # Save the actual times at which the simulation spatial distributions are \
    # saved (cdp, 01/2021)
    file_name = "Time_sd_actual.tsv"
    path_save = os.path.join(f_path, file_name)
    if conductor.i_save == 0:
        with open(path_save, "w") as writer:
            np.savetxt(
                writer,
                np.array([conductor.cond_time[-1]]),
                header="time (s)",
                comments="",
                delimiter="\t",
            )
    else:
        with open(path_save, "a") as writer:
            np.savetxt(
                writer, np.array([conductor.cond_time[-1]]), comments="", delimiter="\t"
            )
    # end if conductor.i_save (cdp, 01/2021)
    # update conductor.i_save (cdp, 01/2021)
    conductor.i_save = conductor.i_save + 1


# end function Save_simulation_space (cdp, 10/2020)


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
    list_sol_key_full = ["temperature", "T_cur_sharing"]
    list_sol_key_reduced = ["temperature"]
    list_sol_key_gauss = ["current_along", "voltage_drop_along", "P_along"]
    list_jk = ["temperature"]
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
                        list_sol_key = list_sol_key_full
                    else:
                        list_sol_key = list_sol_key_reduced
                elif s_comp.KIND == "StrandStabilizerComponent":
                    list_sol_key = list_sol_key_reduced
                else:  # Jacket
                    list_sol_key = list_jk

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
                for jj, _ in enumerate(list_sol_key_gauss):
                    prop = list_sol_key_gauss[jj]
                    # decompose the data frame in four dataframes (cdp, 11/2020)
                    dict_df_new[prop] = dict_df[file_name].filter(
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
                        prop = list_sol_key_gauss[jj]
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


def save_te_on_file(conductor, val, file_name, tend, ind_zcoord):
    """Function that saves the time evolution of selectet variables at given saptial coordinates.

    Args:
        conductor ([type]): [description]
        df ([type]): [description]
        file_name ([type]): [description]
        tend ([type]): [description]

    Returns:
        [type]: [description]
    """
    if len(val["time (s)"]) == conductor.CHUNCK_SIZE:
        pd.DataFrame(val, columns=list(val.keys()), dtype=float).to_csv(
            file_name,
            sep="\t",
            mode="a",
            chunksize=conductor.CHUNCK_SIZE,
            index=False,
            header=False,
        )
        val = initialize_dictionaty_te(val, ind_zcoord)
    elif abs(conductor.cond_time[-1] - tend) / tend <= 1e-6:
        pd.DataFrame(val, columns=list(val.keys()), dtype=float).to_csv(
            file_name,
            sep="\t",
            mode="a",
            chunksize=conductor.CHUNCK_SIZE,
            index=False,
            header=False,
        )
    # End if len(df.index).
    return val


# End function save_te_on_file.


def save_convergence_data(cond, f_path, *n_digit_time, space_conv=True):

    """
    Function that saves data for the space convergence: saved information is the solution spatial distribution at TEND for the given number of elements, together with the global mass and energy balance on channels.
    The spatial distributions are stored in a folder whose name is given by the combination of TEND and STPMIN, files name is comp.ID_(NELEMS).tsv.
    Mass and energy balance for every NELEMS are collected in file cond.ID_mass_energy_sc.tsv in each channel folder.
    """

    if space_conv:
        # Save data for the Space convergence analysis (cdp, 12/2020)
        # compute spatial discretization pitch (cdp, 12/2020)
        discr = cond.inputs["ZLENGTH"] / cond.grid_input["NELEMS"]
        folder_path = os.path.join(f_path, cond.identifier)
        # Create the path of the file {cond.identifier}_delta_x.tsv (cdp, 11/2020)
        file_path_name = os.path.join(folder_path, f"{cond.identifier}_delta_x.tsv")
        discr_header = "delta_x (m)"
        # desinence to sictinguisch among space and time convergence (cdp, 12/2020)
        des = "sc"
        # the content of the round brackets in the file name (cdp, 12/2020)
        brackets = cond.grid_input["NELEMS"]
        # convergence on mass and energy balance (cdp, 12/2020)
        AA = np.zeros((1, 4))
        AA[0, 0] = cond.grid_input["NELEMS"]
        AA[0, 1] = discr
        AA[0, 2] = cond.mass_balance
        AA[0, 3] = cond.energy_balance
        # discretization values for file CONDUCTOR_ID_delta_x.tsv
        val = np.zeros((1, 2))
        val[0, 0] = cond.grid_input["NELEMS"]
        val[0, 1] = discr
    elif space_conv == False:
        # Save data for the Time convergence analysis (cdp, 12/2020)
        # compute time step, remember that n_digit_time is a tuple (cdp, 12/2020)
        discr = round(cond.time_step, n_digit_time[0])
        folder_path = os.path.join(f_path, cond.identifier)
        # Create the path of the file {cond.identifier}_delta_t.tsv (cdp, 11/2020)
        file_path_name = os.path.join(folder_path, f"{cond.identifier}_delta_t.tsv")
        discr_header = "delta_t (s)"
        # desinence to sictinguisch among space and time convergence (cdp, 12/2020)
        des = "tc"
        # the content of the round brackets in the file name (cdp, 12/2020)
        brackets = f"{discr}s"
        # convergence on mass and energy balance (cdp, 12/2020)
        AA = np.zeros((1, 3))
        AA[0, 0] = discr
        AA[0, 1] = cond.mass_balance
        AA[0, 2] = cond.energy_balance
        # discretization values for file CONDUCTOR_ID_delta_t.tsv
        val = np.array([discr])
    else:
        raise ValueError(f"Not valid option {space_conv}.\n")
    # end if option (cdp, 12/2020)

    os.makedirs(name=folder_path, exist_ok=True)

    # Create file file_path_name.tsv if don't exist; it stores the number of \
    # elements and discretization pitches used to perform the space convergence \
    # analysis or the time steps used to perform the time convergence analysis \
    # (cdp, 11/2020)
    if not os.path.exists(file_path_name):
        print(
            f"Created file {file_path_name} and wrote first two lines: headers and first row of data,\n"
        )
        if space_conv:
            # build header and array (cdp, 12/2020)
            header = "Nelems\t" + discr_header
        elif space_conv == False:
            # build header and array (cdp, 12/2020)
            header = discr_header
        # end if space_conv (cdp, 12/2020)
        # Write file file_path_name.tsv for the first time: headings \
        # and first row of data (cdp, 11/2020)
        with open(file_path_name, "w") as writer:
            np.savetxt(writer, val, header=header, comments="", delimiter="\t")
    else:
        print(f"Directory {file_path_name} already exists; appended row\n")
        # Append new row of data to file f_name (cdp, 11/2020)
        with open(file_path_name, "a") as writer:
            np.savetxt(writer, val, delimiter="\t")
    # end if (cdp, 11/2020)
    # save mass and energy balance at conductor level (cdp, 12/2020)
    file_path = os.path.join(folder_path, f"{cond.identifier}_mass_energy_{des}.tsv")
    if not os.path.exists(file_path):
        print(
            f"Created file {f_path} and wrote first two lines: headers and first \
    row of data,\n"
        )
        if space_conv:
            # build header for space convergence (cdp, 12/2020)
            mass_energy_header = (
                f"Nelems\t{discr_header}\tmass_bal (kg)\tenergy_bal (J)"
            )
        elif space_conv == False:
            # build header for time convergence (cdp, 12/2020)
            mass_energy_header = discr_header + "\tmass_bal (kg)\tenergy_bal (J)"
        # end if space_conv (cdp, 12/2020)
        # Write file f_name for the first time: headings and first row of data \
        # (cdp, 12/2020)
        with open(file_path, "w") as writer:
            np.savetxt(
                writer, AA, delimiter="\t", header=mass_energy_header, comments=""
            )
    else:
        print(f"Directory {f_path} already exists; appended row\n")
        # Append new row of data to file f_name (cdp, 09/2020)
        with open(file_path, "a") as writer:
            np.savetxt(writer, AA, delimiter="\t")
    # end if not (cdp, 12/2020)
    # Loop on FluidComponent (cdp, 12/2020)
    for fluid_comp in cond.inventory["FluidComponent"].collection:
        # save FluidComponent solution spatial distribution at TEND: velocity, \
        # pressure and temperature only (cdp, 11/2020)
        folder_path = os.path.join(f_path, cond.identifier, fluid_comp.identifier)
        os.makedirs(name=folder_path, exist_ok=True)
        file_path = os.path.join(
            folder_path, f"{fluid_comp.identifier}_({brackets}).tsv"
        )
        A_chan = np.zeros(
            (
                cond.grid_features["N_nod"],
                int(
                    cond.dict_N_equation["FluidComponent"]
                    / cond.inventory["FluidComponent"].number
                ),
            )
        )
        header_chan = "velocity (m/s)\tpressure (Pa)\ttemperature (K)"
        A_chan[:, 0] = fluid_comp.coolant.dict_node_pt["velocity"]
        A_chan[:, 1] = fluid_comp.coolant.dict_node_pt["pressure"]
        A_chan[:, 2] = fluid_comp.coolant.dict_node_pt["temperature"]
        with open(file_path, "w") as writer:
            np.savetxt(writer, A_chan, delimiter="\t", header=header_chan, comments="")
    # end for fluid_comp (cdp, 11/2020)
    # Loop on SolidComponent (cdp, 12/2020)
    for s_comp in cond.inventory["SolidComponent"].collection:
        # save SolidComponent solution spatial distribution at TEND: temperature \
        # only (cdp, 11/2020)
        folder_path = os.path.join(f_path, cond.identifier, s_comp.identifier)
        os.makedirs(name=folder_path, exist_ok=True)
        file_path = os.path.join(folder_path, f"{s_comp.identifier}_({brackets}).tsv")
        headers_s_comp = "temperature (K)"
        with open(file_path, "w") as writer:
            np.savetxt(
                writer,
                s_comp.dict_node_pt["temperature"],
                delimiter="\t",
                header=headers_s_comp,
                comments="",
            )
    # end for s_comp (cdp, 11/2020)


# end function Save_convergence_data (cdp, 12/2020)


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
