import enum
from sys import flags
from properties_of_materials.stainless_steel import density_ss
import numpy as np
import pandas as pd
import os


def save_properties(conductor, f_path):

    """Functions that save .tsv files with suitable file names, FluidComponents and SolidComponents initialization and final solution, together with spatial coordinate discretization. Channels saved variables are: temperature, pressure, density, viscosity, specific heat at constant pressure, thermal conductivity, velocity, Reynolds number and Prandtl number. Strands saved variables are: temperature, density, specific heat at constant pressure, thermal conductivity, magnetic field, electrical resistivity, current sharing temperature; jackets saved variables are temperature, specific heat at constant pressure, thermal conductivity, magnetic field, electrical resistivity."""

    list_prop_chan = list(
        conductor.dict_obj_inventory["FluidComponents"]["Objects"][
            0
        ].coolant.dict_node_pt.keys()
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
    header_chan = "xcoord (m)"
    for jj in range(len(list_prop_chan)):
        header_chan = f"{header_chan}	{list_prop_chan[jj]} {list_units[jj]}"
    header_st = "xcoord (m)	temperature (K)	B_field (T)	T_cur_sharing (K)"
    header_stab = "xcoord (m)	temperature (K)	B_field (T)"
    header_jk = "xcoord (m)	temperature (K)"
    for fluid_comp in conductor.dict_obj_inventory["FluidComponents"]["Objects"]:
        A_chan = np.zeros(
            (
                conductor.dict_discretization["N_nod"],
                len(fluid_comp.coolant.dict_node_pt) + 2,
            )
        )
        file_path = os.path.join(f_path, f"{fluid_comp.ID}.tsv")
        A_chan[:, 0] = conductor.dict_discretization["xcoord"]
        for ii, prop_value in enumerate(fluid_comp.coolant.dict_node_pt.values(), 1):
            A_chan[:, ii] = prop_value
        # Save total friction factor
        A_chan[:, -1] = fluid_comp.channel.dict_friction_factor[True]["total"]
        with open(file_path, "w") as writer:
            np.savetxt(writer, A_chan, delimiter="\t", header=header_chan, comments="")
    for strand in conductor.dict_obj_inventory["Strands"]["Objects"]:
        file_path = os.path.join(f_path, f"{strand.ID}.tsv")
        if strand.NAME != conductor.dict_obj_inventory["Stabilizer"]["Name"]:
            A_strand = np.zeros((conductor.dict_discretization["N_nod"], 4))
            A_strand[:, 3] = strand.dict_node_pt["T_cur_sharing"]
        else:
            A_strand = np.zeros((conductor.dict_discretization["N_nod"], 3))
        A_strand[:, 0] = conductor.dict_discretization["xcoord"]
        A_strand[:, 1] = strand.dict_node_pt["temperature"]
        A_strand[:, 2] = strand.dict_node_pt["B_field"]
        with open(file_path, "w") as writer:
            if strand.NAME != conductor.dict_obj_inventory["Stabilizer"]["Name"]:
                np.savetxt(
                    writer, A_strand, delimiter="\t", header=header_st, comments=""
                )
            else:
                np.savetxt(
                    writer, A_strand, delimiter="\t", header=header_stab, comments=""
                )
    for jacket in conductor.dict_obj_inventory["Jacket"]["Objects"]:
        file_path = os.path.join(f_path, f"{jacket.ID}.tsv")
        A_jacket = np.zeros((conductor.dict_discretization["N_nod"], 2))
        A_jacket[:, 0] = conductor.dict_discretization["xcoord"]
        A_jacket[:, 1] = jacket.dict_node_pt["temperature"]
        with open(file_path, "w") as writer:
            np.savetxt(writer, A_jacket, delimiter="\t", header=header_jk, comments="")


# end function Save_properties


def save_simulation_space(conductor, f_path, n_digit):

    """
    Function that save on files with suitable file names transient solution,
    spatial coordinate discretization, time and time step for each conductor
    (cdp, 08/2020)
    """

    # Save transient solution for each Conductors objects at each \
    # iteration_store iteration. For each conductor component, solution is \
    # stored in a dedicated file. FluidComponents stored variables are velocity, \
    # pressure and temperature together with spatial discretization, while for \
    # SolidComponents only temperature and spatial discretization are stored. \
    # (cdp, 08/2020)
    # Hypothesis: for different input files parameters there are different \
    # simulation names (cdp, 08/2020)

    # get the time at which the saving is made
    time = round(conductor.Space_save[conductor.i_save], n_digit)
    conductor.num_step_save[conductor.i_save] = conductor.cond_num_step
    # end if len(tt[0]) (cdp, 12/2020)
    prop_chan = [
        "xcoord",
        "velocity",
        "pressure",
        "temperature",
        "total_density",
        "friction_factor",
    ]
    header_chan = "xcoord (m)	velocity (m/s)	pressure (Pa)	temperature (K)	total_density (kg/m^3)\tfriction_factor (~)"
    for fluid_comp in conductor.dict_obj_inventory["FluidComponents"]["Objects"]:
        file_path = os.path.join(
            f_path, f"{fluid_comp.ID}_({conductor.cond_num_step})_sd.tsv"
        )
        A_chan = np.zeros((conductor.dict_discretization["N_nod"], len(prop_chan)))
        for ii in range(len(prop_chan)):
            if prop_chan[ii] == "xcoord":
                A_chan[:, ii] = conductor.dict_discretization[prop_chan[ii]]
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
    # headers_s_comp = "xcoord (m)	temperature (K)	density (kg/m^3)	spec_heat_p (J/kg/K)	ther_cond (W/m/K)	EXFTLX (W/m)	JHTFLX (W/m^2)"
    # prop_s_comp = ["xcoord", "temperature", "total_density", "total_isobaric_specific_heat", \
    # 							"total_thermal_conductivity", "EXTFLX", "JHTFLX"]
    headers_full = "xcoord (m)	temperature (K)\tcurrent_sharing_temperature (K)"
    headers_reduced = "xcoord (m)	temperature (K)"
    prop_full = ["xcoord", "temperature", "T_cur_sharing"]
    prop_reduced = ["xcoord", "temperature"]
    # Loop to save strand properties spatial distribution.
    for strand in conductor.dict_obj_inventory["Strands"]["Objects"]:
        file_path = os.path.join(
            f_path, f"{strand.ID}_({conductor.cond_num_step})_sd.tsv"
        )
        if strand.KIND != "Stabilizer":
            # Check if current sharing temperature is evaluated at each
            # thermal time step.
            if strand.dict_operation["TCS_EVALUATION"]:
                headers_strand = headers_full
                prop_strand = prop_full
            else:
                headers_strand = headers_reduced
                prop_strand = prop_reduced
        else: # Stabilizer
            headers_strand = headers_reduced
            prop_strand = prop_reduced
        
        A_strand = np.zeros((conductor.dict_discretization["N_nod"], len(prop_strand)))
        for ii in range(len(prop_strand)):
            if prop_strand[ii] == "xcoord":
                A_strand[:, ii] = conductor.dict_discretization[prop_strand[ii]]
            else:
                A_strand[:, ii] = strand.dict_node_pt[prop_strand[ii]]
            # end if prop_strand[ii] (cdp, 01/2021)
        # end for ii (cdp, 01/2021)
        with open(file_path, "w") as writer:
            np.savetxt(
                writer, A_strand, delimiter="\t", header=headers_strand, comments=""
            )
        
    headers_jk = "xcoord (m)	temperature (K)"
    prop_jk = ["xcoord", "temperature"]
    # Loop to save jacket properties spatial distribution.
    for jk in conductor.dict_obj_inventory["Jacket"]["Objects"]:
        file_path = os.path.join(
            f_path, f"{jk.ID}_({conductor.cond_num_step})_sd.tsv"
        )
        A_jk = np.zeros((conductor.dict_discretization["N_nod"], len(prop_jk)))
        for ii in range(len(prop_jk)):
            if prop_jk[ii] == "xcoord":
                A_jk[:, ii] = conductor.dict_discretization[prop_jk[ii]]
            else:
                A_jk[:, ii] = jk.dict_node_pt[prop_jk[ii]]
            # end if prop_s_comp[ii] (cdp, 01/2021)
        # end for ii (cdp, 01/2021)
        with open(file_path, "w") as writer:
            np.savetxt(
                writer, A_jk, delimiter="\t", header=headers_jk, comments=""
            )
    # end for s_comp (cdp, 10/2020)
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


def reorganize_spatial_distribution(cond, f_path, n_digit):
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
    # 							 "EXTFLX", "JHTFLX"]
    list_sol_key = ["temperature"]
    # lists all the file .tsv in subfolder Spatial_distribution (cdp, 11/2020)
    # Round the time to save to n_digit digits only once
    time = np.around(cond.Space_save, n_digit)
    # loop on FluidComponents (cdp, 11/2020)
    for fluid_comp in cond.dict_obj_inventory["FluidComponents"]["Objects"]:
        # create a list of files that have the fluid_comp.ID and User in the name \
        # exploiting list compreension: these files are the ones that will be \
        # reorganized by this function (cdp, 11/2020)
        # list_ch_file = [ff for ff in list_file if (fluid_comp.ID in ff and "User" in ff)]
        # declare the dictionary of data frame (cdp, 11/2020)
        dict_df = dict()
        dict_df_new = dict()
        if fluid_comp.ID == cond.dict_obj_inventory["FluidComponents"]["Objects"][0].ID:
            # declare dictionary to store the spatial diccretizations only once \
            # (cdp, 01/2021)
            dict_xcoord = dict()
        # end if fluid_comp (cdp, 01/2021)
        for ii in range(len(cond.Space_save)):
            file_name = f"{fluid_comp.ID}_({cond.num_step_save[ii]})_sd.tsv"
            file_load = os.path.join(f_path, file_name)
            # Load file file_name as data frame as a value of dictionary \
            # corresponding to key file_name (cdp, 11/2020)
            dict_df[file_name] = pd.read_csv(
                filepath_or_buffer=file_load, delimiter="\t"
            )
            # Delete the old file format.
            os.remove(file_load)
            if (
                fluid_comp.ID
                == cond.dict_obj_inventory["FluidComponents"]["Objects"][0].ID
            ):
                # store the spatial discretizations at each required time step in file \
                # xcoord.tsv only once (cdp,01/2021)
                dict_xcoord[f"time = {time[ii]} (s)"] = dict_df[file_name]["xcoord (m)"]
            # end if fluid_comp.ID (cdp, 01/2021)
            if ii == 0:
                # get columns names only the first time (cdp, 11/2020)
                header = list(dict_df[file_name].columns.values.tolist())
                for jj in range(len(list_ch_key)):
                    prop = list_ch_key[jj]
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
                for jj in range(len(list_ch_key)):
                    prop = list_ch_key[jj]
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
            file_name = f"{fluid_comp.ID}_{prop}_sd.tsv"
            # build path to save the file (cdp, 11/2020)
            path_save = os.path.join(f_path, file_name)
            # save the data frame, without the row index name (cdp, 11/2020)
            dict_df_new[prop].to_csv(path_save, sep="\t", index=False)
        # end for prop (cdp, 11/2020)
        if fluid_comp.ID == cond.dict_obj_inventory["FluidComponents"]["Objects"][0].ID:
            # convert the dictionary to a DataFrame (cdp,01/2021)
            df_xcoord = pd.DataFrame(dict_xcoord)
            # build file name (cdp, 01/2021)
            file_name = f"xcoord.tsv"
            path_save = os.path.join(f_path, file_name)
            # save the DataFrame as file xcoord.tsv
            df_xcoord.to_csv(path_save, sep="\t", index=False)
            # end if fluid_comp.ID (cdp, 01/2021)
    # end for fluid_comp (cdp, 11/2020)
    # loop on SolidComponents (cdp, 11/2020)
    for s_comp in cond.dict_obj_inventory["SolidComponents"]["Objects"]:
        # declare the dictionary of data frame (cdp, 11/2020)
        dict_df = dict()
        dict_df_new = dict()
        for ii in range(len(cond.Space_save)):
            file_name = f"{s_comp.ID}_({cond.num_step_save[ii]})_sd.tsv"
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
                for jj in range(len(list_sol_key)):
                    prop = list_sol_key[jj]
                    # decompose the data frame in four dataframes (cdp, 11/2020)
                    dict_df_new[prop] = dict_df[file_name].filter(
                        items=[header[jj + 1]]
                    )
                    # rename data frames columns (cdp, 11/2020)
                    dict_df_new[prop].rename(
                        columns={header[jj + 1]: f"time = {time[ii]} (s)"}, inplace=True
                    )
            else:
                for jj in range(len(list_sol_key)):
                    prop = list_sol_key[jj]
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
            # end if ii (cdp, 11/2020)
        # end for ii (cdp, 11/2020)
        # for loop to save the new data frame (cdp, 11/2020)
        for prop in list_sol_key:
            # build file name (cdp, 11/2020)
            file_name = f"{s_comp.ID}_{prop}_sd.tsv"
            # build path to save the file (cdp, 11/2020)
            path_save = os.path.join(f_path, file_name)
            # save the data frame, without the row index name (cdp, 11/2020)
            dict_df_new[prop].to_csv(path_save, sep="\t", index=False)
    # end for s_comp (cdp, 11/2020)

    # Manage files with heat exhanged between inner jackets by radiation.
    reorganize_heat_sd(cond, f_path, "Heat_rad_inner", "Heat_rad", n_digit)
    # Manage files with heat exhanged between outer conductor surface and environment by convection and/or radiation.
    reorganize_heat_sd(cond, f_path, "Heat_exch_env", "Heat_exch", n_digit)

    # Manage files with open heat transfer coefficients between fluid components.
    reorganize_heat_sd(cond, f_path, "HTC_ch_ch_o", "HTC_open", n_digit)
    # Manage files with close heat transfer coefficients between fluid components.
    reorganize_heat_sd(cond, f_path, "HTC_ch_ch_c", "HTC_close", n_digit)
    # Manage files with heat transfer coefficient between fluid and solid components.
    reorganize_heat_sd(cond, f_path, "HTC_ch_sol", "HTC", n_digit)
    # Manage files with conductive heat transfer coefficients between solid components.
    reorganize_heat_sd(cond, f_path, "HTC_sol_sol_cond", "HTC_cond", n_digit)
    # Manage files with radiative heat transfer coefficients between solid components.
    reorganize_heat_sd(cond, f_path, "HTC_sol_sol_rad", "HTC_rad", n_digit)


# end function Reorganize_spatial_distribution (cdp, 11/2020)


def reorganize_heat_sd(cond, f_path, radix_old, radix_new, n_digit):
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
    time = np.around(cond.Space_save, n_digit)
    for ii in range(len(cond.Space_save)):
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

    ind_xcoord = {
        f"xcoord = {conductor.Time_save[ii]} (m)": np.max(
            np.nonzero(
                conductor.dict_discretization["xcoord"] <= conductor.Time_save[ii]
            )
        )
        for ii in range(conductor.Time_save.size)
    }
    # construct file header only once (cdp, 08/2020)
    if simulation.num_step == 0:
        headers = ["time (s)"]
        headers.extend([str(key) for key in ind_xcoord.keys()])
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
        for f_comp in conductor.dict_obj_inventory["FluidComponents"]["Objects"]:
            # Loop on velocity, pressure, temperature and total density.
            for key, value in f_comp.coolant.time_evol.items():
                # Inizialize dictionary corresponding to key to a dictionary of empty lists for the first time.
                f_comp.coolant.time_evol[key] = initialize_dictionaty_te(
                    value, ind_xcoord
                )
                # Save the headings only ones.
                pd.DataFrame(columns=headers).to_csv(
                    os.path.join(
                        simulation.dict_path[
                            f"Output_Time_evolution_{conductor.ID}_dir"
                        ],
                        f"{f_comp.ID}_{key}_te.tsv",
                    ),
                    sep="\t",
                    index=False,
                    header=True,
                )
            # End for key.
            # Inizialize dictionary corresponding to key to a dictionary of empty lists for the first time.
            f_comp.channel.time_evol["friction_factor"] = initialize_dictionaty_te(
                f_comp.channel.time_evol["friction_factor"], ind_xcoord
            )
            # Save the headings only ones.
            pd.DataFrame(columns=headers).to_csv(
                os.path.join(
                    simulation.dict_path[f"Output_Time_evolution_{conductor.ID}_dir"],
                    f"{f_comp.ID}_friction_factor_te.tsv",
                ),
                sep="\t",
                index=False,
                header=True,
            )
            # Save the headings only ones.
            pd.DataFrame(columns=headers_inl_out).to_csv(
                os.path.join(
                    simulation.dict_path[f"Output_Time_evolution_{conductor.ID}_dir"],
                    f"{f_comp.ID}_inlet_outlet_te.tsv",
                ),
                sep="\t",
                index=False,
                header=True,
            )
        # End for f_comp.
        for s_comp in conductor.dict_obj_inventory["SolidComponents"]["Objects"]:
            # Loop on velocity, pressure, temperature and total density.
            for key, value in s_comp.time_evol.items():
                # Inizialize dictionary corresponding to key to a dictionary of empty lists for the first time.
                s_comp.time_evol[key] = initialize_dictionaty_te(value, ind_xcoord)
                # Save the headings only ones.
                pd.DataFrame(columns=headers).to_csv(
                    os.path.join(
                        simulation.dict_path[
                            f"Output_Time_evolution_{conductor.ID}_dir"
                        ],
                        f"{s_comp.ID}_{key}_te.tsv",
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

    # FluidComponents objects (cdp, 08/2020)
    for fluid_comp in conductor.dict_obj_inventory["FluidComponents"]["Objects"]:
        # Loop on velocity, pressure, temperature and total density.
        for key, value in fluid_comp.coolant.time_evol.items():
            # Update the contend of the dictionary of lists with propertiy values at selected xcoord and current time.
            fluid_comp.coolant.time_evol[key] = update_values(
                value, fluid_comp.coolant.dict_node_pt[key], time, ind_xcoord
            )
            # Write the content of the dictionary to file, if conditions are satisfied.
            fluid_comp.coolant.time_evol[key] = save_te_on_file(
                conductor,
                fluid_comp.coolant.time_evol[key],
                os.path.join(
                    simulation.dict_path[f"Output_Time_evolution_{conductor.ID}_dir"],
                    f"{fluid_comp.ID}_{key}_te.tsv",
                ),
                simulation.transient_input["TEND"],
                ind_xcoord,
            )
        # End for key.

        # Save friction factor time evolution.
        # Update the contend of the dictionary of lists with propertiy values at selected xcoord and current time.
        fluid_comp.channel.time_evol["friction_factor"] = update_values(
            fluid_comp.channel.time_evol["friction_factor"],
            fluid_comp.channel.dict_friction_factor[True]["total"],
            time,
            ind_xcoord,
        )
        # Write the content of the dictionary to file, if conditions are satisfied.
        fluid_comp.channel.time_evol["friction_factor"] = save_te_on_file(
            conductor,
            fluid_comp.channel.time_evol["friction_factor"],
            os.path.join(
                simulation.dict_path[f"Output_Time_evolution_{conductor.ID}_dir"],
                f"{fluid_comp.ID}_friction_factor_te.tsv",
            ),
            simulation.transient_input["TEND"],
            ind_xcoord,
        )

        if fluid_comp.channel.flow_dir[0] == "forward":
            index_inl = 0
            index_out = -1
        elif fluid_comp.channel.flow_dir[0] == "backward":
            index_inl = -1
            index_out = 0

        # Inlet and outlet quantities (cdp, 08/2020)
        file_name_io = os.path.join(
            simulation.dict_path[f"Output_Time_evolution_{conductor.ID}_dir"],
            f"{fluid_comp.ID}_inlet_outlet_te.tsv",
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

    # SolidComponents objects (cdp, 08/2020)
    for s_comp in conductor.dict_obj_inventory["SolidComponents"]["Objects"]:
        for key, value in s_comp.time_evol.items():
            # Update the contend of the dictionary of lists with propertiy values at selected xcoord and current time.
            s_comp.time_evol[key] = update_values(
                value, s_comp.dict_node_pt[key], time, ind_xcoord
            )
            # Write the content of the dictionary to file, if conditions are satisfied.
            s_comp.time_evol[key] = save_te_on_file(
                conductor,
                s_comp.time_evol[key],
                os.path.join(
                    simulation.dict_path[f"Output_Time_evolution_{conductor.ID}_dir"],
                    f"{s_comp.ID}_{key}_te.tsv",
                ),
                simulation.transient_input["TEND"],
                ind_xcoord,
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
                simulation.dict_path[f"Output_Time_evolution_{conductor.ID}_dir"],
                "Time.tsv",
            ),
            sep="\t",
            header=True,
            index=False,
        )
    # End if abs.


# end function Save_simulation_time (cdp, 08/2020)


def initialize_dictionaty_te(val, ind_xcoord):

    val = {"time (s)": list()}
    val.update({key: list() for key in ind_xcoord.keys()})
    return val


# End function initialize_dictionaty_te.


def update_values(val, prop, time, ind_xcoord):

    val["time (s)"].append(time)
    # Use dict.update to avoid error (do not understood why with val.update does not work).
    dict.update(
        {
            key: value.append(prop[ind_xcoord[key]])
            for key, value in val.items()
            if "xcoord" in key
        }
    )
    return val


# End function update_values.


def save_te_on_file(conductor, val, file_name, tend, ind_xcoord):
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
        val = initialize_dictionaty_te(val, ind_xcoord)
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


def save_convergence_data(cond, f_path, *n_digit, space_conv=True):

    """
    Function that saves data for the space convergence: saved information is the solution spatial distribution at TEND for the given number of elements, together with the global mass and energy balance on channels.
    The spatial distributions are stored in a folder whose name is given by the combination of TEND and STPMIN, files name is comp.ID_(NELEMS).tsv.
    Mass and energy balance for every NELEMS are collected in file cond.ID_mass_energy_sc.tsv in each channel folder.
    """

    if space_conv:
        # Save data for the Space convergence analysis (cdp, 12/2020)
        # compute spatial discretization pitch (cdp, 12/2020)
        discr = (
            cond.dict_input["XLENGTH"]
            / cond.dict_discretization["Grid_input"]["NELEMS"]
        )
        folder_path = os.path.join(f_path, cond.ID)
        # Create the path of the file {cond.ID}_delta_x.tsv (cdp, 11/2020)
        file_path_name = os.path.join(folder_path, f"{cond.ID}_delta_x.tsv")
        discr_header = "delta_x (m)"
        # desinence to sictinguisch among space and time convergence (cdp, 12/2020)
        des = "sc"
        # the content of the round brackets in the file name (cdp, 12/2020)
        brackets = cond.dict_discretization["Grid_input"]["NELEMS"]
        # convergence on mass and energy balance (cdp, 12/2020)
        AA = np.zeros((1, 4))
        AA[0, 0] = cond.dict_discretization["Grid_input"]["NELEMS"]
        AA[0, 1] = discr
        AA[0, 2] = cond.mass_balance
        AA[0, 3] = cond.energy_balance
        # discretization values for file CONDUCTOR_ID_delta_x.tsv
        val = np.zeros((1, 2))
        val[0, 0] = cond.dict_discretization["Grid_input"]["NELEMS"]
        val[0, 1] = discr
    elif space_conv == False:
        # Save data for the Time convergence analysis (cdp, 12/2020)
        # compute time step, remember that n_digit is a tuple (cdp, 12/2020)
        discr = round(cond.time_step, n_digit[0])
        folder_path = os.path.join(f_path, cond.ID)
        # Create the path of the file {cond.ID}_delta_t.tsv (cdp, 11/2020)
        file_path_name = os.path.join(folder_path, f"{cond.ID}_delta_t.tsv")
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
            header = "Nelems	" + discr_header
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
    file_path = os.path.join(folder_path, f"{cond.ID}_mass_energy_{des}.tsv")
    if not os.path.exists(file_path):
        print(
            f"Created file {f_path} and wrote first two lines: headers and first \
    row of data,\n"
        )
        if space_conv:
            # build header for space convergence (cdp, 12/2020)
            mass_energy_header = f"Nelems	{discr_header}	mass_bal (kg)	energy_bal (J)"
        elif space_conv == False:
            # build header for time convergence (cdp, 12/2020)
            mass_energy_header = discr_header + "	mass_bal (kg)	energy_bal (J)"
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
    # Loop on FluidComponents (cdp, 12/2020)
    for fluid_comp in cond.dict_obj_inventory["FluidComponents"]["Objects"]:
        # save FluidComponents solution spatial distribution at TEND: velocity, \
        # pressure and temperature only (cdp, 11/2020)
        folder_path = os.path.join(f_path, cond.ID, fluid_comp.ID)
        os.makedirs(name=folder_path, exist_ok=True)
        file_path = os.path.join(folder_path, f"{fluid_comp.ID}_({brackets}).tsv")
        A_chan = np.zeros(
            (
                cond.dict_discretization["N_nod"],
                int(
                    cond.dict_N_equation["FluidComponents"]
                    / cond.dict_obj_inventory["FluidComponents"]["Number"]
                ),
            )
        )
        header_chan = "velocity (m/s)	pressure (Pa)	temperature (K)"
        A_chan[:, 0] = fluid_comp.coolant.dict_node_pt["velocity"]
        A_chan[:, 1] = fluid_comp.coolant.dict_node_pt["pressure"]
        A_chan[:, 2] = fluid_comp.coolant.dict_node_pt["temperature"]
        with open(file_path, "w") as writer:
            np.savetxt(writer, A_chan, delimiter="\t", header=header_chan, comments="")
    # end for fluid_comp (cdp, 11/2020)
    # Loop on SolidComponents (cdp, 12/2020)
    for s_comp in cond.dict_obj_inventory["SolidComponents"]["Objects"]:
        # save SolidComponents solution spatial distribution at TEND: temperature \
        # only (cdp, 11/2020)
        folder_path = os.path.join(f_path, cond.ID, s_comp.ID)
        os.makedirs(name=folder_path, exist_ok=True)
        file_path = os.path.join(folder_path, f"{s_comp.ID}_({brackets}).tsv")
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
