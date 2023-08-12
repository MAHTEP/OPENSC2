import matplotlib

matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

import os
import numpy as np
import pandas as pd


def plot_properties(simulation, cond, what="initialization"):

    """
    Function that makes singles plots of the properties of initialization or of the final solution. (cdp, 08/2020)
    """

    flag_chan = False
    if what == "initialization":
        dict_path = dict(
            Load=simulation.dict_path[f"Output_Initialization_{cond.identifier}_dir"],
            Save=simulation.dict_path[f"Figures_Initialization_{cond.identifier}_dir"],
        )
    elif what == "solution":
        dict_path = dict(
            Load=simulation.dict_path[f"Output_Solution_{cond.identifier}_dir"],
            Save=simulation.dict_path[f"Figures_Solution_{cond.identifier}_dir"],
        )
    # end if what (cdp, 12/2020)
    # Loop on FluidComponent (cdp, 12/2020)
    for fluid_comp in cond.inventory["FluidComponent"].collection:
        # load values
        file_load = os.path.join(dict_path["Load"], f"{fluid_comp.identifier}.tsv")
        # Load data in file_load as pandas DataFrame
        load_chan = pd.read_csv(filepath_or_buffer=file_load, delimiter="\t")
        folder_save = os.path.join(dict_path["Save"], fluid_comp.identifier)
        # Create target Directory if do not exist (cdp, 07/2020)
        if not os.path.exists(folder_save):
            os.makedirs(folder_save)
            print(f"Created directory {folder_save}\n")
        else:
            print(f"Directory {folder_save} already exists\n")
        # end if (cdp, 07/2020)
        if flag_chan == False:
            # get the list of saved properties (cdp, 12/2020)
            prop_chan = list(load_chan.columns.values.tolist())
            flag_chan = True
        # end if flag_chan (cdp, 12/2020)
        for p_name in prop_chan[1:]:
            ind = p_name.find(" ")
            y_name, units = plot_nomenclature_y(p_name[:ind])
            # path to save figures (cdp, 07/2020)
            save_fig_path_svg = os.path.join(folder_save, f"{p_name[:ind]}.svg")
            # save_fig_path_pdf = os.path.join(folder_save, f"{p_name[:ind]}.pdf")
            Figure, ax = plt.subplots(
                num=(
                    simulation.transient_input["SIMULATION"]
                    + f" {cond.identifier}, {fluid_comp.identifier}:"
                    f" {p_name[:ind]}"
                ),
                figsize=(7.0, 6.0),
            )
            ax.plot(load_chan[prop_chan[0]], load_chan[p_name], "k-", linewidth=2.0)
            ax.grid(True)
            ax.set(
                xlabel="$x\ (m)$",
                ylabel=f"{y_name} {units}",
                title=f"{cond.identifier}, {fluid_comp.identifier}: {p_name[:ind]}",
            )
            plt.savefig(save_fig_path_svg, orientation="portrait", transparent=False)
            # plt.savefig(save_fig_path_pdf, orientation = "portrait", \
            # 																										transparent = False)
            plt.close()
    # end for fluid_comp (cdp, 12/2020)
    # Loop on SolidComponent (cdp, 12/2020)
    for s_comp in cond.inventory["SolidComponent"].collection:
        # load values
        file_load = os.path.join(dict_path["Load"], f"{s_comp.identifier}.tsv")
        # Load data in file_load as pandas DataFrame
        load_s_comp = pd.read_csv(filepath_or_buffer=file_load, delimiter="\t")
        folder_save = os.path.join(dict_path["Save"], s_comp.identifier)
        # Create target Directory if do not exist (cdp, 07/2020)
        if not os.path.exists(folder_save):
            os.makedirs(folder_save)
            print(f"Created directory {folder_save}\n")
        else:
            print(f"Directory {folder_save} already exists\n")
        # end if (cdp, 07/2020)
        # get the list of saved properties (cdp, 12/2020)
        prop_s_comp = list(load_s_comp.columns.values.tolist())
        # end if flag_s_comp (cdp, 12/2020)
        for p_name in prop_s_comp[1:]:
            ind = p_name.find(" ")
            y_name, units = plot_nomenclature_y(p_name[:ind])
            # path to save figures (cdp, 07/2020)
            save_fig_path_svg = os.path.join(folder_save, f"{p_name[:ind]}.svg")
            # save_fig_path_pdf = os.path.join(folder_save, f"{p_name[:ind]}.pdf")
            Figure, ax = plt.subplots(
                num=(
                    simulation.transient_input["SIMULATION"]
                    + f" {cond.identifier}, {s_comp.identifier}:"
                    f" {p_name[:ind]}"
                ),
                figsize=(7.0, 6.0),
            )
            ax.plot(
                load_s_comp[prop_s_comp[0]], load_s_comp[p_name], "k-", linewidth=2.0
            )
            ax.grid(True)
            ax.set(
                xlabel="$x\ (m)$",
                ylabel=f"{y_name} {units}",
                title=f"{cond.identifier}, {s_comp.identifier}: {p_name[:ind]}",
            )
            plt.savefig(save_fig_path_svg, orientation="portrait", transparent=False)
            # plt.savefig(save_fig_path_pdf, orientation = "portrait", \
            # 																								transparent = False)
            plt.close()
    # end for s_comp (cdp, 12/2020)
    print("Plotted " + what + "\n")


# end function Plot_properties (cdp, 08/2020)


def plot_nomenclature_y(prop):

    """
    Function that gives the correct nomenclature to y axes in plots (cdp, 08/2020)
    """

    dict_nomenclature = dict(
        temperature=("$T$", "$(K)$"),
        pressure=("$P$", "$(Pa)$"),
        total_density=(r"$\rho$", "$(kg/m^3)$"),
        total_dynamic_viscosity=(r"$\mu$", "$(Pa*s)$"),
        total_isobaric_specific_heat=("$cp$", "$(J/kg/K)$"),
        total_isochoric_specific_heat=("$cv$", "$(J/kg/K)$"),
        total_thermal_conductivity=("$k$", "$(W/m/K)$"),
        velocity=("$v$", "$(m/s)$"),
        total_enthalpy=("$w$", "$(J/kg)$"),
        total_entropy=("$s$", "$(J/kg/K)$"),
        total_speed_of_sound=("$c$", "$(m/s)$"),
        Gruneisen=(r"$\Phi$", "$(~)$"),
        Reynolds=("$Re$", "$(~)$"),
        Prandtl=("$Pr$", "$(~)$"),
        B_field=("$B$", "$(T)$"),
        T_cur_sharing=("$T_csh$", "$(K)$"),
        mass_flow_rate=("$mdot$", "$(kg/s)$"),
        Heat_rad=("Q_rad", "(W/m)"),
        Heat_exch=("Q_exch", "(W/m)"),
        isobaric_expansion_coefficient=(r"$\alpha_L$", "(1/K)"),
        isothermal_compressibility=(r"$\beta$", "(1/Pa)"),
        friction_factor=("f", "(~)"),
    )
    return dict_nomenclature[prop]


# end function Plot_nomenclature_y (cdp, 11/2020)


def make_plots(simulation, kind="Space_distr"):
    """
    Function that makes plots of spatial discretization and of time evolution according to the **option argument. (cdp, 11/2020)
    """
    # list of colors (cdp, 11/2020)
    colors = ["k", "r", "b", "g", "c"]
    # maximum number of axes for each figure 4 (cdp, 11/2020)
    N_ax_max = int(4)
    # define the number of lines in each plot (min: 1, max: 5, default: 5) \
    # (cdp, 10/2020)
    N_lines_max = int(5)
    # compute the maximum number of lines for each figure (cdp, 01/2021)
    N_lines_tot_max = N_ax_max * N_lines_max
    # Decide what kind of plot make: spatial distribution or time evolution \
    # (cdp, 11/2020)
    if kind == "Space_distr":
        # make plot of spatial distribution (cdp, 11/2020)
        # list of dictionary keys (cdp, 09/2020)
        prop_chan = ["velocity", "pressure", "temperature", "total_density"]
        # Dictionary to easily distinguish the outcome to be saved according to the name of the strand.
        prop_st = dict(
            STR_MIX=["temperature"], STACK=["temperature"], STR_STAB=["temperature"]
        )
        prop_jk = ["temperature"]
    elif kind == "Time_evol":
        # make plot of time evolution (cdp, 11/2020)
        # list of dictionary keys (cdp, 09/2020)
        prop_chan = [
            "velocity",
            "pressure",
            "temperature",
            "total_density",
            "inlet_outlet",
        ]
        # Dictionary to easily distinguish the outcome to be saved according to the name of the strand.
        prop_st = dict(
            STR_MIX=["B_field", "temperature", "T_cur_sharing"],
            STACK=["B_field", "temperature", "T_cur_sharing"],
            STR_STAB=["B_field", "temperature"],
        )
        prop_jk = ["temperature"]
    # end if kind (cdp, 11/2020)
    dict_values = {}  # dictionary declaration (cdp, 09/2020)
    # loop on conductorts objects (cdp,11/2020)
    for cond in simulation.list_of_Conductors:
        dict_values[cond.identifier] = {}  # dictionary declaration (cdp, 09/2020)
        if kind == "Space_distr":
            # unify the array (cdp, 11/2020)
            kind_save = np.around(cond.Space_save, simulation.n_digit_time)
            root_load_path = simulation.dict_path[
                f"Output_Spatial_distribution_{cond.identifier}_dir"
            ]
            des = "sd"
            root_save_path = simulation.dict_path[
                f"Figures_Spatial_distribution_{cond.identifier}_dir"
            ]
            abscissa = pd.read_csv(
                os.path.join(root_load_path, "zcoord.tsv"), delimiter="\t"
            )
            # plot features (cdp, 01/2021)
            # unify the array (cdp, 11/2020)
            # uni_label = cond.Space_save
            # x axis label (cdp, 11/2020)
            ascissa = "$x\ (m)$"
            # legend title (cdp, 11/2020)
            leg_title = "$t\ (s)$"
            # figure title complemen (cdp, 11/2020)
            title_comp = "s. d."
        elif kind == "Time_evol":
            # unify the array (cdp, 11/2020)
            kind_save = cond.Time_save
            root_load_path = simulation.dict_path[
                f"Output_Time_evolution_{cond.identifier}_dir"
            ]
            des = "te"
            root_save_path = simulation.dict_path[
                f"Figures_Time_evolution_{cond.identifier}_dir"
            ]
            # plot features (cdp, 01/2021)
            # unify the array (cdp, 11/2020)
            # uni_label = cond.Time_save
            # x axis label (cdp, 11/2020)
            ascissa = "$t\ (s)$"
            # legend title (cdp, 11/2020)
            leg_title = "$x\ (m)$"
            title_comp = "t. e."
        # end if kind
        # Evaluate the number of figures to make
        N_figure = int(np.ceil(len(kind_save) / N_lines_tot_max))
        # Get the number of axes (subplots) in each Figure; each axes has no more \
        # than N_lines_max lines (cdp, 10/2020)
        # initialize the array with the number of axes for each figure \
        # (cdp, 11/2020)
        N_axes = np.zeros(N_figure, dtype=int)
        # initialize the actual number of lines for each axes (cdp, 01/2021)
        N_lines = N_lines_max * np.ones(N_figure, dtype=int)
        # initialize the total number of lines for each figure, max \
        # N_lines_tot_max (cdp, 01/2021)
        N_lines_tot = np.zeros(N_figure, dtype=int)
        for nn in range(N_figure):
            # evaluate the number of curves to plot in each figure
            N_lines_tot[nn] = np.minimum(
                len(kind_save) - sum(N_lines_tot[0:nn]), N_lines_tot_max
            )
            # evaluate the number of axes for each figure (cdp, 11/2020)
            N_axes[nn] = np.minimum(
                int(np.ceil(N_lines_tot[nn] / N_lines_max)), N_ax_max
            )
            N_lines[nn] = round(N_lines_tot[nn] / N_axes[nn])
        # end for nn (cdp, 11/2020)
        # Loop on FluidComponent (cdp, 11/2020)
        for fluid_comp in cond.inventory["FluidComponent"].collection:
            # dictionary declaration (cdp, 09/2020)
            dict_values[cond.identifier][fluid_comp.identifier] = {}
            # Load properties value for channels (cdp, 09/2020)
            for prop in prop_chan:
                # load file (cdp, 09/2020)
                file_load = os.path.join(
                    root_load_path, f"{fluid_comp.identifier}_{prop}_{des}.tsv"
                )
                # dictionary declaration (cdp, 09/2020)
                dict_values[cond.identifier][fluid_comp.identifier][prop] = pd.read_csv(
                    filepath_or_buffer=file_load, delimiter="\t"
                )
                folder_save = os.path.join(root_save_path, fluid_comp.identifier)
                # Create target Directory if do not exist (cdp, 09/2020)
                if not os.path.exists(folder_save):
                    os.makedirs(folder_save)
                    print(f"Created directory {folder_save}\n")
                else:
                    print(f"Directory {folder_save} already exists\n")

                # Construct sup title for all the axes of the figure.
                sup_title = (
                    f"{cond.identifier} {fluid_comp.identifier} {prop}: {title_comp}"
                )
                if kind == "Space_distr":
                    # call function Make_plots_sd_actually to make plots of spatial \
                    # distributions (cpd, 10/2020)
                    make_plots_sd_actually(
                        abscissa,
                        dict_values[cond.identifier][fluid_comp.identifier][prop],
                        prop,
                        folder_save,
                        N_lines_tot_max,
                        N_lines_tot,
                        N_lines,
                        N_axes,
                        colors,
                        kind_save,
                        ascissa,
                        sup_title,
                        leg_title,
                        figures=N_figure,
                    )
                elif kind == "Time_evol":
                    # call function Make_plots_te_actually to make plots of time \
                    # evolutions (cpd, 10/2020)
                    if prop == "inlet_outlet":
                        # mass flow rate time evolution is plotted only at inlet and \
                        # outlet coordinates so a single figure with a single axes is \
                        # sufficient in this case (cdp, 10/2020)
                        # Max total number of lines (cdp, 01/2021)
                        N_l_t_m = 2 * np.ones(1, dtype=int)
                        # Total number of lines (cdp, 01/2021)
                        N_l_t = 2 * np.ones(1, dtype=int)
                        # Number of axes (cdp, 01/2021)
                        N_ax = np.ones(1, dtype=int)
                        make_plots_te_actually(
                            dict_values[cond.identifier][fluid_comp.identifier][prop],
                            prop,
                            folder_save,
                            N_l_t_m,
                            N_l_t,
                            N_lines,
                            N_ax,
                            colors,
                            kind_save,
                            ascissa,
                            sup_title,
                            leg_title,
                        )
                    else:
                        make_plots_te_actually(
                            dict_values[cond.identifier][fluid_comp.identifier][prop],
                            prop,
                            folder_save,
                            N_lines_tot_max,
                            N_lines_tot,
                            N_lines,
                            N_axes,
                            colors,
                            kind_save,
                            ascissa,
                            sup_title,
                            leg_title,
                            figures=N_figure,
                        )
                    # end if prop (cdp, 01/2021)
                # end if kind (cdp, 01/2021)
            # end for prop (cdp, 10/2020)
        # end for fluid_comp (cdp, 10/2020)
        # Loop on SolidComponent (cdp, 11/2020)
        for s_comp in cond.inventory["SolidComponent"].collection:
            dict_values[cond.identifier][s_comp.identifier] = {}
            if s_comp.name != cond.inventory["JacketComponent"].name:
                # StrandComponent objects (cdp, 09/2020)
                # Load properties value for strands (cdp, 09/2020)
                prop_s_comp = prop_st[s_comp.name]
            else:
                # JacketComponents objects (cdp, 09/2020)
                # Load properties value for jackets (cdp, 09/2020)
                prop_s_comp = prop_jk
            for prop in prop_s_comp:
                # load file (cdp, 09/2020)
                file_load = os.path.join(
                    root_load_path, f"{s_comp.identifier}_{prop}_{des}.tsv"
                )
                # dictionary declaration (cdp, 09/2020)
                dict_values[cond.identifier][s_comp.identifier][prop] = pd.read_csv(
                    filepath_or_buffer=file_load, delimiter="\t"
                )
                folder_save = os.path.join(root_save_path, s_comp.identifier)
                # Create target Directory if do not exist (cdp, 09/2020)
                if not os.path.exists(folder_save):
                    os.makedirs(folder_save)
                    print(f"Created directory {folder_save}\n")
                else:
                    print(f"Directory {folder_save} already exists\n")
                # Construct sup title for all the axes of the figure.
                sup_title = (
                    f"{cond.identifier} {s_comp.identifier} {prop}: {title_comp}"
                )
                if kind == "Space_distr":
                    # call function Make_plots_sd_actually to make plots of spatial \
                    # distributions (cpd, 10/2020)
                    make_plots_sd_actually(
                        abscissa,
                        dict_values[cond.identifier][s_comp.identifier][prop],
                        prop,
                        folder_save,
                        N_lines_tot_max,
                        N_lines_tot,
                        N_lines,
                        N_axes,
                        colors,
                        kind_save,
                        ascissa,
                        sup_title,
                        leg_title,
                        figures=N_figure,
                    )
                elif kind == "Time_evol":
                    # call function Make_plots_te_actually to make plots of time \
                    # evolutions (cpd, 10/2020)
                    make_plots_te_actually(
                        dict_values[cond.identifier][s_comp.identifier][prop],
                        prop,
                        folder_save,
                        N_lines_tot_max,
                        N_lines_tot,
                        N_lines,
                        N_axes,
                        colors,
                        kind_save,
                        ascissa,
                        sup_title,
                        leg_title,
                        figures=N_figure,
                    )
                # end if kind (cdp, 01/2021)
            # end for prop (cdp, 11/2020)
        # end for s_comp (cdp, 11/2020)

        if kind == "Space_distr":
            # Evaluate the Gauss coordinate to make the plot, conversion to numpy necessaty otherwise nan values arise.
            x_gauss = pd.DataFrame(
                {
                    abscissa.columns[ii]: (
                        abscissa.iloc[:-1, ii].to_numpy()
                        + abscissa.iloc[1:, ii].to_numpy()
                    )
                    / 2.0
                    for ii in range(len(abscissa.columns))
                }
            )
            for rr, jk_r in enumerate(cond.inventory["SolidComponent"].collection):
                for _, jk_c in enumerate(
                    cond.inventory["SolidComponent"].collection[rr + 1 :]
                ):
                    if (
                        abs(
                            cond.dict_df_coupling["HTC_choice"].at[
                                jk_r.identifier, jk_c.identifier
                            ]
                        )
                        == 3
                    ):
                        prop = f"Heat_rad_{jk_r.identifier}_{jk_c.identifier}"
                        # Build file path.
                        file_load = os.path.join(root_load_path, f"{prop}_{des}.tsv")
                        # Load file as dataframe.
                        values = pd.read_csv(file_load, delimiter="\t")
                        folder_save = os.path.join(root_save_path)
                        sup_title = f"{cond.identifier} {jk_r.identifier} {jk_c.identifier} Heat rad: {title_comp}"
                        # Plot the heat exchanged by radiation between jackets.
                        make_plots_sd_actually(
                            x_gauss,
                            values,
                            prop,
                            folder_save,
                            N_lines_tot_max,
                            N_lines_tot,
                            N_lines,
                            N_axes,
                            colors,
                            kind_save,
                            ascissa,
                            sup_title,
                            leg_title,
                            figures=N_figure,
                        )
                    # End if abs().
                # End for cc.
            if (
                abs(cond.dict_df_coupling["contact_perimeter_flag"].at[
                    simulation.environment.KIND, jk_r.identifier
                ]
                ) == 1
            ):
                prop = f"Heat_exch_{simulation.environment.KIND}_{jk_r.identifier}"
                # Build file path.
                file_load = os.path.join(root_load_path, f"{prop}_{des}.tsv")
                # Load file as dataframe.
                values = pd.read_csv(file_load, delimiter="\t")
                folder_save = os.path.join(root_save_path)
                sup_title = f"{cond.identifier} {simulation.environment.KIND} {jk_r.identifier} Heat exch: {title_comp}"
                # Plot the heat exchanged by radiation between jackets.
                make_plots_sd_actually(
                    x_gauss,
                    values,
                    prop,
                    folder_save,
                    N_lines_tot_max,
                    N_lines_tot,
                    N_lines,
                    N_axes,
                    colors,
                    kind_save,
                    ascissa,
                    sup_title,
                    leg_title,
                    figures=N_figure,
                )
            # End if cond.dict_df_coupling["contact_perimeter_flag"].
        # End for rr.

    # end for cond (cdp, 11/2020)


# end function Make_plots (cdp, 11/2020)


def make_plots_sd_actually(
    xvalues,
    yvalues,
    p_name,
    folder_save,
    N_lines_tot_max,
    N_lines_tot,
    N_lines,
    N_axes,
    colors,
    tstep_full,
    ascissa,
    sup_title,
    leg_title,
    figures=1,
):
    """
    Function that actually makes plots of the spatial distribution. (cdp, 11/2020)
    Updated (cdp, 01/2021)
    """
    header_full = list(xvalues.columns.values.tolist())
    # end if p_name (cdp, 10/2020)
    # Loop on figures (cdp, 11/2020)
    for nn in range(figures):
        # Path to save figures (cdp, 10/2020)
        if figures > 1:
            dict_save_fig_path = dict(
                eps=os.path.join(folder_save, f"{p_name}_{nn+1}.svg"),
                pdf=os.path.join(folder_save, f"{p_name}_{nn+1}.pdf"),
            )
            # extract the header from header_full to correctly build the lengend and \
            # do not exceed the number of lines for each figures (cdp, 01/2021)
            if nn < figures - 1:
                header = header_full[nn * N_lines_tot_max : (nn + 1) * N_lines_tot_max]
                tstep = tstep_full[nn * N_lines_tot_max : (nn + 1) * N_lines_tot_max]
            else:
                header = header_full[nn * N_lines_tot_max :]
                tstep = tstep_full[nn * N_lines_tot_max :]
            # end if nn (cdp, 01/2021)
        else:
            dict_save_fig_path = dict(
                eps=os.path.join(folder_save, f"{p_name}.svg"),
                pdf=os.path.join(folder_save, f"{p_name}.pdf"),
            )
            header = header_full
            tstep = tstep_full
        # end if figures (cdp, 11/2020)
        if N_axes[nn] == 1:
            # one subplot -> 1 row, 1 column (cdp, 10/2020)
            Figure, ax = plt.subplots(nrows=N_axes[nn], ncols=1, figsize=(7.0, 6.0))
            ii = 0
            for head in header:
                ax.plot(
                    xvalues[head],
                    yvalues[head],
                    colors[ii],
                    linewidth=2,
                    label=tstep[ii],
                )
                ii += 1
            # end for head (cdp, 11/2020)
            # end if p_name (cdp, 11/2020)
            # Call fuction Add_figure_features to coplete and save the figure \
            # (cdp, 01/2021)
            ax = [ax]  # make list to avoid errors (cdp, 01/2021)
            add_figure_features(
                Figure,
                ax,
                p_name,
                ascissa,
                sup_title,
                leg_title,
                dict_save_fig_path["eps"],
            )
        elif N_axes[nn] == 2:
            # 2 subplot -> 2 row, 1 column (cdp, 10/2020)
            Figure, ax = plt.subplots(nrows=2, ncols=1, figsize=(7.0, 6.0), sharex=True)
            for ii in range(N_lines[nn]):
                ax[0].plot(
                    xvalues[header[ii]],
                    yvalues[header[ii]],
                    colors[ii],
                    linewidth=2,
                    label=tstep[ii],
                )
                if ii + N_lines[nn] < N_lines_tot[nn]:
                    ax[1].plot(
                        xvalues[header[ii]],
                        yvalues[header[ii + N_lines[nn]]],
                        colors[ii],
                        linewidth=2,
                        label=tstep[ii + N_lines[nn]],
                    )
                # end if ii + N_lines[nn] (cdp, 01/2021)
            # end for ii (cdp, 01/2021)
            # Call fuction Add_figure_features to coplete and save the figure \
            # (cdp, 01/2021)
            add_figure_features(
                Figure,
                ax,
                p_name,
                ascissa,
                sup_title,
                leg_title,
                dict_save_fig_path["eps"],
            )
        elif N_axes[nn] == 3:
            # 3 subplot -> 3 row, 1 column (cdp, 10/2020)
            Figure, ax = plt.subplots(
                nrows=3, ncols=1, figsize=(10.0, 8.5), sharex=True
            )
            # loop on saved values (cdp, 10/2020)
            for ii in range(N_lines[nn]):
                ax[0].plot(
                    xvalues[header[ii]],
                    yvalues[header[ii]],
                    colors[ii],
                    linewidth=2,
                    label=tstep[ii],
                )
                ax[1].plot(
                    xvalues[header[ii + N_lines[nn]]],
                    yvalues[header[ii + N_lines[nn]]],
                    colors[ii],
                    linewidth=2,
                    label=tstep[ii + N_lines[nn]],
                )
                if ii + 2 * N_lines[nn] < N_lines_tot[nn]:
                    ax[2].plot(
                        xvalues[header[ii + 2 * N_lines[nn]]],
                        yvalues[header[ii + 2 * N_lines[nn]]],
                        colors[ii],
                        linewidth=2,
                        label=tstep[ii + 2 * N_lines[nn]],
                    )
                # end if ii + 2*N_lines[nn] (cdp, 01/2021)
            # end for ii (cdp, 01/2021)
            if N_lines_tot[nn] == 13:
                # special case (cdp, 01/2021)
                ax[2].plot(
                    xvalues[header[-1]],
                    yvalues[header[-1]],
                    colors[-1],
                    linewidth=2,
                    label=tstep[-1],
                )
            # end if N_lines_tot (cdp, 01/2021)
            # Call fuction Add_figure_features to coplete and save the figure \
            # (cdp, 01/2021)
            add_figure_features(
                Figure,
                ax,
                p_name,
                ascissa,
                sup_title,
                leg_title,
                dict_save_fig_path["eps"],
            )
        elif N_axes[nn] == 4:
            # 4 subplot -> 4 row, 1 column (cdp, 10/2020)
            Figure, ax = plt.subplots(
                nrows=4, ncols=1, figsize=(10.0, 8.5), sharex=True
            )
            if N_lines_tot[nn] != 18:
                for ii in range(N_lines[nn]):
                    ax[0].plot(
                        xvalues[header[ii]],
                        yvalues[header[ii]],
                        colors[ii],
                        linewidth=2,
                        label=tstep[ii],
                    )
                    ax[1].plot(
                        xvalues[header[ii + N_lines[nn]]],
                        yvalues[header[ii + N_lines[nn]]],
                        colors[ii],
                        linewidth=2,
                        label=tstep[ii + N_lines[nn]],
                    )
                    ax[2].plot(
                        xvalues[header[ii + 2 * N_lines[nn]]],
                        yvalues[header[ii + 2 * N_lines[nn]]],
                        colors[ii],
                        linewidth=2,
                        label=tstep[ii + 2 * N_lines[nn]],
                    )
                    if ii + 3 * N_lines[nn] < N_lines_tot[nn]:
                        ax[3].plot(
                            xvalues[header[ii + 3 * N_lines[nn]]],
                            yvalues[header[ii + 3 * N_lines[nn]]],
                            colors[ii],
                            linewidth=2,
                            label=tstep[ii + 3 * N_lines[nn]],
                        )
                    # end if ii + 3*N_lines[nn] (cdp, 01/2021)
                # end for ii (cdp, 01/2021)
                if N_lines_tot[nn] == 17:
                    # special case (cdp, 01/2021)
                    ax[3].plot(
                        xvalues[header[-1]],
                        yvalues[header[-1]],
                        colors[-1],
                        linewidth=2,
                        label=tstep[-1],
                    )
                # end if N_lines_tot (cdp, 01/2021)
            else:
                # N_lines_tot = 18: special case (cdp, 01/2021)
                for ii in range(N_lines[nn] + 1):
                    ax[0].plot(
                        xvalues[header[ii]],
                        yvalues[header[ii]],
                        colors[ii],
                        linewidth=2,
                        label=tstep[ii],
                    )
                    ax[1].plot(
                        xvalues[header[ii + N_lines[nn] + 1]],
                        yvalues[header[ii + N_lines[nn] + 1]],
                        colors[ii],
                        linewidth=2,
                        label=tstep[ii + N_lines[nn] + 1],
                    )
                # end for ii (cdp, 01/2021)
                for ii in range(N_lines[nn]):
                    ax[2].plot(
                        xvalues[header[ii + 11]],
                        yvalues[header[ii + 11]],
                        colors[ii],
                        linewidth=2,
                        label=tstep[ii + 11],
                    )
                    ax[3].plot(
                        xvalues[header[ii + 14]],
                        yvalues[header[ii + 14]],
                        colors[ii],
                        linewidth=2,
                        label=tstep[ii + 14],
                    )
                # end for ii (cdp, 01/2021)
            # end if N_lines_tot (cdp, 01/2021)
            # Call fuction Add_figure_features to coplete and save the figure \
            # (cdp, 01/2021)
            add_figure_features(
                Figure,
                ax,
                p_name,
                ascissa,
                sup_title,
                leg_title,
                dict_save_fig_path["eps"],
            )
        elif N_axes[nn] > 4:
            raise ValueError(
                "Not possible to plot more tha 4 axes for each figure. There is something odd in function Make_plots during evaluation of the number of figures and axes.\n"
            )
        # end if N_axes[nn] (cdp, 10/2020)
    # end for nn (cdp, 11/2020)


# end function Make_plots_sd_actually (cdp, 01/2021)


def make_plots_te_actually(
    values,
    p_name,
    folder_save,
    N_lines_tot_max,
    N_lines_tot,
    N_lines,
    N_axes,
    colors,
    zcoord_full,
    ascissa,
    sup_title,
    leg_title,
    figures=1,
):
    """
    Function that actually makes plots of the time evolution. (cdp, 11/2020)
    Updated (cdp, 01/2021)
    """
    if p_name == "inlet_outlet":
        # when data came from file inlet_outlet, the required property is the mass \
        # flow rate (cdp, 10/2020)
        p_name = "mass_flow_rate"
    else:
        # if data are not from file inlet_outlet.tsv, construct the list of \
        # columns header exploited to make plots (cdp, 11/2020)
        header_full = list(values.columns.values.tolist())
    # end if p_name (cdp, 10/2020)
    # Loop on figures (cdp, 11/2020)
    for nn in range(figures):
        # Path to save figures (cdp, 10/2020)
        if figures > 1:
            dict_save_fig_path = dict(
                eps=os.path.join(folder_save, f"{p_name}_{nn+1}.svg"),
                pdf=os.path.join(folder_save, f"{p_name}_{nn+1}.pdf"),
            )
            # extract the header from header_full to correctly build the lengend and \
            # do not exceed the number of lines for each figures (cdp, 01/2021)
            if nn < figures - 1:
                header = header_full[
                    1 + nn * N_lines_tot_max : (nn + 1) * N_lines_tot_max + 1
                ]
                zcoord = zcoord_full[nn * N_lines_tot_max : (nn + 1) * N_lines_tot_max]
            else:
                header = header_full[nn * N_lines_tot_max + 1 :]
                zcoord = zcoord_full[nn * N_lines_tot_max :]
            # end if nn (cdp, 01/2021)
        else:
            dict_save_fig_path = dict(
                eps=os.path.join(folder_save, f"{p_name}.svg"),
                pdf=os.path.join(folder_save, f"{p_name}.pdf"),
            )
            if p_name != "mass_flow_rate":
                header = header_full[1:]
                zcoord = zcoord_full
        # end if figures (cdp, 11/2020)
        if N_axes[nn] == 1:
            # one subplot -> 1 row, 1 column (cdp, 10/2020)
            Figure, ax = plt.subplots(nrows=N_axes[nn], ncols=1, figsize=(7.0, 6.0))
            if p_name == "mass_flow_rate":
                # plot inlet and outlet mass flow rate (cdp, 10/2020)
                mfr_label = [f"{p_name}_inl (kg/s)", f"{p_name}_out (kg/s)"]
                io_label = ["inlet", "outlet"]
                for ii in range(len(mfr_label)):
                    ax.plot(
                        values["time (s)"],
                        values.loc[:, mfr_label[ii]],
                        colors[ii],
                        linewidth=2,
                        label=io_label[ii],
                    )
                    # end for ii (cdp, 11/2020)
            else:
                # plot other properties at given spatial coordinates (cdp, 11/2020)
                ii = 0
                for head in header:
                    ax.plot(
                        values["time (s)"],
                        values.loc[:, head],
                        colors[ii],
                        linewidth=2,
                        label=zcoord[ii],
                    )
                    ii += 1
                # end for head (cdp, 11/2020)
            # end if p_name (cdp, 11/2020)
            # Call fuction Add_figure_features to coplete and save the figure \
            # (cdp, 01/2021)
            ax = [ax]  # to avoid error (cdp, 01/2021)
            add_figure_features(
                Figure,
                ax,
                p_name,
                ascissa,
                sup_title,
                leg_title,
                dict_save_fig_path["eps"],
            )
        elif N_axes[nn] == 2:
            # 2 subplot -> 2 row, 1 column (cdp, 10/2020)
            Figure, ax = plt.subplots(nrows=2, ncols=1, figsize=(7.0, 6.0), sharex=True)
            for ii in range(N_lines[nn]):
                ax[0].plot(
                    values["time (s)"],
                    values[header[ii]],
                    colors[ii],
                    linewidth=2,
                    label=zcoord[ii],
                )
                if ii + N_lines[nn] < N_lines_tot[nn]:
                    ax[1].plot(
                        values["time (s)"],
                        values[header[ii + N_lines[nn]]],
                        colors[ii],
                        linewidth=2,
                        label=zcoord[ii + N_lines[nn]],
                    )
                # end if ii + N_lines[nn] (cdp, 01/2021)
            # end for ii (cdp, 01/2021)
            # Call fuction Add_figure_features to coplete and save the figure \
            # (cdp, 01/2021)
            add_figure_features(
                Figure,
                ax,
                p_name,
                ascissa,
                sup_title,
                leg_title,
                dict_save_fig_path["eps"],
            )
        elif N_axes[nn] == 3:
            # 3 subplot -> 3 row, 1 column (cdp, 10/2020)
            Figure, ax = plt.subplots(
                nrows=3, ncols=1, figsize=(10.0, 8.5), sharex=True
            )
            # loop on saved values (cdp, 10/2020)
            for ii in range(N_lines[nn]):
                ax[0].plot(
                    values["time (s)"],
                    values[header[ii]],
                    colors[ii],
                    linewidth=2,
                    label=zcoord[ii],
                )
                ax[1].plot(
                    values["time (s)"],
                    values[header[ii + N_lines[nn]]],
                    colors[ii],
                    linewidth=2,
                    label=zcoord[ii + N_lines[nn]],
                )
                if ii + 2 * N_lines[nn] < N_lines_tot[nn]:
                    ax[2].plot(
                        values["time (s)"],
                        values[header[ii + 2 * N_lines[nn]]],
                        colors[ii],
                        linewidth=2,
                        label=zcoord[ii + 2 * N_lines[nn]],
                    )
                # end if ii + 2*N_lines[nn] (cdp, 01/2021)
            # end for ii (cdp, 01/2021)
            if N_lines_tot[nn] == 13:
                # special case (cdp, 01/2021)
                ax[2].plot(
                    values["time (s)"],
                    values[header[-1]],
                    colors[-1],
                    linewidth=2,
                    label=zcoord[-1],
                )
            # end if N_lines_tot (cdp, 01/2021)
            # Call fuction Add_figure_features to coplete and save the figure \
            # (cdp, 01/2021)
            add_figure_features(
                Figure,
                ax,
                p_name,
                ascissa,
                sup_title,
                leg_title,
                dict_save_fig_path["eps"],
            )
        elif N_axes[nn] == 4:
            # 4 subplot -> 4 row, 1 column (cdp, 10/2020)
            Figure, ax = plt.subplots(
                nrows=4, ncols=1, figsize=(10.0, 8.5), sharex=True
            )
            if N_lines_tot[nn] != 18:
                for ii in range(N_lines[nn]):
                    ax[0].plot(
                        values["time (s)"],
                        values[header[ii]],
                        colors[ii],
                        linewidth=2,
                        label=zcoord[ii],
                    )
                    ax[1].plot(
                        values["time (s)"],
                        values[header[ii + N_lines[nn]]],
                        colors[ii],
                        linewidth=2,
                        label=zcoord[ii + N_lines[nn]],
                    )
                    ax[2].plot(
                        values["time (s)"],
                        values[header[ii + 2 * N_lines[nn]]],
                        colors[ii],
                        linewidth=2,
                        label=zcoord[ii + 2 * N_lines[nn]],
                    )
                    if ii + 3 * N_lines[nn] < N_lines_tot[nn]:
                        ax[3].plot(
                            values["time (s)"],
                            values[header[ii + 3 * N_lines[nn]]],
                            colors[ii],
                            linewidth=2,
                            label=zcoord[ii + 3 * N_lines[nn]],
                        )
                    # end if ii + 3*N_lines[nn] (cdp, 01/2021)
                # end for ii (cdp, 01/2021)
                if N_lines_tot[nn] == 17:
                    # special case (cdp, 01/2021)
                    ax[3].plot(
                        values["time (s)"],
                        values[header[-1]],
                        colors[-1],
                        linewidth=2,
                        label=zcoord[-1],
                    )
                # end if N_lines_tot (cdp, 01/2021)
            else:
                # N_lines_tot = 18: special case (cdp, 01/2021)
                for ii in range(N_lines[nn] + 1):
                    ax[0].plot(
                        values["time (s)"],
                        values[header[ii]],
                        colors[ii],
                        linewidth=2,
                        label=zcoord[ii],
                    )
                    ax[1].plot(
                        values["time (s)"],
                        values[header[ii + N_lines[nn] + 1]],
                        colors[ii],
                        linewidth=2,
                        label=zcoord[ii + N_lines[nn] + 1],
                    )
                # end for ii (cdp, 01/2021)
                for ii in range(N_lines[nn]):
                    ax[2].plot(
                        values["time (s)"],
                        values[header[ii + 11]],
                        colors[ii],
                        linewidth=2,
                        label=zcoord[ii + 11],
                    )
                    ax[3].plot(
                        values["time (s)"],
                        values[header[ii + 14]],
                        colors[ii],
                        linewidth=2,
                        label=zcoord[ii + 14],
                    )
                # end for ii (cdp, 01/2021)
            # end if N_lines_tot (cdp, 01/2021)
            # Call fuction Add_figure_features to coplete and save the figure \
            # (cdp, 01/2021)
            add_figure_features(
                Figure,
                ax,
                p_name,
                ascissa,
                sup_title,
                leg_title,
                dict_save_fig_path["eps"],
            )
        elif N_axes[nn] > 4:
            raise ValueError(
                "Not possible to plot more than 4 axes for each figure. There is something odd in function Make_plots during evaluation of the number of figures and axes.\n"
            )
        # end if N_axes[nn] (cdp, 10/2020)
    # end for nn (cdp, 11/2020)


# end function Make_plots_te_actually (cdp, 01/2021)


def add_figure_features(
    Figure, ax, p_name, ascissa, sup_title, leg_title, save_fig_path
):
    """
    Function that adds fetures to the figure like gigure title, grid, legend axes labels and saves the figure. (cdp, 01/2021)
    """
    if "Heat_rad" in p_name:
        p_name = "Heat_rad"
    elif "Heat_exch" in p_name:
        p_name = "Heat_exch"
    # End if "Heat_rad".
    # get symbol and unit of the property p_name calling functions \
    # Plot_nomenclature_y (cdp, 10/2020)
    y_name, units = plot_nomenclature_y(p_name)
    # Single title for all the axes (cdp, 10/2020)
    Figure.suptitle(sup_title)
    # Compile ax (cdp, 10/2020)
    for i_ax in range(len(ax)):
        ax[i_ax].grid(True)
        ax[i_ax].set(xlabel=ascissa, ylabel=f"{y_name} {units}")
        ax[i_ax].legend(
            bbox_to_anchor=(1.135, 1.0),
            loc="upper right",
            fontsize=8,
            framealpha=0.2,
            title=leg_title,
        )
    # end for i_ax (cdp, 10/2020)
    plt.savefig(save_fig_path, orientation="landscape", transparent=False)
    plt.close()


# end function Add_figure_features (cdp, 01/2021)


def create_real_time_plots(simulation, conductor):
    """Function create_real_time_plots allows to create the required attributes to buld the real time plots (rtp).

    Args:
        simulation (object): object simulation instance of class Simulation.
        conductor (object): object conductor instance of class Conductor.
    """
    for f_comp in conductor.inventory["FluidComponent"].collection:
        f_comp.create_rtp_max_temperature = {
            True: create_real_time_plots_max_temperature,
            False: do_nothing,
        }
        f_comp.create_rtp_io_mfr = {
            True: create_real_time_plots_inlet_outlet_mfr,
            False: do_nothing,
        }
        f_comp.update_rtp_max_temperature = {
            True: update_real_time_plots_max_temperature,
            False: do_nothing,
        }
        f_comp.update_rtp_io_mfr = {
            True: update_real_time_plots_inlet_outlet_mfr,
            False: do_nothing,
        }
        f_comp.add_legend_rtp_io_mfr = {True: add_legend_rtp, False: do_nothing}
        # When the GUI works add the default to true.

        # If flag Show_fig is set to TRUE define the required attributes to make the real time plot of maximum f_comp temperature (invoke function create_real_time_plots_max_temperature); else does nothing.
        f_comp.create_rtp_max_temperature[f_comp.coolant.inputs["Show_fig"]](
            simulation, conductor, f_comp
        )
        # If flag Show_fig is set to TRUE define the required attributes to make the real time plot of inlet and outlet f_comp mass flow rate (invoke function create_real_time_plots_inlet_outlet_mfr); else does nothing.
        f_comp.create_rtp_io_mfr[f_comp.coolant.inputs["Show_fig"]](
            simulation, conductor, f_comp
        )
    # End for f_comp.

    for s_comp in conductor.inventory["SolidComponent"].collection:
        s_comp.create_rtp_max_temperature = {
            True: create_real_time_plots_max_temperature,
            False: do_nothing,
        }
        s_comp.update_rtp_max_temperature = {
            True: update_real_time_plots_max_temperature,
            False: do_nothing,
        }
        # If flag Show_fig is set to TRUE define the required attributes to make the real time plot of maximum s_comp temperature (invoke function create_real_time_plots_max_temperature); else does nothing.
        s_comp.create_rtp_max_temperature[s_comp.inputs["Show_fig"]](
            simulation, conductor, s_comp
        )
    # End for s_comp.


# End function _create_real_time_plots.


def create_real_time_plots_max_temperature(simulation, conductor, comp):
    """Function create_real_time_plots_max_temperature creates the required attribures of object comp to make the real time plots of the maximum temperature time evolution.

    Args:
        simulation (object): object simulation instance of class Simulation.
        conductor (object): object conductor instance of class Conductor.
        comp (object): can be any object instance of classes FluidComponent, MixSCStabilizet, StackComponent StrandStabilizerComponent or JacketComponent.
    """
    comp.figure_max_temp, comp.axes_max_temp = plt.subplots(
        num=f"{simulation.transient_input['SIMULATION']} ({conductor.number}): maximum {comp.identifier} temperature",
        figsize=(5, 5),
    )
    # set axes features (cdp, 10/2020)
    comp.axes_max_temp.grid(True)
    comp.axes_max_temp.set(
        xlabel="$t\ (s)$",
        ylabel="$T_{max}\ (K)$",
        title=f"{conductor.identifier} {comp.identifier} max temperature time evol",
    )
    comp.axes_max_temp.set_xlim([0.0, simulation.transient_input["TEND"]])


# End function create_real_time_plots_max_temperature.


def create_real_time_plots_inlet_outlet_mfr(simulation, conductor, f_comp):
    """Function create_real_time_plots_inlet_outlet_mfr creates the required attribures of object fluid_component to make the real time plots of the inlet and outlet mass flow rates time evolution.

    Args:
        simulation (object): object simulation instance of class Simulation.
        conductor (object): object conductor instance of class Conductor.
        f_comp (object): object fluid_component instance of class FluidComponent.
    """
    f_comp.figure_io_mfr, f_comp.axes_io_mfr = plt.subplots(
        num=f"{simulation.transient_input['SIMULATION']} ({conductor.number}): {f_comp.identifier} mass flow rates",
        figsize=(5, 5),
    )

    # set axes features (cdp, 10/2020)
    f_comp.axes_io_mfr.grid(True)
    f_comp.axes_io_mfr.set(
        xlabel="$t\ (s)$",
        ylabel="$mdot\ (kg/s)$",
        title=f"{conductor.identifier} {f_comp.identifier} inlet and outlet mfr time evol",
    )
    f_comp.axes_io_mfr.set_xlim([0.0, simulation.transient_input["TEND"]])


# End function create_real_time_plots_inlet_outlet_mfr.


def do_nothing(simulation, cond=None, comp=None):
    """Function do_nothing does nothing; necessay in the branchless coding.

    Args:
        simulation (object): object simulation instance of class Simulation
        cond (objetc, optional): object conductor instance of class Conductor.. Defaults to None.
        comp (object, optional): can be any object instance of classes FluidComponent, MixSCStabilizet, StackComponent StrandStabilizerComponent or JacketComponent. Defaults to None.
    """
    pass


# End function do_nothing.


def update_real_time_plots(conductor):
    """Function update_real_time_plots is the main function that allows to update the real time plots at each time step.

    Args:
        conductor (object): object conductor instance of class Conductor.
    """
    for f_comp in conductor.inventory["FluidComponent"].collection:

        # When the GUI works add the default to true.

        # If flag Show_fig is set to TRUE update the real time plot of maximum f_comp temperature (invoke function update_real_time_plots_max_temperature); else does nothing.
        f_comp.update_rtp_max_temperature[f_comp.coolant.inputs["Show_fig"]](
            conductor, f_comp
        )
        # If flag Show_fig is set to TRUE update the real time plot of inlet and outlet f_comp mass flow rate (invoke function update_real_time_plots_inlet_outlet_mfr); else does nothing.
        f_comp.update_rtp_io_mfr[f_comp.coolant.inputs["Show_fig"]](conductor, f_comp)
    # End for f_comp.

    for s_comp in conductor.inventory["SolidComponent"].collection:
        # If flag Show_fig is set to TRUE update the real time plot of maximum s_comp temperature (invoke function update_real_time_plots_max_temperature); else does nothing.
        s_comp.update_rtp_max_temperature[s_comp.inputs["Show_fig"]](conductor, s_comp)
        # End for s_comp.

    plt.pause(5e-3)

    # End function update_real_time_plots.


def update_real_time_plots_max_temperature(conductor, comp):
    """Function update_real_time_plots_max_temperature updates the real time plots of the maximum temperature time evolution of each component (comp) object.

    Args:
        conductor (object): object conductor instance of class Conductor.
        comp (object): can be any object instance of classes FluidComponent, MixSCStabilizet, StackComponent StrandStabilizerComponent or JacketComponent.
    """
    if comp.KIND == "Fluid_component":
        comp.axes_max_temp.plot(
            conductor.cond_time[-1],
            comp.coolant.dict_node_pt["temperature"].max(),
            "k.",
        )
    else:
        comp.axes_max_temp.plot(
            conductor.cond_time[-1], comp.dict_node_pt["temperature"].max(), "k."
        )


# End function update_real_time_plots_max_temperature.


def update_real_time_plots_inlet_outlet_mfr(conductor, f_comp):
    """Function update_real_time_plots_inlet_outlet_mfr updates the real time plots of the inlet and outlet time evolution of fluid_component (f_comp) objects.

    Args:
        conductor (object): object conductor instance of class Conductor.
        f_comp (object): object fluid_component instance of class FluidComponent.
    """
    # Aggiustare per tener conto della effettiva direzione del flusso di refrigerante.
    f_comp.axes_io_mfr.plot(
        conductor.cond_time[-1],
        f_comp.coolant.dict_node_pt["mass_flow_rate"][0],
        "co",
        label="Inlet",
    )
    f_comp.axes_io_mfr.plot(
        conductor.cond_time[-1],
        f_comp.coolant.dict_node_pt["mass_flow_rate"][-1],
        "b.",
        label="Outlet",
    )


# End function update_real_time_plots_inlet_outlet_mfr.


def create_legend_rtp(conductor):
    """Function create_legend_rtp is the main fucntion to add legend to the real time plots.

    Args:
        conductor (object): object conductor instance of class Conductor.
    """
    for f_comp in conductor.inventory["FluidComponent"].collection:
        # If flag Show_fig is set to TRUE create the legend in the time plot of inlet and outlet f_comp mass flow rate (invoke function add_legend_mfr); else does nothing.
        f_comp.add_legend_rtp_io_mfr[f_comp.coolant.inputs["Show_fig"]](f_comp)
    # End for create_legend_rtp.


# End function create_legend_rtp.


def add_legend_rtp(f_comp):
    """Function add_legend_rtp adds legend to the real time plot of inlet and outlet mass flow rate time evolutions for fluid_components objects.

    Args:
        f_comp (object): object fluid_component instance of class FluidComponent.
    """
    f_comp.axes_io_mfr.legend(loc="best", fontsize=8, framealpha=0.2)
    # comp.axes_max_temp.legend(loc = 'best', fontsize = 8, framealpha = 0.2)


# End function add_legend_rtp.


def plot_time_animation(simulation, conductor):

    """
    Function that make an animated plot of the time evolution of the maximum temperature of the strands. (cdp, 10/2020)
    """

    list_types = ["FluidComponent", "StrandComponent", "JacketComponent"]
    for l_type in list_types:
        if conductor.cond_time[-1] == 0:
            # Define some immutable features of the plots only once (cdp, 10/2020)
            # Construct Figure and axes object related to maximum temperature plots, \
            # storing them in dictionaries (cdp, 10/2020)
            (
                conductor.dict_Figure_animation["T_max"][l_type],
                conductor.dict_axes_animation["T_max"][l_type],
            ) = plt.subplots(
                num=f"{simulation.transient_input['SIMULATION']} ({conductor.number}): maximum {l_type} temperature",
                figsize=(5, 5),
            )
            if l_type == "FluidComponent":
                # Create a number of Figure and axes objects equal to the number of \
                # FluidComponent objects to plot the inlet and outlet mass flow rate \
                # for each FluidComponent object in a dedicated plot. (cdp, 10/2020)
                # conductor.dict_Figure_animation["mfr"][l_type] = dict()
                # conductor.dict_axes_animation["mfr"][l_type] = dict()
                for fluid_comp in conductor.inventory["FluidComponent"].collection:
                    (
                        conductor.dict_Figure_animation["mfr"][fluid_comp.identifier],
                        conductor.dict_axes_animation["mfr"][fluid_comp.identifier],
                    ) = plt.subplots(
                        num=f"{simulation.transient_input['SIMULATION']} ({conductor.number}): {fluid_comp.identifier} "
                        + f"mass flow rate",
                        figsize=(5, 5),
                    )
                    # set axes features (cdp, 10/2020)
                    conductor.dict_axes_animation["mfr"][fluid_comp.identifier].grid(
                        True
                    )
                    conductor.dict_axes_animation["mfr"][fluid_comp.identifier].set(
                        xlabel="$t\ (s)$",
                        ylabel="$mdot\ kg/s$",
                        title=f"{fluid_comp.identifier} inlet and outlet mfr time evol",
                    )
                    conductor.dict_axes_animation["mfr"][
                        fluid_comp.identifier
                    ].set_xlim([0.0, simulation.transient_input["TEND"]])
                # end for fluid_comp (cdp, 10/2020)
            # end if l_type (cdp, 10/2020)
            # set axes features (cdp, 10/2020)
            conductor.dict_axes_animation["T_max"][l_type].grid(True)
            conductor.dict_axes_animation["T_max"][l_type].set(
                xlabel="$t\ (s)$",
                ylabel="$T_{max}\ K$",
                title=f"{l_type} max temperature time evol",
            )
            conductor.dict_axes_animation["T_max"][l_type].set_xlim(
                [0.0, simulation.transient_input["TEND"]]
            )
        # end if conductor.cond_time[-1] (cdp, 10/2020)
        if l_type == "FluidComponent":
            # loop on FluidComponent (cdp, 10/2020)
            for fluid_comp in conductor.inventory["FluidComponent"].collection:
                # make the plot of channels maximum temperature (cdp, 10/2020)
                # conductor.dict_axes_animation["T_max"][l_type].plot(
                #     conductor.cond_time[-1], fluid_comp.coolant.dict_node_pt["temperature"].max(), conductor.color[ii], label = fluid_comp.identifier) # choose the color
                conductor.dict_axes_animation["T_max"][l_type].plot(
                    conductor.cond_time[-1],
                    fluid_comp.coolant.dict_node_pt["temperature"].max(),
                    label=fluid_comp.identifier,
                )
                # make plot of inlet and outlet mass flow rate; each channel has its \
                # own figures (cdp, 10/2020)
                conductor.dict_axes_animation["mfr"][fluid_comp.identifier].plot(
                    conductor.cond_time[-1],
                    fluid_comp.coolant.dict_node_pt["mass_flow_rate"][0],
                    conductor.color[0],
                    label="Inlet",
                )
                conductor.dict_axes_animation["mfr"][fluid_comp.identifier].plot(
                    conductor.cond_time[-1],
                    fluid_comp.coolant.dict_node_pt["mass_flow_rate"][-1],
                    conductor.color[1],
                    label="Outlet",
                )
                if conductor.cond_time[-1] == 0:
                    # define legend only once (cdp, 10/2020)
                    conductor.dict_axes_animation["T_max"][l_type].legend(
                        loc="best", fontsize=8, framealpha=0.2
                    )
                    # conductor.dict_axes_animation["T_max"][l_type].tick_params(\
                    #   labelsize = 8, colors = "r", length = 4)
                    conductor.dict_axes_animation["mfr"][fluid_comp.identifier].legend(
                        loc="best", fontsize=8, framealpha=0.2
                    )
                # end if conductor.cond_time[-1] (cdp, 10/2020)
            # end for ii (cdp, 10/2020)
        elif l_type == "StrandComponent":
            # loop on StrandComponent (cdp, 10/2020)
            for strand in conductor.inventory["StrandComponent"].collection:
                # plot the maximum strand temperature (cdp, 10/2020)
                # conductor.dict_axes_animation["T_max"][l_type].plot(
                #     conductor.cond_time[-1], strand.dict_node_pt["temperature"].max(), conductor.color[ii], label = strand.identifier) # choose the color.
                conductor.dict_axes_animation["T_max"][l_type].plot(
                    conductor.cond_time[-1],
                    strand.dict_node_pt["temperature"].max(),
                    label=strand.identifier,
                )
            # end for strand (cdp, 10/2020)
            if conductor.cond_time[-1] == 0:
                # define legend only once (cdp, 10/2020)
                conductor.dict_axes_animation["T_max"][l_type].legend(
                    loc="best", fontsize=8, framealpha=0.2
                )
        elif l_type == "JacketComponent":
            # for the time being do not make figures for JacketComponents objects \
            # (cdp, 10/2020)
            pass
        # end if l_type (cdp, 10/2020)
    # end for l_type (cdp, 10/2020)

    if conductor.cond_time[-1] == 0:
        # plt.show(block = False)
        # close figure related to JacketComponent objects (cdp, 10/2020)
        plt.close(conductor.dict_Figure_animation["T_max"]["JacketComponent"])
        # gui.plots_window.deiconify()
    plt.pause(0.01)


# end function Plot_time_animation (cdp, 10/2020)
