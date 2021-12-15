# Moule that manages the code Graphycal User Interface (GUI) (cdp, 12/2020).

import tkinter as tk
from tkinter import messagebox
from tkinter import filedialog
from tkinter import ttk
from PIL import ImageTk, Image
import os
import subprocess

from simulations import Simulations


class OPENSC2_GUI:
    """
    Class that manages the code GUI (cdp, 12/2020).
    """

    # GUI attributes
    root_window = ""
    main_window = ""
    # cascades
    input_menu = ""
    output_path_menu = ""
    control_panel_menu = ""
    interaction_mode_menu = ""
    help_menu = ""
    # icons
    img_gears_run = ""
    dict_input = ""
    dict_save_res = ""
    dict_control_panel = ""
    dict_drivers = ""
    dict_help = ""

    simulation = ""
    current_dir = ""

    def __init__(self):
        """
        Constructor method of class OPENSC2_GUI (cdp, 12/2020).
        """

        self.current_dir = os.path.abspath("")
        icons_dir = os.path.join(self.current_dir, "GUI_icons")
        self.root_window = tk.Tk()
        self.root_window.option_add("*tearOff", False)
        # Does not show the root window of the GUI (cdp, 12/2020)
        self.root_window.withdraw()
        # create a main window
        self.main_window = tk.Toplevel(master=self.root_window)
        self.main_window.wm_iconbitmap(os.path.join(icons_dir, "MAHTEP_LOGO.ico"))
        self.main_window.title("OPENSC2 GUI Main Window")
        # Start: Create images objects (cdp, 08/2020)
        # Simulation input cascade images (cdp, 08/2020)
        self.dict_input = dict(
            Load_data=(
                "Load input data",
                ImageTk.PhotoImage(
                    Image.open(os.path.join(icons_dir, "super_conductor_1_r.png"))
                ),
            ),
            Check_topology=(
                "Check conductor topology",
                ImageTk.PhotoImage(
                    Image.open(os.path.join(icons_dir, "super_conductor_2_r.png"))
                ),
            ),
            Initialization=(
                "Initialize variables",
                ImageTk.PhotoImage(
                    Image.open(os.path.join(icons_dir, "chart_scatter_plot.png"))
                ),
            ),
        )
        # Save Simulation results cascade images (cdp, 12/2020)
        self.dict_save_res = dict(
            Create=(
                "Create new directory",
                ImageTk.PhotoImage(
                    Image.open(os.path.join(icons_dir, "folder_plus.png"))
                ),
            ),
            Open=(
                "Open existing directory",
                ImageTk.PhotoImage(
                    Image.open(os.path.join(icons_dir, "folder_open.png"))
                ),
            ),
        )
        # Simulation control panel images (cdp, 08/2020)
        self.dict_control_panel = dict(
            Run=(
                "Run",
                ImageTk.PhotoImage(
                    Image.open(os.path.join(icons_dir, "play_circle_outline.png"))
                ),
            ),
            Pause=(
                "Pause",
                ImageTk.PhotoImage(
                    Image.open(os.path.join(icons_dir, "pause_circle_outline.png"))
                ),
            ),
            Continue=(
                "Continue",
                ImageTk.PhotoImage(
                    Image.open(os.path.join(icons_dir, "play_box_outline.png"))
                ),
            ),
            Stop=(
                "Stop",
                ImageTk.PhotoImage(
                    Image.open(os.path.join(icons_dir, "stop_circle_outline.png"))
                ),
            ),
            Save_status=(
                "Save status",
                ImageTk.PhotoImage(
                    Image.open(os.path.join(icons_dir, "content_save_outline.png"))
                ),
            ),
            Close=(
                "Close",
                ImageTk.PhotoImage(
                    Image.open(os.path.join(icons_dir, "location_exit.png"))
                ),
            ),
        )
        # Change simulation drivers images (cdp, 08/2020)
        self.dict_drivers = dict(
            Plot_sd=(
                "Plot spatial distributions",
                ImageTk.PhotoImage(
                    Image.open(os.path.join(icons_dir, "chart_scatter_plot.png"))
                ),
            ),
            Plot_te=(
                "Plot temperature evolutions",
                ImageTk.PhotoImage(
                    Image.open(os.path.join(icons_dir, "chart_scatter_plot.png"))
                ),
            ),
            Change=(
                "Change drivers",
                ImageTk.PhotoImage(
                    Image.open(os.path.join(icons_dir, "help_circle_outline.png"))
                ),
            ),
            Back=(
                "Go back",
                ImageTk.PhotoImage(
                    Image.open(
                        os.path.join(icons_dir, "arrow_left_thin_circle_outline.png")
                    )
                ),
            ),
        )
        # Help images (cdp, 08/2020)
        self.dict_help = dict(
            Input_files=(
                "User guide: input files",
                ImageTk.PhotoImage(
                    Image.open(os.path.join(icons_dir, "help_circle_outline.png"))
                ),
            )
        )
        # End: Create images objects (cdp, 08/2020)
        # menubar widget (cdp, 12/2020)
        menu_bar = tk.Menu(master=self.main_window, font="TkMenuFont")
        chosen_font = "Helvetica 12"
        # Create Simulation input cascade (cdp, 08/2020)
        self.input_menu = tk.Menu(master=menu_bar, font=chosen_font)
        self.input_menu.add_command(
            label=self.dict_input["Load_data"][0],
            image=self.dict_input["Load_data"][1],
            compound="left",
            command=self.load_input_file,
        )
        self.input_menu.add_command(
            label=self.dict_input["Check_topology"][0],
            image=self.dict_input["Check_topology"][1],
            compound="left",
            state="disabled",
        )
        self.input_menu.add_command(
            label=self.dict_input["Initialization"][0],
            image=self.dict_input["Initialization"][1],
            compound="left",
            state="disabled",
        )
        # Create Save Simulation results cascade (cdp, 12/2020)
        self.output_path_menu = tk.Menu(master=menu_bar, font=chosen_font)
        self.output_path_menu.add_command(
            label=self.dict_save_res["Create"][0],
            state="disabled",
            image=self.dict_save_res["Create"][1],
            compound="left",
            command=self.create_directories,
        )
        self.output_path_menu.add_command(
            label=self.dict_save_res["Open"][0],
            image=self.dict_save_res["Open"][1],
            compound="left",
            command=self.open_existing_directories,
            state="disabled",
        )
        # Create Simulation control panel cascade (cdp, 08/2020)
        self.control_panel_menu = tk.Menu(master=menu_bar, font=chosen_font)
        self.control_panel_menu.add_command(
            label=self.dict_control_panel["Run"][0],
            state="disabled",
            image=self.dict_control_panel["Run"][1],
            compound="left",
            command=self.run_simulation,
        )
        self.control_panel_menu.add_command(
            label=self.dict_control_panel["Pause"][0],
            state="disabled",
            image=self.dict_control_panel["Pause"][1],
            compound="left",
        )
        self.control_panel_menu.add_command(
            label=self.dict_control_panel["Continue"][0],
            state="disabled",
            image=self.dict_control_panel["Continue"][1],
            compound="left",
        )
        self.control_panel_menu.add_command(
            label=self.dict_control_panel["Stop"][0],
            state="disabled",
            image=self.dict_control_panel["Stop"][1],
            compound="left",
        )
        self.control_panel_menu.add_separator()
        self.control_panel_menu.add_command(
            label=self.dict_control_panel["Save_status"][0],
            state="disabled",
            image=self.dict_control_panel["Save_status"][1],
            compound="left",
        )
        self.control_panel_menu.add_separator()
        self.control_panel_menu.add_command(
            label=self.dict_control_panel["Close"][0],
            state="disabled",
            image=self.dict_control_panel["Close"][1],
            compound="left",
            command=self.close_simulation,
        )
        # Create Change simulation drivers cascade (cdp, 08/2020)
        self.interaction_mode_menu = tk.Menu(master=menu_bar, font=chosen_font)
        self.interaction_mode_menu.add_command(
            label=self.dict_drivers["Plot_sd"][0],
            image=self.dict_drivers["Plot_sd"][1],
            compound="left",
            state="disabled",
        )
        self.interaction_mode_menu.add_command(
            label=self.dict_drivers["Plot_te"][0],
            image=self.dict_drivers["Plot_te"][1],
            compound="left",
            state="disabled",
        )
        self.interaction_mode_menu.add_separator()
        self.interaction_mode_menu.add_command(
            label=self.dict_drivers["Change"][0],
            image=self.dict_drivers["Change"][1],
            compound="left",
            state="disabled",
        )
        self.interaction_mode_menu.add_separator()
        self.interaction_mode_menu.add_command(
            label=self.dict_drivers["Back"][0],
            image=self.dict_drivers["Back"][1],
            compound="left",
            state="disabled",
        )
        # Create Help? cascade (cdp, 08/2020)
        self.help_menu = tk.Menu(master=menu_bar, font=chosen_font)
        self.help_menu.add_command(
            label=self.dict_help["Input_files"][0],
            image=self.dict_help["Input_files"][1],
            compound="left",
            command=self.read_user_guide,
        )
        # list of tuple to create the menu (cdp, 12/2020)
        menu_bar_items = [
            ("Simulation input", self.input_menu),
            (
                "Assign output path",
                self.output_path_menu,
            ),  # old save simulation results
            ("Simulation control panel", self.control_panel_menu),
            (
                "Interaction mode                                   ",
                self.interaction_mode_menu,
            ),  # old Change simulation driver
            ("Help?", self.help_menu),
        ]
        # Create the menu (cdp, 08/2020)
        for lbl, cascade in menu_bar_items:
            menu_bar.add_cascade(label=lbl, menu=cascade)
        self.main_window["menu"] = menu_bar

        # self.plots_window = tk.Toplevel(master = self.root_window)
        # self.plots_window.title("Real time plots window")
        # self.plots_window.iconify()
        ## START: NOTEBOOK (cdp, 12/2020)
        ## create the Notebook objects (cdp, 12/2020)
        # self.notebook = tk.ttk.Notebook(master = self.plots_window)
        # self.notebook.pack()
        ## create the nested Notebook objects (cdp, 12/2020)
        # self.np_plots = tk.ttk.Notebook(master = self.plots_window)
        # self.np_plots.pack()
        ## FluidComponents sheet (cdp, 12/2020)
        # self.chan_sheet = tk.Frame(master = self.np_plots)
        ## Strands sheet (cdp, 12/2020)
        # self.str_sheet = tk.Frame(master = self.np_plots)
        ## Jackets sheet (cdp, 12/2020)
        # self.jk_sheet = tk.Frame(master = self.np_plots)
        ## add sheets to the notebook (cdp, 12/2020)
        # self.np_plots.add(self.chan_sheet, text = "FluidComponents")
        # self.np_plots.add(self.str_sheet, text = "Strands")
        # self.np_plots.add(self.jk_sheet, text = "Jackets")
        ## END: NOTEBOOK (cdp, 12/2020)

    # end method __init__ (cdp, 12/2020)

    def load_input_file(self):
        """
        Method that allows to select the folder with the input files and makes an instance of class Simulation that reads the first main input file, transitory_input.xlsx (cdp, 12/2020).
        """

        # Main directory with all the input files (cdp, 10/2020)
        main_input = "Description_of_Components"
        sub_input = tk.filedialog.askdirectory(
            parent=self.main_window,
            title="Select input files directory",
            initialdir=main_input,
            mustexist=True,
        )
        base_path = os.path.join(main_input, sub_input)

        # simulation instance (cdp, 08/2020)
        self.simulation = Simulations(base_path)
        # Enables the Save Simulation results cascade options (cdp, 12/2020)
        for key in list(self.dict_save_res.keys()):
            self.output_path_menu.entryconfigure(
                self.dict_save_res[key][0],
                image=self.dict_save_res[key][1],
                state="normal",
            )
        # end for key (cdp, 12/2020)
        # Disables Load data in Simulation input cascade (cdp, 12/2020)
        self.input_menu.entryconfigure(
            self.dict_input["Load_data"][0],
            image=self.dict_input["Load_data"][1],
            state="disable",
        )

    # end method Load_input_file

    def create_directories(self):

        """
        Method that asks user to provide the name to the main simulation directory (cdp, 10/2020). Updated (cdp, 12/2020)
        """
        Title = "Give a name to the set of simulation main directory"
        # Get the new folder name (cdp, 10/2020)
        self.simulation.dict_path["Main_dir"] = tk.filedialog.asksaveasfilename(
            parent=self.main_window,
            title=Title,
            initialdir=self.simulation.dict_path["Results_dir"],
        )
        # Disables the Save Simulation results cascade options (cdp, 12/2020)
        for key in list(self.dict_save_res.keys()):
            self.output_path_menu.entryconfigure(
                self.dict_save_res[key][0],
                image=self.dict_save_res[key][1],
                state="disabled",
            )
        # end for key (cdp, 12/2020)
        # Enable option Run simulation in Simulation control panel cascade \
        # (cdp, 12/2020)
        self.control_panel_menu.entryconfigure(
            self.dict_control_panel["Run"][0],
            state="normal",
            image=self.dict_control_panel["Run"][1],
            compound="left",
        )

    # end method Create_directories (cdp, 10/2020)

    def open_existing_directories(self):

        """
        Methoid that asks to select the existing simulation main directory to save the simulation results (cdp, 10/2020). Updated (cdp, 12/2020)
        """
        Title = "Open existing simulation main directory"
        # Get the existing folder name (cdp, 10/2020)
        self.simulation.dict_path["Main_dir"] = tk.filedialog.askdirectory(
            parent=self.main_window,
            mustexist=True,
            title=Title,
            initialdir=self.simulation.dict_path["Results_dir"],
        )
        # Disables the Save Simulation results cascade options (cdp, 12/2020)
        for key in list(self.dict_save_res.keys()):
            self.output_path_menu.entryconfigure(
                self.dict_save_res[key][0],
                image=self.dict_save_res[key][1],
                state="disabled",
            )
        # end for key (cdp, 12/2020)
        # Enable option Run simulation in Simulation control panel cascade \
        # (cdp, 12/2020)
        self.control_panel_menu.entryconfigure(
            self.dict_control_panel["Run"][0],
            state="normal",
            image=self.dict_control_panel["Run"][1],
            compound="left",
        )

    # end Method Open_existing_directories (cdp, 10/2020)

    def run_simulation(self):

        """
        Method that allows user to execute the simulation. (cdp, 08/2020)
        """
        list_cp_keys = list(self.dict_control_panel.keys())
        # Loop to enalble all the options except Run that is desabled (cdp, 12/2020)
        for key in list_cp_keys:
            if key == "Run":
                stato = "disabled"
            else:
                stato = "normal"
            # end if key (cdp, 12/2020)
            self.control_panel_menu.entryconfigure(
                self.dict_control_panel[key][0],
                state=stato,
                image=self.dict_control_panel[key][1],
                compound="left",
            )
        # end for key (cdp, 12/2020)
        self.simulation.flag_start = True
        # Start simulation message (cdp, 12/2020)
        messaggio = (
            "Launched simulation called "
            + self.simulation.transient_input["SIMULATION"]
            + "\n"
        )
        tk.messagebox.showinfo(
            parent=self.main_window, title=f"Start Simulation", message=messaggio
        )
        # Create and instance of class Conductor for each user defined conductor \
        # (cdp, 08/2020)
        self.simulation.conductor_instance()
        # Create the whole tree of folders to store the simulation data invoking \
        # method Simulation_result_manager (cdp, 10/2020)
        self.simulation.simulation_folders_manager()
        # Initialize each user defined conductor (cdp, 08/2020)
        self.simulation.conductor_initialization(self)
        # Solve the linear system of equations at each time steps (cdp, 08/2020)
        self.simulation.conductor_solution(self)
        # Create plots of time evlustions and spatial distributions according to \
        # user requirements (cdp, 08/2020)
        self.simulation.conductor_post_processing()
        # End simulation message (cdp, 12/2020)
        messaggio = (
            "Simulation called "
            + self.simulation.transient_input["SIMULATION"]
            + " ends.\n"
            + "End of data processing and saving of figures.\n"
        )
        tk.messagebox.showinfo(
            parent=self.main_window, title=f"End Simulation", message=messaggio
        )
        # Loop to disable all the commands except Close, still available \
        # (cdp, 12/2020)
        for key in list_cp_keys[1:]:
            if key != "Close":
                self.control_panel_menu.entryconfigure(
                    self.dict_control_panel[key][0],
                    state="disabled",
                    image=self.dict_control_panel[key][1],
                    compound="left",
                )
            # end if key (cdp, 12/2020)
        # end for key (cdp, 12/2020)

    # end method Run_simulation (cdp, 08/2020)

    def close_simulation(self):
        """
        Method that closes the simulation session killing the simulation if in progress (cdp, 12/2020)
        """
        messaggio = (
            "Do you want to close the simulation called "
            + self.simulation.transient_input["SIMULATION"]
            + " and kill it?\n"
        )
        res = tk.messagebox.askyesno(title="Exit SC2", message=messaggio)
        if res == True:
            # destroy the root window: kill everything (cdp, 12/2020)
            self.root_window.destroy()
        else:
            # return to the code (cdp, 12/2020)
            tk.messagebox.showinfo(title="Return", message="Returning to SC2")
        # end if (cdp, 12/2020)

    # end method Close_simulation (cdp, 12/2020)

    def read_user_guide(self):
        """
        Method that allows to open and read the file of the User Guide. (cdp, 12/2020)
        """
        # Path of the user guide folder (cdp, 12/2020)
        user_guide = os.path.join(self.current_dir, "User_guide")
        # get the name of the file to be opened (cdp, 12/2020)
        f_name = tk.filedialog.askopenfilename(
            parent=self.main_window, title="User guide", initialdir=user_guide
        )
        # Open the selected file in pdf (cdp, 12/2020)
        subprocess.Popen(os.path.join(user_guide, f_name), shell=True)

    # end method Read_user_guide (cdp, 12/2020)
