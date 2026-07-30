import argparse
import logging
import logging.config
import time
from pathlib import Path

import yaml

from opensc2_gui import OPENSC2_GUI
from simulation import Simulation


class HeadlessGUI:
    """Minimal GUI-like object used to run OPENSC2 without opening the Tk GUI."""

    def __init__(self, simulation):
        self.simulation = simulation
        self.main_window = None


def parse_command_line_arguments():
    parser = argparse.ArgumentParser(
        description="Run OPENSC2 with or without the graphical user interface."
    )
    parser.add_argument(
        "--no-head",
        action="store_true",
        help="Run OPENSC2 without opening the GUI.",
    )
    parser.add_argument(
        "--io-path",
        default="io_path.yaml",
        help="Path to the YAML file used in --no-head mode.",
    )
    return parser.parse_args()


def load_io_path_from_yaml(io_path):
    io_path = Path(io_path)

    if not io_path.is_file():
        raise FileNotFoundError(f"YAML IO file not found: {io_path}")

    with io_path.open("r", encoding="utf-8") as stream:
        data = yaml.safe_load(stream)

    if data is None:
        raise ValueError(f"YAML IO file is empty: {io_path}")

    required_keys = ("input_dir", "output_dir")
    for key in required_keys:
        if key not in data:
            raise KeyError(f"Missing required key '{key}' in YAML IO file: {io_path}")

    input_dir = Path(data["input_dir"]).expanduser()
    output_dir = Path(data["output_dir"]).expanduser()

    if not input_dir.is_dir():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")

    return input_dir, output_dir


def format_elapsed_time(elapsed_time):
    if elapsed_time < 60:
        return elapsed_time, "s"
    if elapsed_time < 3600:
        return elapsed_time / 60, "min"
    if elapsed_time < 86400:
        return elapsed_time / 3600, "h"
    return elapsed_time / 86400, "days"


def run_headless_simulation(io_path):
    input_dir, output_dir = load_io_path_from_yaml(io_path)

    simulation = Simulation(str(input_dir))
    simulation.dict_path["Main_dir"] = str(output_dir)
    simulation.flag_start = True

    headless_gui = HeadlessGUI(simulation)

    print("Running OPENSC2 in headless mode.")
    print(f"Input directory: {input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Simulation: {simulation.transient_input['SIMULATION']}")

    tt = time.time()

    simulation.conductor_instance()
    simulation.simulation_folders_manager()
    simulation.save_input_files()
    simulation.conductor_initialization(headless_gui)
    simulation.conductor_solution(headless_gui)
    simulation.conductor_post_processing()

    elapsed_time = time.time() - tt
    elapsed_time_print, unit = format_elapsed_time(elapsed_time)

    print(f"\nSimulation completed in {elapsed_time_print:.2f} {unit}")


def run_gui():
    gui = OPENSC2_GUI()
    gui.main_window.mainloop()


logging.config.fileConfig(fname="logging_opensc2.conf", disable_existing_loggers=True)

logger = logging.getLogger("opensc2Logger")

args = parse_command_line_arguments()

if args.no_head:
    run_headless_simulation(args.io_path)
else:
    run_gui()