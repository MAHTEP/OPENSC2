import yaml
import logging
import logging.config

from simulation import Simulation

if __name__ == "__main__":

    logging.config.fileConfig(fname="logging_opensc2.conf", disable_existing_loggers=True)

    # Get the logger specified in the file, this will be the parent logger.
    logger = logging.getLogger("opensc2Logger")

    # Open .json file with paths to the input files of the simulation and the 
    # path to the folder where to save the output of the simulation.
    with open("io_path.yaml","r") as read_file:
        io_path = yaml.safe_load(read_file)

    # Create an instance of class Simulation
    simulation = Simulation(io_path)

    # Save the path to the folder where to save the output of the simulation 
    # in an attribute of class Simulation
    simulation.dict_path["Main_dir"] = io_path["output"]

    # Create and instance of class Conductor for each user defined conductor.
    simulation.conductor_instance()
    # Create the whole tree of folders to store the simulation data invoking 
    # method Simulation_result_manager.
    simulation.simulation_folders_manager()
    # Save the input files in read only as metadata of the simulation outcome.
    simulation.save_input_files()
    # Initialize each user defined conductor.
    simulation.conductor_initialization()
    # Solve the linear system of equations at each time steps.
    simulation.conductor_solution()
    # Create plots of time evlustions and spatial distributions according to 
    # user requirements.
    simulation.conductor_post_processing()
    # End simulation message.
    messaggio = (
        "Simulation called "
        + simulation.transient_input["SIMULATION"]
        + " ends.\n"
        + "End of data processing and saving of figures.\n"
    )
    
    print(messaggio)