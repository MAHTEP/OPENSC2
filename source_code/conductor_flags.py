from collections import namedtuple

# Flags for the current definition

# User does not define a current: do not use the electric module
IOP_NOT_DEFINED = None
# Constant value read from input file conductor_definition.xlsx
IOP_CONSTANT = 0
# Initial current spatial distribution load form auxiliary input file
IOP_FROM_FILE = -1
# Current beavior as function of space and time from user defined function
IOP_FROM_EXT_FUNCTION = -2

# Flags for electric conductance

# Electric conductance is defined per unit length
ELECTRIC_CONDUCTANCE_UNIT_LENGTH = 1
# Electric conductance is not defined per unit length
ELECTRIC_CONDUCTANCE_NOT_UNIT_LENGTH = 2

# Analytical self inductance is evaluated according to mode 1.
SELF_INDUCTANCE_MODE_1 = 1
# Analytical self inductance is evaluated according to mode 2.
SELF_INDUCTANCE_MODE_2 = 2

# Flag to evaluate inductance analytically
ANALYTICAL_INDUCTANCE = 0
# Flag to evaluate inductance using an approximation
APPROXIMATE_INDUCTANCE = 1

# Flag to solve the electric problem in steady state conditions.
STATIC_ELECTRIC_SOLVER = 0

# Default number for electric time step
ELECTRIC_TIME_STEP_NUMBER = 10

# Flags for contact perimeter
# Variable contact perimeter (from auxiliary input file)
VARIABLE_CONTACT_PERIMETER = -1
# Constant contact perimeter (from sheet contact_perimeter in file 
# conductor_coupling.xlsx)
CONSTANT_CONTACT_PERIMETER = 1

# Component sheet names
Comp_sheet_name = namedtuple("Comp_sheet_name",
    (
        "fluid_comp",
        "stack",
        "str_mix",
        "stab",
        "jacket",
    )
)

# Dictionary with all the valid sheet names of each input file.
SHEET_NAME = dict(
    conductor_coupling = {
            "contact_perimeter_flag",
            "contact_perimeter",
            "HTC_choice",
            "contact_HTC",
            "thermal_contact_resistance",
            "HTC_multiplier",
            "electric_conductance_mode",
            "electric_conductance",
            "open_perimeter_fract",
            "interf_thickness",
            "trans_transp_multiplier",
            "view_factors",
            },
    conductor_definition = {
        "CONDUCTOR_files",
        "CONDUCTOR_input",
        "CONDUCTOR_operation",
        "CONDUCTOR_coupling",
    },
    conductor_diagnostic = {
        "Spatial_distribution",
        "Time_evolution",
        "Voltage_tap_name",
        "Voltage_tap_coordinate",
    },
    conductor_grid = {"GRID"},
    # Key conductor_input is a namedtuple and not a set because I want to 
    # exploit access by field in method conductor_component_instance.
    conductor_input = Comp_sheet_name(
        fluid_comp = "CHAN",
        stack = "STACK",
        str_mix = "STR_MIX",
        stab = "STR_STAB",
        jacket = "Z_JACKET",
        ),
    # Key conductor_input is a namedtuple and not a set because I want to 
    # exploit access by field in method conductor_component_instance.
    conductor_operation = Comp_sheet_name(
        fluid_comp = "CHAN",
        stack = "STACK",
        str_mix = "STR_MIX",
        stab = "STR_STAB",
        jacket = "Z_JACKET",
        ),
    environment_input = {"ENVIRONMENT"},
    transitory_intput = {"TRANSIENT"},
)