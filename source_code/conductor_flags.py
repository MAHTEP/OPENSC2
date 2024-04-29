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

# Self inductance readed from input file conductor_definition.xlsx
SELF_INDUCTANCE_MODE_0 = 0
# Analytical self inductance is evaluated according to mode 1.
SELF_INDUCTANCE_MODE_1 = 1
# Analytical self inductance is evaluated according to mode 2.
SELF_INDUCTANCE_MODE_2 = 2

# Constant mutual inductance readed from input file conductor_definition.xlsx
CONSTANT_INDUCTANCE = 0
# Flag to evaluate inductance analytically
ANALYTICAL_INDUCTANCE = 1
# Flag to evaluate inductance using an approximation
APPROXIMATE_INDUCTANCE = 2

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

# Namedtuple constructor for conductor sheet names
Cond_sheet_name = namedtuple("Cond_sheet_name",
    (
        "files",
        "inputs",
        "operation",
        "coupling",
    )
)

# Namedtuple constructor for diagnostic sheet names
Diagno_sheet_name = namedtuple("Cond_sheet_name",
    (
        "space_distr",
        "time_evol",
    )
)

# Namedtuple constructor for component sheet names
Comp_sheet_name = namedtuple("Comp_sheet_name",
    (
        "FluidComponent",
        "StackComponent",
        "StrandMixedComponent",
        "StrandStabilizerComponent",
        "JacketComponent",
    )
)

# Namedtuple constructor for environment sheet name
Env_sheet_name = namedtuple("Env_sheet_name",("environment"))

# Namedtuple constructor for transient sheet name
Trans_sheet_name = namedtuple("Trans_sheet_name",("transient"))

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
    conductor_definition = Cond_sheet_name(
        files = "CONDUCTOR_files",
        inputs = "CONDUCTOR_input",
        operation = "CONDUCTOR_operation",
        coupling = "CONDUCTOR_coupling",
        ),
    conductor_diagnostic = Diagno_sheet_name(
        space_distr = "Spatial_distribution",
        time_evol = "Time_evolution",
        ),
    conductor_grid = {"GRID"},
    # Key conductor_input is a namedtuple and not a set because I want to 
    # exploit access by field in method conductor_component_instance.
    conductor_input = Comp_sheet_name(
        FluidComponent = "CHAN",
        StackComponent = "STACK",
        StrandMixedComponent = "STR_MIX",
        StrandStabilizerComponent = "STR_STAB",
        JacketComponent = "Z_JACKET",
        ),
    # Key conductor_input is a namedtuple and not a set because I want to 
    # exploit access by field in method conductor_component_instance.
    conductor_operation = Comp_sheet_name(
        FluidComponent = "CHAN",
        StackComponent = "STACK",
        StrandMixedComponent = "STR_MIX",
        StrandStabilizerComponent = "STR_STAB",
        JacketComponent = "Z_JACKET",
        ),
    environment_input = Env_sheet_name(environment = "ENVIRONMENT"),
    transitory_input = Trans_sheet_name(transient = "TRANSIENT"),
)

# Flag for mesh definition.
UNIFORM_MESH = 0
REFINED_MESH = 1
# Adaptive mesh from initial uniform mesh.
ADAPTIVE_UNIFORM_MESH = 2
# Adaptive mesh from initial refined mesh.
ADAPTIVE_REFINED_MESH = 3
MESH_FROM_FILE = -1