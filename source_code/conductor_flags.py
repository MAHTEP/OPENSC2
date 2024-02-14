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