# Flags for the current definition

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
# Flag to solve the electric problem in transient conditions.
TRANSIENT_ELECTRIC_SOLVER = 1