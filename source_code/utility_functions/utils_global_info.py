# Valid values for flag IADAPTIME (to set thermal-hydraulic time step behaviour)
# Please, specify new valid values for the flag as shown below.
IADAPTIME_VALUES = (
    -2, # adaptive time step from user defiend function
    -1, # adaptive time step from user defined auxiliary input file
    0, # no adaptive time step (time_step = t_step_min)
    1, # adaptive time step accounting for variation in the whole thermal-hydraulic solution
    2, # adaptive time step accounting only for variation in temperature solution
)

# Dictionary of tuple with valid flag values, used in function 
# check_flag_value. Dictionary keys are the name of the flag as specified in 
# the input file.
VALID_FLAG_VALUES = dict(
    IADAPTIME = IADAPTIME_VALUES, # time step beaviour
)