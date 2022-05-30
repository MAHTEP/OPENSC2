import numpy as np
from typing import Union

def custom_current_function(time: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """User defined custom function for the current beavior in time (and maybe in space).

    Args:
        time (Union[float, np.ndarray]): time at which evaluate the current

    Returns:
        Union[float, np.ndarray]: current value at given time
    """
    # current customizable by the user.
    CURRENT_AMPLITUDE = 1.0  # [A] current amplitude
    FREQUENCY = 50.0  # [Hz] frequency
    return CURRENT_AMPLITUDE * np.cos(2 * np.pi * FREQUENCY * time)