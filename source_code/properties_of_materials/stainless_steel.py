# Import python libraries and other functions
import numpy as np

# Function condss.py starts here
def thermal_conductivity_ss(t):
    """
    # Source: Cryosoft
    # Translation from MatLab to Python: S.Poccia UniTo & D.Placido PoliTo 03/2020
    # Tested against corresponding MatLab function in temperature range [4,400] K
    """

    t = np.array(t)
    kk = np.zeros(t.shape)  # variable initialization

    AA = 21224.0426
    BB = 103.007623
    CC = 135.870511
    DD = 0.07942859
    a = 10844822.8
    b = 6800262.44
    c = 34.4411414
    d = 0.23194142
    na = 1.07801291
    nb = 1.88806818
    nc = 3.44233222
    nd = 3.20629457

    TT = t  # why?? comment Placido Daniele

    kk = (
        AA * TT / (a + TT) ** na
        + BB * TT**2 / (b + TT) ** nb
        + CC * TT**3 / (c + TT) ** nc
        + DD * TT**4 / (d + TT) ** nd
    )

    return kk  # end of the function


# Function cpss.py starts here
def isobaric_specific_heat_ss(T):
    """
    # Source: Cryosoft
    # Translation from MatLab to Python: S.Poccia UniTo & D.Placido PoliTo 03/2020
    # Correction after test: D.Placido PoliTo 03/04/2020
    # Tested against corresponding MatLab function in temperature range [4,400] K
    """

    T = np.array(T)
    cp = np.zeros(T.shape)  # variable initialization

    T0 = 12.8973011

    A1 = 0.46437387
    A2 = 0.00035157
    AA = -468.772387
    BB = -81.2353649
    CC = 827.209072
    DD = 29.9811823

    a = 422.117262
    b = 1235598.51
    c = 51.2732976
    d = 0.44591988
    na = 0.57149148
    nb = 2.0006633
    nc = 2.79384889
    nd = 3.25864797

    TMIN = 1.0
    TMAX = 1253.0

    if max(T) <= TMIN:
        cp = 0.498 * T + 3.71e-4 * T**3
    elif max(T) < T0:
        T = np.minimum(T, TMAX)
        cp = A1 * T + A2 * T**3
    else:
        T = np.minimum(T, TMAX)
        cp = (
            AA * T / (a + T) ** na
            + BB * T**2 / (b + T) ** nb
            + CC * T**3 / (c + T) ** nc
            + DD * T**4 / (d + T) ** nd
        )

    return cp  # end of the function


# Function rhoess.py starts here
def electrical_resistivity_ss(t):
    """
    # Source: Cryosoft
    # Translation from MatLab to Python: D.Placido PoliTo 03/2020
    # Tested against corresponding MatLab function in temperature range [4,400] K
    """

    t = np.array(t)
    RHOESS = np.array(t.shape)  # variable initialization

    A = 4.7589144
    B = 0.0099003
    NB = 0.95536674
    RHOESS = (A + B * t**NB) * 1.0e-7
    return RHOESS


# Function rho_ss starts here
def density_ss(nn: int) -> np.nodarrayF:

    """
    Function that evaluates stainless steel density, assumed constant.

    Args:
        nn (int): number of elements of the array.

    Returns:
        np.ndarray: stainless steel density array in kg/m^3.
    """
    return 7800.0 * np.ones(nn)


# end function rho_ss (cdp, 01/2021)
