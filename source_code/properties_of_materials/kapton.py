# Import python libraries and other functions
import numpy as np

# Function condpl starts here
def thermal_conductivity_kapton(t):
    """
    # Source: Cryosoft
    # Translation from MatLab to Python: S.Poccia UniTo & D.Placido PoliTo 03/2020
    # Correction after test: D.Placido PoliTo 02/04/2020
    # Tested against corresponding MatLab function in temperature range [4,400] K
    """

    t = np.array(t)
    kk = np.zeros(t.shape)  # variable initialization

    aa = 9701.61616
    bb = -0.05967072
    cc = -0.0148391
    dd = 0.00945874
    a = 228.205261
    b = 5.61102451
    c = 0.54581157
    d = 0.66347987
    na = 2.62128508
    nb = 1.90212806
    nc = 2.75642208
    nd = 3.42669585
    tmin = 1.0
    tmax = 833.0
    t = np.minimum(t, tmax)
    t = np.maximum(t, tmin)

    kk = (
        aa * t / (a + t) ** na
        + bb * t ** 2 / (b + t) ** nb
        + cc * t ** 3 / (c + t) ** nc
        + dd * t ** 4 / (d + t) ** nd
    )

    return kk  # end of the function


# Function cppl starts here
def isobaric_specific_heat_kapton(t):
    """
    # Source: Cryosoft
    # Translation from MatLab to Python: S.Poccia UniTo & D.Placido PoliTo 03/2020
    # Tested against corresponding MatLab function in temperature range [4,400] K
    """

    t = np.array(t)
    cp = np.zeros(t.shape)  # variable initialization

    AA = 0.01650911
    BB = 0.00011155
    CC = 2.60292523
    DD = -3.87742905
    a = 247.518453
    b = 8.23106376
    c = 1.65332503
    d = 1.47379679

    na = 2.26547298
    nb = 7.60272577
    nc = 6.18461619
    nd = 6.44133927
    TMIN = 1.0
    TMAX = 300.0

    t = np.minimum(t, TMAX)
    t = np.maximum(t, TMIN)
    cp = (
        AA * t ** na / (1 + t / a) ** na
        + BB * t ** nb / (1 + t / b) ** nb
        + CC * t ** nc / (1 + t / c) ** nc
        + DD * t ** nd / (1 + t / d) ** nd
    )

    return cp  # end of the function

def density_kapton(nn:int)->np.ndarray:
    """
    Function that evaluates kapton density, assumed constant.

    Reference: CRYOSOFT
    Args:
        nn (int): number of elements of the array.

    Returns:
        np.ndarray: kapton density array in kg/m^3.
    """
    return 1.380e3 * np.ones(nn)

def electrical_resistivity_kapton(nn:int)->np.ndarray:
    """
    Function that evaluates kapton electrical resistivity, assumed constant.

    N.B: no reference for this value.

    Args:
        nn (int): number of elements of the array.

    Returns:
        np.ndarray: kapton electrical resistivity array in Ohm*m.
    """
    return 1.5e15 * np.ones(nn)