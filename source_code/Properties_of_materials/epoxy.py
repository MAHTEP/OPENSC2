# Import python libraries and other functions
import numpy as np

# Function condep starts here
def thermal_conductivity_ep(t):
    """
    # Source: Cryosoft
    # Translation from MatLab to Python: S.Poccia UniTo & D.Placido PoliTo 03/2020
    # Correction after test: D.Placido PoliTo 02/04/2020
    # Tested against corresponding MatLab function in temperature range [4,400] K
    """

    t = np.array(t)
    kk = np.zeros(t.shape)  # variable initialization

    aa = -0.37683848
    bb = 100.999824
    cc = 0.01286509
    dd = 0.09211400
    a = 4483.59742
    b = 6977179.28
    c = 9.45111414
    d = 0.43241200
    na = 0.48088594
    nb = 2.64828981
    nc = 2.10623669
    nd = 3.89389714
    tmin = 1.0
    tmax = 833.0

    t = np.minimum(t, tmax)
    t = np.maximum(t, tmin)

    kk = (
        aa * t / (a + t) ** na
        + (bb * t ** 2) / (b + t) ** nb
        + cc * t ** 3 / (c + t) ** nc
        + dd * t ** 4 / (d + t) ** nd
    )

    return kk


# Function cpep starts here
def isobaric_specific_heat_ep(t):
    """
    # Source: Cryosoft
    # Translation from MatLab to Python: S.Poccia UniTo & D.Placido PoliTo 03/2020
    # Correction after test: D.Placido PoliTo 02/04/2020
    # Tested against corresponding MatLab function in temperature range [4,400] K
    """

    # %epoxy (check tmax-tmin)
    # % bf coefficients of the polynomial expression of cpep

    t = np.array(t)
    cp = np.zeros(t.shape)  # variable initialization

    AA = 10.1978618
    BB = 0.00021728
    CC = 0.00760138
    DD = 1.8199e-06
    A = 0.68320779
    B = 14.8392381
    C = 1.92507754
    D = 96.6186651
    NA = 7.03975410
    NB = 5.72616511
    NC = 14.7550871
    ND = 4.62647203
    TMIN = 1.0
    TMAX = 833.0

    t = np.minimum(t, TMAX)
    t = np.maximum(t, TMIN)

    cp = (
        AA * (t / (1 + t / A)) ** NA
        + BB * (t / (1 + t / B)) ** NB
        + CC * (t / (1 + t / C)) ** NC
        + DD * (t / (1 + t / D)) ** ND
    )

    return cp  # end of the function
