# Import python libraries and other functions
import numpy as np

# Function condin starts here
def thermal_conductivity_incoloy908(t):
    """
    # Source: Cryosoft
    # Translation from MatLab to Python: S.Poccia UniTo & D.Placido PoliTo 03/2020
    # Correction after test: D.Placido PoliTo 02/04/2020
    # Tested against corresponding MatLab function in temperature range [4,400] K
    """

    t = np.array(t)
    kk = np.zeros(t.shape)  # variable initialization
    aa = 644.74407
    bb = -141.184107
    cc = -0.85851077
    dd = 0
    a = 5154.07845
    b = 439.708273
    c = 1.93083487
    d = 0
    na = 1
    nb = 2
    nc = 3
    nd = 4

    tmin = 1.0
    tmax = 1000.0
    t = np.minimum(t, tmax)
    t = np.maximum(t, tmin)

    kk = (
        aa * t / (a + t) ** na
        + bb * t ** 2 / (b + t) ** nb
        + cc * t ** 3 / (c + t) ** nc
        + dd * t ** 4 / (d + t) ** nd
    )

    return kk  # end of the function


# Function cpin starts here
def isobaric_specific_heat_incoloy908(t):
    """
    # function cp = cpin(t);
    #
    ######################################################################
    # NOTE: This routine has been vectorised for MATLAB!!
    # The input t MUST be a row vector t = [t1 t2 ... tn]
    # A. Portone @ EFDA, Garching, 5.12.2000
    ######################################################################
    # Specific Heat of INCOLOY 908, in J/Kg K as a function of temperature
    # For 4.27 <= T <= 1423 K.
    #
    #            References
    #            ----------
    # L.S.Toma, M.M.Steeves, R.P.Reed, Incoloy Alloy 908 Data Handbook,
    # PFC/RR-94-2, 1994
    #
    # variable  I/O         meaning            units
    # --------------------------------------------------------------------
    #   t     x      absolute temperature          K
    #   cp      x      specific heat             J/Kg K
    #
    #
    # Author : L.Bottura & C. Rosso at Cryosoft
    # Version: 3.0   March 1997
    #
    # Translation from MatLab to Python: S.Poccia UniTo & D.Placido PoliTo 03/2020
    # Tested against corresponding MatLab function in temperature range [4,400] K
    ######################################################################
    """

    t = np.array(t)
    cp = np.zeros(t.shape)  # variable initialization

    aa = 0.00219167
    bb = 0.1333227
    cc = 0.00018546
    dd = -4.8176e-05
    a = 4813.28489
    b = 10.8112775
    c = 58.178214
    d = 136.765771
    na = -0.22447346
    nb = 1.41030517
    nc = 3.84780106
    nd = 3.30329979
    tmin = 4.27
    tmax = 1423.0

    t = np.minimum(t, tmax)
    t = np.maximum(t, tmin)

    num1 = aa * np.power(t, na)  # aa * t. ^ na
    den1 = np.power((1 + t / a), na)  # (1 + t / a). ^ na
    num2 = bb * np.power(t, nb)  # bb * t. ^ nb
    den2 = np.power((1 + t / b), nb)  # (1 + t / b). ^ nb
    num3 = cc * np.power(t, nc)  # cc * t. ^ nc
    den3 = np.power((1 + t / c), nc)  # (1 + t / c). ^ nc
    num4 = dd * np.power(t, nd)  # dd * t. ^ nd
    den4 = np.power((1 + t / d), nd)  # (1 + t / d). ^ nd

    cp = (
        np.divide(num1, den1)
        + np.divide(num2, den2)
        + np.divide(num3, den3)
        + np.divide(num4, den4)
    )

    return cp  # end of the function


# Function rhoein starts here
def electrical_resistivity_incoloy908(t):
    """
    #
    # function rho=rhoein(t)
    #
    ######################################################################
    # NOTE: This routine has been vectorised for MATLAB!!
    # The input t MUST be a row vector t = [t1 t2 ... tn]
    # A. Portone @ EFDA, Garching, 5.12.2000
    ######################################################################
    #
    # Electrical resistivity of INCOLOY 908, in Ohm m as a function of T
    # NOTE: the resistivity has been assumed to be equal to that of
    # inconel. Accuracy unknown. For 1 <= T <= 300 K.
    #
    #
    #                        References
    #                        ----------
    #
    # variable    I/O               meaning                        units
    # --------------------------------------------------------------------
    #   t         x            absolute temperature                  K
    #   rho         x          Electrical resistivity                Ohm m
    #
    #
    # Author : L.Bottura & C. Rosso at Cryosoft
    # Version: 3.0   March 1997
    #
    # Translation from MatLab to Python: D.Placido PoliTo 03/2020
    # Correction after test: D.Placido PoliTo 02/04/2020
    # Tested against corresponding MatLab function in temperature range [4,400] K
    ######################################################################
    """

    t = np.array(t)
    rho = np.array(t.shape)

    a = 11.8
    b = 0.00056648
    nb = 1.18930588
    tmin = 1.0
    tmax = 300.0
    t = np.minimum(t, tmax)
    t = np.maximum(t, tmin)

    rho = (a + b * t ** nb) * 1e-7

    return rho  # end of the function
