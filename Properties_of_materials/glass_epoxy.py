# Import python libraries and other functions
import numpy as np

# Function condge starts here
def thermal_conductivity_ge(t):
    """
    # Source: Cryosoft
    # Translation from MatLab to Python: S.Poccia UniTo & D.Placido PoliTo 03/2020
    # Correction after test: D.Placido PoliTo 02/04/2020
    # Tested against corresponding MatLab function in temperature range [4,400] K
    """

    t = np.array(t)
    kk = np.zeros(t.shape)  # variable initialization

    aa = 0.106635834
    bb = 4.254937095
    cc = 0.369710568
    dd = -0.220612963
    a = 4.374557267
    b = 78.36483015
    c = 0.896107679
    d = 0.261916360
    na = 0.551227580
    nb = 2.274402853
    nc = 2.979875183
    nd = 3.644905343

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

    return kk


# Function cpge starts here
def isobaric_specific_heat_ge(t):
    """
    # Source: Cryosoft
    # Translation from MatLab to Python: S.Poccia UniTo & D.Placido PoliTo 03/2020
    # Correction after test: D.Placido PoliTo 02/04/2020
    # Tested against corresponding MatLab function in temperature range [4,400] K
    """

    # %glass-epoxy (check tmax-tmin)

    t = np.array(t)
    cp = np.zeros(t.shape)  # variable initialization

    AA = 0.309399374
    BB = -1387.509141
    CC = 283.0320589
    DD = 1.822153805
    A = 113.1687618
    B = 80.88333136
    C = 26.27075203
    D = 1.104229335
    NA = 0.254475819
    NB = 1.938523579
    NC = 2.603466143
    ND = 3.522809815
    TMIN = 1.0
    TMAX = 1000.0

    t = np.minimum(t, TMAX)
    t = np.maximum(t, TMIN)

    cp = (
        AA * t / (A + t) ** NA
        + (BB * t ** 2) / (B + t) ** NB
        + CC * t ** 3 / (C + t) ** NC
        + DD * t ** 4.0 / (D + t) ** ND
    )

    return cp  # end of the function


# Function rho_ge starts here
def density_ge():
    """
    Glass-epoxy density kg/m^3. It is assumed constant.
    Autor: D. Placido Polito 21/01/2021
    """
    return 2000.0


# end function rho_ge (cdp, 01/2021)
