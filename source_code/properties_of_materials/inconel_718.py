import numpy as np

# Function CONDNE starts here
def thermal_conductivity_inconel718(TT):

    """
    ##############################################################################
    #                  FUNCTION CONDNE(TT)
    ##############################################################################
    #
    # Thermal conductivity of Inconel 718, in W/m K as a function of TT
    # for 1 <= TT <= 500 K
    #
    #                        References
    #                        ----------
    # Cryocomp version 3.0
    #
    # variable    I/O               meaning                        units
    # --------------------------------------------------------------------
    #   TT         x            absolute temperature                  K
    #   CONDNE      x          Thermal conductivity                  W/m K
    #
    #
    # Author : Cryosoft
    # Version: 3.0   March 1997
    #
    ##############################################################################
    # Translation from Fortran to Python: D.Placido PoliTo 10/07/2020
    # Tested against temperature in range [1,500] K: D.Placido PoliTo 11/07/2020
    ##############################################################################
    """

    AA = 15.7823796
    BB = -7.11845509
    CC = -0.49416263
    DD = 0.0017008
    a = 48.6907438
    b = 9.65510385
    c = 0.06978251
    d = -0.54474549
    na = 0.86455662
    nb = 2.14212483
    nc = 2.29137737
    nd = 2.44590736
    TMIN = 1.0
    TMAX = 500.0

    TT = np.array(TT)
    CONDNE = np.zeros(TT.shape)

    TT = np.minimum(TT, TMAX)
    TT = np.maximum(TT, TMIN)

    CONDNE = (
        AA * TT / (a + TT) ** na
        + BB * TT ** 2 / (b + TT) ** nb
        + CC * TT ** 3 / (c + TT) ** nc
        + DD * TT ** 4 / (d + TT) ** nd
    )

    return CONDNE


# Function CPNE starts here
def isobaric_specific_heat_inconel718(TT):

    """
    ##############################################################################
    #                        FUNCTION CPNE(TT)
    ##############################################################################
    #
    # Specific Heat of Inconel 718, in J/Kg K as a function of temperature
    # for 1 <= T <= 1000 K
    #
    #                        References
    #                        ----------
    # Cryocomp version 3.0
    #
    # variable    I/O               meaning                        units
    # --------------------------------------------------------------------
    #   T         x            absolute temperature                  K
    #   CPNE        x          specific heat                       J/Kg K
    #
    #
    # Author : Cryosoft
    # Version: 3.0   March 1997
    #
    ##############################################################################
    # Translation from Fortran to Python: D.Placido PoliTo 10/07/2020
    # Tested against temperature in range [2,1000] K: D.Placido PoliTo 11/07/2020
    ##############################################################################
    """

    T0 = 22.1488673
    a1 = 0.41309481
    a2 = -0.01300475
    a3 = 0.00150658
    a4 = -0.00010202
    a5 = 3.5377e-06
    AA = 33.2554045
    BB = 0.00033495
    CC = -5.80523108
    DD = 53.5383062
    a = 52.1477536
    b = 4.56140518
    c = 7.80826623
    d = 0.32466724
    na = 0.99906362
    nb = 1.99999986
    nc = 2.54414531
    nd = 5.9794021
    TMIN = 1.0
    TMAX = 1000.0

    TT = np.array(TT)
    CPNE = np.zeros(TT.shape)

    TT = np.minimum(TT, TMAX)
    TT = np.maximum(TT, TMIN)

    intervals = [(TT <= T0), (TT > T0)]
    behavior = [
        lambda TT: a1 * TT + a2 * TT ** 2 + a3 * TT ** 3 + a4 * TT ** 4 + a5 * TT ** 5,
        lambda TT: AA * TT ** na / (1 + TT / a) ** na
        + BB * TT ** nb / (1 + TT / b) ** nb
        + CC * TT ** nc / (1 + TT / c) ** nc
        + DD * TT ** nd / (1 + TT / d) ** nd,
    ]
    CPNE = np.piecewise(TT, intervals, behavior)

    return CPNE


# Function RHOENE starts here
def electrical_resistivity_inconel718(TT):

    """
    ##############################################################################
    #             FUNCTION RHOENE(TT)
    ##############################################################################
    #
    # Electrical resistivity of Inconel 718, in Ohm m as a function of T
    # for 1 <= T <= 300 K
    #
    #                        References
    #                        ----------
    # Cryocomp version 3.0
    #
    # variable    I/O               meaning                        units
    # --------------------------------------------------------------------
    #   T         x            absolute temperature                  K
    #   RHOENE      x          Electrical resistivity                Ohm m
    #
    #
    # Author : Cryosoft
    # Version: 3.0   March 1997
    #
    ##############################################################################
    # Translation from Fortran to Python: D.Placido PoliTo 10/07/2020
    # Tested against temperature in range [4,300] K: D.Placido PoliTo 11/07/2020
    ##############################################################################
    """

    A = 11.8
    B = 0.00056648
    nb = 1.18930588
    TMIN = 1.0
    TMAX = 300.0

    TT = np.array(TT)
    RHOENE = np.zeros(TT.shape)

    TT = np.minimum(TT, TMAX)
    TT = np.maximum(TT, TMIN)

    RHOENE = (A + B * TT ** nb) * 1.0e-7

    return RHOENE

def density_inconel718(temperature: np.ndarray)->np.ndarray:
    """
    Function that evaluates inconel 718 density, assumed constant.

    Args:
        temperature (np.ndarray): temperature array, used to get the shape of density array.

    Returns:
        np.ndarray: inconel 718 density array in kg/m^3.
    """
    return 8.17e3 * np.ones(temperature.shape)
