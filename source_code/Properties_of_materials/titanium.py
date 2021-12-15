# Import python libraries and other functions
import numpy as np

# Function CONDTI.py starts here
def thermal_conductivity_ti(TT):

    """
    ##############################################################################
    #             FUNCTION CONDTI(TT)
    ##############################################################################
    #
    # Thermal conductivity of titanium, in W/m K as a function of T
    # for 1.5 <= T <= 300 K
    #
    #                        References
    #                        ----------
    # Cryocomp version 3.0
    #
    # variable    I/O               meaning                        units
    # --------------------------------------------------------------------
    #   T         x            absolute temperature                  K
    #   CONDTI      x          Thermal conductivity                  W/m K
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

    AA = 76813749.1
    BB = -399771.242
    CC = 86.5283574
    DD = -48.3988419
    a = 46.7934973
    b = 17.2222627
    c = 4.56620444
    d = 0.15360579
    na = 3.77004775
    nb = 3.84520961
    nc = 2.663728
    nd = 3.57572912
    TMIN = 1.5
    TMAX = 300.0

    TT = np.array(TT)
    CONDTI = np.zeros(TT.shape)

    TT = np.minimum(TT, TMAX)
    TT = np.maximum(TT, TMIN)

    CONDTI = (
        AA * TT / (a + TT) ** na
        + BB * TT ** 2 / (b + TT) ** nb
        + CC * TT ** 3 / (c + TT) ** nc
        + DD * TT ** 4 / (d + TT) ** nd
    )

    return CONDTI


# Function CPTI.py starts here
def isobaric_specific_heat_ti(TT):

    """
    ##############################################################################
    #                     FUNCTION CPTI(T)
    ##############################################################################
    #
    # Specific Heat of titanium, in J/Kg K as a function of temperature
    # for 1 <= T <= 1000 K
    #
    #                        References
    #                        ----------
    # Cryocomp version 3.0
    #
    # variable    I/O               meaning                        units
    # --------------------------------------------------------------------
    #   T         x            absolute temperature                  K
    #   CPTI        x          specific heat                       J/Kg K
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

    T0 = 20.5360584
    A0 = 0.04795678
    A1 = 0.00238309
    A2 = 0.02173058
    A3 = -0.00239383
    A4 = 0.00016843
    A5 = -2.9575e-06
    AA = -20.2555084
    BB = 0.29234465
    CC = 103.66234
    DD = 19.2186704
    a = 307.498376
    b = 59.346198
    c = 1.15298171
    d = 2.40601976
    na = 8.61469954
    nb = 2.62878227
    nc = 2.88408921
    nd = 8.7579746
    TMIN = 1.0
    TMAX = 1000.0

    TT = np.array(TT)
    CPTI = np.zeros(TT.shape)

    TT = np.minimum(TT, TMAX)
    TT = np.maximum(TT, TMIN)

    intervals = [(TT <= T0), (TT > T0)]
    behavior = [
        lambda TT: A0
        + A1 * TT
        + A2 * TT ** 2
        + A3 * TT ** 3
        + A4 * TT ** 4
        + A5 * TT ** 5,
        lambda TT: AA * TT / (1 + TT / a) ** na
        + BB * TT ** 2 / (1 + TT / b) ** nb
        + CC * TT ** 3 / (1 + TT / c) ** nc
        + DD * TT ** 4 / (1 + TT / d) ** nd,
    ]

    CPTI = np.piecewise(TT, intervals, behavior)

    return CPTI


# Function RHOETI.py starts here
def electrical_resistivity_ti(TT):

    """
    ##############################################################################
    #                     FUNCTION RHOETI(T)
    ##############################################################################
    #
    # Electrical resistivity of titanium, in Ohm m as a function of T
    # for 3 <= T <= 300 K
    #
    #                        References
    #                        ----------
    # Cryocomp version 3.0
    #
    # variable    I/O               meaning                        units
    # --------------------------------------------------------------------
    #   T         x            absolute temperature                  K
    #   RHOETI      x          Electrical resistivity                Ohm m
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

    A = 0.00370359
    B = 0.00169855
    C = -0.00186073
    D = -645.552135
    E = 0.00134204
    n = 1.43164096
    TMIN = 3.0
    TMAX = 300.0

    TT = np.array(TT)
    RHOETI = np.zeros(TT.shape)

    TT = np.minimum(TT, TMAX)
    TT = np.maximum(TT, TMIN)

    intervals = [(TT <= 8.6) | (TT >= 8.7), (TT > 8.6) & (TT < 8.7)]
    behavior = [
        lambda TT: (
            A * TT ** 3 / (1 + B * TT ** 4) + C * TT ** 3 / (D + TT ** 3) + E * TT ** n
        )
        * 1.0e-7,
        lambda TT: 3.828579e-8 + (TT - 8.6) * 2.203794e-8 / 0.1,
    ]
    RHOETI = np.piecewise(TT, intervals, behavior)

    return RHOETI
