# Import python libraries and other functions
import numpy as np

# Function condal starts here
def thermal_conductivity_al(TT):

    """
    ######################################################################
    #
    # Thermal conductivity of Al alloy (RRR=5.5) in W/m K, as a function
    # of temperature T, for 2 <= T <= 1000 K.
    #
    #            References
    #            ----------
    # P.Reed and A.F.Clark, Materials at low Temperature, ASM, 1983 (hand
    # picked-up data points, but good accuracy)
    #
    # variable  I/O         meaning            units
    # --------------------------------------------------------------------
    #   T     x      absolute temperature          K
    #   CONDAL    x      thermal conductivity        W/m K
    #
    # Author : R.Heller @ ITP
    # Version: 2.0   10.10.1991
    #
    # Translation from MatLab to Python: S.Poccia UniTo & D.Placido PoliTo 03/2020
    # Correction after test: D.Placido PoliTo 01/04/2020
    # Tested against corresponding MatLab function in temperature range [4,400] K
    ######################################################################
    """

    TT = np.array(TT)
    TMIN = 2.0
    TMAX = 1000.0
    k_Al = np.zeros(TT.shape)
    if (min(TT) < TMIN) or (max(TT) > TMAX):

        print(
            f"ERROR>Temperature outside valid range: T >= {TMIN} K and T <= {TMAX} K.\n"
        )

    else:

        intervals = [
            (TT >= 2) & (TT <= 4.0),
            (TT > 4.0) & (TT <= 6.0),
            (TT > 6.0) & (TT <= 8.0),
            (TT > 8.0) & (TT <= 10.0),
            (TT > 10.0) & (TT <= 20.0),
            (TT > 20.0) & (TT <= 40.0),
            (TT > 40.0) & (TT <= 60.0),
            (TT > 60.0) & (TT <= 80.0),
            (TT > 80.0) & (TT <= 100.0),
            (TT > 100.0) & (TT <= 200.0),
            (TT > 200.0) & (TT <= 300.0),
            (TT > 300.0) & (TT <= 1000.0),
        ]

        behavior = [
            lambda TT: 4.2345 * (TT - 2.0) + 6.6510,
            lambda TT: 4.2000 * (TT - 4.0) + 15.1200,
            lambda TT: 4.3650 * (TT - 6.0) + 23.5200,
            lambda TT: 4.8550 * (TT - 8.0) + 32.2500,
            lambda TT: 4.3430 * (TT - 10.0) + 41.9600,
            lambda TT: 3.8805 * (TT - 20.0) + 85.3900,
            lambda TT: 2.0375 * (TT - 40.0) + 163.000,
            lambda TT: 0.5360 * (TT - 60.0) + 203.75,
            lambda TT: 0.1145 * (TT - 80.0) + 214.47,
            lambda TT: -0.0452 * (TT - 100.0) + 216.76,
            lambda TT: -0.0291 * (TT - 200.0) + 212.24,
            lambda TT: 209.3300,
        ]

        k_Al = np.piecewise(TT, intervals, behavior)

    return k_Al


# Function cpal starts here
def isobaric_specific_heat_al(TT):
    """
    ######################################################################
    #
    # Specific Heat of Al, in J/Kg K as a function of temperature
    # T, for 2 <= T <= 1000 K.
    #
    #            References
    #            ----------
    # Cryocomp version 3.0
    # V.J. Johnson, Properties of materials at low temperatures (phase 1)
    # Pergamon Press, 1961
    #
    # variable  I/O         meaning            units
    # --------------------------------------------------------------------
    #   T     x      absolute temperature          K
    #   CPAL    x      specific heat             J/Kg K
    #
    #
    # Author : R.Heller @ ITP
    # Version: 2.0   10.10.1991
    #
    # Translation from MatLab to Python: S.Poccia UniTo & D.Placido PoliTo 03/2020
    # Correction after test: D.Placido PoliTo 01/04/2020
    # Tested against corresponding MatLab function in temperature range [4,400] K

    #######################################################################
    """

    TT = np.array(TT)
    cp_Al = np.zeros(TT.shape)  # variables initialization

    T0 = 11.0578562
    T1 = 53.2685055
    a1 = 0.04257809
    a2 = 0.00340313
    a3 = 0.00063472
    b0 = -2.50413265
    b1 = 0.71747923
    b2 = -0.06270275
    b3 = 0.0031683
    b4 = -2.0154e-05
    AA = 11458.0127
    BB = -438.603857
    CC = 331.057118
    DD = -10001.7006
    a = -7.12302329
    b = -34.6373149
    c = -29.4649313
    d = -3.30673113
    na = 0.9531781
    nb = 2.00135991
    nc = 2.974805
    nd = 3.94384164
    TMIN = 2.0
    TMAX = 1000.0

    if (min(TT) < TMIN) or (max(TT) > TMAX):
        print(
            f"ERROR>Temperature outside valid range: T >= {TMIN} K and T <= {TMAX} K.\n"
        )

    else:

        intervals = [
            (TT >= TMIN) & (TT <= T0),
            (TT > T0) & (TT <= T1),
            (TT > T1) & (TT <= TMAX),
        ]

        behaviour = [
            lambda TT: (a1 * TT + a2 * TT ** 2 + a3 * TT ** 3),
            lambda TT: (b0 + b1 * TT + b2 * TT ** 2 + b3 * TT ** 3 + b4 * TT ** 4),
            lambda TT: (
                AA * TT / (a + TT) ** na
                + BB * TT ** 2 / (b + TT) ** nb
                + CC * TT ** 3 / (c + TT) ** nc
                + DD * TT ** 4 / (d + TT) ** nd
            ),
        ]

        cp_Al = np.piecewise(TT, intervals, behaviour)

    return cp_Al  # end of the function


# Function RHOEAL starts here
def electrical_resistivity_al(TT):

    """
    ##############################################################################
    #                 FUNCTION RHOEAL(TT)
    ##############################################################################
    #
    # Electrical resistivity of Al in Ohm m, as a function of
    # temperature T, for 2 <= T <= 1000 K.
    #
    #                        References
    #                        ----------
    # P.Reed & A.F.Clark, Materials at low Temperature, ASM, 1983 (hand
    # picked-up data points, but good accuracy)
    #
    #
    # variable    I/O               meaning                        units
    # --------------------------------------------------------------------
    #   T         x            absolute temperature                  K
    #   RHOEAL      x          resistivity                         Ohm m
    #
    # Author : R.Heller @ ITP
    # Version: 2.0   10.10.1991
    #
    ##############################################################################
    # Translation from Fortran to Python: D.Placido PoliTo 10/07/2020
    # Tested against temperature in range [2,1000] K: D.Placido PoliTo 11/07/2020
    ##############################################################################
    """

    TT = np.array(TT)
    RHOEAL = np.zeros(TT.shape)
    X = np.array(
        [
            0.200000e1,
            0.400000e1,
            0.600000e1,
            0.800000e1,
            0.100000e2,
            0.200000e2,
            0.400000e2,
            0.600000e2,
            0.800000e2,
            0.100000e3,
            0.200000e3,
            0.300000e3,
            0.500000e3,
            0.100000e4,
        ]
    )
    Y = np.array(
        [
            0.576000e-8,
            0.576000e-8,
            0.576000e-8,
            0.576000e-8,
            0.576000e-8,
            0.576000e-8,
            0.600000e-8,
            0.720000e-8,
            0.912000e-8,
            0.112800e-7,
            0.230400e-7,
            0.350400e-7,
            0.584800e-7,
            0.116800e-6,
        ]
    )
    C = np.array(
        [
            [
                -0.341170e-14,
                0.170585e-14,
                -0.341170e-14,
                0.119409e-13,
                -0.443521e-13,
                0.472520e-12,
                0.332536e-10,
                0.825132e-10,
                0.104684e-09,
                0.110712e-09,
                0.120792e-09,
                0.118922e-09,
                0.115284e-09,
            ],
            [
                0.255877e-14,
                0.123260e-30,
                -0.255877e-14,
                0.102351e-13,
                -0.383816e-13,
                0.900689e-13,
                0.154898e-11,
                0.913994e-12,
                0.195038e-12,
                0.105852e-12,
                -0.505247e-14,
                -0.136421e-13,
                -0.454738e-14,
            ],
            [
                -0.426462e-15,
                -0.426462e-15,
                0.213231e-14,
                -0.810279e-14,
                0.428168e-14,
                0.243153e-13,
                -0.105832e-13,
                -0.119826e-13,
                -0.148644e-14,
                -0.369682e-15,
                -0.286322e-15,
                0.151579e-15,
                0.151579e-15,
            ],
        ]
    )

    if (min(TT) < X[0]) or (max(TT) > X[-1]):
        raise ValueError(
            f"ERROR>Temperature outside valid range (T>= {X[0]} K & T <= {X[-1]} K).\n"
        )
    else:
        # no need of function NBISEC in this way

        intervals = [
            ((TT >= X[0]) & (TT <= X[1])),
            ((TT > X[1]) & (TT <= X[2])),
            ((TT > X[2]) & (TT <= X[3])),
            ((TT > X[3]) & (TT <= X[4])),
            ((TT > X[4]) & (TT <= X[5])),
            ((TT > X[5]) & (TT <= X[6])),
            ((TT > X[6]) & (TT <= X[7])),
            ((TT > X[7]) & (TT <= X[8])),
            ((TT > X[8]) & (TT <= X[9])),
            ((TT > X[9]) & (TT <= X[10])),
            ((TT > X[10]) & (TT <= X[11])),
            ((TT > X[11]) & (TT <= X[12])),
            ((TT > X[12]) & (TT <= X[13])),
        ]

        behavior = [
            lambda TT: (np.sum(C[:, 0]) * (TT - X[0]) + Y[0]),
            lambda TT: (np.sum(C[:, 1]) * (TT - X[1]) + Y[1]),
            lambda TT: (np.sum(C[:, 2]) * (TT - X[2]) + Y[2]),
            lambda TT: (np.sum(C[:, 3]) * (TT - X[3]) + Y[3]),
            lambda TT: (np.sum(C[:, 4]) * (TT - X[4]) + Y[4]),
            lambda TT: (np.sum(C[:, 5]) * (TT - X[5]) + Y[5]),
            lambda TT: (np.sum(C[:, 6]) * (TT - X[6]) + Y[6]),
            lambda TT: (np.sum(C[:, 7]) * (TT - X[7]) + Y[7]),
            lambda TT: (np.sum(C[:, 8]) * (TT - X[8]) + Y[8]),
            lambda TT: (np.sum(C[:, 9]) * (TT - X[9]) + Y[9]),
            lambda TT: (np.sum(C[:, 10]) * (TT - X[10]) + Y[10]),
            lambda TT: (np.sum(C[:, 11]) * (TT - X[11]) + Y[11]),
            lambda TT: (np.sum(C[:, 12]) * (TT - X[12]) + Y[12]),
        ]

        RHOEAL = np.piecewise(TT, intervals, behavior)

    return RHOEAL


# Function rho_Al starts here
def density_al():
    """
    Aluminium density kg/m^3. It is assumed constant.
    Autor: D. Placido Polito 21/01/2021
    """
    return 2700.0


# end function rho_Al (cdp, 01/2021)
