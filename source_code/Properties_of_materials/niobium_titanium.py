# Importing python libraries and other funcions
import numpy as np

# Function bcnbti starts here
def bcnbti(T, TC0, BC20):
    """
    ######################################################################
    #
    # Upper critical field, in T, for NbTi as a function of field.
    #
    #            References
    #            ----------
    # M.S.Lubell, Scaling formulas for critical current and critical field
    # for commercial NbTi, IEEE Trans. Mag. ,19, (1983).
    #
    # variable  I/O         meaning            units
    # --------------------------------------------------------------------
    #   T     x      magnetic field            K
    #   TC0     x      upper critical temperature (B=0)    K
    #   BC20    x      upper critical field   (T=0)      T
    #   BCNBTI_old    x      critical field            T
    #
    # Other functions called: NONE
    #
    # Author : L.Bottura @ CryoSoft
    # Version: 3.0  16.3.1996
    #
    # Translation from MatLab to Python: S.Poccia UniTo & D.Placido PoliTo 03/2020
    ######################################################################
    #    IMPLICIT NONE
    #    REAL   T,TC0,BC20
    # *
    #    REAL   TLIM
    # *
    #    REAL   N,ALPHA,BETA,GAMMA,TLOW,BLOW
    #    COMMON   /NBTI/N,ALPHA,BETA,GAMMA,TLOW,BLOW
    #
    # * SET THE LOWER LIMIT FOR THE TEMPERATURE
    #    TLIM=AMAX1(T,TLOW)

    # %crb PFCI (2008/2009)
    # % N=parameter

    # %crb DEMO (2014)
    """

    N = 1.7
    T = np.array(T)
    BCNBTI = np.zeros(T.shape)
    index = np.nonzero(T <= TC0)
    if index[0].size == 0:  # empty index array (cdp, 06/2020)
        return BCNBTI

    BCNBTI[index[0]] = BC20 * (
        1.0 - (T[index[0]] / TC0) ** N
    )  # always > 0 (cdp, 06/2020)

    return BCNBTI  # end of the function


# Function CONDNBTI starts here
def thermal_conductivity_nbti(TT):

    """
    ######################################################################
    #                     FUNCTION CONDNBTI(TT)
    ######################################################################
    #
    # Thermal conductivity of NbTi. As the data is scarcely available, a
    # crude fit from White and Woods, (1957) of the K of Nb only has been
    # used. Errors bounds unknown. For 1 <= TT <= 750.
    #
    #                        References
    #                        ----------
    # G.K.White, S.B.Woods, Can. J. Physics, 35, 892-900, (1957)
    #
    # variable    I/O               meaning                        units
    # --------------------------------------------------------------------
    #   TT         x            absolute temperature                  K
    #   CONDNBTI      x          thermal conductivity                W/m K
    #
    # Author : Cryosoft
    # Version: 3.0   March 1997
    #
    ##############################################################################
    # Translation from Fortran to Python: D.Placido PoliTo 10/07/2020
    # Tested against temperature in range [4,750] K: D.Placido PoliTo 11/07/2020
    ##############################################################################
    """

    A = 6.5
    Tm = 25.0
    Tnb = 50.0
    K0 = 21.5
    TMIN = 1.0
    TMAX = 750.0

    TT = np.array(TT)
    CONDNBTI = np.zeros(TT.shape)

    TT = np.minimum(TT, TMAX)
    TT = np.maximum(TT, TMIN)

    intervals = [(TT < Tnb), (TT >= Tnb)]
    behavior = [
        lambda TT: A * TT / (1.0 + (TT / Tm) ** 3),
        lambda TT: A * TT / (1.0 + (TT / Tm) ** 3)
        + K0 * (1.0 - np.exp(-(TT - Tnb) / Tnb)),
    ]

    CONDNBTI = np.piecewise(TT, intervals, behavior)

    return CONDNBTI


# Function CPNBTI starts here
def isobaric_specific_heat_nbti(T, B, TCS, TC):
    """
    ######################################################################
    #
    # Specific Heat of NbTi in J/Kg K as a function of temperature and
    # field for T <= 1000 K and 0 < B < 7 T.  A combination of the
    # superconducting and normal Cp's is used in the current sharing re-
    # gime ( Tcs < T < Tc ). For B > 7 T the coefficients of B = 7 T
    # are used (approximation), while for T > 20 K a fit is used based
    # on predicted data from the Debye function (Debye T = 222 K).
    # For T > 1000 K a constant Cp is assumed.
    # NOTE : in the range 0...20 K the fit from V.Arp, (1980) is used
    # directly. This fit takes into account the superconducting transition
    #
    #            References
    #            ----------
    # V.D. Arp, Stability and Thermal Quenches in Force-Cooled Supercon-
    # ducting Cables, Superconducting MHD Magnet Design Conference, MIT,
    # pp 142-157, 1980.
    # G.K.White, S. B. Woods, Can.J.Phys, 35, 892-900, 1957
    #
    # variable    I/O      meaning            units
    # --------------------------------------------------------------------
    #   T      x      absolute temperature         K
    #   TCS    x      current sharing temperature     K
    #   TC     x      critical temperature         K
    #   CPNBTI   x    specific heat           J/Kg K
    #
    # Author : Cryosoft
    # Version: 3.0   March 1997
    #
    # Translation from MatLab to Python: S.Poccia UniTo & D.Placido PoliTo 03/2020
    ######################################################################
    """

    T = np.array(T)
    B = np.array(B)

    # T and B should have the same shape
    if len(T) == 1:
        T = T * np.ones(B.shape)
    if len(B) == 1:
        B = B * np.ones(T.shape)

    # variable initialization (cdp, 06/2020)
    CPN = np.zeros(T.shape)
    CPS = np.zeros(T.shape)
    F = np.ones(T.shape)
    CP = np.zeros(T.shape)

    AA = 17.971311
    BB = -3215.18087
    CC = 4325.5737
    DD = -711.296119
    a = 88.768525
    b = 3.66034399
    c = 38.0453266
    d = 14.67806
    na = 1.0456476
    nb = 51.7204943
    nc = 4.20661339
    nd = 6.01104574
    TMIN = 1.0
    TMAX = 1000.0

    X = np.array(
        [
            0.000000e01,
            0.200000e01,
            0.300000e01,
            0.400000e01,
            0.500000e01,
            0.600000e01,
            0.700000e01,
        ]
    )
    AB = np.array(
        [
            0.341000e-01,
            0.282900e-01,
            0.219200e-01,
            0.158600e-01,
            0.993000e-02,
            0.491000e-02,
            0.152000e-02,
        ]
    )
    E = np.array(
        [
            0.233300e01,
            0.244600e01,
            0.259700e01,
            0.280600e01,
            0.310700e01,
            0.356800e01,
            0.434800e01,
        ]
    )

    I = np.zeros(B.shape, dtype=int)
    # For each value of B find index for which abs(B) >= X and take the maximum one into account for following evaluation
    for ii in range(len(B)):
        JJ = np.nonzero(abs(B[ii]) >= X[1 : len(X)])
        if len(JJ[0]) > 0:
            I[ii] = (
                max(JJ[0]) + 1
            )  # since we extract the index for a portion of X of length len(X) - 1

    # Find index in T such that T <= TMIN
    ind_TMINa = np.nonzero(T <= TMIN)
    # Find index T such that T[ind_TMINa[0]] > TC
    IND = np.nonzero(T[ind_TMINa[0]] > TC[ind_TMINa[0]])
    ind_TCa = ind_TMINa[0][IND[0]]
    # CP evaluation
    CP[ind_TCa] = (0.161) * T[ind_TCa] + (0.00279) * (T[ind_TCa] ** 3)
    # Find index T such that T[ind_TMINa[0]] < TCS
    IND = np.nonzero(T[ind_TMINa[0]] < TCS[ind_TMINa[0]])
    ind_TCSa = ind_TMINa[0][IND[0]]
    # CP evaluation
    CP[ind_TCSa] = AB[I[ind_TCSa]] * (abs(T[ind_TCSa]) ** E[I[ind_TCSa]])
    # Find index T such that T[ind_TMINa[0]] <= TC and T[ind_TMINa[0]] >= TCS
    IND = np.nonzero(
        (T[ind_TMINa[0]] <= TC[ind_TMINa[0]]) & (T[ind_TMINa[0]] >= TCS[ind_TMINa[0]])
    )
    ind = ind_TMINa[0][IND[0]]
    # Evaluate CPN and CPS for these case
    CPS[ind] = AB[I[ind]] * (abs(T[ind]) ** E[I[ind]])
    CPN[ind] = (0.161) * T[ind] + (0.00279) * (T[ind] ** 3)
    # find index in ind such that TCS < TC
    IND = np.nonzero(TCS[ind] < TC[ind])
    ind_TCS = ind[IND[0]]
    # Evaluate F for these index; otherwise it is set to 1.0 by initialization
    F[ind_TCS] = (T[ind_TCS] - TCS[ind_TCS]) / (TC[ind_TCS] - TCS[ind_TCS])
    CP[ind] = F[ind] * CPN[ind] + (1.0 - F[ind]) * CPS[ind]
    # Find index in T such that T > TMIN
    ind_TMINb = np.nonzero(T > TMIN)
    T[ind_TMINb[0]] = np.minimum(T[ind_TMINb[0]], TMAX)
    CP[ind_TMINb[0]] = (
        AA * (T[ind_TMINb[0]] / (a + T[ind_TMINb[0]])) ** na
        + BB * (T[ind_TMINb[0]] / (b + T[ind_TMINb[0]])) ** nb
        + CC * (T[ind_TMINb[0]] / (c + T[ind_TMINb[0]])) ** nc
        + DD * (T[ind_TMINb[0]] / (d + T[ind_TMINb[0]])) ** nd
    )

    return CP  # end of the function


# Function jcnbti starts here
def critical_current_density_nbti(T, B, TC0, BC20, C0):
    """
    ######################################################################
    #
    # Critical (non-copper) current density, in A/m**2, for NbTi as a
    # function of temperature and field. A combined B and T dependence
    # has been fitted to the data from the source below
    #
    #                        References for data
    #                        -------------------
    # M.A. Green, Calculating the Jc, B, T Surface for Niobium Titanium
    # Using a Reduced State Model, IEEE Trans. Mag., 25, 2, (1989).
    #
    # G. Morgan, A Comparison of Two Analytic Forms for the Jc(B,T)
    # surface, SSC Magnet Division Notes, 310-1 (SSC-MD-218), (1989)
    #
    # variable    I/O               meaning                        units
    # --------------------------------------------------------------------
    #   T         x            absolute temperature                  K
    #   B         x            magnetic field                        T
    #   TC0       x            upper critical temperature (B = 0)      K
    #   BC20      x            upper critical field   (T = 0)          T
    #   C0        x            Jc normalization constant          A T/m**2
    #   JCNBTI      x          critical current density            A/m**2
    #
    # Other functions called: BCNBTI
    #
    # Author : L.Bottura @ CryoSoft
    # Version: 3.0  9.7.1996
    #
    # Translation from MatLab to Python: D.Placido PoliTo 03/2020
    ######################################################################
    # IMPLICIT NONE
    # REAL     T,B,TC0,BC20,C0
    #
    # REAL     BC2,BCNBTI,TLCASE,BLCASE,TLIM,BLIM
    #
    # REAL     N,ALPHA,BETA,GAMMA,TLOW,BLOW
    # COMMON   /NBTI/N,ALPHA,BETA,GAMMA,TLOW,BLOW
    """

    # crb PFCI (2008/2009)
    # ALPHA = 1.69
    # BETA = 1.91
    # GAMMA = 2.13
    # N = parameter

    # crb DEMO (2014)
    ALPHA = 1.0
    BETA = 1.54
    GAMMA = 2.1
    N = 1.7

    T = np.array(T)  # is always a vector (cdp, 06/2020)
    B = np.array(B)  # is always a vector (cdp, 06/2020)

    if len(T) == 1:
        T = T * np.ones(B.shape)
    if len(B) == 1:
        B = B * np.ones(T.shape)

    JCNBTI = np.zeros(B.shape)  # JCNBTI initialization
    BC2 = -1 * np.ones(
        B.shape
    )  # initialization to negative value since it will be always > 0
    TLCASE = np.zeros(B.shape)
    BLCASE = np.zeros(B.shape)

    if TC0 < 1.0e-6:
        raise ValueError("ERROR> From JcNbTi\nTc0 = 0.0\nSTOP jcnbti")

    if BC20 < 1.0e-6:
        raise ValueError("ERROR> From JcNbTi\nBc0 = 0.0\nSTOP jcnbti")

    if C0 < 1.0e-6:
        raise ValueError("ERROR> From JcNbTi\nc0 = 0.0\nSTOP jcnbti")

    B = np.abs(B)

    # Find element index in T < TC0 (cdp, 06/2020)
    T_ind = np.nonzero(T < TC0)  # this is a tuple (cdp, 06/2020)
    if T_ind[0].size == 0:  # empty index array (cdp, 06/2020)
        return JCNBTI
    # Evaluate BC2 only for those index, elsewere its value is -1. by initialization (cdp, 06/2020)
    BC2[T_ind[0]] = bcnbti(T[T_ind[0]], TC0, BC20)
    # Find element index in B[T_ind] such that B[T_ind] < BC2[T_ind], only for those index JCNBTI will be evaluated, elsewere JCNBTI = 0.0 by initialization (cdp, 06/2020)
    ind = np.nonzero(B[T_ind[0]] < BC2[T_ind[0]])  # this is a tuple (cdp, 06/2020)
    B_ind = T_ind[0][ind[0]]  # this is an array (cdp, 06/2020)
    if B_ind.size == 0:  # empty index array (cdp, 06/2020)
        return JCNBTI
    # JCNBTI evaluation (cdp, 06/2020)
    TLCASE[B_ind] = T[B_ind] / TC0
    BLCASE[B_ind] = B[B_ind] / BC2[B_ind]
    JCNBTI[B_ind] = (
        C0
        / B[B_ind]
        * BLCASE[B_ind] ** ALPHA
        * (1.0 - BLCASE[B_ind]) ** BETA
        * (1.0 - TLCASE[B_ind] ** N) ** GAMMA
    )

    return JCNBTI  # end of the function


# Function TCNBTI starts here
def critical_temperature_nbti(B, TC0, BC20):
    """
    ######################################################################
    #    REAL FUNCTION TCNBTI(B     ,TC0   ,BC20  )
    ######################################################################
    #
    # Critical temperature, in K, for NbTi as a function of field.
    #
    #                        References
    #                        ----------
    # M.S.Lubell, Scaling formulas for critical current and critical field
    # for commercial NbTi, IEEE Trans. Mag. ,19, (1983).
    #
    # variable    I/O               meaning                        units
    # --------------------------------------------------------------------
    #   B         x            magnetic field                        T
    #   TC0       x            upper critical temperature (B=0)      K
    #   BC20      x            upper critical field   (T=0)          T
    #   TCNBTI      x          critical temperature                  K
    #
    # Other functions called: NONE
    #
    # Author : CryoSoft
    # Version: 3.0  16.3.1996
    #
    # Translation from MatLab to Python: D.Placido PoliTo 03/2020
    ######################################################################
    """
    # USE SCALING !crb (April 28, 2015)
    # ALPHA = 1.0
    # BETA = 1.54
    # GAMMA = 2.1
    N_NbTi = 1.7
    BLOW = 0.01

    B = np.array(B)
    TC = np.zeros(B.shape)  # variable initialization

    # !crb Add check and error message (January 15, 2015)
    if TC0 < 1e-6:
        raise ValueError("ERROR> From TcNbTi:\nTc0 = 0.0!\nSTOP TCNBTI")

    # !crb Add check and error message (January 15, 2015)
    if BC20 < 1e-6:
        raise ValueError("ERROR> From TcNbTi:\nBc0 = 0.0!\nSTOP TCNBTI")

    BLIM = np.maximum(abs(B), BLOW)
    index = np.nonzero(BLIM <= BC20)  # this is a tuple (cdp, 06/2020)
    if index[0].size == 0:  # empty index array (cdp, 06/2020)
        return TC
    TC[index[0]] = TC0 * (1.0 - BLIM[index[0]] / BC20) ** (1.0 / N_NbTi)

    return TC  # end of the function


# Function TCSNBTI starts here
def current_sharing_temperature_nbti(B, JOP, TC0, BC20, C0):
    """
    ######################################################################
    #    REAL FUNCTION TCSNBTI(B     ,JOP   ,TC0   ,BC20  ,C0    ,ICOND)
    ######################################################################
    #
    # Current sharing temperature, in K, for NbTi as a function of field
    # and operating current density. This function is the inverse of the
    # Jc function of JCNBTI
    #
    # variable    I/O               meaning                        units
    # --------------------------------------------------------------------
    #   B         x            magnetic field                        T
    #   JOP       x            operating current density           A/m**2
    #   TC0       x            upper critical temperature (B = 0)      K
    #   BC20      x            upper critical field   (T = 0)          T
    #   C0        x            Jc normalization constant         A T/m**2
    #   TCSNBTI      x          current sharing temperature           K
    #
    # Other functions called: JCNBTI,TCNBTI
    #
    # Author : CryoSoft
    # Version: 3.0  15.7.1996
    #
    # Translation from MatLab to Python: D.Placido PoliTo 03/2020
    ######################################################################
    """

    B = np.array(B)

    # variable initialization
    Tcs = np.zeros(B.shape)
    JC = np.zeros(B.shape)
    TC = np.zeros(B.shape)
    BLIM = np.zeros(B.shape)
    BLCASE = np.zeros(B.shape)
    JLCASE = np.zeros(B.shape)
    TLCASE = np.zeros(B.shape)
    TLUPPR = np.zeros(B.shape)
    TLLOWR = np.zeros(B.shape)
    R = np.zeros(B.shape)

    # USE SCALING !crb (April 28, 2015)

    ALPHA = 1.0
    BETA = 1.54
    GAMMA = 2.1
    N_NbTi = 1.7
    # TLOW = 0.01
    BLOW = 0.01
    TOLRNC = 1.0e-5

    # C * FUNCTIONS DEFINITION
    def BN(TT, BB, nnn):

        # nnn is N_NbTi
        bn = BB / (1 - TT ** nnn)
        return bn

    def RSDL(TT, BB, JJ, aa, bb, cc, nnn):

        # aa is ALPHA
        # bb is BETA
        # cc is GAMMA
        # nnn is N_NbTi

        rsdl = JJ - (BN(TT, BB, nnn) ** aa) * ((1 - BN(TT, BB, nnn)) ** bb) * (
            (1 - TT ** nnn) ** cc
        )
        return rsdl

    # !crb Add check and error message (January 15, 2015)
    if TC0 < 1e-6:
        raise ValueError("ERROR > From TcsNTi: Tc0 = 0!")

    # !crb Add check and error message (January 15, 2015)
    if BC20 < 1e-6:
        raise ValueError("ERROR > From TcsNTi: Bc20 = 0!")

    # !crb Add check and error message (January 15, 2015)
    if C0 < 1e-6:
        raise ValueError("ERROR > From TcsNTi: C0 = 0!")

    # Find element index in B such that B < BC20 (cdp, 06/2020)
    B_ind = np.nonzero(B < BC20)  # this is a tuple (cdp, 06/2020)
    if B_ind[0].size == 0:
        return Tcs
    # Evaluate JC @ T = 0 only for elements corresponding to B_ind, elsewere JC = 0.0 by initialization (cdp, 06/2020)
    JC[B_ind[0]] = critical_current_density_nbti(
        np.zeros(B[B_ind[0]].shape), B[B_ind[0]], TC0, BC20, C0
    )
    # Find element index in JC[B_ind] such that JC[B_ind] < JOP (cdp, 06/2020)
    ind = np.nonzero(JC[B_ind[0]] > JOP)  # this is a tuple (cdp, 06/2020)
    JC_ind = B_ind[0][ind[0]]  # this is an array (cdp, 06/2020)
    if JC_ind.size == 0:
        return Tcs
    # Evaluate TC only for B[JC_ind]
    TC[JC_ind] = critical_temperature_nbti(B[JC_ind], TC0, BC20)
    if JOP <= 0.0:
        Tcs = TC
        return Tcs

    # Following evaluation is only for survivor index, i.e JC_ind; elsewere Tcs = 0.0 (cdp, 06/2020)
    # C * SET THE LOWER LIMIT FOR THE FIELD
    BLIM[JC_ind] = np.maximum(abs(B[JC_ind]), BLOW)

    # C * NORMALISED FIELD
    BLCASE[JC_ind] = BLIM[JC_ind] / BC20

    # C * NORMALISED CURRENT DENSITY
    JLCASE[JC_ind] = JOP / C0 * BLIM[JC_ind]

    # C * UPPER AND LOWER LIMITS FOR THE ITERATIVE ROOT SEARCH
    TLUPPR[JC_ind] = TC[JC_ind] / TC0

    # C * FIND THE NORMALISED TEMPERATURE TCS/TC0 BY ITERATION
    for ii in range(len(JC_ind)):
        CONVER = np.array([False])
        TLLOWR = 0.0
        while CONVER == False:
            # C * ROOT SEARCH BY STRAIGHT BISECTION (MOST STABLE WAY)
            TLCASE[JC_ind[ii]] = 0.5 * (TLLOWR + TLUPPR[JC_ind[ii]])

            # C * COMPUTE RESIDUAL
            R = RSDL(
                TLCASE[JC_ind[ii]],
                BLCASE[JC_ind[ii]],
                JLCASE[JC_ind[ii]],
                ALPHA,
                BETA,
                GAMMA,
                N_NbTi,
            )

            # C * DECISION PROCESS
            if R > TOLRNC:
                TLUPPR[JC_ind[ii]] = TLCASE[JC_ind[ii]]
            elif R < -TOLRNC:
                TLLOWR = TLCASE[JC_ind[ii]]

            # C * CHECK CONVERGENCE
            CONVER = (abs(R) <= TOLRNC) + ((TLUPPR[JC_ind[ii]] - TLLOWR) <= TOLRNC)
        # end while
    # end for

    # C * CONVERT NORMALISED TO ABSOLUTE TEMPERATURE
    Tcs[JC_ind] = TLCASE[JC_ind] * TC0

    return Tcs  # end function


# Function rho_NbTi starts here
def density_nbti():
    """
    NbTi density kg/m^3. It is assumed constant.
    Autor: D. Placido Polito 21/01/2021
    """
    return 6160.0


# end function rho_NbTi (cdp, 01/2021)
