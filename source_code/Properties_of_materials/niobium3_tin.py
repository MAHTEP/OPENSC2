# Importing python libraries and other funcions
import numpy as np
from scipy import optimize
import warnings

# Function BCNBSN starts here
def BCNBSN(T, EPSLON, TC0M, BC20M):

    """
    ######################################################################
    #
    # [BC]=BCNBSN(T,EPSLON,TC0M,BC20M)
    #
    # Critical field, in T, for Nb3Sn as a function of temperature and
    # strain. The Nb3Sn material is characterized by the parameters
    # Tc0m and Bc20m.
    #
    #                        References
    #                        ----------
    # L. Bottura, Jc(B,T,epsilon) Parametrizations for the ITER Nb3Sn Production,# April 2, 2008
    #
    #
    # variable    I/O               meaning                        units
    # --------------------------------------------------------------------
    #   T         x            absolute temperature                  K
    #   EPSLON    x            strain                                -
    #   TC0M      x            critical temperature (B=0)            K
    #   BC20M     x            upper critical field (T=0)            T
    #   BCNBSN      x          Upper critical field                  T
    #
    # Other functions called: NONE
    #
    # Author : L.Bottura @ CERN
    # Version: 2  2.4.2008
    #
    # Translation from MatLab to Python: S.Poccia UniTo & D.Placido PoliTo 03/2020
    ######################################################################
    """

    T = np.array(T)
    EPSLON = np.array(EPSLON)
    if T.size == 1:
        T = T * np.ones(EPSLON.shape)
    if EPSLON.size == 1:
        EPSLON = EPSLON * np.ones(T.shape)

    BC = np.zeros(T.shape)  # variable initialization

    # %*NORMALISED TEMPERATURE T/TC0
    TLCASE = T / (TC0M * (SNBSN(EPSLON)) ** (1.0e0 / 3.0e0))
    index = np.nonzero(TLCASE < 1.0)[0]  # this is a numpy array (cdp, 06/2020)
    if index.size == 0:  # empty index array (cdp, 06/2020)
        return BC

    # CRITICAL FIELD BC2(T,EPSLON)
    BC[index] = (
        BC20M * SNBSN(EPSLON[index]) * (1.0e0 - TLCASE[index] ** 1.52e0)
    )

    return BC


# Function CONDNBSN starts here
def thermal_conductivity_nb3sn(TT):

    """
    ##############################################################################
    #                 FUNCTTION CONDNBSN(TT)
    ##############################################################################
    #
    # Thermal conductivity of Nb3Sn as a function of temperature TT, for
    # 4 <= TT <= 750 K. NOTE: the conductivity is assumed to be continuous at
    # the superconducting to normal state transition, mainly because of
    # lack of directly measured data. The data has been picked-up by hand,
    # but is accurate below 100 K. Over 100 K, they are only indicative
    # (the values have been extrapolated from the t < 100 K curve). Above
    # 300 K K is assumed constant, although I have no evidence of this.
    #
    #                        References
    #                        ----------
    # Handbook on Materials for S.C. Machinery, NBS Boulder (yellow book),
    # 1977
    # V.D.Arp, Stability and Thermal Quenches in Force-Cooled
    # Superconducting Cables, Superconductiong MHD Magnet Design
    # Conference, MIT, 142-157, 1980
    # G.S.Knapp, S.D.Bader, Z.Fisk, Phonon Properties of A-15
    # Superconductors Obtained from Heat Capacity Measurements, Phys. rev.
    # B, 13, No. 9, 3783-3789, 1976
    #
    # variable    I/O               meaning                        units
    # --------------------------------------------------------------------
    #   TT         x            absolute temperature                  K
    #   CONDNBSN      x          thermal conductivity                W/m K
    #
    # Author : Cryosoft
    # Version: 3.0   March 1997
    #
    ##############################################################################
    # Translation from Fortran to Python: D.Placido PoliTo 10/07/2020
    # Tested against temperature in range [4,750] K: D.Placido PoliTo 11/07/2020
    ##############################################################################
    """

    T0 = 21.9620369
    A = 1.3905e14
    B = 76.6127237
    C = -7.1614565e-1
    D = 5.40229689
    n = 3.66343961
    m = 9.3220843

    TT = np.array(TT)
    CONDNBSN = np.zeros(TT.shape)

    TMIN = 4.0
    TMAX = 750.0

    TT = np.minimum(TT, TMAX)
    TT = np.maximum(TT, TMIN)

    intervals = [(TT <= T0), (TT > T0)]
    behavior = [lambda TT: A * TT ** n / (B + TT) ** m, lambda TT: C * np.log(TT) + D]
    CONDNBSN = np.piecewise(TT, intervals, behavior)

    return CONDNBSN


# Function CPNBSN starts here
def isobaric_specific_heat_nb3sn(T, TCS, TC, TC0):
    """
    # REAL FUNCTION CPNBSN(T,TCS,TC,TC0)
    ######################################################################
    #
    # Specific Heat of binary Nb3Sn in J/Kg K as a function of temperature
    # field and strain for  T <= 1000 K.  A combination of
    # the superconducting and normal Cp's is used in the current sharing
    # regime ( Tcs < T < Tc ). An analytic expression is used below 20 K
    # while a spline fit is used over 20 K. For T < 20 K in the normal
    # state two different expressions are used below and over 10 K
    # For T > 400 K a constant Cp is assumed (Dulong-Petit rule)
    #
    #            References
    #            ----------
    #  0.. 20 K: V.D. Arp, Stability and Thermal Quenches in Force-Cooled
    #      Superconducting Cables, Superconducting MHD Magnet Design
    #      Conference, MIT, pp 142-157, 1980.
    # 20..400 K: G.S.Knapp,S.D.Bader,Z.Fisk, Phonon properties of A-15
    #      superconductors obtained from heat capacity measurements,
    #      Phys.Rev B, 13, no.9, pp 3783-3789, 1976.
    #
    # variable      I/O         meaning                              units
    # --------------------------------------------------------------------
    # T             x           absolute temperature                 K
    # TCS           x           current sharing temperature          K
    # TC            x           critical temperature Tc(B)           K
    # TC0           x           critical temperature at B=0 Tc(0)    K
    # CPNBSN        x           specific heat                        J/Kg K
    #
    # Other functions called: NONE
    #
    # Author : L.Bottura @ NET
    # Version: 2.0  23.11.1994
    #
    # Translation from MatLab to Python: S.Poccia UniTo & D.Placido PoliTo 03/2020
    ######################################################################
    """

    T = np.array(T)
    TCS = np.array(TCS)
    TC = np.array(TC)

    if T.size == 1:
        T = T * np.ones(TC.shape)

    CPN = np.zeros(T.shape)
    CPS = np.zeros(T.shape)
    DBC2DT = np.zeros(T.shape)
    DELCP = np.zeros(T.shape)
    CPNTC = np.zeros(T.shape)
    TC_norm = np.zeros(T.shape)
    F = np.ones(T.shape)
    CPNBSN = np.zeros(T.shape)

    AA = 38.2226877
    BB = -848.364226
    CC = 1415.13808
    DD = -346.837966
    a = 6.80458608
    b = 59.9209182
    c = 25.8286334
    d = 8.77918335
    na = 1
    nb = 2
    nc = 3
    nd = 4
    TMIN = 20.0
    TMAX = 1000.0

    T = np.minimum(T, TMAX)
    intervals_CPN = [(T <= 10.0), (T > 10.0) & (T <= 20.0), (T > 20.0)]
    behavior_CPN = [
        lambda T: 7.5475e-3 * T ** 2,
        lambda T: (-0.3 + 0.00375 * T ** 2) / 0.09937,
        lambda T: AA * T / (a + T) ** na
        + BB * T ** 2 / (b + T) ** nb
        + CC * T ** 3 / (c + T) ** nc
        + DD * T ** 4 / (d + T) ** nd,
    ]

    # NORMAL COMPONENT OF CP
    # piecewise fuction CPN
    CPN = np.piecewise(T, intervals_CPN, behavior_CPN)

    #### added on 28/03/2020 by Placido Daniele: translation from CPNBSN.m

    intervals_CPNTC = [(TC <= 10.0), (TC > 10.0) & (TC <= 20.0)]
    behavior_CPNTC = [
        lambda TC: 7.5475e-3 * TC ** 2,
        lambda TC: (-0.3 + 0.00375 * TC ** 2) / 0.09937,
    ]
    CPNTC = np.piecewise(TC, intervals_CPNTC, behavior_CPNTC)

    # SUPERCONDUCTING COMPONENT OF CP
    # Find index in T such that T <= TC
    ind_TCa = np.nonzero(T <= TC)[0]
    # For those index evaluate TC/TC0 ratio, elsewere it is = 0 by initialization
    TC_norm[ind_TCa] = TC[ind_TCa] / TC0
    # For those index evaluate DBC2DT and DELCP, elsewere they are = 0 by initialization

    DBC2DT[ind_TCa] = -0.46306 - 0.067830 * TC[ind_TCa]
    DELCP[ind_TCa] = (
        1500.0
        * (DBC2DT[ind_TCa] ** 2)
        / (2.0 * (27.2 / (1.0 + (0.34 * TC_norm[ind_TCa]))) ** 2 - 1.0)
    )
    # For those index evaluate piecewise fuction CPNTC, elsewere it is = 0 by initialization
    # CPNTC[ind_TCa] = np.piecewise(TC[ind_TCa], intervals_CPNTC, behavior_CPNTC)
    # For those index evaluate CPS, elsewere it is = 0 by initialization
    CPS[ind_TCa] = (CPNTC[ind_TCa] + DELCP[ind_TCa]) * (
        T[ind_TCa] / TC[ind_TCa]
    ) ** 3

    # COMPOSE THE COMPONENTS
    # Find index in T such that T <= TCS
    ind_TCSa = np.nonzero(T <= TCS)[0]  # this is a numpy array (cdp, 06/2020)
    # For those index evaluate CPNBSN
    CPNBSN[ind_TCSa] = CPS[ind_TCSa]
    # Find index in T such that T > TCS and t <= TC
    ind_inter = np.nonzero((T > TCS) & (T <= TC))[0]
    # Find index in TCS such that TCS < TC
    IND = np.nonzero(TCS[ind_inter] < TC[ind_inter])  # this is a numpy array (cdp, 06/2020)
    ind = ind_inter[IND]  # this is an array (cdp, 06/2020)
    # For those index evaluate F, elsewere it is = 1 by initialization
    F[ind] = (T[ind] - TCS[ind]) / (TC[ind] - TCS[ind])
    # For index in ind_inter evaluate CPNBSN
    CPNBSN[ind_inter] = (
        F[ind_inter] * CPN[ind_inter]
        + (1.0 - F[ind_inter]) * CPS[ind_inter]
    )
    # Find index in T such that T > TC
    ind_TCb = np.nonzero(T > TC)[0]  # this is a numpy array (cdp, 06/2020)
    CPNBSN[ind_TCb] = CPN[ind_TCb]

    return CPNBSN  # end of the function


# Function JCNBSN starts here
def critical_current_density_nb3sn(T, B, EPSLON, TC0M, BC20M, C0):
    """
    ######################################################################
    #
    # [JC] = JCNBSN(T,B,EPSLON,TC0M,BC20M,C0)
    #
    # Critical (non-copper) current density, in A/m**2, for Nb3Sn as a
    # function of temperature, field and strain. The Nb3Sn material is
    # characterized by the parameters Tc0m and Bc20m. The constant C0
    # determines the overall scaling of the Jc curve. At the moment the
    # low field value is obtained by using a lower boundary on B/Bc2.
    # Note that a lower limit is used on field and temperature to avoid
    # the 0-singularity.
    #
    #                        References
    #                        ----------
    # L. Bottura, Jc(B,T,epsilon) Parametrizations for the ITER Nb3Sn Production, # April 2, 2008
    #
    #
    # variable    I/O               meaning                        units
    # --------------------------------------------------------------------
    #   T         x            absolute temperature                  K
    #   B         x            magnetic field                        T
    #   EPSLON    x            strain                                -
    #   TC0M      x            critical temperature (B = 0)            K
    #   BC20M     x            upper critical field (T = 0)          T
    #   C0        x            normalization constant             A T/m**2
    #   JCNBSN      x          critical current density            A/m**2
    #
    # Other functions called: NONE
    #
    # Author : L.Bottura @ CERN
    # Version: 2  2.4.2008
    #
    # Translation from MatLab to Python: D.Placido PoliTo 03/2020
    ######################################################################
    """

    # ppp = 0.56
    # qqq = 1.75
    # Actual used values
    ppp = 0.63
    qqq = 2.1

    # from file scaling_input.dat (cdp, 10/2020)
    # ppp = 0.84 # [] real - low field exponent of the pinning force (ppp~0.5)
    # qqq = 2.57 # [] real - high field exponent of the pinning force (qqq~2)

    BLOW = 0.01

    T = np.array(T)
    B = np.array(B)
    EPSLON = np.array(EPSLON)

    if T.size == 1:
        T = T * np.ones(B.shape)
    if B.size == 1:
        B = B * np.ones(T.shape)

    JC = np.zeros(T.shape)  # JC initialization
    BLCASE = np.zeros(T.shape)

    # SET THE LOWER LIMIT FOR THE FIELD
    BLIM = np.maximum(B, BLOW)

    # NORMALISED TEMPERATURE T/TC0
    TLCASE = T / (TC0M * (SNBSN(EPSLON)) ** (1.0e0 / 3.0e0))

    # Find element index in TLCASE such that TLCASE < 1.0 (cdp, 06/2020)
    TLCASE_ind = np.nonzero(TLCASE < 1.0)[0]  # this is a numpy array (cdp, 06/2020)
    if TLCASE_ind.size == 0:  # empty index array (cdp, 06/2020)
        return JC

    # NORMALISED FIELD B/BC0
    # Evaluate BLCASE only for those index, elsewere its value is 0.0 by initialization (cdp, 06/2020)
    BLCASE[TLCASE_ind] = BLIM[TLCASE_ind] / (
        BCNBSN(T[TLCASE_ind], EPSLON[TLCASE_ind], TC0M, BC20M)
    )
    # Find element index such that BLCASE[TLCASE_ind[0]] < 1.0, only for those index JC will be evaluated, elsewere JC = 0.0 by initialization (cdp, 06/2020)
    ind = np.nonzero(BLCASE[TLCASE_ind] < 1.0)[0]  # this is a numpy array (cdp, 06/2020)
    BLCASE_ind = TLCASE_ind[ind]  # this is an array (cdp, 06/2020)
    if BLCASE_ind.size == 0:  # empty index array (cdp, 06/2020)
        return JC

    # JCNBTI evaluation (cdp, 06/2020)
    # JC(T, B, EPSLON)
    JC[BLCASE_ind] = (
        C0
        / BLIM[BLCASE_ind]
        * SNBSN(EPSLON[BLCASE_ind])
        * (1.0e0 - TLCASE[BLCASE_ind] ** 1.52e0)
        * (1.0e0 - TLCASE[BLCASE_ind] ** 2e0)
        * BLCASE[BLCASE_ind] ** ppp
        * (1 - BLCASE[BLCASE_ind]) ** qqq
    )
    # Find element index such that JC[BLCASE_ind] < 0, for these index JC = 0.0 (cdp, 2020)
    ind = np.nonzero(JC[BLCASE_ind] < 0.0)[0]  # this is a numpy array (cdp, 06/2020)
    JC_ind = BLCASE_ind[ind]  # this is an array (cdp, 06/2020)
    JC[JC_ind] = 0.0

    return JC  # end of the function


# Function SNBSN starts here
def SNBSN(EPSLON):

    """
    ######################################################################
    #
    # [SS]=SNBSN(EPSLON)
    #
    # Strain function for Nb3Sn as a function of strain.
    #
    #            References
    #            ----------
    # L. Bottura, Jc(B,T,epsilon) Parametrizations for the ITER Nb3Sn Production, # April 2, 2008
    #
    #
    # variable  I/O         meaning            units
    # --------------------------------------------------------------------
    #   EPSLON  x      strain                -
    #   SNBSN     x      strain function             -
    #
    # Other functions called: NONE
    #
    # Author : L.Bottura @ CERN
    # Version: 2  2.4.2008
    #
    # Translation from MatLab to Python: D.Placido PoliTo 03/2020
    ######################################################################
    """

    # Ca1 = 45.46
    # Ca2 = 6.52
    # EPS0a = 0.00244
    # Ca1 = 44.48 #[] real - strain fitting constant
    # Ca2 = 0.0 #[] real - strain fitting constant
    # EPS0a = 0.256e-2 #[] real - residual strain component

    EPSLON = np.array(EPSLON)
    SS = np.zeros(EPSLON.shape)  # variable initialization

    # Actual used values
    Ca1 = 45.74  # [] real - strain fitting constant
    Ca2 = 4.431  # [] real - strain fitting constant
    EPS0a = 0.232e-2  # [] real - residual strain component

    # from file scaling_input.dat (cdp, 10/2020)
    # Ca1 = 47.02 # [] real - strain fitting constant
    # Ca2 = 11.76 # [] real - strain fitting constant
    # EPS0a = 2.31e-3 #[] real - residual strain component
    # EPSm = 3.97e-3 #[] real - tensile strain at which the maximum critical

    # %*ESSE(EPSsh) !crb (July 24, 2011)
    EPSsh = Ca2 * EPS0a / np.sqrt(Ca1 ** 2.0e0 - Ca2 ** 2.0e0)

    # %*ESSE(EPSLON) !crb (July 24, 2011)
    SS = 1.0e0 + 1.0e0 / (1.0e0 - Ca1 * EPS0a) * (
        Ca1
        * (
            np.sqrt(EPSsh ** 2.0e0 + EPS0a ** 2.0e0)
            - np.sqrt((EPSLON - EPSsh) ** 2.0e0 + EPS0a ** 2.0e0)
        )
        - Ca2 * EPSLON
    )

    return SS


# Function TCNBSN starts here
def critical_temperature_nb3sn(B, EPSLON, TC0M, BC20M):
    """
    ######################################################################
    #
    # [TC]=TCNBSN(B,EPSLON,TC0M,BC20M)
    #
    # Critical temperature, in K, for Nb3Sn as a function of field and
    # strain. The Nb3Sn material is characterized by the parameters Tc0m
    # and Bc20m.
    #
    #                        References
    #                        ----------
    # L. Bottura, Jc(B,T,epsilon) Parametrizations for the ITER Nb3Sn Production, # April 2, 2008
    #
    #
    # variable    I/O               meaning                        units
    # --------------------------------------------------------------------
    #   B         x            magnetic field                        T
    #   EPSLON    x            strain                                -
    #   TC0M      x            critical temperature (B=0)            K
    #   BC20M     x            upper critical field (T=0)            T
    #   TCNBSN      x          critical temperature                  K
    #
    # Other functions called: NONE
    #
    # Author : L.Bottura @ CERN
    # Version: 2  2.4.2008
    #
    # Translation from MatLab to Python: D.Placido PoliTo 03/2020
    ######################################################################
    """

    B = np.array(B)
    EPSLON = np.array(EPSLON)
    TC = np.zeros(max(B.shape, EPSLON.shape))  # initialization
    BLOW = 0.01

    # *SET THE LOWER LIMIT FOR THE FIELD
    BLIM = np.maximum(B, BLOW)

    # *NORMALISED FIELD
    BLCASE = BLIM / (BC20M * SNBSN(EPSLON))

    index = np.nonzero(BLCASE < 1.0)[0]  # this is a numpy array (cdp, 06/2020)
    if index.size == 0:  # empty index array (cdp, 06/2020)
        return TC
    TC[index] = (
        TC0M
        * (SNBSN(EPSLON[index])) ** (1.0e0 / 3.0e0)
        * (1.0e0 - BLCASE[index]) ** (1.0 / 1.52)
    )

    return TC  # end of the function


# Function TCSNSN starts here
def current_sharing_temperature_nb3sn(B, EPSLON, JOP, TC0M, BC20M, C):
    """
    ######################################################################
    #
    # [TCS] = TCSNSN(B,EPSLON,JOP,TC0M,BC20M,C)
    #
    # Current sharing temperature, in K, for Nb3Sn as a function of field
    # and strain. The Nb3Sn material is characterized by the parameters
    # Tc0m and Bc20m. The critical temperature is computed by iterative
    # inversion of the Jc(T,B,epslon) relation. Note that a minimum value
    # is used for the field and temperature to avoid the 0-singularity
    #
    #                        References
    #                        ----------
    # L. Bottura, Jc(B,T,epsilon) Parametrizations for the ITER Nb3Sn Production, 	April 2, 2008
    #
    #
    # variable    I/O               meaning                        units
    # --------------------------------------------------------------------
    #   B         x            magnetic field                        T
    #   EPSLON    x            strain                                -
    #   JOP       x            operating current density           A/m**2
    #   TC0M      x            critical temperature (B = 0)            K
    #   BC20M     x            upper critical field (T = 0)            T
    #   C0        x            normalization constant             A T/m**2
    #   TCSNSN      x          current sharing temperature           K
    #
    # Other functions called: JCNBSN
    #
    # Author : L.Bottura @ CERN
    # Version: 2  2.4.2008
    #
    # Translation from MatLab to Python: D.Placido PoliTo 03/2020
    ######################################################################
    """

    def critical_current_density_bisection_nb3sn(TT, BB, EPSLON, JOP, TC0M, BC20M, C):
        return critical_current_density_nb3sn([TT], [BB], [EPSLON], TC0M, BC20M, C) - JOP

    # ppp = 0.56
    # qqq = 1.75
    # ppp = 0.63
    # qqq = 2.1
    # Actual used values
    # ppp = 0.556
    # qqq = 1.698

    # from file scaling_input.dat (cdp, 10/2020)
    # ppp = 0.84 # [] real - low field exponent of the pinning force (ppp~0.5)
    # qqq = 2.57 # [] real - high field exponent of the pinning force (qqq~2)
    BLOW = 0.01

    B = np.array(B)
    if B.size == 1:
        B = B * np.ones(EPSLON.shape)

    # variable initialization
    TCS = np.zeros(B.shape)
    JC = np.zeros(B.shape)
    # * SET THE LOWER LIMIT FOR THE FIELD
    BLIM = np.maximum(B, BLOW)

    Bstar = BLIM / BC20M / SNBSN(EPSLON)
    # * CHECK THAT THE FIELD IS BELOW THE UPPER CRITICAL VALUE
    # Find element index in such that Bstar < 1.0 (cdp, 06/2020)
    Bstar_ind = np.nonzero(Bstar < 1.0)[0]  # this is a numpy array (cdp, 06/2020)
    if Bstar_ind.size == 0:  # empty index array (cdp, 06/2020)
        return TCS
    # Evaluate JC @ T = 0 only for elements corresponding to Bstar_ind, elsewere JC = 0.0 by initialization (cdp, 06/2020)
    JC[Bstar_ind] = critical_current_density_nb3sn(
        np.zeros(B[Bstar_ind].shape),
        B[Bstar_ind],
        EPSLON[Bstar_ind],
        TC0M,
        BC20M,
        C,
    )
    # * CHECK THAT JOP IS BELOW THE UPPER CRITICAL VALUE
    # Find element index such that 0 < JOP < JC[Bstar_ind] (cdp, 06/2020)
    ind = np.nonzero((JOP[Bstar_ind] < JC[Bstar_ind]) & (JOP[Bstar_ind] > 0.0))[0]  # this is a numpy array (cdp, 06/2020)
    JC_ind = Bstar_ind[ind]  # this is an array (cdp, 06/2020)
    if JC_ind.size == 0:
        return TCS

    # FIND CURRENT SHARING TEMPERATURE WITH BISECTION (cdp, 09/2020)

    T_upper = critical_temperature_nb3sn(B[JC_ind], EPSLON[JC_ind], TC0M, BC20M)

    for ii, vv in enumerate(ind):

        T_lower = np.array([3.5])

        ex_args = (
            B[vv], EPSLON[vv], JOP[vv], TC0M, BC20M, C
        )
        # Evaluate current sharing temperature with bisection method.
        TCS[vv] = optimize.bisect(
            critical_current_density_bisection_nb3sn,
            T_lower,
            T_upper[ii],
            ex_args,
            xtol=1e-5,
        )
    # End for ii.

    return TCS  # end of the function


# Function rho_Nb3Sn starts here
def density_nb3sn():
    """
    Nb3Sn density kg/m^3. It is assumed constant.
    Autor: D. Placido Polito 21/01/2021
    """
    return 8950.0


# end function rho_Nb3Sn (cdp, 01/2021)
