# Importing python libraries and other funcions
import numpy as np
import pandas as pd
from scipy import optimize
import warnings


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


# Function rho_NbTi starts here
def density_nbti(nn: int) -> np.ndarray:
    """
    Function that evaluates NbTi density, assumed constant.

    Args:
        nn (int): number of elements of the array.

    Returns:
        np.ndarray: NbTi density array in kg/m^3.
    """
    return 6160.0 * np.ones(nn)


# end function rho_NbTi (cdp, 01/2021)


def reduced_temperature_nbti(temperature, T_c0):

    return temperature / T_c0


# End function reduced_temperature_nbti


def _convert_to_nparray(value):

    if isinstance(value, (np.ndarray, pd.DataFrame, pd.Series)):
        return value
    else:
        # Necessary to deal with scalar: conversion to numpy array
        return np.array([value])


# End function convert_to_nparray


def critical_magnetic_field_nbti(temperature, B_c20, T_c0, nn=1.7):

    temperature = _convert_to_nparray(temperature)
    critical_mag_field = np.zeros(temperature.shape)
    # Find index such that temperature <= maximum critical temperature
    ind = temperature <= T_c0
    # Evaluate critical magnetic field only for those values of temperature, elsewhere critical magnetic field is 0 by initialization.
    critical_mag_field[ind] = B_c20 * (
        1.0 - reduced_temperature_nbti(temperature[ind], T_c0) ** nn
    )
    return critical_mag_field


# End function critical_magnetic_field_nbti


def reduced_magnetic_field_nbti(magnetic_field, temperature, B_c20, T_c0, nn=1.7):

    return magnetic_field / critical_magnetic_field_nbti(
        temperature, B_c20, T_c0, nn=nn
    )


# End function reduced_magnetic_field_nbti


def critical_temperature_nbti(magnetic_field, B_c20, T_c0, nn=1.7):

    magnetic_field = _convert_to_nparray(magnetic_field)

    critical_temp = np.zeros(magnetic_field.shape)
    # Find index such that magnetic field <= maximum critical magnetic field
    ind = magnetic_field <= B_c20
    # Evaluate critical temperature only for those values of magnetic field, elsewhere critical temperature is 0 by initialization.
    critical_temp[ind] = T_c0 * (1 - magnetic_field[ind] / B_c20) ** (1 / nn)
    return critical_temp


# End function critical_temperature_nbti


def g_func(alpha, beta):

    return ((alpha / (alpha + beta)) ** alpha) * ((beta / (alpha + beta)) ** beta)


# End function g_func


def critical_current_density_nbti(
    temperature,
    magnetic_field,
    B_c20,
    C_0,
    T_c0,
    alpha=[3.2, 0.65],
    beta=[2.43, 2],
    gamma=1.8,
    delta=0.45,
    nn=1.7,
):

    """
    Critical current density scaling of niobium titanium.

    For reference:
        Muzzi, L., De Marzi, G., Zignani, C.F., Affinito, L., Napolitano, M., Viola, R., Domínguez, C.O., Bottura, L., Le Naour, S., Richter, D. and Della Corte, A., 2010. Pinning properties of commercial Nb-Ti wires described by a 2-components model. IEEE Transactions on Applied Superconductivity, 20(3), pp.1496-1499.

    Returns:
        _type_: _description_
    """

    bb = reduced_magnetic_field_nbti(magnetic_field, temperature, B_c20, T_c0, nn=nn)
    tau = reduced_temperature_nbti(temperature, T_c0)

    # Compute critical current density:
    return (
        C_0
        / magnetic_field
        * (1.0 - tau**nn) ** gamma
        * (
            delta * bb ** alpha[0] * (1.0 - bb) ** beta[0] / g_func(alpha[0], beta[0])
            + (1.0 - delta)
            * bb ** alpha[1]
            * (1.0 - bb) ** beta[1]
            / g_func(alpha[1], beta[1])
        )
    )


# End function critical_current_density_nbti


def critical_current_density_bisection_nbti(
    temperature,
    magnetic_field,
    operative_current_density,
    B_c20,
    C_0,
    T_c0,
    alpha=[3.2, 0.65],
    beta=[2.43, 2],
    gamma=1.8,
    delta=0.45,
    nn=1.7,
):

    return operative_current_density - critical_current_density_nbti(
        temperature, magnetic_field, B_c20, C_0, T_c0, alpha, beta, gamma, delta, nn=nn
    )


# End function critical_current_density_bisection_nbti


def current_sharing_temperature_nbti(
    magnetic_field,
    op_current_density,
    B_c20,
    C_0,
    T_c0,
    alpha=[3.2, 0.65],
    beta=[2.43, 2],
    gamma=1.8,
    delta=0.45,
    nn=1.7,
    temp_lb=3.5,
):

    magnetic_field = _convert_to_nparray(magnetic_field)
    op_current_density = _convert_to_nparray(op_current_density[: magnetic_field.size])

    # Initialize with the original shape of the magnetic field.
    curr_shar_temp = np.zeros(magnetic_field.shape)
    current_density = np.zeros(magnetic_field.shape)
    critical_temp = np.zeros(magnetic_field.shape)

    # Find element index in magnetic_field such that magnetic_field < B_c20
    B_ind = magnetic_field < B_c20  # this is a boolean array
    # Check that the magnetic field is below the maximum critical magnetic field.
    if all(B_ind) == False:
        warnings.warn(
            "Magnetic field must be lower than maximum critical magnetic field!"
        )
        return curr_shar_temp

    # Evaluate current_density @ T = 0 only for elements corresponding to B_ind, elsewere current_density = 0.0 by initialization.
    current_density[B_ind] = critical_current_density_nbti(
        np.zeros(magnetic_field[B_ind].shape),
        magnetic_field[B_ind],
        B_c20,
        C_0,
        T_c0,
        alpha,
        beta,
        gamma,
        delta,
        nn=nn,
    )
    # Find element index in JC[B_ind] such that op_current_density[B_ind] < JC[B_ind] (boolean array)
    op_ind = op_current_density[B_ind] < current_density[B_ind]
    # Check that the operating current density is below the critical current density.
    if all(op_ind) == False:
        warnings.warn(
            "Operating current density must be below the critical current density!"
        )
        return curr_shar_temp

    critical_temp = critical_temperature_nbti(
        magnetic_field[op_ind], B_c20, T_c0, nn=nn
    )
    # Find index in op_current_density[op_ind] such that op_current_density[op_ind] > 0.0 (boolean array).
    op_ind_0 = op_current_density[op_ind] > 0.0
    if all(op_ind_0) == False:
        return critical_temp[np.invert(op_ind_0)]

    magnetic_field[op_ind_0] = np.maximum(magnetic_field[op_ind_0], 0.01)
    temp_ub = critical_temp[op_ind_0]
    root_result_obj = list()

    # Convert boolean array to an array of index.
    ind = np.nonzero(op_ind_0 == True)[0]
    for ii in range(ind.size):

        ex_args = (
            magnetic_field[ind[ii]],
            op_current_density[ind[ii]],
            B_c20,
            C_0,
            T_c0,
            alpha,
            beta,
            gamma,
            delta,
            nn,
        )
        # Evaluate current sharing temperature with bisection method.
        curr_shar_temp[ind[ii]], rr_obj = optimize.bisect(
            critical_current_density_bisection_nbti,
            temp_lb,
            temp_ub[ii],
            ex_args,
            xtol=1e-5,
            full_output=True,
        )
        root_result_obj.append(rr_obj)
    # End for ii.

    return curr_shar_temp


# End function current_sharing_temperature_nbti.


def electrical_resistivity_nbti(
    curr_dens: np.ndarray, crit_curr_dens: np.ndarray, E0: float = 1e-5, nn: int = 20
) -> np.ndarray:
    """Function that evaluathe the electrical resistivity of NbTi in Ohm*m.

    Args:
        curr_dens (np.ndarray): current density of the strand.
        crit_curr_dens (np.ndarray): critical current density of the strand.
        E0 (float, optional): electric field. Defaults to 1e-5.
        nn (int, optional): exponent of the correlation. Defaults to 20.

    Returns:
        np.ndarray: electrical resistivity of NbTi in Ohm*m.
    """
    # Is this valid in general or it is valid only in steady state (static) conditions?
    return E0 / crit_curr_dens * (curr_dens / crit_curr_dens) ** (nn - 1)
