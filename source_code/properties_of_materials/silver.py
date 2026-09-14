# Import python libraries and other functions
import numpy as np

# Function CONDAG starts here
def thermal_conductivity_ag(TT):

    """
    ##############################################################################
    #								FUNCTION CONDAG(TT)
    ##############################################################################
    #
    # FUNCTION FOR SILVER - CONDAG.M
    # Thermophysical Properties of Matter, v1, Y.S. Touloukian, R.W. Powell,
    # C.Y. Ho & P.G. Klemens, 1970, IFI/Plenum, NY, NY
    # Purity 99.999#  well-annealed with residual resistivity of 0.000620
    # uohm-cm
    # error is 2#  near RT, 2-5#  at others
    #
    # INPUT
    # TT [K]
    # OUTPUT
    # k [W/mK]
    #
    ##############################################################################
    # Translation from Fortran to Python: D.Placido PoliTo 10/07/2020
    # Tested against temperature in range [1,1234] K: D.Placido PoliTo 11/07/2020
    ##############################################################################
    """

    TMIN = 1.0
    TMAX = 1234.0

    TT = np.array(TT)
    CONDAG = np.zeros(TT.shape)

    TT = np.minimum(TT, TMAX)
    TT = np.maximum(TT, TMIN)

    intervals = [
        (TT >= 1.0) & (TT < 13.0),
        (TT >= 13.0) & (TT < 35.0),
        (TT >= 35.0) & (TT < 100.0),
        (TT >= 100.0) & (TT < 273.0),
        (TT >= 273.0) & (TT < 1234.0),
    ]
    behavior = [
        lambda TT: 3.5389 * TT**4
        - 8.2624e1 * TT**3
        + 2.7138e2 * TT**2
        + 3.6897e3 * TT,
        lambda TT: -1.3578 * TT**3 + 1.2791e2 * TT**2 - 4.1285e3 * TT + 4.7377e4,
        lambda TT: 1.8407e-4 * TT**4
        - 5.9684e-2 * TT**3
        + 7.241 * TT**2
        - 3.9169e2 * TT
        + 8.4877e3,
        lambda TT: 1.4878e-7 * TT**4
        - 1.2548e-4 * TT**3
        + 3.9209e-2 * TT**2
        - 5.4107 * TT
        + 7.0959e2,
        lambda TT: 9.3528e-9 * TT**3
        - 2.6629e-5 * TT**2
        - 5.4314e-2 * TT
        + 4.4529e2,
    ]

    CONDAG = np.piecewise(TT, intervals, behavior)

    return CONDAG


# Function CPAG starts here
def isobaric_specific_heat_ag(TT):

    """
    ##############################################################################
    #                   FUNCTION CPAG(TT)
    ##############################################################################
    #
    # SPECIFIC HEAT FOR SILVER - CPAG.M
    # From "Low-Temperature Properties of Silver", D.R. Smith and F.R. Fickett,
    # Journal of Research of the National Institute of Standards and
    # Technology, Volume 100, Number 2, March?April 1995.
    #
    # INPUT
    # Input [K]
    # OUTPUT
    # cp [J/kgK]
    #
    #
    ##############################################################################
    # Translation from Fortran to Python: D.Placido PoliTo 10/07/2020
    # Tested against temperature in range [4,300] K: D.Placido PoliTo 11/07/2020
    ##############################################################################
    """

    TMIN = 1.0
    TMAX = 300.0

    TT = np.array(TT)
    CPAG = np.zeros(TT.shape)

    TT = np.minimum(TT, TMAX)
    TT = np.maximum(TT, TMIN)

    intervals = [(TT <= 50.0), (TT > 50.0) & (TT < 285.0), (TT >= 285.0)]
    behavior = [
        lambda TT: 0.047220784283425 * TT**2
        - 0.018765676876719 * TT
        - 1.946839691056091,
        lambda TT: -0.000000104620955 * TT**4
        + 0.000091884831207 * TT**3
        - 0.030295300516950 * TT**2
        + 4.598581892940776 * TT
        - 52.331891977926269,
        lambda TT: (2.343447902350777e2 - 2.343414867962295e2)
        / (285 - 284.9)
        * (TT - 295)
        + 2.343447902350777e2,
    ]

    CPAG = np.piecewise(TT, intervals, behavior)

    CPAG = np.maximum(CPAG, 1.0)

    return CPAG


# Function RHOEAG starts here
def electrical_resistivity_ag(TT):

    """
    ##############################################################################
    #                  FUNCTION RHOEAG(TT)
    ##############################################################################
    #
    # Electrical resistivity SILVER
    # R.A. Matula, J. Phys. Chem. Ref. Data, vol 8, no. 4, p 1147 (1979) purity
    # 99.995% or higher data below 40K is for Ag with a residual resistivity of
    # 0.001 x 10E-8 ohm-m (RRR=1450)
    #
    # INPUT
    # T [K]
    # OUTPUT
    # rho_el [Ohm x m]
    #
    ##############################################################################
    # Translation from Fortran to Python: D.Placido PoliTo 10/07/2020
    # Tested against temperature in range [1,1234] K: D.Placido PoliTo 11/07/2020
    ##############################################################################
    """

    TMIN = 1.0
    TMAX = 1235.0

    TT = np.array(TT)
    RHOEAG = np.zeros(TT.shape)

    TT = np.minimum(TT, TMAX)
    TT = np.maximum(TT, TMIN)

    intervals = [
        (TT >= 1.0) & (TT < 15.0),
        (TT >= 15.0) & (TT < 30.0),
        (TT >= 30.0) & (TT < 60.0),
        (TT >= 60.0) & (TT < 200.0),
        (TT >= 200.0) & (TT < 1235.0),
    ]
    behavior = [
        lambda TT: 6.144183e-015 * TT**3
        - 6.690094e-014 * TT**2
        + 2.259567e-013 * TT
        + 9.822048e-012,
        lambda TT: 2.026667e-014 * TT**3
        - 6.160000e-013 * TT**2
        + 7.473333e-012 * TT
        - 2.330000e-011,
        lambda TT: -1.200000e-014 * TT**3
        + 2.210952e-012 * TT**2
        - 7.586429e-011 * TT
        + 8.015476e-010,
        lambda TT: 6.841184e-017 * TT**3
        - 4.447028e-014 * TT**2
        + 6.974508e-011 * TT
        - 2.428741e-009,
        lambda TT: 8.269045e-018 * TT**3
        - 3.077059e-015 * TT**2
        + 6.074742e-011 * TT
        - 1.812752e-009,
    ]

    RHOEAG = np.piecewise(TT, intervals, behavior)

    return RHOEAG


# Function rho_ag starts here
def density_ag(temperature: np.ndarray) -> np.ndarray:
    """
    Function that evaluates silver density, assumed constant.

    Args:
        temperature (np.ndarray): temperature array, used to get the shape of density array.

    Returns:
        np.ndarray: silver density array in kg/m^3.
    """
    return 10630.0 * np.ones(temperature.shape)


# end function rho_Ag (cdp, 01/2021)
# ---------------------------------------------------------------------------
# CryoSoft silver properties (opt-in)
# ---------------------------------------------------------------------------
# Direct Python translation of the CryoSoft Ag.f routines supplied for the
# H4C benchmark audit. These functions intentionally coexist with the legacy
# OPENSC2 functions above so input files using material key "ag" remain on
# the historical property path.

def _broadcast_ag_cryosoft_inputs(temperature, magnetic_field, rrr):
    """Broadcast and clamp CryoSoft silver inputs to the Fortran ranges."""
    temperature, magnetic_field, rrr = np.broadcast_arrays(
        np.asarray(temperature, dtype=float),
        np.asarray(magnetic_field, dtype=float),
        np.asarray(rrr, dtype=float),
    )
    temperature = np.clip(temperature, 0.1, 1000.0)
    magnetic_field = np.clip(magnetic_field, 0.0, 100.0)
    rrr = np.clip(rrr, 1.5, 10000.0)
    return temperature, magnetic_field, rrr


def magnetoresistance_ag_cryosoft(temperature, magnetic_field, rrr):
    """CryoSoft transverse magnetoresistivity factor for pure silver."""
    tt, bb, rr = _broadcast_ag_cryosoft_inputs(
        temperature, magnetic_field, rrr
    )
    bb = np.clip(bb, 0.0, 30.0)

    rho273 = 1.48e-8
    p1 = 1.474e-17
    p2 = 4.82
    p3 = 1.16e11
    p4 = -1.33
    p5 = 10.0
    p6 = 1.0
    p7 = 0.333
    a1 = 6.5e-05
    a2 = 1.8
    a3 = 3.0e-3
    a4 = 1.1

    rhozero = rho273 / (rr - 1.0)
    arg = np.minimum((p5 / tt) ** p6, 30.0)
    rhoi = p1 * tt**p2 / (
        1.0 + p1 * p3 * tt ** (p2 + p4) * np.exp(-arg)
    )
    rhoi0 = p7 * rhoi * rhozero / (rhoi + rhozero)
    rho0 = rhozero + rhoi + rhoi0

    brr = np.clip(bb * rho273 / rho0, 0.0, 10.0e3)
    increase = np.zeros_like(tt, dtype=float)
    mask = brr > 1.0
    increase[mask] = (
        a1 * brr[mask] ** a2 / (1.0 + a3 * brr[mask] ** a4)
    )
    return increase + 1.0


def thermal_conductivity_ag_cryosoft(temperature, magnetic_field, rrr):
    """CryoSoft silver thermal conductivity k(T, B, RRR), W/(m K)."""
    tt, bb, rr = _broadcast_ag_cryosoft_inputs(
        temperature, magnetic_field, rrr
    )

    rho273 = 1.48e-8
    lorenz = 2.443e-8
    alpha2 = 7.87561879e-8
    m = 2.75
    n = 2.3
    ell = -2.5
    t0 = 35.0
    k0 = 408.59712205

    rhozero = rho273 / (rr - 1.0)
    beta = rhozero / lorenz
    alpha = alpha2 * (beta / n / alpha2) ** ((m - n) * (m + ell))

    w0 = beta / tt
    wi = alpha * tt**n
    wi0 = wi + w0

    wt = wi0.copy()
    mask = tt > t0
    wt[mask] = wi0[mask] / (
        1.0
        + wi0[mask] * k0 * (1.0 - np.exp(-(tt[mask] - t0) / t0))
    )

    return 1.0 / (wt * magnetoresistance_ag_cryosoft(tt, bb, rr))


def electrical_resistivity_ag_cryosoft(temperature, magnetic_field, rrr):
    """CryoSoft silver electrical resistivity rho(T, B, RRR), ohm m."""
    tt, bb, rr = _broadcast_ag_cryosoft_inputs(
        temperature, magnetic_field, rrr
    )

    rho273 = 1.48e-8
    p1 = 1.474e-17
    p2 = 4.82
    p3 = 1.16e11
    p4 = -1.33
    p5 = 10.0
    p6 = 1.0
    p7 = 0.333

    rhozero = rho273 / (rr - 1.0)
    arg = np.minimum((p5 / tt) ** p6, 30.0)
    rhoi = p1 * tt**p2 / (
        1.0 + p1 * p3 * tt ** (p2 + p4) * np.exp(-arg)
    )
    rhoi0 = p7 * rhoi * rhozero / (rhoi + rhozero)
    rho0 = rhozero + rhoi + rhoi0

    return magnetoresistance_ag_cryosoft(tt, bb, rr) * rho0


def isobaric_specific_heat_ag_cryosoft(temperature):
    """CryoSoft silver specific heat cp(T), J/(kg K)."""
    tt = np.clip(np.asarray(temperature, dtype=float), 1.0, 300.0)

    t0 = 9.785133189
    t1 = 34.5580572
    a1 = 6.14e-03
    a2 = -9.09e-04
    a3 = 1.79e-03
    b0 = -3.515882063
    b1 = 0.768177822
    b2 = -0.072031
    b3 = 5.59e-03
    b4 = -7.55e-05
    aa = -14331.2453
    cc = 1415.969135
    a = 217.2435636
    c = 23.862111
    na = 1.475372864
    nc = 3.091515729

    cp = np.empty_like(tt, dtype=float)
    low = tt <= t0
    mid = (tt > t0) & (tt <= t1)
    high = tt > t1

    cp[low] = a1 * tt[low] + a2 * tt[low] ** 2 + a3 * tt[low] ** 3
    cp[mid] = (
        b0
        + b1 * tt[mid]
        + b2 * tt[mid] ** 2
        + b3 * tt[mid] ** 3
        + b4 * tt[mid] ** 4
    )
    cp[high] = (
        aa * tt[high] / (a + tt[high]) ** na
        + cc * tt[high] ** 3 / (c + tt[high]) ** nc
    )
    return cp


def density_ag_cryosoft(temperature):
    """CryoSoft silver density, kg/m3 (constant)."""
    temperature = np.asarray(temperature)
    return 10490.0 * np.ones(temperature.shape)
