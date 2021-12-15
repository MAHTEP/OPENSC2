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
        lambda TT: 3.5389 * TT ** 4
        - 8.2624e1 * TT ** 3
        + 2.7138e2 * TT ** 2
        + 3.6897e3 * TT,
        lambda TT: -1.3578 * TT ** 3 + 1.2791e2 * TT ** 2 - 4.1285e3 * TT + 4.7377e4,
        lambda TT: 1.8407e-4 * TT ** 4
        - 5.9684e-2 * TT ** 3
        + 7.241 * TT ** 2
        - 3.9169e2 * TT
        + 8.4877e3,
        lambda TT: 1.4878e-7 * TT ** 4
        - 1.2548e-4 * TT ** 3
        + 3.9209e-2 * TT ** 2
        - 5.4107 * TT
        + 7.0959e2,
        lambda TT: 9.3528e-9 * TT ** 3
        - 2.6629e-5 * TT ** 2
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
        lambda TT: 0.047220784283425 * TT ** 2
        - 0.018765676876719 * TT
        - 1.946839691056091,
        lambda TT: -0.000000104620955 * TT ** 4
        + 0.000091884831207 * TT ** 3
        - 0.030295300516950 * TT ** 2
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
        lambda TT: 6.144183e-015 * TT ** 3
        - 6.690094e-014 * TT ** 2
        + 2.259567e-013 * TT
        + 9.822048e-012,
        lambda TT: 2.026667e-014 * TT ** 3
        - 6.160000e-013 * TT ** 2
        + 7.473333e-012 * TT
        - 2.330000e-011,
        lambda TT: -1.200000e-014 * TT ** 3
        + 2.210952e-012 * TT ** 2
        - 7.586429e-011 * TT
        + 8.015476e-010,
        lambda TT: 6.841184e-017 * TT ** 3
        - 4.447028e-014 * TT ** 2
        + 6.974508e-011 * TT
        - 2.428741e-009,
        lambda TT: 8.269045e-018 * TT ** 3
        - 3.077059e-015 * TT ** 2
        + 6.074742e-011 * TT
        - 1.812752e-009,
    ]

    RHOEAG = np.piecewise(TT, intervals, behavior)

    return RHOEAG


# Function rho_ag starts here
def density_ag():
    """
    Silver density kg/m^3. It is assumed constant.
    Autor: D. Placido Polito 21/01/2021
    """
    return 10630.0


# end function rho_Ag (cdp, 01/2021)
