# Import python libraries and other functions
import numpy as np

# Function condcu_nist starts here
def thermal_conductivity_cu_nist(t, b, rrr):
    """
    #
    # function kk=condcu(t,b,rrr);
    #
    ######################################################################
    # NOTE: This routine has been vectorised for MATLAB!!
    # The input t MUST be a row vector t = [t1 t2 ... tn]
    # The input b and rrr can be:
    # EITHER a scalar (>=0 and >1 respectively) => same values for all t
    # OR a row vector of the same length as t	  => 1-1 correspondence with t
    # A. Portone @ EFDA, Garching, 5.12.2000
    ######################################################################
    # Thermal conductivity of Copper in W/m K, as a function of T, B and
    # RRR for T > 0, B >= 0 and RRR >= 1.
    #
    #            References
    #            ----------
    # Luehning, Heller, Private Communication, 1989.
    #
    # variable  I/O         meaning            units
    # --------------------------------------------------------------------
    #   t     x      absolute temperature          K
    #   b     x      magnetic field            T
    #   rrr     x      residual resistivity ratio      -
    #   kk      x      thermal conductivity        W/m K
    #
    # Author : L.Bottura @ NET
    # Version: 1.1   24.7.1990
    #
    # Translation from MatLab to Python: S.Poccia UniTo & D.Placido PoliTo 03/2020
    ######################################################################
    """

    t = np.array(t)
    b = np.array(b)
    rrr = np.array(rrr)

    if len(b) == 1:
        b = b * np.ones(t.shape)
    if len(t) == 1:
        t = t * np.ones(t.shape)

    kk = np.zeros(t.shape)  # variable initialization

    beta = 0.634 / rrr
    betat = beta / 0.0003

    P1 = 1.754e-8
    P2 = 2.763
    P3 = 1102
    P4 = -0.165
    P5 = 70
    P6 = 1.756
    P7 = 0.838 / (betat**0.1661)

    W0 = beta / t
    Wi = (P1 * t**P2) / (1 + (P1 * P3 * t ** (P2 + P4) * np.exp(-((P5 / t) ** P6))))
    Wi0 = P7 * Wi * W0 / (Wi + W0)
    kk = (
        1.0
        / (Wi + W0 + Wi0)
        * rhoecu0_nist(t, rrr)
        / electrical_resistivity_cu_nist(t, b, rrr)
    )

    return kk


# Function cpcu_nist starts here
def isobaric_specific_heat_cu_nist(t):
    """
    ######################################################################
    # FUNCTION CPCU(T)
    ######################################################################
    #
    # Specific heat of Copper at constant pressure in J/kg/K, as a
    # function of T.
    #
    # Polynomial interpolation for cp.
    #
    #                References
    #                ----------
    #    E S Drexler, R P Reed, and N J Simon. Properties of copper
    #    and copper alloys at cryogenic temperatures. NIST Mono. NIST,
    #    Gaithersburg, MD, 1992
    #
    # 4 K <= T <= 300 K
    #
    # variable  I/O      meaning            units
    # --------------------------------------------------------------------
    #   T     I      temperature          K
    #   CPCU     O       Specific heat        J/kg/K
    #
    # Author : Zappatore A., Energy Department, Politecnico di Torino.
    # Version: 1.0 (August 2015).
    # Translation from MatLab to Python: S.Poccia UniTo & D.Placido PoliTo 03/2020
    # Correction after test: D.Placido PoliTo 07/04/2020
    # Tested against corresponding MatLab function in temperature range [4,300] K
    ######################################################################
    """

    t = np.array(t)
    cp = np.zeros(t.shape)  # variable initialization

    A0 = -1.91844
    A1 = -0.15973
    A2 = 8.61013
    A3 = -18.996
    A4 = 21.9661
    A5 = -12.7328
    A6 = 3.54322
    A7 = -0.3797

    intervals = [(t >= 4) & (t <= 300)]
    behavior = [
        lambda t: 10
        ** (
            A0
            + A1 * np.log10(t)
            + A2 * np.log10(t) ** 2
            + A3 * np.log10(t) ** 3
            + A4 * np.log10(t) ** 4
            + A5 * np.log10(t) ** 5
            + A6 * np.log10(t) ** 6
            + A7 * np.log10(t) ** 7
        )
    ]

    cp = np.piecewise(t, intervals, behavior)

    return cp  # end of the function


# Function rhoecu_nist starts here
def electrical_resistivity_cu_nist(t, b, rrr):
    """
    # function rho=rhoecu(t,b,rrr)
    #
    ######################################################################
    # NOTE: This routine has been vectorised for MATLAB!!
    # The input t MUST be a row vector t = [t1 t2 ... tn]
    # The input b and rrr can be:
    # EITHER a scalar (>=0 and >1 respectively) => same values for all t
    # OR a row vector of the same length as t	  => 1-1 correspondence with t
    # A. Portone @ EFDA, Garching, 5.12.2000
    ######################################################################
    # Electrical resistivity of Copper in Ohm m, as a function of T, B
    # and RRR for T > 0 , B >= 0 and RRR >= 1
    #
    #            References
    #            ----------
    # RHO(T,RRR): Luehning, Heller, Private Communication, 1989
    # RHO(B)  : F.R.Fickett, Int.Copper Res. Rep. 186, 1972 (recovered
    #       through a publication by V. Arp on the quench code)
    #
    # variable  I/O         meaning            units
    # --------------------------------------------------------------------
    #   t     x      absolute temperature          K
    #   b     x      magnetic field            T
    #   rrr     x      residual resistivity ratio      -
    #   rho     x      resistivity             Ohm m
    #
    # Author : L.Bottura @ NET
    # Version: 1.1   24.7.1990
    #
    # Translation from MatLab to Python: S.Poccia UniTo & D.Placido PoliTo 03/2020
    ######################################################################
    """

    a = np.array([[-2.662], [0.3168], [0.6229], [-0.1839], [0.01827]])
    t = np.array(t)
    b = np.array(b)
    rrr = np.array(rrr)

    T0 = 273.0  # [K]

    if len(b) == 1:
        b = b * np.ones(t.shape)
    if len(t) == 1:
        t = t * np.ones(b.shape)

    rho = np.zeros(t.shape)
    rhoN = np.zeros(t.shape)

    jok = np.nonzero((rrr > 1.0) & (b >= 0) & (t > 0))  # modified by Placido Daniele
    jng = np.setdiff1d(np.r_[:: len(t)], jok)  # modified by Placido Daniele

    jokb = np.nonzero((rrr > 1.0) & (b > 0) & (t > 0))  # modified by Placido Daniele

    if jok[0].size > 0:  # jok is not empty, modified by Placido Daniele
        rhoN[jok] = rhoecu0_nist(t[jok], rrr)

        # Add magnetoresistivity effect with respect to rhoecu0

        aaa = np.zeros(t.shape)
        xxx = np.zeros(t.shape)
        ccc = np.zeros(t.shape)

        xxx[jokb] = rhoecu0_nist(T0, rrr) / rhoecu0_nist(t[jokb], rrr) * b[jokb]

        for ij in range(5):  # changed range (1,5) -- > (5), modified by Placido Daniele
            aaa[jokb] = aaa[jokb] + a[ij] * (np.log10(xxx[jokb])) ** (ij)
            # end for loop

        ccc[jokb] = 10 ** (aaa[jokb])

        rho[jok] = rhoN[jok] * (1.0 + ccc[jok])

    if jng.size > 0:  # jng is not empty, modified by Placido Daniele
        print(" *** Warning: some data out of range in rhoecu")

    return rho  # end of the function


# Function rhoecu0_nist starts here
def rhoecu0_nist(t, rrr):
    """
    # function rho=rhoecu0(t,rrr);
    #
    ######################################################################
    # NOTE: This routine has been vectorised for MATLAB!!
    # The input t MUST be a row vector t = [t1 t2 ... tn]
    # The input b and rrr can be:
    # EITHER a scalar (>=0 and >1 respectively) => same values for all t
    # OR a row vector of the same length as t	  => 1-1 correspondence with t
    # A. Portone @ EFDA, Garching, 5.12.2000
    ######################################################################
    # Electrical resistivity of Copper in Ohm m, as a function of T, B
    # and RRR for T > 0 , B >= 0 and RRR >= 1
    #
    #            References
    #            ----------
    # RHO(T,RRR): Luehning, Heller, Private Communication, 1989
    # RHO(B)  : F.R.Fickett, Int.Copper Res. Rep. 186, 1972 (recovered
    #       through a publication by V. Arp on the quench code)
    #
    # variable  I/O         meaning            units
    # --------------------------------------------------------------------
    #   t     x      absolute temperature          K
    #   b     x      magnetic field            T
    #   rrr     x      residual resistivity ratio      -
    #   rho     x      resistivity             Ohm m
    #
    # Author : L.Bottura @ NET
    # Version: 1.1   24.7.1990
    #
    # Translation from MatLab to Python: S.Poccia Unito & D.Placido PoliTo 03/2020
    # Tested against corresponding MatLab function in temperature range [4,400] K
    ######################################################################
    """

    P0 = 1.553e-8
    P1 = 1.171e-17
    P2 = 4.49
    P3 = 3.841e10
    P4 = 1.14
    P5 = 50
    P6 = 6.428
    P7 = 0.4531

    t = np.array(t)
    rrr = np.array(rrr)
    rho = np.zeros(t.shape)  # variables initialization

    rho0 = P0 / rrr
    rhoi = (P1 * t**P2) / (1.0 + P1 * P3 * t ** (P2 - P4) * np.exp(-((P5 / t) ** P6)))
    rhoi0 = P7 * rhoi * rho0 / (rhoi + rho0)
    rho = rho0 + rhoi + rhoi0

    return rho  # end of the function


# Function rho_cu starts here
def density_cu(temperature: np.ndarray) -> np.ndarray:
    """
    Function that evaluates copper density, assumed constant.

    Args:
        temperature (np.ndarray): temperature array, used to get the shape of density array.
    Returns:
        np.ndarray: copper density array in kg/m^3.
    """
    return 8900.0 * np.ones(temperature.shape)


# end function rho_Cu (cdp, 01/2021)
