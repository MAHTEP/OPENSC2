# Import python libraries and other functions
import numpy as np
from scipy import optimize
import warnings

# Is the same function in module copper.py: at the time being mgb2 thermal 
# conductivity is assumed to be the same of copper with RRR = 100
def thermal_conductivity_mgb2(t, b, rrr=100):
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


# Is the same function in module copper.py: at the time being mgb2 isobaric 
# specific heat is assumed to be the same of copper.
def isobaric_specific_heat_mgb2(t):
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

def electrical_resistivity_mgb2(
    curr_dens: np.ndarray, crit_curr_dens: np.ndarray, E0: float = 1e-4, nn: int = 20
) -> np.ndarray:
    """Function that evaluathe the electrical resistivity of magnesium diboride (MgB2) in Ohm*m.

    Args:
        curr_dens (np.ndarray): current density of the strand.
        crit_curr_dens (np.ndarray): critical current density of the strand.
        E0 (float, optional): electric field. Defaults to 1e-4.
        nn (int, optional): exponent of the correlation. Defaults to 20.

    Returns:
        np.ndarray: electrical resistivity of magnesium diboride (MgB2) in Ohm*m.
    """
    # Is this valid in general or it is valid only in steady state (static) 
    # conditions?
    return E0 / crit_curr_dens * (curr_dens / crit_curr_dens) ** (nn - 1)

def density_mgb2(temperature: np.ndarray) -> np.ndarray:
    """
    Function that evaluates magnesium diboride (MgB2) density, assumed constant.

    Args:
        temperature (np.ndarray): temperature array, used to get the shape of density array.
    Returns:
        np.ndarray: magnesium diboride density array in kg/m^3.
    """
    return 2570.0 * np.ones(temperature.shape)


def critical_magnetic_field_mgb2(temp, Bc20, Tc0):
    """Function that evaluates the critical magnetic field of magnesium diboride.
    Ref: G. Giunchi et al."High performance new MgB2 superconducting hollow wires,” Supercond. Sci. Technol., vol. 16, pp. 285-291, 2003

    Args:
        temp (numpy array of float): temperature in K
        Bc20 (float): maximum critical magnetic field at 0 K in T
        Tc0 (float): maximum critical temperature at 0 T in K

    Returns:
        numpy array of float: critical magnetic field in T
    """
    alpha = 1.2
    return Bc20 * (1.0 - (temp / Tc0) ** alpha)


# End function critical_magnetic_field


def critical_current_density_mgb2(
    temp,
    magnetic_field,
    Bc20,
    C0,
    Tc0,
    b_low = 0.01,
    ):
    """Function that evaluates the critical current density of magnesium diboride.

    Args:
        temp (numpy array of float): temperature in K
        magnetic_field (numpy array of float): magnetic field in T
        Bc20 (float): maximum critical magnetic field at 0 K in T
        C0 (float): scaling coefficient in A*T/m^2
        Tc0 (float): maximum critical temperature at 0 T in K

    Returns:
        numpy array of float: critical current density in A/m^2
    """
    pp = 0.5
    qq = 5.0
    beta = 1.55
    gamma = 1.89

    magnetic_field = np.maximum(magnetic_field,b_low)
    tau = temp / Tc0
    bb = magnetic_field / critical_magnetic_field_mgb2(temp, Bc20, Tc0)

    return (
        (C0 / magnetic_field) * ((1.0 - tau ** beta) ** gamma) * (bb ** pp) 
        * ((1.0 - bb) ** qq)
    )


# End function critical_current_density


def _critical_current_density_residual(
    temp, magnetic_field, Bc20, C0, Tc0, operative_current_density
):
    # Function used to evaluate the current sharign temperaure with the bisection method
    return (
        critical_current_density_mgb2(temp, magnetic_field, Bc20, C0, Tc0)
        - operative_current_density
    )


# End function _critical_current_density_residual


def _d_critical_current_density(temp, magnetic_field, Bc20, C0, Tc0, dummy):
    # First derivative of the critical current density wrt temperature at constant magnetic field.
    pp = 0.5
    qq = 5.0
    alpha = 1.2
    beta = 1.55
    gamma = 1.89

    tau = temp / Tc0
    bb = magnetic_field / critical_magnetic_field_mgb2(temp, Bc20, Tc0)
    # la derivata è corretta!
    return (
        critical_current_density_mgb2(temp, magnetic_field, Bc20, C0, Tc0)
        / temp
        * (
            alpha * tau ** alpha / (1.0 - tau ** alpha) * (pp - qq * bb / (1.0 - bb))
            - beta * gamma * tau ** beta / (1.0 - tau ** beta)
        )
    )


# End function _d_critical_current_density


def current_sharing_temperature_mgb2(
    magnetic_field, op_current_density, Bc20, C0, Tc0, temp_lb=3.5
):

    """Function that evaluate the current sharing temperature inverting the equation of the critical current.

    Returns:
        np.ndarray: current sharing temperature
    """

    # magnetic_field = _convert_to_nparray(magnetic_field.dropna())
    # op_current_density = _convert_to_nparray(op_current_density[:magnetic_field.size])

    # Initialize with the original shape of the magnetic field.
    curr_shar_temp = np.zeros(magnetic_field.shape)
    current_density = np.zeros(magnetic_field.shape)

    # Find element index in magnetic_field such that magnetic_field < Bc20
    B_ind = magnetic_field < Bc20  # this is a boolean array
    # Check that the magnetic field is below the maximum critical magnetic field.
    if all(B_ind) == False:
        warnings.warn(
            "Magnetic field must be lower than maximum critical magnetic field!"
        )
        return curr_shar_temp

    # Evaluate current_density @ T = 0 only for elements corresponding to B_ind, elsewere current_density = 0.0 by initialization.
    current_density[B_ind] = critical_current_density_mgb2(
        np.zeros(magnetic_field[B_ind].shape), magnetic_field[B_ind], Bc20, C0, Tc0
    )
    # Find element index in JC[B_ind] such that op_current_density[B_ind] < JC[B_ind] (boolean array)
    op_ind = op_current_density[B_ind] < current_density[B_ind]
    # Check that the operating current density is below the critical current density.
    if all(op_ind) == False:
        warnings.warn(
            "Operating current density must be below the critical current density!"
        )
        return curr_shar_temp

    # Find index in op_current_density[op_ind] such that op_current_density[op_ind] > 0.0 (boolean array).
    op_ind_0 = op_current_density[op_ind] > 0.0
    if all(op_ind_0) == False:
        return Tc0 * np.ones_like(magnetic_field)

    magnetic_field[op_ind_0] = np.maximum(magnetic_field[op_ind_0], 0.01)
    temp_ub = Tc0

    # Convert boolean array to an array of index.
    ind = np.nonzero(op_ind_0 == True)[0]
    for ii in range(ind.size):

        ex_args = (magnetic_field[ii], Bc20, C0, Tc0, op_current_density[ii])
        # Evaluate current sharing temperature with bisection method.
        curr_shar_temp[ind[ii]] = optimize.bisect(
            _critical_current_density_residual,
            temp_lb,
            temp_ub,
            ex_args,
            xtol=1e-5,
            maxiter=10,
            disp=False
        )
    # End for ii.

    ex_args = (magnetic_field, Bc20, C0, Tc0, op_current_density)

    # Evaluate current sharing temperature with secant method (working!)
    # curr_shar_temp = optimize.newton(
    #     _critical_current_density_residual,
    #     curr_shar_temp,
    #     args = ex_args,
    #     full_output = True,
    # )

    # Evaluate current sharing temperature with Newthon-Rampson method
    curr_shar_temp = optimize.newton(
        _critical_current_density_residual,
        curr_shar_temp,
        args=ex_args,
        fprime=_d_critical_current_density,
    )

    return curr_shar_temp


# End function current_sharing_temperature.


# Needed for the thermal conductivity.
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
