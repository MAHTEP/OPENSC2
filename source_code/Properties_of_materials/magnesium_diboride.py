import numpy as np
from scipy import optimize
import warnings


def thermal_conductivity_mgb2(temp):
    """Very ugly fit of thermal conductivity: cut at 5 K (for temperatures below 5 K thermal conductivity is 0 W/m/k) and cut at 39 K (the critical temperature) (for temperature above 39 K thermal conductivity is assumed 120.0 W/m/k). For temperature between 5 K and 39 K thermal conductivity is evaluated with the line equation.

    Args:
        temp (numpy array of float): temperature array
    """

    Ta = 5.0  # K
    Tb = 39.0  # K
    ka = 0.0  # W/m/K
    kb = 120.0  # W/m/K

    temp = np.maximum(temp, Ta)
    temp = np.minimum(temp, Tb)
    # Line equation:
    # k = ka + (kb - ka) * (temp - Ta) / (Tb - Ta)
    return ka + (kb - ka) * (temp - Ta) / (Tb - Ta)


# End function thermal_conductivity_mgb2


def density_mgb2():
    """Magnesium diboride density. It is assumed constant.

    Returns:
        [float]: density value in kg/m^3.
    """
    return 2.57e3  # kg/m^3


# End function density_mgb2


def isobaric_specific_heat_mgb2(temp):
    """Isobaric specific heat of magnesium diboride (MgB2). For the time being the same correlation valid for copper is used. To be fixed as soon as possible.

    Args:
        temp (numpy array of float): temperature

    Returns:
        numpy array of float: isobaric specific heat in J/kg/K
    """

    temp = np.array(temp)
    cp = np.zeros(temp.shape)  # variable initialization

    A0 = -1.91844
    A1 = -0.15973
    A2 = 8.61013
    A3 = -18.996
    A4 = 21.9661
    A5 = -12.7328
    A6 = 3.54322
    A7 = -0.3797

    intervals = [(temp >= 4) & (temp <= 300)]
    behavior = [
        lambda temp: 10
        ** (
            A0
            + A1 * np.log10(temp)
            + A2 * np.log10(temp) ** 2
            + A3 * np.log10(temp) ** 3
            + A4 * np.log10(temp) ** 4
            + A5 * np.log10(temp) ** 5
            + A6 * np.log10(temp) ** 6
            + A7 * np.log10(temp) ** 7
        )
    ]

    cp = np.piecewise(temp, intervals, behavior)

    return cp  # J/Kg/K
    # End function isobaric_specific_heat_mgb2


def critical_magnetic_field_mgb2(temp, Bc20, Tc0):
    """Function that evaluates the critical magnetic field of magnesium diboride.
    Ref: G. Giunchi et al."High performance new MgB2 superconducting hollow wires,” Supercond. Sci. Technol., vol. 16, pp. 285–291, 2003

    Args:
        temp (numpy array of float): temperature in K
        Bc20 (float): maximum critical magnetic field at 0 K in T
        Tc0 (float): maximum critical temperature at 0 T in K

    Returns:
        numpy array of float: critical magnetic field in T
    """
    alpha = 1.2
    return Bc20 * (1.0 - (temp / Tc0) ** alpha)


# End function critical_magnetic_field_mgb2


def critical_current_density_mgb2(temp, magnetic_field, Bc20, C0, Tc0):
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

    # Added to avoid division by 0
    magnetic_field = np.maximum(magnetic_field, 0.01)

    tau = temp / Tc0
    bb = magnetic_field / critical_magnetic_field_mgb2(temp, Bc20, Tc0)
    return (
        C0 / magnetic_field * (1.0 - tau ** beta) ** gamma * bb ** pp * (1.0 - bb) ** qq
    )


# End function critical_current_density_mgb2


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
        return Tc0

    magnetic_field[op_ind_0] = np.maximum(magnetic_field[op_ind_0], 0.01)
    temp_ub = Tc0
    root_result_obj = list()

    # Convert boolean array to an array of index.
    ind = np.nonzero(op_ind_0 == True)[0]
    for ii in range(ind.size):

        ex_args = (magnetic_field[ii], Bc20, C0, Tc0, op_current_density[ii])
        # Evaluate current sharing temperature with bisection method.
        curr_shar_temp[ind[ii]], rr_obj = optimize.bisect(
            _critical_current_density_residual,
            temp_lb,
            temp_ub,
            ex_args,
            xtol=1e-5,
            full_output=True,
        )
        root_result_obj.append(rr_obj)
    # End for ii.

    ex_args = (magnetic_field, Bc20, C0, Tc0, op_current_density)


    # # Evaluate current sharing temperature with secant method (not working)
    # curr_shar_temp, converged, _ = optimize.newton(
    #     _critical_current_density_residual, curr_shar_temp, args = ex_args, full_output = True)


    # # Evaluate current sharing temperature with Newthon-Rampson method (not working)
    # curr_shar_temp, converged, _ = optimize.newton(
    #     _critical_current_density_residual,
    #     curr_shar_temp,
    #     args=ex_args,
    #     fprime=_d_critical_current_density,
    #     full_output=True,
    # )

    return curr_shar_temp


# End function current_sharing_temperature_mgb2
