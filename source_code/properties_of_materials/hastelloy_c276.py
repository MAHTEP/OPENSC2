# Import python libraries and other functions
import numpy as np

# Function CONDHC276 starts here
def thermal_conductivity_hc276(TT):

    """
    ##############################################################################
    #       FUNCTION CONDHC276(T)
    ##############################################################################
    #
    # Thermal conductivity HASTELLOY - C276
    # From "Experimental Analysis and Numerical Simulation
    # of Quench in Superconducting HTS Tapes and Coils", Ph. D. Thesis
    # by Marco Casali, 2014.
    #
    # INPUT
    # T [K]
    # OUTPUT
    # k [W/mK]
    #
    ##############################################################################
    # Translation from Fortran to Python: D.Placido PoliTo 27/07/2020
    # Tested against temperature in range [1,350] K: D.Placido PoliTo 27/07/2020
    ##############################################################################
    """

    TMIN = 1.0
    TMAX = 350.0

    TT = np.array(TT)
    CONDHC276 = np.zeros(TT.shape)

    TT = np.minimum(TT, TMAX)
    TT = np.maximum(TT, TMIN)

    intervals = [(TT <= 60.0), (TT > 60.0) & (TT < 300.0), (TT >= 300.0)]
    behavior = [
        lambda T: 0.000028975912976 * T ** 3
        - 0.004655231435231 * T ** 2
        + 0.293992451992451 * T
        + 0.166939393939405,
        lambda T: 0.000000124577333 * T ** 3
        - 0.000079993515434 * T ** 2
        + 0.035243632195403 * T
        + 5.503358179018721,
        lambda T: (12.240619439579623 - 12.219768038241552) / (300 - 299) * (T - 300)
        + 12.240619439579623,
    ]

    CONDHC276 = np.piecewise(TT, intervals, behavior)

    return CONDHC276


# Function CPHC276 starts here
def isobaric_specific_heat_hc276(TT):

    """
    ##############################################################################
    #      FUNCTION CPHC276(TT)
    ##############################################################################
    #
    # From "Experimental Analysis and Numerical Simulation
    # of Quench in Superconducting HTS Tapes and Coils", Ph. D. Thesis
    # by Marco Casali, 2014.
    #
    # INPUT
    # T  [K]
    #
    # OUTPUT
    # cp [J/kgK]
    #
    ##############################################################################
    # Translation from Fortran to Python: D.Placido PoliTo 27/07/2020
    # Tested against temperature in range [1,350] K: D.Placido PoliTo 27/07/2020
    ##############################################################################
    """

    TMIN = 1.0
    TMAX = 350.0

    TT = np.array(TT)
    CPHC276 = np.zeros(TT.shape)

    TT = np.minimum(TT, TMAX)
    TT = np.maximum(TT, TMIN)

    intervals = [(TT <= 60.0), (TT > 60.0) & (TT < 300.0), (TT >= 300.0)]
    behavior = [
        lambda T: 0.000092969696970 * T ** 3
        + 0.046381858141858 * T ** 2
        - 0.738681152181147 * T
        + 4.091909090909040,
        lambda T: 0.0000302086983 * T ** 3
        - 0.0226456716880 * T ** 2
        + 6.0225562313795 * T
        - 179.3891242698739,
        lambda T: (4.049021473239762e2 - 4.043151575124803e002)
        / (300 - 299)
        * (T - 300)
        + 4.049021473239762e2,
    ]

    CPHC276 = np.piecewise(TT, intervals, behavior)

    return CPHC276


# Function RHOEHC276 starts here
def electrical_resistivity_hc276(TT):

    """
    ##############################################################################
    #          FUNCTION RHOEHC276(TT)
    ##############################################################################
    #
    # Electrical resistivity HASTELLOY - C276
    # From "Experimental Analysis and Numerical Simulation
    # of Quench in Superconducting HTS Tapes and Coils", Ph. D. Thesis
    # by Marco Casali, 2014.
    #
    # INPUT
    # T [K]
    # OUTPUT
    # rho_el [Ohm x m]
    #
    ##############################################################################
    # Translation from Fortran to Python: D.Placido PoliTo 27/07/2020
    # Tested against temperature in range [4,300] K: D.Placido PoliTo 27/07/2020
    ##############################################################################
    """

    TMIN = 4.0
    TMAX = 300.0

    TT = np.array(TT)
    RHOEHC276 = np.zeros(TT.shape)

    TT = np.minimum(TT, TMAX)
    TT = np.maximum(TT, TMIN)

    RHOEHC276 = 1.309449e-10 * TT + 1.220104e-06

    return RHOEHC276


# Function rho_HC276 starts here
def density_hc276(temperature: np.ndarray) -> np.ndarray:
    """
    Function that evaluates hastelloy HC276 density, assumed constant.

    Args:
        temperature (np.ndarray): temperature array, used to get the shape of density array.

    Returns:
        np.ndarray: hastelloy HC276 density array in kg/m^3.
    """
    return 8890.0 * np.ones(temperature.shape)


# end function rho_HC276 (cdp, 01/2021)
