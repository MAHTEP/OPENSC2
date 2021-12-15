# Import python libraries and other functions
import numpy as np

# Function CONDSN60PB40.py starts here
def thermal_conductivity_sn60pb40(TT):

    """
  ##############################################################################
  #            FUNCTION CONDSN60PB40(TT)
  ##############################################################################
  #
  # Thermal conductivity of solder Sn60Pb40
  #
  # N. Bagrets, C. Barth, and K.-P. Weiss,
  # Low Temperature Thermal and Thermo-Mechanical Properties of Soft Solders for
  # Superconducting Applications
  # IEEE Transaction on Applied Superconductivity, Vol. 24, No. 3, Jun. 2014, \
  # Art. ID. 7800203
  #
  # INPUT
  # T [K]
  # OUTPUT
  # k [W/mK]
  #
  ##############################################################################
  # Translation from Fortran to Python: D.Placido PoliTo 27/07/2020
  # Tested against temperature in range [4,350] K: D.Placido PoliTo 27/07/2020
  ##############################################################################
  """

    TMIN = 4.0
    TMAX = 350.0

    TT = np.array(TT)
    CONDSN60PB40 = np.zeros(TT.shape)

    TT = np.minimum(TT, TMAX)
    TT = np.maximum(TT, TMIN)

    intervals = [(TT <= 300.0), (TT >= 300.0)]
    behavior = [
        lambda T: 2.667e-010 * T ** 5
        - 2.483e-7 * T ** 4
        + 8.769e-005 * T ** 3
        - 0.01453 * T ** 2
        + 1.147 * T
        + 10.22,
        lambda T: 0.0902 * (T - 300.0) + 51.1010,
    ]

    CONDSN60PB40 = np.piecewise(TT, intervals, behavior)

    return CONDSN60PB40


# Function CPSN60PB40.py starts here
def isobaric_specific_heat_sn60pb40(TT):

    """
  ##############################################################################
  #            FUNCTION CPSN60PB40(TT)
  ##############################################################################
  #
  # Specific heat of solder Sn60Pb40.
  #
  # N. Bagrets, C. Barth, and K.-P. Weiss,
  # Low Temperature Thermal and Thermo-Mechanical Properties of Soft Solders for
  # Superconducting Applications
  # IEEE Transaction on Applied Superconductivity, Vol. 24, No. 3, Jun. 2014, \
  # Art. ID. 7800203
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
    CPSN60PB40 = np.zeros(TT.shape)

    TT = np.minimum(TT, TMAX)
    TT = np.maximum(TT, TMIN)

    intervals = [(TT <= 300.0), (TT >= 300.0)]
    behavior = [
        lambda T: 6.856e-010 * T ** 5
        - 6.764e-7 * T ** 4
        + 0.0002604 * T ** 3
        - 0.0492 * T ** 2
        + 4.728 * T
        - 27.32,
        lambda T: 0.2316 * (T - 300.0) + 181.0480,
    ]
    CPSN60PB40 = np.maximum(np.piecewise(TT, intervals, behavior), 1.0)
    # CPSN60PB40 = np.piecewise(TT, intervals, behavior)

    return CPSN60PB40


# Function RHOESN60PB40.py starts here
def electrical_resistivity_sn60pb40(TT):

    """
  ##############################################################################
  #            FUNCTION RHOESN60PB40(TT)
  ##############################################################################
  #
  # Electrical resistivity of solder Sn60Pb40.
  #
  # N. Bagrets, C. Barth, and K.-P. Weiss,
  # Low Temperature Thermal and Thermo-Mechanical Properties of Soft Solders for
  # Superconducting Applications
  # IEEE Transaction on Applied Superconductivity, Vol. 24, No. 3, Jun. 2014, \
  # Art. ID. 7800203
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
    RHOESN60PB40 = np.zeros(TT.shape)

    TT = np.minimum(TT, TMAX)
    TT = np.maximum(TT, TMIN)

    RHOESN60PB40 = np.maximum(5.146e-010 * TT - 7.346e-009, 0.01e-9)
    return RHOESN60PB40


# Function rho_Sn60Pb40 starts here
def density_sn60pb40():
    """
    Sn60Pb40 density kg/m^3. It is assumed constant.
    Autor: D. Placido Polito 21/01/2021
    """
    return 7310.0 * 0.6 + 11340.0 * 0.4


# end function rho_Sn60Pb40 (cdp, 01/2021)
