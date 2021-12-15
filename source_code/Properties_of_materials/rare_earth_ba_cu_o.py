# Import python libraries and other functions
import numpy as np

# Function CONDREBCO starts here
def CONDREthermal_conductivity_rebco(TT):

    """
    ##############################################################################
    #                     FUNCTION CONDREBCO(TT)
    ##############################################################################
    #
    # Thermal conductivity of REBCO tape.
    #
    # from R.Heller data (measurement of a SuperOx tape for HTS-CL application
    #
    # INPUT
    # T [K]
    # OUTPUT
    # cond [W / m / K]
    #
    ##############################################################################
    # Translation from Fortran to Python: D.Placido PoliTo 10/07/2020
    # Tested against temperature in range [4,300] K: D.Placido PoliTo 11/07/2020
    ##############################################################################
    """

    TT = np.array(TT)

    Temp = np.array(
        [
            4.10023,
            6.05824,
            8.06484,
            10.06678,
            15.12171,
            20.12922,
            25.13238,
            30.22605,
            35.25125,
            40.29514,
            45.3187,
            50.3598,
            55.3828,
            60.41351,
            65.42644,
            71.14248,
            90.4347,
            120.78317,
            140.77557,
            160.73988,
            180.76685,
            211.27165,
            231.58216,
            250.70677,
            276.13281,
            291.47819,
            301.42616,
        ]
    )

    cond = np.array(
        [
            3.21495,
            4.84998,
            6.66456,
            8.53395,
            13.26289,
            17.41152,
            21.02355,
            24.1332,
            26.74567,
            28.93452,
            30.98772,
            33.01002,
            34.71014,
            36.51584,
            38.21827,
            40.41303,
            46.70009,
            55.98624,
            60.97617,
            66.83899,
            71.92989,
            80.74836,
            86.19238,
            90.10422,
            96.58573,
            100.56546,
            102.91945,
        ]
    )

    CONDREBCO = np.interp(TT, Temp, cond)

    return CONDREBCO


# Function CPREBCO starts here
def isobaric_specific_heat_rebco(TT):

    """
    ##############################################################################
    #                         FUNCTION CPREBCO(TT)
    ##############################################################################
    #
    # Specific heat of REBCO tape.
    #
    # from R.Heller data (measurement of a SuperOx tape for HTS-CL application
    #
    # INPUT
    # T [K]
    # OUTPUT
    # cp [J / kg / K]
    #
    ##############################################################################
    # Translation from Fortran to Python: D.Placido PoliTo 10/07/2020
    # Tested against temperature in range [4,300] K: D.Placido PoliTo 11/07/2020
    ##############################################################################
    """

    TT = np.array(TT)
    Temp = np.array(
        [
            4.0,
            6.0,
            8.0,
            10.0,
            15.0,
            20.0,
            25.0,
            30.0,
            40.0,
            50.0,
            60.0,
            70.0,
            80.0,
            90.0,
            100.0,
            120.0,
            140.0,
            160.0,
            180.0,
            200.0,
            220.0,
            240.0,
            260.0,
            280.0,
            300.0,
        ]
    )

    cp = np.array(
        [
            0.32473,
            0.72767,
            1.61794,
            2.65951,
            7.4452,
            15.33016,
            27.71197,
            43.23714,
            80.63019,
            117.34127,
            153.17621,
            182.20992,
            209.76962,
            233.62775,
            252.57224,
            277.24257,
            301.65733,
            319.30137,
            330.43024,
            341.40578,
            347.43717,
            353.46857,
            358.12862,
            361.46844,
            364.75715,
        ]
    )

    CPREBCO = np.interp(TT, Temp, cp)

    return CPREBCO
