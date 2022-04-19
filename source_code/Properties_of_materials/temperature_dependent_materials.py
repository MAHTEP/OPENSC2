import numpy as np
import warnings

CRYOSOFT = 0
KSTAR = 1
NIST = 1


class TemperatureDependentMaterial:
    """Class TemperatureDependentMaterial is the parent class of materials whose properties are temperature dependent only."""

    property_name = "Prop"

    def _convert_to_nparray(self, value):
        """Private method that covert value to numpy array if it is not a numpy array.

        Args:
            value (scalar, list or np.ndarray): quantiy to be converted to numpy arraty if it is not yet a numpy array.

        Returns:
            [np.ndarray]: value converted to numpy array. If value is already a numpy array it returns value.
        """
        if isinstance(value, np.ndarray):
            return value
        else:
            # Necessary to deal with scalar: conversion to numpy array
            return np.array(value)

    # End method __convert_to_nparray

    def _fix_temperature_lb(self, temp, method_name):
        """Private method that checks the temperature values. If the minimum temperature is below the lower boundary, warns that the properties at temperature below the lower boundary are evaluated at lb.

        Args:
            temp (np.ndarray): temperature array at which evaluate material propetries in K.
            method_name (string): name of the method for which the temperature must be checked in order to select the sutable temperature range, which changes with the material properties to be evaluated.

        Returns:
            [np.ndarray]: modified temperature array with a cut off the temperature below the lower boundary.
        """

        if temp.min() < self.temperature_ranges[method_name][0]:
            warnings.warn(
                f"Temperaure below the validity range. {self.__class__.__name__} {self.property_name} evaluated at the lower bound temperature {self.temperature_ranges[method_name][0] = } K for those temperature values. \n"
            )

            temp[
                temp < self.temperature_ranges[method_name][0]
            ] = self.temperature_ranges[method_name][0]
        return temp

    # End method _fix_temperature_lb

    def _fix_temperature_ub(self, temp, method_name):
        """Private method that checks the temperature values. If the maximum temperature is above the upper boundary, warns that the properties at temperature above the upper boundary are evaluated at ub.

        Args:
            temp (np.ndarray): temperature array at which evaluate material propetries in K.
            method_name (string): name of the method for which the temperature must be checked in order to select the sutable temperature range, which changes with the material properties to be evaluated.

        Returns:
            [np.ndarray]: modified temperature array with a cut off the temperature above the upper boundary.
        """

        if temp.max() > self.temperature_ranges[method_name][1]:
            warnings.warn(
                f"Temperaure above the validity range. {self.__class__.__name__} {self.property_name} evaluated at the upper bound temperature {self.temperature_ranges[method_name][1] = } K for those temperature values. \n"
            )

            temp[
                temp > self.temperature_ranges[method_name][1]
            ] = self.temperature_ranges[method_name][1]
        return temp

    # End method _fix_temperature_ub

    def _preprocessing(self, temp, method_name):
        f"""Priovate method that preproces temp applying private methods {self._convert_to_nparray.__name__}, {self._fix_temperature_lb.__name__} and {self._fix_temperature_ub.__name__}. The result of the pre processing is stored in attribute self.temperature that is used in the method called to evaluate material properties.

        Args:
            temp (np.ndarray): temperature array at which evaluate material propetries in K.
            method_name (string): name of the method for which the temperature must be checked in order to select the sutable temperature range, which changes with the material properties to be evaluated.
        """
        # make a copy of temp to presere from changes the array passed to the function.
        temp = self._convert_to_nparray(temp)
        temp = temp.copy()
        temp = self._fix_temperature_lb(temp, method_name)
        temp = self._fix_temperature_ub(temp, method_name)

        self.temperature = temp

    # End method _preprocessing

    def aluminium(self):
        """Method Object to build instance of class Aluminium.

        Returns:
            [Aluminium]: instance of class Aluminium.
        """
        return Aluminium()

    def epoxy(self):
        """Method Object to build instance of class Epoxy.

        Returns:
            [Epoxy]: instance of class Epoxy.
        """
        return Epoxy()

    # End methdo epoxy

    def glass_epoxy(self, source=0):
        """Method Object to build instance of class GlassEpoxy.

        Args:
            source ([type], optional): Flag to select the formulation used to evaluate properties: 0 for Cryosoft formulation, 1 for NIST formulation. Defaults to 0.

        Returns:
            [GlassEpoxy]: instance of class GlassEpoxy.
        """
        return GlassEpoxy(source=source)

    # End method glass_epoxy

    def hastelloy(self):
        """Method Object to build instance of class Hastelloy.

        Returns:
            [Hastelloy]: instance of class Hastelloy.
        """
        return Hastelloy()

    # End method hastelloy

    def high_mn_austenitic_steinless_steel_jk2lb(self):
        """Method Object to build instance of class HighMnAusteniticSteinlessSteelJK2LB.

        Returns:
            [HighMnAusteniticSteinlessSteelJK2LB]: instance of class HighMnAusteniticSteinlessSteelJK2LB.
        """
        return HighMnAusteniticSteinlessSteelJK2LB()

    def incoloy(self, source=0):
        """Method Object to build instance of class Incoloy.

        Args:
            source ([type], optional): Flag to select the formulation used to evaluate properties: 0 for Cryosoft formulation, 1 for KSTAR formulation. Defaults to 0.
        Returns:
            [Incoloy]: instance of class Incoloy.
        """
        return Incoloy(source)

    # End method incoloy

    def inconel(self):
        """Method Object to build instance of class Iconel.

        Returns:
            [Inconel]: instance of class Inconel.
        """
        return Inconel()

    # End method inconel

    def kapton(self, source=0):
        """Method Object to build instance of class Kapton.

        Returns:
            [Kapton]: instance of class Kapton.
        """
        return Kapton(source)

    # End method inconel


# End class TemperatureDependentMaterial


class Aluminium(TemperatureDependentMaterial):
    """Class to evaluate Aluminium properties: density, electrical resistivity, isobaric specific heat and thermal conductivity.

    Args:
        TemperatureDependentMaterial ([class]): parent class of Class Aluminium with some useful private method.
    """

    def __init__(self):
        """Make an instance of class Aluminium."""
        self.temperature_ranges = dict(
            electrical_resistivity=(2.0, 1000.0),  # K
            isobaric_specific_heat=(2.0, 1000.0),  # K
            thermal_conductivity=(2.0, 1000.0),  # K
        )
        self.material = self.__class__.__name__

    # End method __init__

    def _build_interval(self, points):
        """Private methods used to build the interval to be passed to np.picewise in order to evaluate the aluminium properties of interest.

        Args:
            points (np.ndarray): array of points used to indentify the intervals

        Returns:
            [list]: list of the intervals obtained starting from the points.
        """
        # Build the intervals.
        intervals = [
            (self.temperature > val) & (self.temperature <= points[ii + 1])
            for ii, val in enumerate(points)
            if ii < len(points) - 1
        ]
        # Overwrite the first interval for which the correct writing is:
        # (temp >= p[0]) & (temp <= p[1])
        return [
            (self.temperature >= points[0]) & (self.temperature <= points[1]),
            *intervals[1:],
        ]

    # End method _build_interval

    def thermal_conductivity(self, temp):
        """Method to evaluate the thermal conductivity of aluminium alloy (RRR=5.5) in W/m/K, as a function of temperature T, for 2 <= T <= 1000 K.
        References: P.Reed and A.F.Clark, Materials at low Temperature, ASM, 1983 (hand picked-up data points, but good accuracy).

        Author: R.Heller @ ITP, Version: 2.0   10.10.1991

        Translation to Python: Daniele Placido.

        Args:
            temp (scalar, list or np.ndarray): temperature at which evaluate the thermal conductivity in K; converted to a np.ndarray if it is not a np.ndarray.

        Returns:
            [np.ndarray]: aluminium thermal conductivity in W/m/K; array as the same shape of temp.
        """

        self.property_name = self.thermal_conductivity.__name__.replace(
            "_", " "
        ).capitalize()

        self._preprocessing(temp, self.thermal_conductivity.__name__)

        piecewise_points = np.array(
            [
                2.0,
                4.0,
                6.0,
                8.0,
                10.0,
                20.0,
                40.0,
                60.0,
                80.0,
                100.0,
                200.0,
                300.0,
                1000.0,
            ]
        )

        intervals = self._build_interval(piecewise_points)

        behavior = [
            lambda TT: 4.2345 * (TT - 2.0) + 6.6510,
            lambda TT: 4.2000 * (TT - 4.0) + 15.1200,
            lambda TT: 4.3650 * (TT - 6.0) + 23.5200,
            lambda TT: 4.8550 * (TT - 8.0) + 32.2500,
            lambda TT: 4.3430 * (TT - 10.0) + 41.9600,
            lambda TT: 3.8805 * (TT - 20.0) + 85.3900,
            lambda TT: 2.0375 * (TT - 40.0) + 163.000,
            lambda TT: 0.5360 * (TT - 60.0) + 203.75,
            lambda TT: 0.1145 * (TT - 80.0) + 214.47,
            lambda TT: -0.0452 * (TT - 100.0) + 216.76,
            lambda TT: -0.0291 * (TT - 200.0) + 212.24,
            lambda TT: 209.3300,
        ]
        return np.piecewise(self.temperature, intervals, behavior)

    # End method thermal_conductivity

    def isobaric_specific_heat(self, temp):
        """Method to evaluate the isobaric specific heat of alumimium, in J/Kg K as a function of temperature T, for 2 <= T <= 1000 K.
        Reference: Cryocomp version 3.0 V.J. Johnson, Properties of materials at low temperatures (phase 1) Pergamon Press, 1961.

        Author: R.Heller @ ITP, Version: 2.0   10.10.1991

        Translation to Python: Daniele Placido.

        Args:
            temp (scalar, list or np.ndarray): temperature at which evaluate the isobaric specific heat in K; converted to a np.ndarray if it is not a np.ndarray.

        Returns:
            [np.ndarray]: aluminium isobaric specific heat in J/kg/K; array as the same shape of temp.
        """

        self.property_name = self.isobaric_specific_heat.__name__.replace(
            "_", " "
        ).capitalize()
        self._preprocessing(temp, self.isobaric_specific_heat.__name__)

        T0 = 11.0578562
        T1 = 53.2685055
        a1 = 0.04257809
        a2 = 0.00340313
        a3 = 0.00063472
        b0 = -2.50413265
        b1 = 0.71747923
        b2 = -0.06270275
        b3 = 0.0031683
        b4 = -2.0154e-05
        AA = 11458.0127
        BB = -438.603857
        CC = 331.057118
        DD = -10001.7006
        a = -7.12302329
        b = -34.6373149
        c = -29.4649313
        d = -3.30673113
        na = 0.9531781
        nb = 2.00135991
        nc = 2.974805
        nd = 3.94384164

        piecewise_points = [
            self.temperature_ranges[self.isobaric_specific_heat.__name__][0],
            T0,
            T1,
            self.temperature_ranges[self.isobaric_specific_heat.__name__][1],
        ]

        intervals = self._build_interval(piecewise_points)

        behaviour = [
            lambda TT: (a1 * TT + a2 * TT ** 2 + a3 * TT ** 3),
            lambda TT: (b0 + b1 * TT + b2 * TT ** 2 + b3 * TT ** 3 + b4 * TT ** 4),
            lambda TT: (
                AA * TT / (a + TT) ** na
                + BB * TT ** 2 / (b + TT) ** nb
                + CC * TT ** 3 / (c + TT) ** nc
                + DD * TT ** 4 / (d + TT) ** nd
            ),
        ]

        return np.piecewise(self.temperature, intervals, behaviour)

    # End method isobaric_specific_heat

    def electrical_resistivity(self, temp):
        """Method to evaluate the electrical resistivity of aluminium in Ohm*m, as a function of temperature T, for 2 <= T <= 1000 K.
        References: P.Reed & A.F.Clark, Materials at low Temperature, ASM, 1983 (hand picked-up data points, but good accuracy).

        Author: R.Heller @ ITP, Version: 2.0   10.10.1991

        Translation to Python: Daniele Placido.

        Args:
            temp (scalar, list or np.ndarray): temperature at which evaluate the electrical resistivity in K; converted to a np.ndarray if it is not a np.ndarray.

        Returns:
            [np.ndarray]: aluminium electrical resistivity in Ohm*m; array as the same shape of temp.
        """

        self.property_name = self.electrical_resistivity.__name__.replace(
            "_", " "
        ).capitalize()

        self._preprocessing(temp, self.electrical_resistivity.__name__)

        piecewise_points = np.array(
            [
                0.200000e1,
                0.400000e1,
                0.600000e1,
                0.800000e1,
                0.100000e2,
                0.200000e2,
                0.400000e2,
                0.600000e2,
                0.800000e2,
                0.100000e3,
                0.200000e3,
                0.300000e3,
                0.500000e3,
                0.100000e4,
            ]
        )
        rhoel = np.array(
            [
                0.576000e-8,
                0.576000e-8,
                0.576000e-8,
                0.576000e-8,
                0.576000e-8,
                0.576000e-8,
                0.600000e-8,
                0.720000e-8,
                0.912000e-8,
                0.112800e-7,
                0.230400e-7,
                0.350400e-7,
                0.584800e-7,
                0.116800e-6,
            ]
        )
        coeff = np.array(
            [
                [
                    -0.341170e-14,
                    0.170585e-14,
                    -0.341170e-14,
                    0.119409e-13,
                    -0.443521e-13,
                    0.472520e-12,
                    0.332536e-10,
                    0.825132e-10,
                    0.104684e-09,
                    0.110712e-09,
                    0.120792e-09,
                    0.118922e-09,
                    0.115284e-09,
                ],
                [
                    0.255877e-14,
                    0.123260e-30,
                    -0.255877e-14,
                    0.102351e-13,
                    -0.383816e-13,
                    0.900689e-13,
                    0.154898e-11,
                    0.913994e-12,
                    0.195038e-12,
                    0.105852e-12,
                    -0.505247e-14,
                    -0.136421e-13,
                    -0.454738e-14,
                ],
                [
                    -0.426462e-15,
                    -0.426462e-15,
                    0.213231e-14,
                    -0.810279e-14,
                    0.428168e-14,
                    0.243153e-13,
                    -0.105832e-13,
                    -0.119826e-13,
                    -0.148644e-14,
                    -0.369682e-15,
                    -0.286322e-15,
                    0.151579e-15,
                    0.151579e-15,
                ],
            ]
        )

        intervals = self._build_interval(piecewise_points)

        # List compreension to build behavior does not work for some unclear reason.
        # behavior = [lambda TT: (np.sum(coeff[:, ii]) * (TT - piecewise_points[ii]) + rhoel[ii]) for ii, _ in enumerate(intervals)]

        behavior = [
            lambda TT: (np.sum(coeff[:, 0]) * (TT - piecewise_points[0]) + rhoel[0]),
            lambda TT: (np.sum(coeff[:, 1]) * (TT - piecewise_points[1]) + rhoel[1]),
            lambda TT: (np.sum(coeff[:, 2]) * (TT - piecewise_points[2]) + rhoel[2]),
            lambda TT: (np.sum(coeff[:, 3]) * (TT - piecewise_points[3]) + rhoel[3]),
            lambda TT: (np.sum(coeff[:, 4]) * (TT - piecewise_points[4]) + rhoel[4]),
            lambda TT: (np.sum(coeff[:, 5]) * (TT - piecewise_points[5]) + rhoel[5]),
            lambda TT: (np.sum(coeff[:, 6]) * (TT - piecewise_points[6]) + rhoel[6]),
            lambda TT: (np.sum(coeff[:, 7]) * (TT - piecewise_points[7]) + rhoel[7]),
            lambda TT: (np.sum(coeff[:, 8]) * (TT - piecewise_points[8]) + rhoel[8]),
            lambda TT: (np.sum(coeff[:, 9]) * (TT - piecewise_points[9]) + rhoel[9]),
            lambda TT: (np.sum(coeff[:, 10]) * (TT - piecewise_points[10]) + rhoel[10]),
            lambda TT: (np.sum(coeff[:, 11]) * (TT - piecewise_points[11]) + rhoel[11]),
            lambda TT: (np.sum(coeff[:, 12]) * (TT - piecewise_points[12]) + rhoel[12]),
        ]

        return np.piecewise(self.temperature, intervals, behavior)

    # End method electrical_resistivity

    def density(self):
        """Method to evaluate aluminium density in kg/m^3. It is assumed independent from the temperature.

        Returns:
            [scalar]: aluminium density in kg/m^3.
        """
        self.property_name = self.density.__name__.capitalize()
        return 2700.0

    # End method density


# End class Aluminium


class Epoxy(TemperatureDependentMaterial):
    """Class to evaluate Epoxy properties: density, isobaric specific heat and thermal conductivity.

    Args:
        TemperatureDependentMaterial ([class]): parent class of Class Epoxy with some useful private method.
    """

    def __init__(self):
        """Make an instance of class Epoxy."""
        self.temperature_ranges = dict(
            isobaric_specific_heat=(1.0, 833.0),  # K
            thermal_conductivity=(1.0, 833.0),  # K
        )
        self.material = self.__class__.__name__

    # End method __init__

    def thermal_conductivity(self, temp):
        """Method to evaluate the thermal conductivity of epoxy in W/m/K, as a function of temperature T, for 1 <= T <= 833 K.

        Source: Cryosoft.

        Translation to Python: Daniele Placido.

        Args:
            temp (scalar, list or np.ndarray): temperature at which evaluate the thermal conductivity in K; converted to a np.ndarray if it is not a np.ndarray.

        Returns:
            [np.ndarray]: epoxy thermal conductivity in W/m/K; array as the same shape of temp.
        """

        self.property_name = self.thermal_conductivity.__name__.replace(
            "_", " "
        ).capitalize()

        self._preprocessing(temp, self.thermal_conductivity.__name__)

        AA = -0.37683848
        BB = 100.999824
        CC = 0.01286509
        DD = 0.09211400
        ALPHA = 4483.59742
        BETA = 6977179.28
        GAMMA = 9.45111414
        DELTA = 0.43241200
        NA = 0.48088594
        NB = 2.64828981
        NC = 2.10623669
        ND = 3.89389714

        return (
            AA * self.temperature / (ALPHA + self.temperature) ** NA
            + (BB * self.temperature ** 2) / (BETA + self.temperature) ** NB
            + CC * self.temperature ** 3 / (GAMMA + self.temperature) ** NC
            + DD * self.temperature ** 4 / (DELTA + self.temperature) ** ND
        )

    # End method thermal_conductivity

    def isobaric_specific_heat(self, temp):
        """Method to evaluate the isobric specific heat of epoxy in J/kg/K, as a function of temperature T, for 1 <= T <= 833 K.

        Source: Cryosoft.

        Translation to Python: Daniele Placido.

        Args:
            temp (scalar, list or np.ndarray): temperature at which evaluate the isobaric specific heat in K; converted to a np.ndarray if it is not a np.ndarray.

        Returns:
            [np.ndarray]: epoxy ispbaric specific heat in J/kg/K; array as the same shape of temp.
        """

        # epoxy (check tmax-tmin)
        # bf coefficients of the polynomial expression of cpep

        self.property_name = self.isobaric_specific_heat.__name__.replace(
            "_", " "
        ).capitalize()

        self._preprocessing(temp, self.isobaric_specific_heat.__name__)

        AA = 10.1978618
        BB = 0.00021728
        CC = 0.00760138
        DD = 1.8199e-06
        ALPHA = 0.68320779
        BETA = 14.8392381
        GAMMA = 1.92507754
        DELTA = 96.6186651
        NA = 7.03975410
        NB = 5.72616511
        NC = 14.7550871
        ND = 4.62647203

        return (
            AA * (self.temperature / (1.0 + self.temperature / ALPHA)) ** NA
            + BB * (self.temperature / (1.0 + self.temperature / BETA)) ** NB
            + CC * (self.temperature / (1.0 + self.temperature / GAMMA)) ** NC
            + DD * (self.temperature / (1.0 + self.temperature / DELTA)) ** ND
        )

    # End method isobaric_specific_heat

    def density(self):
        """Method to evaluate epoxy density in kg/m^3. It is assumed independent from the temperature.

        Reference: doi = 10.1088/1742-6596/1240/1/012080, 07/2019, IOP Publishing, 1240, 1, 012080, Sunny Bhatia and Surjit Angra and Sabah Khan, Mechanical and wear properties of epoxy matrix composite reinforced with varying ratios of solid glass microspheres, Journal of Physics: Conference Series

        Returns:
            [scalar]: epoxy density in kg/m^3.
        """
        self.property_name = self.density.__name__.capitalize()

        return 1.15e3

    # End method density


# End class Epoxy


class GlassEpoxy(TemperatureDependentMaterial):
    """Class to evaluate GlassEpoxyCryosoft properties: density, isobaric specific heat and thermal conductivity.

    Args:
        TemperatureDependentMaterial ([class]): parent class of Class GlassEpoxyCryosoft with some useful private method.
    """

    def __init__(self, source=CRYOSOFT):
        """Make an instance of class GlassEpoxyCryosoft.

        Args:
            source ([type], optional): Flag to select the formulation with which evaluate the glass epoxy properties: CRYOSOFT for Cryosoft formulation and NIST for NIST formulation. Defaults to CRYOSOFT.
        """

        self.material = self.__class__.__name__
        sources = {CRYOSOFT: "Cryosoft", NIST: "NIST"}
        default_ranges = {CRYOSOFT: (1.0, 1000.0), NIST: (4.0, 300.0)}
        self.source_name = sources[source]
        self.source = source
        self.temperature_ranges = dict(
            isobaric_specific_heat=default_ranges[source],  # K
            thermal_conductivity=default_ranges[source],  # K
        )
        self._curve_fit_equation()

    # End method __init__

    def _curve_fit_equation(self):
        """Private method to switch the fit curve equation for the evaluation of isobaric specific heat and thermal conductivity according to the selected source; assigned to private attribute _fit_equation."""

        switch = {CRYOSOFT: self._fit_cryosoft, NIST: self._fit_nist}

        self._fit_equation = switch[self.source]

    def _fit_cryosoft(self, **kwargs):
        """Private method whit the correlation used to evaluate the thermal conductivity and the isobaric specific heat of glass epoxy according to Cryosoft formulation in temperature range [1.0, 1000.0] K.

        Source: Cryosoft.

        Translation to Python: Daniele Placido

        Returns:
            [np.ndarray]: the array of the thermal conductivity or of the isobaric specific heat accordig to the coefficients passed to **kwargs.
        """
        # Fit curve equation:
        # aa * t / (a + t) ** na + bb * t ** 2 / (b + t) ** nb
        # + cc * t ** 3 / (c + t) ** nc + dd * t ** 4 / (d + t) ** nd
        return (
            kwargs["AA"]
            * self.temperature
            / (kwargs["ALPHA"] + self.temperature) ** kwargs["NA"]
            + kwargs["BB"]
            * self.temperature ** 2
            / (kwargs["BETA"] + self.temperature) ** kwargs["NB"]
            + kwargs["CC"]
            * self.temperature ** 3
            / (kwargs["GAMMA"] + self.temperature) ** kwargs["NC"]
            + kwargs["DD"]
            * self.temperature ** 4
            / (kwargs["DELTA"] + self.temperature) ** kwargs["ND"]
        )

    def _fit_nist(self, **kwargs):
        """Private method whit the correlation used to evaluate the thermal conductivity and the isobaric specific heat of glass epoxy according to NIST formulation in temperature range [4.0, 300.0] K.

        Source: NIST.

        Translation to Python: Daniele Placido

        Returns:
            [np.ndarray]: the array of the thermal conductivity or of the isobaric specific heat accordig to the coefficients passed to **kwargs.
        """
        # Fit curve equation:
        # 10 ** (aa + bb * log10(T) + cc * log10(T)**2 + dd * log10(T) ** 3 + ee * log10(T) ** 4 + ff * log10(T) ** 5 + gg * log10(T) ** 6 + hh * log10(T) ** 7 + ii * log10(T) ** 8)

        exponent = np.zeros(self.temperature.shape)

        for ii, value in enumerate(kwargs.values()):
            exponent += value * np.log10(self.temperature) ** ii

        return 10.0 ** exponent

    def thermal_conductivity(self, temp):
        """Method to evaluate the thermal conductivity according to the selected source (Cryosoft or NIST) in W/m/K.

        Args:
            temp (scalar, list or np.ndarray): temperature at which evaluate the thermal conductivity in K; converted to a np.ndarray if it is not a np.ndarray.

        Returns:
            [np.ndarray]: glass epoxy thermal conductivity in W/m/K; array as the same shape of temp.
        """

        self.property_name = self.thermal_conductivity.__name__.replace(
            "_", " "
        ).capitalize()

        self._preprocessing(temp, self.thermal_conductivity.__name__)

        coefficients = {
            CRYOSOFT: dict(
                AA=0.106635834,
                BB=4.254937095,
                CC=0.369710568,
                DD=-0.220612963,
                ALPHA=4.374557267,
                BETA=78.36483015,
                GAMMA=0.896107679,
                DELTA=0.261916360,
                NA=0.551227580,
                NB=2.274402853,
                NC=2.979875183,
                ND=3.644905343,
            ),
            NIST: dict(
                AA=-2.64827,  # W/m/K
                BB=8.80228,  # W/m/K
                CC=-24.8998,  # W/m/K
                DD=41.1625,  # W/m/K
                EE=-39.8754,  # W/m/K
                FF=23.1778,  # W/m/K
                GG=-7.95635,  # W/m/K
                HH=1.48806,  # W/m/K
                II=-0.11701,  # W/m/K
            ),
        }

        # self._fit_equation is either _fit_cryosoft or _fit_nist according to
        # self.source; pass to it the method the unpacked dictionary with the
        # coefficients for the fit equation; the method returns the thermal
        # conductivity in W/m/K.
        return self._fit_equation(**coefficients[self.source])

    # End method thermal_conductivity

    def isobaric_specific_heat(self, temp):
        """Method to evaluate the isobaric specific heat according to the selected source (Cryosoft or NIST) in J/kg/K.

        Args:
            temp (scalar, list or np.ndarray): temperature at which evaluate the isobaric specific heat in K; converted to a np.ndarray if it is not a np.ndarray.

        Returns:
            [np.ndarray]: glass epoxy isobaric specific heat in J/kg/K; array as the same shape of temp.
        """

        # %glass-epoxy (check tmax-tmin)

        self.property_name = self.isobaric_specific_heat.__name__.replace(
            "_", " "
        ).capitalize()

        self._preprocessing(temp, self.isobaric_specific_heat.__name__)

        coefficients = {
            CRYOSOFT: dict(
                AA=0.309399374,
                BB=-1387.509141,
                CC=283.0320589,
                DD=1.822153805,
                ALPHA=113.1687618,
                BETA=80.88333136,
                GAMMA=26.27075203,
                DELTA=1.104229335,
                NA=0.254475819,
                NB=1.938523579,
                NC=2.603466143,
                ND=3.522809815,
            ),
            NIST: dict(
                AA=-2.4083,  # J/kg/K
                BB=7.6006,  # J/kg/K
                CC=-8.2982,  # J/kg/K
                DD=7.3301,  # J/kg/K
                EE=-4.2386,  # J/kg/K
                FF=1.4294,  # J/kg/K
                GG=-0.24396,  # J/kg/K
                HH=0.015236,  # J/kg/K
                II=0.0,  # J/kg/K
            ),
        }

        # self._fit_equation is either _fit_cryosoft or _fit_nist according to
        # self.source; pass to it the method the unpacked dictionary with the
        # coefficients for the fit equation; ; the method returns the isobaric specific heat in J/kg/K.
        return self._fit_equation(**coefficients[self.source])

    # End method isobaric_specific_heat

    def density(self):
        """Method to evaluate epoxy density according to the selected source (Cryosoft or NIST) in kg/m^3. It is assumed independent from the temperature.

        Returns:
            [scalar]: glass epoxy density in kg/m^3.
        """
        self.property_name = self.density.__name__.capitalize()

        rho = {CRYOSOFT: 2e3, NIST: 1.8e3}
        return rho[self.source]

    # End method density


# End class GlassEpoxy


class Hastelloy(TemperatureDependentMaterial):
    """Class to evaluate Hastelloy properties: density, electrical resitivity, isobaric specific heat and thermal conductivity.

    Args:
        TemperatureDependentMaterial ([class]): parent class of Class Hastelloy with some useful private method.
    """

    def __init__(self):
        """Make an instance of class Hatelloy."""

        self.temperature_ranges = dict(
            electrical_resistivity=(4.0, 300.0),  # K
            isobaric_specific_heat=(1.0, 350.0),  # K
            thermal_conductivity=(1.0, 350.0),  # K
        )
        self.material = self.__class__.__name__

    # End method __init__

    def thermal_conductivity(self, temp):
        """Method to evaluate the thermal conductivity of hastelloy-C276 in W/m/K, as a function of temperature T, for 1 <= T <= 350 K.

        Source: "Experimental Analysis and Numerical Simulation of Quench in Superconducting HTS Tapes and Coils", Ph. D. Thesis by Marco Casali, 2014.

        Translation to Python: Daniele Placido

        Args:
            temp (scalar, list or np.ndarray): temperature at which evaluate the thermal conductivity in K; converted to a np.ndarray if it is not a np.ndarray.

        Returns:
            [np.ndarray]: hastelloy-C276 thermal conductivity in W/m/K; array as the same shape of temp.
        """

        self.property_name = self.thermal_conductivity.__name__.replace(
            "_", " "
        ).capitalize()

        self._preprocessing(temp, self.thermal_conductivity.__name__)

        intervals = [
            (self.temperature <= 60.0),
            (self.temperature > 60.0) & (self.temperature < 300.0),
            (self.temperature >= 300.0),
        ]
        behavior = [
            lambda T: 0.000028975912976 * T ** 3
            - 0.004655231435231 * T ** 2
            + 0.293992451992451 * T
            + 0.166939393939405,
            lambda T: 0.000000124577333 * T ** 3
            - 0.000079993515434 * T ** 2
            + 0.035243632195403 * T
            + 5.503358179018721,
            lambda T: (12.240619439579623 - 12.219768038241552)
            / (300 - 299)
            * (T - 300)
            + 12.240619439579623,
        ]

        return np.piecewise(self.temperature, intervals, behavior)

    # End method thermal_conductivity

    def isobaric_specific_heat(self, temp):
        """Method to evaluate the isobaric specific heat of hastelloy-C276 in J/kg/K, as a function of temperature T, for 1 <= T <= 350 K.

        Source: "Experimental Analysis and Numerical Simulation of Quench in Superconducting HTS Tapes and Coils", Ph. D. Thesis by Marco Casali, 2014.

        Translation to Python: Daniele Placido

        Args:
            temp (scalar, list or np.ndarray): temperature at which evaluate the isobaric specific heat in K; converted to a np.ndarray if it is not a np.ndarray.

        Returns:
            [np.ndarray]: hastelloy-C276 ispbaric specific heat in J/kg/K; array as the same shape of temp.
        """

        self.property_name = self.isobaric_specific_heat.__name__.replace(
            "_", " "
        ).capitalize()

        self._preprocessing(temp, self.isobaric_specific_heat.__name__)

        intervals = [
            (self.temperature <= 60.0),
            (self.temperature > 60.0) & (self.temperature < 300.0),
            (self.temperature >= 300.0),
        ]
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

        return np.piecewise(self.temperature, intervals, behavior)

    # End methdo isobaric_specific_heat

    def electrical_resistivity(self, temp):
        """Method to evaluate the electrical resistivity of hastelloy-C276 in Ohm*m, as a function of temperature T, for 4 <= T <= 300 K.

        Source: "Experimental Analysis and Numerical Simulation of Quench in Superconducting HTS Tapes and Coils", Ph. D. Thesis by Marco Casali, 2014.

        Translation to Python: Daniele Placido

        Args:
            temp (scalar, list or np.ndarray): temperature at which evaluate the electrical resistivity in K; converted to a np.ndarray if it is not a np.ndarray.

        Returns:
            [np.ndarray]: hastelloy-C276 electrical resistivity in Ohm*m; array as the same shape of temp.
        """

        self.property_name = self.electrical_resistivity.__name__.replace(
            "_", " "
        ).capitalize()

        self._preprocessing(temp, self.electrical_resistivity.__name__)

        return 1.309449e-10 * self.temperature + 1.220104e-06

    # End mehod electrical_resistivity

    def density(self):
        """Method to evaluate hastelloy-c276 density in kg/m^3. It is assumed independent from the temperature.

        Returns:
            [scalar]: hastelloy-c276 density in kg/m^3.
        """
        self.property_name = self.density.__name__.capitalize()
        return 8.89e3

    # End method density


# End class Hastelloy


class HighMnAusteniticSteinlessSteelJK2LB(TemperatureDependentMaterial):
    """Class to evaluate JK2LB (high Mn austenitic steinless steel) properties: density, electrical resitivity, isobaric specific heat and thermal conductivity.

    Args:
        TemperatureDependentMaterial ([class]): parent class of Class HighMnAusteniticSteinlessSteelJK2LB with some useful private method.
    """

    def __init__(self):
        """Make an instance of class HighMnAusteniticSteinlessSteelJK2LB"""
        self.temperature_ranges = dict(
            electrical_resistivity=(1.0, 1000.0),  # K
            isobaric_specific_heat=(1.0, 1000.0),  # K
            thermal_conductivity=(1.0, 1000.0),  # K
        )
        self.material = self.__class__.__name__

    # End method __init__

    def density(self):
        """Method to evaluate high Mn austenitic steinless steel (JK2LB) density in kg/m^3. It is assumed independent from the temperature.

        Returns:
            [scalar]: high Mn austenitic steinless steel (JK2LB) density in kg/m^3.
        """
        self.property_name = self.density.__name__.capitalize()
        return 7.86e3

    # End method density

    def thermal_conductivity(self, temp):
        """Method to evaluate the thermal conductivity of high Mn austenitic steinless steel (JK2LB) in W/m/K, as a function of temperature T, for 1 <= T <= 1000 K.

        Polinomial of degree 3 for 1.0 K<= T <= 20.0 K

        Polinomial of degree 6 for 20.0 K < T >= 300.0 K

        Linear interpolation for 3000.0 K < T >= 1000.0 K

        Reference: J. Lu, R.P. Walsh, K. Han Low temperature physical properties of a high Mn austenitic steel JK2LB, Cryogenics 49, pp. 133-137, 2009.

        Data extracted from the graphs in the paper

        Source: Cryosoft

        Translation to Python: Daniele Placido

        Args:
            temp (scalar, list or np.ndarray): temperature at which evaluate the thermal conductivity in K; converted to a np.ndarray if it is not a np.ndarray.

        Returns:
            [np.ndarray]: high Mn austenitic steinless steel (JK2LB) thermal conductivity in W/m/K; array as the same shape of temp.
        """

        self.property_name = self.thermal_conductivity.__name__.replace(
            "_", " "
        ).capitalize()
        self._preprocessing(temp, self.thermal_conductivity.__name__)

        T0 = 20.0
        T1 = 300.0

        # Coefficients (with 95 % confidence bounds):
        p1a = -3.495e-05  # (-4.704e-05, -2.286e-05)
        p2a = 0.001794  # (0.001323, 0.002265)
        p3a = 0.06738  # (0.06205, 0.07272)
        p4a = -0.03533  # (-0.05162, -0.01904)

        p1b = 2.174e-14  # (-1.541e-15, 4.503e-14)
        p2b = -1.451e-11  # (-3.663e-11, 7.621e-12)
        p3b = -2.417e-10  # (-8.404e-09, 7.92e-09)
        p4b = 2.426e-06  # (9.55e-07, 3.898e-06)
        p5b = -0.0007596  # (-0.0008933, -0.0006259)
        p6b = 0.1175  # (0.1118, 0.1231)
        p7b = -0.2904  # (-0.3743, -0.2065)

        cond250 = 9.7147
        cond300 = 10.7446
        dy = -0.0243

        intervals = [
            (self.temperature <= T0),
            (self.temperature > T0) & (self.temperature <= T1),
            (self.temperature > T1),
        ]

        behavior = [
            lambda TT: p1a * TT ** 3 + p2a * TT ** 2 + p3a * TT + p4a,
            lambda TT: p1b * TT ** 6
            + p2b * TT ** 5
            + p3b * TT ** 4
            + p4b * TT ** 3
            + p5b * TT ** 2
            + p6b * TT ** 1
            + p7b
            + dy,
            lambda TT: cond250 + dy + (cond300 - cond250) / (300 - 250) * (TT - 250),
        ]

        return np.piecewise(self.temperature, intervals, behavior)

    # End method thermal conductivity

    def isobaric_specific_heat(self, temp):
        """Method to evaluate the isobaric specific heat of high Mn austenitic steinless steel (JK2LB) in J/kg/K, as a function of temperature T, for 1 <= T <= 1000 K.

        Polinomial of degree 6 for 1.0 K <= T <= 50.0 K

        Polinomial of degree 3 for 50.0 K < T <= 300.0 K

        Linear interpolation for 300.0 K < T <= 1000.0 K

        Properties from: J. Lu, R.P. Walsh, K. Han Low temperature physical properties of a high Mn austenitic steel JK2LB Cryogenics 49, pp. 133-137, 2009.

        Data extracted from the graphs in the paper

        Translation to Python: Daniele Placido

        Args:
            temp (scalar, list or np.ndarray): temperature at which evaluate the isobaric specific heat in K; converted to a np.ndarray if it is not a np.ndarray.

        Returns:
            [np.ndarray]: hastelloy-C276 ispbaric specific heat in J/kg/K; array as the same shape of temp.
        """

        self.property_name = self.isobaric_specific_heat.__name__.replace(
            "_", " "
        ).capitalize()
        self._preprocessing(temp, self.isobaric_specific_heat.__name__)

        T0 = 50.0
        T1 = 300.0

        # Coefficients (with 95% confidence bounds):

        p1a = -4.18e-09  # (-1.016e-08, 1.796e-09)
        p2a = 3.04e-07  # (-6.271e-07, 1.235e-06)
        p3a = -1.039e-05  # (-6.585e-05, 4.508e-05)
        p4a = 0.001128  # (-0.000453, 0.002709)
        p5a = -0.01538  # (-0.03742, 0.006657)
        p6a = 0.3118  # (0.1757, 0.4478)
        p7a = -0.1974  # (-0.4711, 0.07628)
        p1b = 2.092e-05  # (1.806e-05, 2.378e-05)
        p2b = -0.01945  # (-0.02092, -0.01797)
        p3b = 6.149  # (5.925, 6.373)
        p4b = -179.3  # (-189, -169.5)

        cp299 = 479.7580
        cp300 = 479.8868
        dy = 0.5230

        intervals = [
            (self.temperature <= T0),
            (self.temperature > T0) & (self.temperature <= T1),
            (self.temperature > T1),
        ]

        behavior = [
            lambda TT: p1a * TT ** 6
            + p2a * TT ** 5
            + p3a * TT ** 4
            + p4a * TT ** 3
            + p5a * TT ** 2
            + p6a * TT
            + p7a,
            lambda TT: p1b * TT ** 3 + p2b * TT ** 2 + p3b * TT + p4b + dy,
            lambda TT: cp299 + dy + (cp300 - cp299) / (300 - 299) * (TT - 299),
        ]

        return np.piecewise(self.temperature, intervals, behavior)

    # End method isobaric_specific_heat

    def electrical_resistivity(self, temp):
        """Method to evaluate the electrical resistivity of high Mn austenitic steinless steel (JK2LB) in Ohm*m, as a function of temperature T, for 1 <= T <= 1000 K.

        Polinomial of degree 6 for 1.0 K <= T <= 300.0 K

        Linear interpolation for 300.0 K < T <= 1000.0 K

        Properties from: J. Lu, R.P. Walsh, K. Han, Low temperature physical properties of a high Mn austenitic steel JK2LB, Cryogenics 49, pp. 133-137, 2009.

        Data extracted from the graphs in the paper

        Translation to Python: Daniele Placido

        Args:
            temp (scalar, list or np.ndarray): temperature at which evaluate the electrical resistivity in K; converted to a np.ndarray if it is not a np.ndarray.

        Returns:
            [np.ndarray]: hastelloy-C276 electrical resistivity in Ohm*m; array as the same shape of temp.
        """

        self.property_name = self.electrical_resistivity.__name__.replace(
            "_", " "
        ).capitalize()
        self._preprocessing(temp, self.electrical_resistivity.__name__)

        T1 = 300.0
        # Coefficients (with 95% confidence bounds):
        p1a = -2.506e-21  # (-3.319e-21, -1.693e-21)
        p2a = 2.214e-18  # (1.489e-18, 2.939e-18)
        p3a = -6.662e-16  # (-9.124e-16, -4.201e-16)
        p4a = 7.45e-14  # (3.48e-14, 1.142e-13)
        p5a = -5.208e-13  # (-3.603e-12, 2.562e-12)
        p6a = -1.466e-10  # (-2.487e-10, -4.46e-11)
        p7a = 9.109e-07  # (9.098e-07, 9.119e-07)

        cond280 = 9.7261e-07
        cond300 = 9.8883e-07

        intervals = [(self.temperature <= T1), (self.temperature > T1)]

        behavior = [
            lambda TT: p1a * TT ** 6
            + p2a * TT ** 5
            + p3a * TT ** 4
            + p4a * TT ** 3
            + p5a * TT ** 2
            + p6a * TT ** 1
            + p7a,
            lambda TT: cond280 + (cond300 - cond280) / (300.0 - 280.0) * (TT - 280.0),
        ]

        return np.piecewise(self.temperature, intervals, behavior)

    # End method electrical_resistivity


# End class HighMnAusteniticSteinlessSteelJK2LB


class Incoloy(TemperatureDependentMaterial):
    """Class to evaluate incoloy (908) properties: density, electrical resitivity, isobaric specific heat and thermal conductivity.

    Args:
        TemperatureDependentMaterial ([class]): parent class of Class Incoloy with some useful private method.
    """

    def __init__(self, source=CRYOSOFT):
        self.material = self.__class__.__name__
        sources = {CRYOSOFT: "Cryosoft", KSTAR: "KSTAR"}
        default_ranges = dict(
            thermal_conductivity={CRYOSOFT: (1.0, 1000.0), KSTAR: (4.0, 300.0)},  # K
            isobaric_specific_heat={CRYOSOFT: (4.27, 1423.0), KSTAR: (4.2, 300.0)},  # K
            electrical_resistivity={CRYOSOFT: (1.0, 300.0), KSTAR: (4.0, 300.0)},  # K
        )
        self.source_name = sources[source]
        self.source = source
        # Assign temperature ranges according to the source.
        self.temperature_ranges = dict(
            electrical_resistivity=default_ranges["electrical_resistivity"][
                source
            ],  # K
            isobaric_specific_heat=default_ranges["isobaric_specific_heat"][
                source
            ],  # K
            thermal_conductivity=default_ranges["thermal_conductivity"][source],  # K
        )

        # to switch curve fit equation used to evaluate incoloy properties according to the chosen formulation.
        self.__curve_fit = dict(
            electrical_resistivity={
                CRYOSOFT: self.__electrical_resistivity_cryosoft,
                KSTAR: self.__fit_kstar,
            },
            isobaric_specific_heat={
                CRYOSOFT: self.__isobaric_specific_heat_cryosoft,
                KSTAR: self.__fit_kstar,
            },
            thermal_conductivity={
                CRYOSOFT: self.__thermal_conductivity_cryosoft,
                KSTAR: self.__fit_kstar,
            },
        )
        self.material = self.__class__.__name__

    # End method __init__

    def __electrical_resistivity_cryosoft(self, **kwargs):
        """Private method to evaluate Incoloy electrical resistivity according to the Cryosoft formulation in Ohm * m as function of temperature in range 1.0 K <= T <= 300.0 K.

        Note: the resistivity has been assumed to be equal to that of inconel. Accuracy unknown.


        Source: Cryosoft. A. Portone @ EFDA, Garching, 5.12.2000

        Translation to Python: Daniele Placido.


        Returns:
            [np.array]: Incoloy 908 electrical resistivity in Ohm * m according to Cryosoft formulation; array has the same shape of self.temperature.
        """
        # (a + b * t ** nb) * 1e-7
        return (kwargs["AA"] + kwargs["BB"] * self.temperature ** kwargs["NB"]) * 1e-7

    def __isobaric_specific_heat_cryosoft(self, **kwargs):
        """Private method to evaluate Incoloy isobaric specific hear according to the Cryosoft formulation in J/kg/k as function of temperature in range 4.27 K <= T <= 1423.0 K.

        Source: Cryosoft.
        Reference: L.S.Toma, M.M.Steeves, R.P.Reed, Incoloy Alloy 908 Data Handbook, PFC/RR-94-2, 1994

        Translation to Python: Daniele Placido.

        Returns:
            [np.array]: Incoloy 908 isobaric specific hear in J/kg/k according to Cryosoft formulation; array has the same shape of self.temperature.
        """

        # aa * t. ^ na / (1 + t / a). ^ na + bb * t. ^ nb / (1 + t / b). ^ nb
        # + cc * t. ^ nc / (1 + t / c). ^ nc + dd * t. ^ nd / (1 + t / d). ^ nd
        return (
            kwargs["AA"]
            * self.temperature ** kwargs["NA"]
            / (1.0 + self.temperature / kwargs["ALPHA"]) ** kwargs["NA"]
            + kwargs["BB"]
            * self.temperature ** kwargs["NB"]
            / (1.0 + self.temperature / kwargs["BETA"]) ** kwargs["NB"]
            + kwargs["CC"]
            * self.temperature ** kwargs["NC"]
            / (1.0 + self.temperature / kwargs["GAMMA"]) ** kwargs["NC"]
            + kwargs["DD"]
            * self.temperature ** kwargs["ND"]
            / (1.0 + self.temperature / kwargs["DELTA"]) ** kwargs["ND"]
        )

    # End method __isobaric_specific_heat_cryosoft

    def __thermal_conductivity_cryosoft(self, **kwargs):
        """Private method to evaluate Incoloy thermal conductivity according to the Cryosoft formulation in W/m/K as function of temperature in range 1.0 K <= T <= 1000.0 K.

        Source: Cryosoft.

        Translation to Python: Daniele Placido.

        Returns:
            [np.array]: Incoloy 908 thermal conductivity in W/m/K according to Cryosoft formulation; array has the same shape of self.temperature.
        """

        # aa * T / (alpha + T) ** na + bb * T ** 2 / (beta + T) ** nb
        # + cc * T ** 3 / (gamma + T) ** nc + dd * T ** 4 / (delta + T) ** nd
        return (
            kwargs["AA"]
            * self.temperature
            / (kwargs["ALPHA"] + self.temperature) ** kwargs["NA"]
            + kwargs["BB"]
            * self.temperature ** 2
            / (kwargs["BETA"] + self.temperature) ** kwargs["NB"]
            + kwargs["CC"]
            * self.temperature ** 3
            / (kwargs["GAMMA"] + self.temperature) ** kwargs["NC"]
            + kwargs["DD"]
            * self.temperature ** 4
            / (kwargs["DELTA"] + self.temperature) ** kwargs["ND"]
        )

    # End method __thermal_coductivity_cryosoft

    def __fit_kstar(self, **kwargs):
        """Private method to evaluate Incoloy thermal conductivity in W/m/K or isobaric specific heat in J/kg/K or electrical resistivity in Ohm * m according to the KSTAR formulation function of temperature. Temperature ranges:

        thermal conductivity 4.0 K <= T <= 300.0 K;

        isobaric specific heat 4.2 K <= T <= 300.0 K;

        electrical resistivity 4.0 K <= T <= 300.0 K;

        Source: Material properties of the KSTAR CICC Insulators

        Translation to Python: Daniele Placido


        Returns:
            [np.array]: Incoloy thermal conductivity in W/m/K or isobaric specific heat in J/kg/K or electrical resistivity in Ohm * m according to the KSTAR formulation; array has the same shape of self.temperature.
        """

        kk = np.zeros(self.temperature.shape)

        # aa + bb * T + cc * T ** 2 + dd * T ** 3 + ee * T ** 4 + ff * T ** 5
        for ii, value in enumerate(kwargs.values()):
            kk += value * self.temperature ** ii
        # end for

        return kk

    # End method __fit_kstar

    def density(self):
        """Method to evaluate incoloy 908 density in kg/m^3. It is assumed independent from the temperature.

        Returns:
            [scalar]: incoloy 908 density in kg/m^3.
        """
        self.property_name = self.density.__name__.capitalize()
        rho = {CRYOSOFT: 8.17e3, KSTAR: 8.17e3}
        return rho[self.source]

    # End method density

    def __evaluate_property(self, temp, method_name, fit, coefficient):
        """Private method that actually evaluates the electrical resistivty or the isobaric specific heat or the thermal conductivity of Incoloy 908.

        Args:
            temp (_scalar, list or np.ndarray_): _temperature at which evaluate the thermal conductivity in K; converted to a np.ndarray if it is not a np.ndarray._
            method_name (_string_): _name of the method/properties to be evaluated_
            fit (_method_): _private method with the fit curve equation to evaluate the property_
            coefficient (_dict_): _coefficients to be used with the selected curve fit equation_

        Returns:
           _np.ndarray_: _incoloy 908 electrical resitivity in Ohm*m or isobaric specific heat in J/kg/K or thermal conductivity in W/m/K; array as the same shape of temp._
        """

        self.property_name = method_name.replace("_", " ").capitalize()

        self._preprocessing(temp, method_name)

        return fit(**coefficient[self.source])

    # Ene method __evaluate_properties

    def thermal_conductivity(self, temp):
        """Method to evaluate the thermal conductivity of Incoloy 908 in W/m/K, as a function of temperature T, according to the selected formulation (Cryosoft or KSTAR).

        Translation to Python: Daniele Placido

        Args:
            temp (scalar, list or np.ndarray): temperature at which evaluate the thermal conductivity in K; converted to a np.ndarray if it is not a np.ndarray.

        Returns:
            [np.ndarray]: incoloy 908 thermal conductivity in W/m/K; array as the same shape of temp.
        """

        coefficients = {
            CRYOSOFT: dict(
                AA=644.74407,  # W/m/K
                BB=-141.184107,  # W/m/K
                CC=-0.85851077,  # W/m/K
                DD=0,  # W/m/K
                ALPHA=5154.07845,  # K
                BETA=439.708273,  # K
                GAMMA=1.93083487,  # K
                DELTA=0,  # K
                NA=1,  # -
                NB=2,  # -
                NC=3,  # -
                ND=4,  # -
            ),
            KSTAR: dict(
                AA=0.08469,
                BB=0.10489,
                CC=-4.08513e-4,
                DD=8.12423e-7,
                EE=-7.17421e-10,
                FF=2.30711e-13,
            ),
        }
        return self.__evaluate_property(
            temp,
            self.thermal_conductivity.__name__,
            self.__curve_fit[self.thermal_conductivity.__name__][self.source],
            coefficients,
        )

    # End method thermal_conductivity

    def isobaric_specific_heat(self, temp):
        """Method to evaluate the isobaric specific heat of Incoloy 908 in J/kg/K, as a function of temperature T, according to the selected formulation (Cryosoft or KSTAR).

        Translation to Python: Daniele Placido

        Args:
            temp (scalar, list or np.ndarray): temperature at which evaluate the isobaric specific heat in K; converted to a np.ndarray if it is not a np.ndarray.

        Returns:
            [np.ndarray]: incoloy 908 isobaric specific heat in J/kg/K; array as the same shape of temp.
        """

        coefficients = {
            CRYOSOFT: dict(
                AA=0.00219167,
                BB=0.1333227,
                CC=0.00018546,
                DD=-4.8176e-05,
                ALPHA=4813.28489,
                BETA=10.8112775,
                GAMMA=58.178214,
                DELTA=136.765771,
                NA=-0.22447346,  # -
                NB=1.41030517,  # -
                NC=3.84780106,  # -
                ND=3.30329979,  # -
            ),
            KSTAR: dict(
                AA=4.84789,
                BB=-0.15139,
                CC=0.019,
                DD=-4.27478e-5,
                EE=-6.32387e-8,
                FF=1.84889e-10,
            ),
        }

        return self.__evaluate_property(
            temp,
            self.isobaric_specific_heat.__name__,
            self.__curve_fit[self.isobaric_specific_heat.__name__][self.source],
            coefficients,
        )

    # End method isobaric_specific_heat

    def electrical_resistivity(self, temp):
        """Method to evaluate the electrical resistivity of Incoloy 908 in Ohm * m, as a function of temperature T, according to the selected formulation (Cryosoft or KSTAR).

        Translation to Python: Daniele Placido

        Args:
            temp (scalar, list or np.ndarray): temperature at which evaluate the electrical resistivity in K; converted to a np.ndarray if it is not a np.ndarray.

        Returns:
            [np.ndarray]: incoloy 908 electrical resistivity in Ohm * m; array as the same shape of temp.
        """

        coefficients = {
            CRYOSOFT: dict(
                AA=11.8,  # Ohm*m
                BB=0.00056648,  # Ohm*m/K**nb
                NB=1.18930588,  # -
            ),
            KSTAR: dict(
                AA=1.17806e-6,  # Ohm*m
                BB=1.99188e-10,  # Ohm*m/K
                CC=2.3553e-13,  # Ohm*m/K**2
                DD=0.0,
                EE=0.0,
                FF=0.0,
            ),
        }

        return self.__evaluate_property(
            temp,
            self.electrical_resistivity.__name__,
            self.__curve_fit[self.electrical_resistivity.__name__][self.source],
            coefficients,
        )

    # End method electrical_resistivity


# End class Incoloy


class Inconel(TemperatureDependentMaterial):
    """Class to evaluate Inconel 718 properties: density, electrical resitivity, isobaric specific heat and thermal conductivity.

    Args:
        TemperatureDependentMaterial ([class]): parent class of Class Inconel with some useful private method.
    """

    def __init__(self):
        """Make an instance of class Inconel"""
        self.temperature_ranges = dict(
            electrical_resistivity=(1.0, 300.0),  # K
            isobaric_specific_heat=(1.0, 1000.0),  # K
            thermal_conductivity=(1.0, 500.0),  # K
        )
        self.material = self.__class__.__name__

    # End method __init__

    def density(self):
        """Method to evaluate inconel 718 density in kg/m^3. It is assumed independent from the temperature.

        Returns:
            [scalar]: inconel 718 density in kg/m^3.
        """
        self.property_name = self.density.__name__.capitalize()
        return 8.17e3

    # End method density

    def thermal_conductivity(self, temp):
        """Method to evaluate the thermal conductivity of inconel 718 in W/m/K, as a function of temperature T, for 1 <= T <= 500 K.

        Reference: Cryocomp version 3.0, March 1997.

        Source: Cryosoft

        Translation to Python: Daniele Placido

        Args:
            temp (scalar, list or np.ndarray): temperature at which evaluate the thermal conductivity in K; converted to a np.ndarray if it is not a np.ndarray.

        Returns:
            [np.ndarray]: inconel 718 thermal conductivity in W/m/K; array as the same shape of temp.
        """
        self.property_name = self.thermal_conductivity.__name__.replace(
            "_", " "
        ).capitalize()
        self._preprocessing(temp, self.thermal_conductivity.__name__)

        AA = 15.7823796
        BB = -7.11845509
        CC = -0.49416263
        DD = 0.0017008
        ALPHA = 48.6907438
        BETA = 9.65510385
        GAMMA = 0.06978251
        DELTA = -0.54474549
        NA = 0.86455662
        NB = 2.14212483
        NC = 2.29137737
        ND = 2.44590736

        return (
            AA * self.temperature / (ALPHA + self.temperature) ** NA
            + BB * self.temperature ** 2 / (BETA + self.temperature) ** NB
            + CC * self.temperature ** 3 / (GAMMA + self.temperature) ** NC
            + DD * self.temperature ** 4 / (DELTA + self.temperature) ** ND
        )

    # End method thermal_conductivity

    def isobaric_specific_heat(self, temp):
        """Method to evaluate the isobaric specific heat of inconel 718 in J/kg/K, as a function of temperature T, for 1 <= T <= 1000 K.

        Reference: Cryocomp version 3.0, March 1997.

        Source: Cryosoft

        Translation to Python: Daniele Placido

        Args:
            temp (scalar, list or np.ndarray): temperature at which evaluate the isobaric specific heat in K; converted to a np.ndarray if it is not a np.ndarray.

        Returns:
            [np.ndarray]: inconel 718 isobaric specific heat in J/kg/K; array as the same shape of temp.
        """
        self.property_name = self.isobaric_specific_heat.__name__.replace(
            "_", " "
        ).capitalize()
        self._preprocessing(temp, self.isobaric_specific_heat.__name__)

        T0 = 22.1488673
        A1 = 0.41309481
        A2 = -0.01300475
        A3 = 0.00150658
        A4 = -0.00010202
        A5 = 3.5377e-06
        AA = 33.2554045
        BB = 0.00033495
        CC = -5.80523108
        DD = 53.5383062
        APLHA = 52.1477536
        BETA = 4.56140518
        GAMMA = 7.80826623
        DELTA = 0.32466724
        NA = 0.99906362
        NB = 1.99999986
        NC = 2.54414531
        ND = 5.9794021

        intervals = [(self.temperature <= T0), (self.temperature > T0)]
        behavior = [
            lambda TT: A1 * TT
            + A2 * TT ** 2
            + A3 * TT ** 3
            + A4 * TT ** 4
            + A5 * TT ** 5,
            lambda TT: AA * TT ** NA / (1.0 + TT / APLHA) ** NA
            + BB * TT ** NB / (1.0 + TT / BETA) ** NB
            + CC * TT ** NC / (1.0 + TT / GAMMA) ** NC
            + DD * TT ** ND / (1.0 + TT / DELTA) ** ND,
        ]
        return np.piecewise(self.temperature, intervals, behavior)

    # End method isobaric_specific_heat

    def electrical_resistivity(self, temp):
        """Method to evaluate the electrical resistivity of inconel 718 in Ohm*m, as a function of temperature T, for 1 <= T <= 300 K.

        Reference: Cryocomp version 3.0, March 1997.

        Source: Cryosoft

        Translation to Python: Daniele Placido

        Args:
            temp (scalar, list or np.ndarray): temperature at which evaluate the electrical resistivity in K; converted to a np.ndarray if it is not a np.ndarray.

        Returns:
            [np.ndarray]: inconel 718 electrical resistivity in Ohm*m; array as the same shape of temp.
        """
        self.property_name = self.electrical_resistivity.__name__.replace(
            "_", " "
        ).capitalize()
        self._preprocessing(temp, self.electrical_resistivity.__name__)
        AA = 11.8
        BB = 0.00056648
        NB = 1.18930588

        return (AA + BB * self.temperature ** NB) * 1.0e-7

    # End method electrical_resistivity


# End class Inconel


class Kapton(TemperatureDependentMaterial):
    """Class to evaluate Kapton properties: density, isobaric specific heat and thermal conductivity.

    Args:
        TemperatureDependentMaterial ([class]): parent class of Class Kapton with some useful private method.
    """

    def __init__(self, source=CRYOSOFT):
        """Make an instance of class Kapton.

        Args:
            source ([type], optional): Flag to select the formulation with which evaluate the kapton properties: CRYOSOFT for Cryosoft formulation and NIST for NIST formulation. Defaults to CRYOSOFT.
        """
        sources = {CRYOSOFT: "Cryosoft", NIST: "NIST"}
        self.source_name = sources[source]
        self.source = source

        default_ranges = dict(
            isobaric_specific_heat={CRYOSOFT: (1.0, 300.0), NIST: (4.0, 300.0)},  # K
            thermal_conductivity={CRYOSOFT: (1.0, 833.0), NIST: (4.0, 300.0)},  # K
        )
        self.temperature_ranges = dict(
            isobaric_specific_heat=default_ranges["isobaric_specific_heat"][
                source
            ],  # K
            thermal_conductivity=default_ranges["thermal_conductivity"][source],  # K
        )
        self.material = self.__class__.__name__

    # End method __init__

    def density(self):
        """Method to evaluate kapton density according to the selected source (Cryosoft or NIST) in kg/m^3. It is assumed independent from the temperature.

        Returns:
            [scalar]: glass epoxy density in kg/m^3.
        """
        self.property_name = self.density.__name__.capitalize()
        rho = {CRYOSOFT: 1.380e3, NIST: 1.420e3}
        return rho[self.source]

    # End method density