# Import library packages.
import warnings
import numpy as np
import pandas as pd
from scipy import constants
from CoolProp.CoolProp import PropsSI


class Environment:
    """docstring for Environment."""

    KIND = "Environment"

    def __init__(self, f_path):
        """[summary]

        Args:
            f_path ([type]): [description]
        """
        # Dictionary declaration (cdp, 11/2020)
        self.inputs = dict()
        # self.dict_node_pt = dict()
        # self.dict_Gauss_pt = dict()
        # Dictionary initialization: inputs.
        self.inputs = pd.read_excel(
            f_path,
            sheet_name="ENVIRONMENT",
            header=0,
            index_col=0,
            usecols=["Variable name", "Value"],
        )["Value"].to_dict()
        self.type = self.inputs["Medium"].lower()
        # Declare the dictionary with methods used to evaluate nusselt number.
        self.dict_nusselt_correlations = dict(
            vertical_plate=self._vertical_plate,
            vertical_plate_churchill_chu=self._vertical_plate_churchill_chu,
            vertical_plate_churchill_chu_accurate=self._vertical_plate_churchill_chu_accurate,
            long_horziontal_cylinder_morgan=self._long_horziontal_cylinder_morgan,
            long_horziontal_cylinder_churchill_chu=(
                self._long_horziontal_cylinder_churchill_chu
            ),
        )
        # dry air
        self.fluid_prop_aliases = dict(
            density="Dmass",
            thermal_conductivity="conductivity",
            volumetric_thermal_expansion_coefficient="isobaric_expansion_coefficient",
            dynamic_viscosity="viscosity",
            prandtl="Prandtl",
        )

    # End method __init__.

    def __repr__(self):
        """[summary]

        Returns:
            [type]: [description]
        """
        return f"{self.__class__.__name__}(Type: {self.KIND}"

    def __str__(self):
        """[summary]"""
        pass

    def eval_heat_transfer_coefficient(self, conductor, T_s):
        """[summary]

        Args:
            conductor ([type]): [description]
            T_s ([type]): [description]

        Returns:
            [type]: [description]
        """
        # Method that evaluates the external free convection heat transfer coefficient.

        # Declare dictionary with the characterisctic length.
        dict_characterisctic_length = dict(
            vertical_plate=conductor.grid_features["delta_z"],  # to be checked
            vertical_plate_churchill_chu=conductor.grid_features[
                "delta_z"
            ],  # to be checked
            vertical_plate_churchill_chu_accurate=conductor.grid_features[
                "delta_z"
            ],  # to be checked
            long_horziontal_cylinder_morgan=conductor.inputs["Diameter"],
            long_horziontal_cylinder_churchill_chu=conductor.inputs["Diameter"],
        )

        # Define the film temperature.
        film_temperature = (self.inputs["Temperature"] + T_s) / 2.0  # K
        # Evaluate air propreties.
        dict_air_properties = self.eval_prop(film_temperature)
        # Evaluate Nusselt dimensionless number.
        if conductor.inputs["Is_rectangular"]:
            # Evaluate Grashof dimensionless number for vertical side.
            grashof_side = self.grashof_number(
                dict_air_properties, T_s, conductor.inputs["Height"]
            )
            # Evaluate Rayleigh dimensionless number for vertical side.
            rayleigh_side = self.rayleigh_number(
                grashof_side, dict_air_properties["prandtl"]
            )
            # Evaluate Nussel dimensionless number for vertical side.
            nusselt_side = self._vertical_plate_churchill_chu_accurate(
                rayleigh_side, dict_air_properties["prandtl"], grashof_side, conductor
            )
            # L = A_s/P
            characteristic_length = (
                conductor.inputs["XLENGHT"]
                * conductor.inputs["Width"]
                / (2 * (conductor.inputs["XLENGHT"] + conductor.inputs["Width"]))
            )
            # Evaluate Grashof dimensionless number lower/upper cold plate.
            grashof_lu = self.grashof_number(
                dict_air_properties, T_s, characteristic_length
            )
            # Evaluate Rayleigh dimensionless number for lower/upper cold plate.
            rayleigh_lu = self.rayleigh_number(
                grashof_lu, dict_air_properties["prandtl"]
            )
            # Evaluate Nussel dimensionless number for lower cold plate.
            nusselt_bottom = self._horiziontal_lower_surface_cold_plate(
                rayleigh_lu, dict_air_properties["prandtl"]
            )
            # Evaluate Nussel dimensionless number for upper cold plate.
            nusselt_top = self._horiziontal_upper_surface_cold_plate(
                rayleigh_lu, dict_air_properties["prandtl"]
            )
            # Evaluate external free convection heat transfer coefficients: vertical sides, bottom and upper horiziontal sides.
            return (
                nusselt_side
                * dict_air_properties["thermal_conductivity"]
                / conductor.inputs["Height"],
                nusselt_bottom
                * dict_air_properties["thermal_conductivity"]
                / characteristic_length,
                nusselt_top
                * dict_air_properties["thermal_conductivity"]
                / characteristic_length,
            )
        else:
            # Non rectangular duct (cylinder).
            # Get the characterisctic length needed to evaluare Grashof and Rayleigh dimensionless numbers according to the selected external free convection correlation.
            characteristic_length = dict_characterisctic_length[
                conductor.inputs["external_free_convection_correlation"]
            ]
            # Evaluate Grashof dimensionless numbers.
            grashof = self.grashof_number(
                dict_air_properties, T_s, characteristic_length
            )
            # Evaluate Rayleigh number.
            rayleigh = self.rayleigh_number(grashof, dict_air_properties["prandtl"])
            nusselt = self.dict_nusselt_correlations[
                conductor.inputs["external_free_convection_correlation"]
            ](rayleigh, dict_air_properties["prandtl"], grashof, conductor)
            # Evaluate external free convection heat transfer coefficient.
            return (
                nusselt
                * dict_air_properties["thermal_conductivity"]
                / characteristic_length
            )
        # End if conductor.inputs["Is_rectangular"]

    # End method eval_heat_transfer_coefficient.

    def eval_prop(self, film_temperature):
        """Function eval_prop evaluates the envirmonment medium (eg. air) properties at the environment pressure and at the film temperature exploiting cool prop library. Properties are: density, thermal conductivity, volumetric thermal expansion coefficient, dynamic viscosity and Prandtl number.

        Args:
            film_temperature (np.array): average of the envirmonment medium and outer surface of the jacket temperatures.

        Returns:
            dict: dictionaty of environment properties.
        """
        return {
            prop_name: PropsSI(
                alias,
                "T",
                film_temperature,
                "P",
                self.inputs["Pressure"],
                self.type,
            )
            for prop_name, alias in self.fluid_prop_aliases.items()
        }

    # End method eval_air_prop.

    def grashof_number(self, dict_air_prop, T_s, characteristic_length):
        """[summary]

        Args:
            dict_air_prop ([type]): [description]
            T_s ([type]): [description]
            characteristic_length ([type]): [description]

        Returns:
            [type]: [description]
        """
        # Evaluate Grashof dimensionless number:
        # Gr_L = g*rho^2*betha*(T_s - T_inf)*L^3/mu^2
        return (
            constants.g
            * dict_air_prop["density"] ** 2
            * dict_air_prop["volumetric_thermal_expansion_coefficient"]
            * (T_s - self.inputs["Temperature"])
            * characteristic_length ** 3
            * np.reciprocal(dict_air_prop["dynamic_viscosity"] ** 2)
        )

    # End method grashof_number.

    def rayleigh_number(self, grashof, prandtl):
        """[summary]

        Args:
            grashof ([type]): [description]
            prandtl ([type]): [description]

        Returns:
            [type]: [description]
        """
        # Evaluate Rayleigh dimensionless number:
        # Ra_L = Gr_L*Pr
        return grashof * prandtl

    # End method rayleigh_number.

    def _vertical_plate(self, rayleigh, prandtl, grashof, conductor):
        """[summary]

        Args:
            rayleigh ([type]): [description]
            prandtl ([type]): [description]
            grashof ([type]): [description]
            conductor ([type]): [description]

        Returns:
            [type]: [description]
        """
        # Check if the correlation for the vertical plate can be applyed also to the case of vertical cylinder.
        self._check_validity_vertical_cylinder(conductor, grashof)

        nusselt = np.zeros(rayleigh.shape)
        nusselt[rayleigh >= 1e4 and rayleigh <= 1e9] = self._vertical_plate_laminar(
            rayleigh
        )
        nusselt[rayleigh >= 1e9 and rayleigh <= 1e13] = self._vertical_plate_turbulent(
            rayleigh
        )
        # Evaluate Nusselt dimensionless number.
        return nusselt

    # End method vertical_plate.

    def _vertical_plate_laminar(self, rayleigh):
        """[summary]

        Args:
            rayleigh ([type]): [description]

        Returns:
            [type]: [description]
        """
        # Incropera Foundations of heat transfer, sixth edition chapter 9.6.1 pag 573 eq. 9.24.
        return 0.59 * rayleigh ** (1.0 / 4.0)

    def _vertical_plate_turbulent(self, rayleigh):
        """[summary]

        Args:
            rayleigh ([type]): [description]

        Returns:
            [type]: [description]
        """
        # Incropera Foundations of heat transfer, sixth edition chapter 9.6.1 pag 573 eq. 9.24.
        return 0.10 * rayleigh ** (1.0 / 3.0)

    def _vertical_plate_churchill_chu(self, rayleigh, prandtl, grashof, conductor):
        """[summary]

        Args:
            rayleigh ([type]): [description]
            prandtl ([type]): [description]
            grashof ([type]): [description]
            conductor ([type]): [description]

        Returns:
            [type]: [description]
        """
        # Curchill and Chu correlation Incropera Foundations of heat transfer, sixth edition chapter 9.6.1 pag 573 eq. 9.26.

        # Check if the correlation for the vertical plate can be applyed also to the case of vertical cylinder.
        self._check_validity_vertical_cylinder(conductor, grashof)
        # Evaluate Nusselt dimensionless number.
        return (
            0.825
            + 0.387
            * rayleigh ** (1.0 / 6.0)
            / (1.0 + (0.492 / prandtl) ** (9.0 / 16.0)) ** (8.0 / 27.0)
        ) ** 2.0

    def _vertical_plate_laminar_accurate(self, rayleigh, prandtl, rayleigh_ub=1e9):
        """[summary]

        Args:
            rayleigh ([type]): [description]
            prandtl ([type]): [description]
            rayleigh_ub ([type], optional): [description]. Defaults to 1e9.

        Returns:
            [type]: [description]
        """
        # Incropera Foundations of heat transfer, sixth edition chapter 9.6.1 pag 573 eq. 9.27.
        self._check_validity_rayleigh(rayleigh, rayleigh_ub)
        return 0.68 + 0.670 * rayleigh ** (1.0 / 4.0) / (
            1.0 + (0.492 / prandtl) ** (9.0 / 16.0)
        ) ** (4.0 / 9.0)

    # End method vertical_plate_laminar_accurate.

    def _vertical_plate_churchill_chu_accurate(
        self, rayleigh, prandtl, grashof, conductor, rayleigh_ub=1e9
    ):
        """[summary]

        Args:
            rayleigh ([type]): [description]
            prandtl ([type]): [description]
            grashof ([type]): [description]
            conductor ([type]): [description]
            rayleigh_ub ([type], optional): [description]. Defaults to 1e9.

        Returns:
            [type]: [description]
        """
        # Check if the correlation for the vertical plate can be applyed also to the case of vertical cylinder.
        self._check_validity_vertical_cylinder(conductor, grashof)

        nusselt = np.zeros(rayleigh.shape)
        nusselt[rayleigh <= rayleigh_ub] = self._vertical_plate_laminar_accurate(
            rayleigh[rayleigh <= rayleigh_ub], prandtl[rayleigh <= rayleigh_ub]
        )
        nusselt[rayleigh > rayleigh_ub] = self._vertical_plate_churchill_chu(
            rayleigh[rayleigh > rayleigh_ub],
            prandtl[rayleigh > rayleigh_ub],
            grashof[rayleigh > rayleigh_ub],
            conductor,
        )
        return nusselt

    # End method vertical_plate_churchill_chu_accurate.

    def _check_validity_vertical_cylinder(self, conductor, grashof):
        """[summary]

        Args:
            grashof ([type]): [description]
            conductor ([type]): [description]
        """
        dict_check = {True: self._do_nothing, False: warnings.warn}
        check = any(
            conductor.inputs["Diameter"] / conductor.grid_features["delta_z"]
            > 35.0 / grashof ** (1.0 / 4.0)
        )
        dict_check[check](
            f"External free convection heat transfer coefficient may be inaccurate since the selected correlation for its evaluation {conductor.inputs['external_free_convection_correlation']} can not be applied to the case of a vertical cylinder!\n"
        )

    # End method check_validity_vertical_cylinder.

    def _check_validity_rayleigh(self, rayleigh, ub):
        """[summary]

        Args:
            rayleigh ([type]): [description]
            ub ([type]): [description]
        """
        dict_check = {True: self._do_nothing, False: warnings.warn}
        # Print warning message if Rayleigh dimensionless number > 1e12
        dict_check[rayleigh.max() <= ub](
            f"External free convection heat transfer coefficient may be inaccurate since Rayleigh dimensionless number > {ub}."
        )

    # End method check_validity_rayleigh.

    def _horiziontal_lower_surface_cold_plate(self, rayleigh, prandtl):
        """[summary]

        Args:
            rayleigh ([type]): [description]
            prandtl ([type]): [description]

        Returns:
            [type]: [description]
        """
        # Upper surface of hot plate or lower surface of cold plate (used in the latter way).
        cc = np.zeros(rayleigh.shape)
        nn = np.zeros(rayleigh.shape)
        cc[rayleigh >= 1e4 and rayleigh <= 1e7 and prandtl >= 0.7] = 0.54
        cc[rayleigh >= 1e7 and rayleigh <= 1e11] = 0.15

        nn[rayleigh >= 1e4 and rayleigh <= 1e7 and prandtl >= 0.7] = 1.0 / 4.0
        nn[rayleigh >= 1e7 and rayleigh <= 1e11] = 1.0 / 3.0

        return cc * rayleigh ** nn

    # End method _horiziontal_lower_surface_cold_plate.

    def _horiziontal_upper_surface_cold_plate(self, rayleigh, prandtl):
        """[summary]

        Args:
            rayleigh ([type]): [description]
            prandtl ([type]): [description]

        Returns:
            [type]: [description]
        """
        # Lower surface of hot plate or upper surface of cold plate (used in the latter way).
        # Validity 1e4 <= Ra <= 1e9 Pr >= 0.7
        return 0.52 * rayleigh ** (1.0 / 5.0)

    # End method _horiziontal_upper_surface_cold_plate.

    def _long_horziontal_cylinder_morgan(self, rayleigh, prandtl):
        """[summary]

        Args:
            rayleigh ([type]): [description]
            prandtl ([type]): [description]

        Returns:
            [type]: [description]
        """
        # Morgan correlation Incropera Foundations of heat transfer, sixth edition chapter 9.6.3 pag 581 eq. 9.33.
        cc = np.zeros(rayleigh.shape)
        nn = np.zeros(rayleigh.shape)

        cc[rayleigh >= 1e-10 and rayleigh <= 1e-2] = 0.675
        cc[rayleigh >= 1e-2 and rayleigh <= 1e2] = 1.02
        cc[rayleigh >= 1e2 and rayleigh <= 1e4] = 0.850
        cc[rayleigh >= 1e4 and rayleigh <= 1e7] = 0.480
        cc[rayleigh >= 1e7 and rayleigh <= 1e12] = 0.125

        nn[rayleigh >= 1e-10 and rayleigh <= 1e-2] = 0.058
        nn[rayleigh >= 1e-2 and rayleigh <= 1e2] = 0.148
        nn[rayleigh >= 1e2 and rayleigh <= 1e4] = 0.188
        nn[rayleigh >= 1e4 and rayleigh <= 1e7] = 0.250
        nn[rayleigh >= 1e7 and rayleigh <= 1e12] = 0.333

        # Evaluate nusselt dimensionless number by Morgan.
        return cc * rayleigh ** nn

    def _long_horziontal_cylinder_churchill_chu(
        self, rayleigh, prandtl, rayleigh_ub=1e12
    ):
        """[summary]

        Args:
            rayleigh ([type]): [description]
            prandtl ([type]): [description]
            rayleigh_ub ([type], optional): [description]. Defaults to 1e12.

        Returns:
            [type]: [description]
        """
        # Curchill and Chu correlation Incropera Foundations of heat transfer, sixth edition chapter 9.6.3 pag 581 eq. 9.34.
        self._check_validity_rayleigh(self, rayleigh, rayleigh_ub)
        # Evaluate nusselt dimensionless number by Churchill and Chu.
        return (
            0.60
            + 0.387
            * rayleigh ** (1.0 / 6.0)
            / (1.0 + (0.559 / prandtl) ** (9.0 / 16.0)) ** (8.0 / 27.0)
        ) ** 2.0

    # End method long_horziontal_cylinder_churchill_chu.

    def _do_nothing(self, *args, **kwargs):
        pass

    # End class Environment

    # -------------------------------------------------------------------------------
    # old code

    def eval_air_prop(self, film_temperature):
        """[summary]

        Args:
            film_temperature ([type]): [description]

        Returns:
            [type]: [description]
        """
        # Declare dictionary to decide if evaluate dry or humid properties according to flag self.inputs["Use_humidity"].
        dict_eval_air_properties = {
            True: self._eval_humid_air_properties,
            False: self._eval_dry_air_properties,
        }

        # Evaluate air volumetric thermal expansion coefficient, kinematic viscosity and Prandtl number at film temperature.
        return dict_eval_air_properties[self.inputs["Use_humidity"]](film_temperature)

    # End method eval_air_properties

    def _eval_dry_air_properties(self, temperature):
        """[summary]

        Args:
            temperature ([type]): [description]

        Returns:
            [type]: [description]
        """
        # Evaluate dry air thermal conductivity, volumetric thermal expansion coefficient, kinematic viscosity and Prandtl number.
        return eval_air_properties(
            temperature,
            self.inputs["Pressure"],
            thermal_conductivity=True,
            volumetric_thermal_expansion_coefficient=True,
            kinematic_viscosity=True,
            prandtl=True,
        )

    # End method dry_air_properties.

    def _eval_humid_air_properties(self, temperature):
        """[summary]

        Args:
            temperature ([type]): [description]

        Returns:
            [type]: [description]
        """
        # Evaluate humid air thermal conductivity, volumetric thermal expansion coefficient kinematic viscosity and Prandtl number.
        return eval_air_properties(
            temperature,
            self.inputs["Pressure"],
            self.ict_input["Relative_humidity"],
            thermal_conductivity=True,
            volumetric_thermal_expansion_coefficient=True,
            kinematic_viscosity=True,
            prandtl=True,
        )

    # End method eval_humid_air_properties.
