import numpy as np
import pandas as pd

from solid_component import SolidComponent
from strand_component import StrandComponent

# Aluminium properties
from properties_of_materials.aluminium import (
    thermal_conductivity_al,
    isobaric_specific_heat_al,
    density_al,
    electrical_resistivity_al,
)

# Cu properties
from properties_of_materials.copper import (
    thermal_conductivity_cu_nist,
    isobaric_specific_heat_cu_nist,
    density_cu,
    electrical_resistivity_cu_nist,
)

# NbTi properties
from properties_of_materials.niobium_titanium import (
    thermal_conductivity_nbti,
    isobaric_specific_heat_nbti,
    density_nbti,
)

# Nb3Sn properties
from properties_of_materials.niobium3_tin import (
    thermal_conductivity_nb3sn,
    isobaric_specific_heat_nb3sn,
    density_nb3sn,
)

DENSITY_FUNC = dict(
    al=density_al,
    cu=density_cu,
    nb3sn=density_nb3sn,
    nbti=density_nbti,
)

THERMAL_CONDUCTIVITY_FUNC = dict(
    al=thermal_conductivity_al,
    cu=thermal_conductivity_cu_nist,
    nb3sn=thermal_conductivity_nb3sn,
    nbti=thermal_conductivity_nbti,
)

ISOBARIC_SPECIFIC_HEAT_FUNC = dict(
    al=isobaric_specific_heat_al,
    cu=isobaric_specific_heat_cu_nist,
    nb3sn=isobaric_specific_heat_nb3sn,
    nbti=isobaric_specific_heat_nbti,
)

ELECTRICAL_RESISTIVITY_FUNC = dict(
    al=electrical_resistivity_al,
    cu=electrical_resistivity_cu_nist,
)


class StrandMixedComponent(StrandComponent):

    # Class for mixed strands objects

    ### INPUT PARAMETERS
    # some are inherited form the parent classes StrandComponent and SolidComponent

    ### THERMOPHYSICAL PROPERTIES
    # inherited from class SolidComponent

    ### OPERATIONAL PARAMETERS
    # inherited from parent classes StrandComponent and SolidComponent

    ### COMPUTED IN INITIALIZATION
    # inherited from class SolidComponent

    ### COMPUTED VECTOR FOR MAGNETIC FIELD
    # inherited from class SolidComponent

    KIND = "Mixed_sc_stab"

    def __init__(self, simulation, sheet, icomp, name, dict_file_path, conductor):

        self.NAME = name
        # get channels ID consistently with user definition (cdp, 09/2020)
        self.identifier = sheet.cell(row=3, column=4 + icomp).value

        # dictionary declaration (cdp, 11/2020)
        self.inputs = dict()
        self.operations = dict()
        self.dict_node_pt = dict()
        self.dict_Gauss_pt = dict()
        self.dict_num_step = dict()
        self.coordinate = dict()
        # Empty dictionary of list to save variable time evolutions at selected spatial coordinates.
        self.time_evol = dict(temperature=dict(), B_field=dict(), T_cur_sharing=dict())
        self.dict_scaling_input = dict()
        # Dictionary initialization: inputs.
        self.inputs = pd.read_excel(
            dict_file_path["input"],
            sheet_name=sheet.title,
            skiprows=2,
            header=0,
            index_col=0,
            usecols=["Variable name", self.identifier],
        )[self.identifier].to_dict()
        # Dictionary initialization: operations.
        self.operations = pd.read_excel(
            dict_file_path["operation"],
            sheet_name=sheet.title,
            skiprows=2,
            header=0,
            index_col=0,
            usecols=["Variable name", self.identifier],
        )[self.identifier].to_dict()

        # Superconductor cross section in m^2
        self.sc_cross_section = self.inputs["CROSSECTION"] / (
            1.0 + self.inputs["STAB_NON_STAB"]
        )
        # Stabilizer cross section in m^2
        self.stabilizer_cross_section = (
            self.inputs["CROSSECTION"] - self.sc_cross_section
        )
        # Call SolidComponent class constructor to deal with StrandMixedComponent time \
        # steps for current, external heating and so on (cdp, 11/2020)
        SolidComponent(simulation, self)
        if self.inputs["ISTABILIZER"] != "Cu":
            # remove key RRR from inputs if stabilizer is not Cu (cdp, 07/2020)
            self.inputs.pop("RRR")
        if self.operations["IBIFUN"] != -1:
            # Remove key B_field_units.
            del self.operations["B_field_units"]
        # end if (cdp, 07/2020)

        self.__strand_density_flag = False
        self.__reorganize_input()
        self.__check_consistecy(conductor)

    def __repr__(self):
        return f"{self.__class__.__name__}(Type: {self.NAME}, identifier: {self.identifier})"

    def __str__(self):
        pass

    def __reorganize_input(self):
        """Private method that reorganizes input data stored in dictionary self.inputs to simplify the procedure of properties homogenization."""

        # Create numpy array of string with the identifier of strand mixed material:
        # [stabilizer_material, superconductor_material]
        self.strand_material = np.array(
            [
                value.lower()
                for key, value in self.inputs.items()
                if key.endswith("material")
            ],
            dtype=str,
        )

        # Create numpy array of string with the identifier of strand mixed materials
        # that are not superconducting.
        self.strand_material_not_sc = np.array(
            [
                value.lower()
                for key, value in self.inputs.items()
                if key.endswith("material") and key != "Superconducting_material"
            ],
            dtype=str,
        )

        # Create numpy array of float with coefficient values used in
        # homogenization; order is consistent with values in
        # self.strand_material:
        # [stab_non_stab, 1.0]
        self.homogenization_cefficients = np.array(
            [self.inputs["STAB_NON_STAB"], 1.0],
            dtype=float,
        )

        # Total value of homogenization coefficients.
        self.homogenization_cefficients_sum = self.homogenization_cefficients.sum()

        # Create numpy array with density functions according to the strand mixed
        # material; order is consistent with values in self.strand_material.
        self.density_function = np.array(
            [DENSITY_FUNC[key] for key in self.strand_material]
        )

        # Create numpy array with electrical resistivity functions according to
        # the strand mixed material; order is consistent with values in
        # self.strand_material.
        self.electrical_resistivity_function_not_sc = np.array(
            [ELECTRICAL_RESISTIVITY_FUNC[key] for key in self.strand_material_not_sc]
        )

        # Create numpy array with isobaric specific heat functions according to
        # the strand mixed material; order is consistent with values in
        # self.strand_material.
        self.isobaric_specific_heat_function = np.array(
            [ISOBARIC_SPECIFIC_HEAT_FUNC[key] for key in self.strand_material]
        )

        # Create numpy array with thermal conductivity functions according to
        # the strand mixed material; order is consistent with values in
        # self.strand_material.
        self.thermal_conductivity_function = np.array(
            [THERMAL_CONDUCTIVITY_FUNC[key] for key in self.strand_material]
        )

    def __check_consistecy(self, conductor):
        """Private method that checks consistency of strand mixed user definition.

        Args:
            conductor (Conductor): instance of class Conductor.

        Raises:
            ValueError: if number of strand midex materials given in input is not consistent with user declared materials.

        """
        # Check that number of type materials given in input is consistent with
        # user declared materials.
        if self.strand_material.size != self.inputs["NUM_MATERIAL_TYPES"]:
            # Message to be completed!
            raise ValueError(
                f"{conductor.identifier = } -> {self.identifier = }\nThe number of material constituting the strand mixed ({self.inputs['NUM_MATERIAL_TYPES'] = }) is inconsistent with the number of defined materials ({self.tape_material.size = }).\nPlease check..."
            )

    def strand_density(self, property: dict) -> np.ndarray:
        """Method that evaluates the homogenized denstiy of the strand mixed, in the case it is made by two materials (stabilizer and superconductor). Homogenization is based on material cross sections.

        Args:
            property (dict): dictionary with material properties in nodal points or Gauss points according to the value of flag nodal in method eval_sol_comp_properties of class SolidComponent.

        Returns:
            np.ndarray: array with homogenized density of strand mixed in kg/m^3.
        """
        # Set fleag to true to allow evaluation of homogenized isobaric
        # specific heat.
        self.__strand_density_flag = True
        density = np.array(
            [func(property["temperature"].size) for func in self.density_function]
        )
        # Evaluate homogenized density of the strand mixed:
        # rho_eq = (rho_sc + stab_non_stab * rho_stab)/(1 + stab_non_stab)
        self.__density_numerator = np.array(
            list(map(np.multiply, density, self.homogenization_cefficients))
        )
        self.__density_numerator_sum = self.__density_numerator.sum(axis=0)
        return self.__density_numerator_sum / self.homogenization_cefficients_sum

    def strand_isobaric_specific_heat(self, property: dict) -> np.ndarray:
        """Method that evaluates homogenized isobaric specific heat of the sstrand mixed, in the case it is made by two materials (stabilizer and superconductor). Homogenization is based on material mass.

        Args:
            property (dict): dictionary with material properties in nodal points or Gauss points according to the value of flag nodal in method eval_sol_comp_properties of class SolidComponent.

        Returns:
            np.ndarray: array with homogenized isobaric specific heat of the strand mixed in J/kg/K.
        """
        # Check on homogenized density evaluation before homogenized isobaric
        # specific heat, since some therms are in common and are not evaluated
        # twices.
        if self.__strand_density_flag == False:
            raise ValueError(
                f"Call method {self.strand_density.__name__} before evaluation of homogenized strand mixed isobaric specific heat.\n"
            )

        # Set flag to false to trigger error in the next homogenized isobaric
        # specific heat evaluation if not done properly.
        self.__strand_density_flag = False
        isobaric_specific_heat = np.array(
            [
                func(property["temperature"])
                for func in self.isobaric_specific_heat_function
            ]
        )
        # Evaluate homogenized isobaric specific heat of the strand mixed:
        # cp_eq = (cp_sc*rho_sc + stab_non_stab*cp_stab*rho_stab)/(rho_sc + stab_non_stab * rho_stab)
        return (
            np.array(
                list(map(np.multiply, isobaric_specific_heat, self.__density_numerator))
            ).sum(axis=0)
            / self.__density_numerator_sum
        )

    def strand_thermal_conductivity(self, property: dict) -> np.ndarray:
        """Method that evaluates the homogenized thermal conductivity of the strand mixed, in the case it is made by two materials (stabilizer and superconductor). Homogenization is based on material cross sections.

        Args:
            property (dict): dictionary with material properties in nodal points or Gauss points according to the value of flag nodal in method eval_sol_comp_properties of class SolidComponent.

        Returns:
            np.ndarray: array with homogenized thermal conductivity of the strand mixed in W/m/K.
        """
        thermal_conductivity = np.zeros(
            (self.inputs["NUM_MATERIAL_TYPES"], property["temperature"].size)
        )
        for ii, func in enumerate(self.thermal_conductivity_function):
            if "cu" in func.__name__:
                thermal_conductivity[ii, :] = func(
                    property["temperature"],
                    property["B_field"],
                    self.inputs["RRR"],
                )
            else:
                thermal_conductivity[ii, :] = func(property["temperature"])
        # Evaluate homogenized thermal conductivity of the strand mixed:
        # k_eq = (K_sc + stab_non_stab * K_stab)/(1 + stab_non_stab)
        return (
            np.array(
                list(
                    map(
                        np.multiply,
                        thermal_conductivity,
                        self.homogenization_cefficients,
                    )
                )
            ).sum(axis=0)
            / self.homogenization_cefficients_sum
        )

    def strand_electrical_resistivity_not_sc(self, property: dict) -> np.ndarray:
        """Method that evaluates electrical resisitivity of the stabilizer material of the strand mixed.

        Args:
            property (dict): dictionary with material properties in nodal points or Gauss points according to the value of flag nodal in method eval_sol_comp_properties of class SolidComponent.

        Returns:
            np.ndarray: array with electrical resisitivity of the stabilizer of the strand mixed in Ohm*m.
        """
        if self.strand_material_not_sc == "cu":
            return ELECTRICAL_RESISTIVITY_FUNC[self.strand_material_not_sc](
                property["temperature"],
                property["B_field"],
                self.inputs["RRR"],
            )
        else:
            return ELECTRICAL_RESISTIVITY_FUNC[self.strand_material_not_sc](
                property["temperature"]
            )
