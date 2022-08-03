from pickle import DICT
from solid_component import SolidComponent
import pandas as pd
import numpy as np


# Stainless steel properties
from properties_of_materials.stainless_steel import (
    thermal_conductivity_ss,
    isobaric_specific_heat_ss,
    density_ss,
    electrical_resistivity_ss,
)

# Glass-epoxy properties
from properties_of_materials.glass_epoxy import (
    thermal_conductivity_ge,
    isobaric_specific_heat_ge,
    density_ge,
    electrical_resistivity_ge,
)

DENSITY_FUNC = dict(ge=density_ge, ss=density_ss)
ISOBARIC_SPECIFIC_HEAT_FUNC = dict(
    ge=isobaric_specific_heat_ge, ss=isobaric_specific_heat_ss
)
THERMAL_CONDUCTIVITY_FUNC = dict(ge=thermal_conductivity_ge, ss=thermal_conductivity_ss)
ELECTRICAL_RESISTIVITY_FUNC = dict(
    ge=electrical_resistivity_ge, ss=electrical_resistivity_ss
)


class JacketComponent(SolidComponent):

    # Class for jacket objects

    ### INPUT PARAMETERS
    # some are inherited from class SolidComponent

    ### THERMOPHYSICAL PROPERTIES
    # inherited from class SolidComponent

    ### OPERATIONAL PARAMETERS
    # inherited from class SolidComponent

    ### COMPUTED IN INITIALIZATION
    # inherited from class SolidComponent

    ### COMPUTED VECTOR FOR MAGNETIC FIELD
    # inherited from class SolidComponent

    KIND = "JacketComponent"

    def __init__(self, simulation, sheet, icomp, name, dict_file_path, conductor):

        self.name = name
        # get channels ID consistently with user definition (cdp, 09/2020)
        self.identifier = sheet.cell(row=3, column=4 + icomp).value

        # dictionary declaration (cdp, 11/2020)
        self.inputs = dict()
        self.operations = dict()
        self.dict_node_pt = dict()
        self.dict_Gauss_pt = dict()
        self.dict_num_step = dict()
        self.radiative_heat_env = ""
        self.radiative_heat_inn = dict()
        self.coordinate = dict()
        # Empty dictionary of list to save variable time evolutions at selected spatial coordinates.
        self.time_evol = dict(temperature=dict())
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

        if self.operations["IBIFUN"] != -1:
            # Remove key B_field_units.
            del self.operations["B_field_units"]

        self.__reorganize_input()
        self.__check_consistency(conductor)

        # Flag to check if evaluation of homogenized isobaric specific heat can
        # be done or not (depends on homogenized density evaluation).
        self.__jacket_density_flag = False

        # Call SolidComponent class constructor to deal with JacketComponent time \
        # steps for current, external heating and so on (cdp, 11/2020)
        SolidComponent(simulation, self)

    def __repr__(self):
        return f"{self.__class__.__name__}(Type: {self.name}, identifier: {self.identifier})"

    def __str__(self):
        pass

    def _radiative_source_therm_env(self, conductor, environment):
        """Method that evaluates the heat transferred by radiation with the external environment.

        Args:
            conductor ([type]): [description]
            environment ([type]): [description]
        """
        key = f"{environment.KIND}_{self.identifier}"
        if conductor.inputs["METHOD"] == "BE" or conductor.inputs["METHOD"] == "CN":
            # Backward Euler or Crank-Nicolson.
            if conductor.cond_time[-1] == 0:
                # Initialization.
                self.radiative_heat_env = np.zeros(
                    (conductor.grid_features["N_nod"], 2)
                )
            elif conductor.cond_time[-1] > 0:
                if conductor.cond_num_step == 1:
                    # Store the old values only immediately after the initializzation, \
                    # since after that the whole SYSLOD array is saved and there is no \
                    # need to compute twice the same values.
                    self.radiative_heat_env[:, 1] = self.radiative_heat_env[:, 0].copy()
                # Update value at the current time step.
                self.radiative_heat_env[:, 0] = (
                    conductor.dict_interf_peri["env_sol"][key]
                    * conductor.dict_node_pt["HTC"]["env_sol"][key]["rad"]
                    * (
                        environment.inputs["Temperature"]
                        - self.dict_node_pt["temperature"]
                    )
                )
            # end if conductor.cond_time[-1].
        elif conductor.inputs["METHOD"] == "AM4":
            # Adams-Moulton 4.
            if conductor.cond_time[-1] == 0:
                # Initialization.
                self.radiative_heat_env = np.zeros(
                    (conductor.grid_features["N_nod"], 4)
                )
            elif conductor.cond_time[-1] > 0:
                self.radiative_heat_env[:, 1:4] = self.radiative_heat_env[:, 0:3].copy()
                # Update value at the current time step.
                self.radiative_heat_env[:, 0] = (
                    conductor.dict_interf_peri["env_sol"][key]
                    * conductor.dict_node_pt["HTC"]["env_sol"][key]["rad"]
                    * (
                        environment.inputs["Temperature"]
                        - self.dict_node_pt["temperature"]
                    )
                )
            # end if conductor.cond_time[-1].
        # end if conductor.inputs["METHOD"].

    # End method _radiative_source_therm.

    def _radiative_heat_exc_inner(self, conductor, jk_inner):
        """Method that evaluates the heat transferred by radiation with the inner surface of the enclosure and the inner jackets.

        Args:
            conductor ([type]): [description]
            environment ([type]): [description]
        """
        if self.identifier < jk_inner.identifier:
            key = f"{self.identifier}_{jk_inner.identifier}"
        else:
            key = f"{jk_inner.identifier}_{self.identifier}"
        # End if self.identifier.
        if conductor.inputs["METHOD"] == "BE" or conductor.inputs["METHOD"] == "CN":
            # Backward Euler or Crank-Nicolson.
            if conductor.cond_time[-1] == 0:
                # Initialization.
                self.radiative_heat_inn[key] = np.zeros(
                    (conductor.grid_features["N_nod"], 2)
                )
            elif conductor.cond_time[-1] > 0:
                if conductor.cond_num_step == 1:
                    # Store the old values only immediately after the initializzation, \
                    # since after that the whole SYSLOD array is saved and there is no \
                    # need to compute twice the same values.
                    self.radiative_heat_inn[key][:, 1] = self.radiative_heat_inn[key][
                        :, 0
                    ].copy()
                # Update value at the current time step.
                self.radiative_heat_inn[key][:, 0] = (
                    conductor.dict_interf_peri["sol_sol"][key]
                    * conductor.dict_node_pt["HTC"]["sol_sol"][key]["rad"]
                    * (
                        jk_inner.dict_node_pt["temperature"]
                        - self.dict_node_pt["temperature"]
                    )
                )
            # end if conductor.cond_time[-1].
        elif conductor.inputs["METHOD"] == "AM4":
            # Adams-Moulton 4.
            if conductor.cond_time[-1] == 0:
                # Initialization.
                self.radiative_heat_inn[key] = np.zeros(
                    (conductor.grid_features["N_nod"], 4)
                )
            elif conductor.cond_time[-1] > 0:
                self.radiative_heat_inn[key][:, 1:4] = self.radiative_heat_inn[key][
                    :, 0:3
                ].copy()
                # Update value at the current time step.
                self.radiative_heat_inn[key][:, 0] = (
                    conductor.dict_interf_peri["sol_sol"][key]
                    * conductor.dict_node_pt["HTC"]["sol_sol"][key]["rad"]
                    * (
                        jk_inner.dict_node_pt["temperature"]
                        - self.dict_node_pt["temperature"]
                    )
                )
            # end if conductor.cond_time[-1].
        # end if conductor.inputs["METHOD"].

    # End method _radiative_source_therm.

    def __reorganize_input(self):
        """Private method that reorganizes input data stored in dictionary self.inputs to simplify the procedure of properties homogenization."""

        # Create numpy array of string with the identifier of jacket material:
        # [jacket_material, insulation_material]
        self.materials = np.array(
            [
                self.inputs["jacket_material"].lower(),
                self.inputs["insulation_material"].lower(),
            ],
            dtype=str,
        )

        # Get the indexes corresponding to "none" used for consistency check.
        self.__index_material_none = np.nonzero(self.materials == "none")[0]
        # Remove "none" items, used to identify not used layer
        self.materials = self.materials[np.nonzero(self.materials != "none")[0]]

        # Create numpy array of float with coefficient values used in
        # homogenization; order is consistent with values in
        # self.materials:
        # [jacket_cross_section, insulation_cross_section]
        self.cross_sections = np.array(
            [
                self.inputs["jacket_cross_section"],
                self.inputs["insulation_cross_section"],
            ],
            dtype=float,
        )

        # Get the indexes corresponding to 0 used for consistency check.
        self.__index_cross_section_0 = np.nonzero(self.cross_sections == 0)[0]
        self.cross_sections = self.cross_sections[np.nonzero(self.cross_sections)[0]]

        # Total value of homogenization coefficients.
        self.__cross_section = self.cross_sections.sum()

        # Create numpy array with density functions according to jacket
        # materials; order is consistent with values in self.materials.
        self.density_function = np.array([DENSITY_FUNC[key] for key in self.materials])

        # Create numpy array with electrical resistivity functions according to
        # jacket materials; order is consistent with values in
        # self.materials.
        self.electrical_resistivity_function = np.array(
            [ELECTRICAL_RESISTIVITY_FUNC[key] for key in self.materials]
        )

        # Create numpy array with isobaric specific heat functions according to
        # jacket materials; order is consistent with values in
        # self.materials.
        self.isobaric_specific_heat_function = np.array(
            [ISOBARIC_SPECIFIC_HEAT_FUNC[key] for key in self.materials]
        )

        # Create numpy array with thermal conductivity functions according to
        # jacket material; order is consistent with values in
        # self.materials.
        self.thermal_conductivity_function = np.array(
            [THERMAL_CONDUCTIVITY_FUNC[key] for key in self.materials]
        )

    def __check_consistency(self, conductor):
        """Private method that checks consistency of jacket user definition.

        Args:
            conductor (Conductor): instance of class Conductor.

        Raises:
            ValueError: if number of jacket materials given in input is not consistent with user declared materials.
            ValueError: if number of jacket materials given in input is not consistent with not zero user defined material thicknes.
            ValueError: if the indexes of "none" material are not equal to the indexes of thickness equal to 0.
            ValueError: if jacket cross section given in input is not consistent with the evaluated one.
        """
        # Check that number of jacket materials given in input is consistent
        # with user declared materials.
        if self.materials.size != self.inputs["NUM_MATERIAL_TYPES"]:
            # Message to be completed!
            raise ValueError(
                f"{conductor.identifier = } -> {self.identifier = }\nThe number of material constituting the jacket ({self.inputs['NUM_MATERIAL_TYPES'] = }) is inconsistent with the number of defined materials ({self.materials.size = }).\nPlease check..."
            )

        # Check that number of jacket materials given in input is consistent
        # with not zero user defined material cross sections.
        if self.cross_sections.size != self.inputs["NUM_MATERIAL_TYPES"]:
            # Message to be completed!
            raise ValueError(
                f"{conductor.identifier = } -> {self.identifier = }\nThe number of material constituting the jacket ({self.inputs['NUM_MATERIAL_TYPES'] = }) is inconsistent with the number of defined cross sections ({self.cross_sections.size = }).\nPlease check..."
            )

        # Check that the indexes of "none" material are equal to the indexes of
        # cross sections equal to 0.
        if any(self.__index_material_none != self.__index_cross_section_0):
            # Message to be completed!
            raise ValueError(
                f"{conductor.identifier = } -> {self.identifier = }\nDefined materials and defined cross sections must be consistent.\nPlease check..."
            )

        # Check that jacket cross section given in input is consistent with the
        # evaluated one.
        tol = 1e-3
        self.inputs["CROSSECTION"] = self.inputs["jacket_cross_section"] + self.inputs["insulation_cross_section"]
        if (
            abs(self.__cross_section - self.inputs["CROSSECTION"])
            / self.inputs["CROSSECTION"]
            > tol
        ):
            # Message to be completed!
            raise ValueError(
                f"{conductor.identifier = } -> {self.identifier = }\nInconsistent cross section value: user defines {self.inputs['CROSSECTION'] = } while computed one is {self.__cross_section = }.\nPlease check..."
            )

        # Delete no longer useful attributes.
        del (
            self.__index_material_none,
            self.__index_cross_section_0,
            self.__cross_section,
        )

    def jacket_density(self, property: dict) -> np.ndarray:
        """Method that evaluates the homogenized denstiy of the jacket, in the case it is made by at most by two materials (jacket and insulation). Homogenization is based on material cross sections.

        Args:
            property (dict): dictionary with material properties in nodal points or Gauss points according to the value of flag nodal in method eval_sol_comp_properties of class SolidComponent.

        Returns:
            np.ndarray: array with homogenized density of jacket in kg/m^3.
        """
        # Set fleag to true to allow evaluation of homogenized isobaric
        # specific heat.
        self.__jacket_density_flag = True
        density = np.array(
            [func(property["temperature"].size) for func in self.density_function]
        )
        density = density.T
        # Evaluate homogenized density of the strand mixed:
        # rho_eq = (A_jk*rho_jk + A_in*rho_in)/(A_jk + A_in)
        self.__density_numerator = density * self.cross_sections
        self.__density_numerator_sum = self.__density_numerator.sum(axis=1)
        return self.__density_numerator_sum / self.inputs["CROSSECTION"]

    def jacket_isobaric_specific_heat(self, property: dict) -> np.ndarray:
        """Method that evaluates homogenized isobaric specific heat of the jacket, in the case it is made at most by two materials (jacket and insulation). Homogenization is based on material mass.

        Args:
            property (dict): dictionary with material properties in nodal points or Gauss points according to the value of flag nodal in method eval_sol_comp_properties of class SolidComponent.

        Returns:
            np.ndarray: array with homogenized isobaric specific heat of the jacket in J/kg/K.
        """
        # Check on homogenized density evaluation before homogenized isobaric
        # specific heat, since some therms are in common and are not evaluated
        # twices.
        if self.__jacket_density_flag == False:
            raise ValueError(
                f"Call method {self.jacket_density.__name__} before evaluation of homogenized jacket isobaric specific heat.\n"
            )

        # Set flag to false to trigger error in the next homogenized isobaric
        # specific heat evaluation if not done properly.
        self.__jacket_density_flag = False
        isobaric_specific_heat = np.array(
            [
                func(property["temperature"])
                for func in self.isobaric_specific_heat_function
            ]
        )
        isobaric_specific_heat = isobaric_specific_heat.T
        # Evaluate homogenized isobaric specific heat of the strand mixed:
        # cp_eq = (cp_jk*A_jk*rho_jk + cp_in*A_in*rho_in)/(A_jk*rho_jk +
        # A_in*rho_in)
        return (isobaric_specific_heat * self.__density_numerator).sum(axis=1)/self.__density_numerator_sum

    def jacket_thermal_conductivity(self, property: dict) -> np.ndarray:
        """Method that evaluates the homogenized thermal conductivity of the jacket, in the case it is made by at most by two materials (jacket and insulation). Homogenization is based on material cross sections.

        Args:
            property (dict): dictionary with material properties in nodal points or Gauss points according to the value of flag nodal in method eval_sol_comp_properties of class SolidComponent.

        Returns:
            np.ndarray: array with homogenized thermal conductivity of the jacket in W/m/K.
        """
        thermal_conductivity = np.array(
            [
                func(property["temperature"].size)
                for func in self.thermal_conductivity_function
            ]
        )
        thermal_conductivity = thermal_conductivity.T
        # Evaluate homogenized thermal conductivity of the strand mixed:
        # k_eq = (A_jk*k_jk + A_in*k_in)/(A_jk + A_in)
        return (thermal_conductivity * self.cross_sections).sum(axis=1)/ self.inputs["CROSSECTION"]

    def jaket_electrical_resistivity(self, property: dict) -> np.ndarray:
        """Method that evaluates the homogenized electrical resistivity otf the jacket, in the case it is made by at most by two materials (jacket and insulation). Homogenization is based on material cross sections.

        Args:
            property (dict): dictionary with material properties in nodal points or Gauss points according to the value of flag nodal in method eval_sol_comp_properties of class SolidComponent.

        Returns:
            np.ndarray: array with homogenized electrical resistivity of the jacket in Ohm*m.
        """

        electrical_resistivity = np.array(
            [
                func(property["temperature"].size)
                for func in self.electrical_resistivity_function
            ]
        )
        electrical_resistivity = electrical_resistivity.T
        # Evaluate homogenized electrical resistivity of the jacket:
        # rho_el_eq = (A_jk + A_in) * ((A_jk/rho_el_jk + A_in/rho_el_in))^-1
        return self.inputs["CROSSECTION"]*np.reciprocal((self.cross_sections / electrical_resistivity).sum(axis=1))