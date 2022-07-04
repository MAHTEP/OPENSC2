# Import libraries
import numpy as np
import pandas as pd

# Import classes
from solid_component import SolidComponent
from strand_component import StrandComponent

# Cu properties
from properties_of_materials.copper import (
    thermal_conductivity_cu_nist,
    isobaric_specific_heat_cu_nist,
    density_cu,
    electrical_resistivity_cu_nist,
)

# RE123 properties
from properties_of_materials.rare_earth_123 import (
    thermal_conductivity_re123,
    isobaric_specific_heat_re123,
    density_re123,
)

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
)

# Silver properties
from properties_of_materials.silver import (
    thermal_conductivity_ag,
    isobaric_specific_heat_ag,
    density_ag,
    electrical_resistivity_ag,
)

# HASTELLOY - C276 properties
from properties_of_materials.hastelloy_c276 import (
    thermal_conductivity_hc276,
    isobaric_specific_heat_hc276,
    density_hc276,
    electrical_resistivity_hc276,
)

# Solder Sn60Pb40 properties
from properties_of_materials.solder_sn60_pb40 import (
    thermal_conductivity_sn60pb40,
    isobaric_specific_heat_sn60pb40,
    density_sn60pb40,
    electrical_resistivity_sn60pb40,
)

# Aluminium properties
from properties_of_materials.aluminium import (
    thermal_conductivity_al,
    isobaric_specific_heat_al,
    density_al,
    electrical_resistivity_al,
)

DENSITY_FUNC = dict(
    ag=density_ag,
    al=density_al,
    cu=density_cu,
    ge=density_ge,
    hc276=density_hc276,
    re123=density_re123,
    sn60pb40=density_sn60pb40,
    ss=density_ss,
)

THERMAL_CONDUCTIVITY_FUNC = dict(
    ag=thermal_conductivity_ag,
    al=thermal_conductivity_al,
    cu=thermal_conductivity_cu_nist,
    ge=thermal_conductivity_ge,
    hc276=thermal_conductivity_hc276,
    re123=thermal_conductivity_re123,
    sn60pb40=thermal_conductivity_sn60pb40,
    ss=thermal_conductivity_ss,
)

ISOBARIC_SPECIFIC_HEAT_FUNC = dict(
    ag=isobaric_specific_heat_ag,
    al=isobaric_specific_heat_al,
    cu=isobaric_specific_heat_cu_nist,
    ge=isobaric_specific_heat_ge,
    hc276=isobaric_specific_heat_hc276,
    re123=isobaric_specific_heat_re123,
    sn60pb40=isobaric_specific_heat_sn60pb40,
    ss=isobaric_specific_heat_ss,
)

ELECTRICAL_RESISTIVITY_FUNC = dict(
    ag=electrical_resistivity_ag,
    al=electrical_resistivity_al,
    cu=electrical_resistivity_cu_nist,
    # ge=electrical_resistivity_ge, not defined
    hc276=electrical_resistivity_hc276,
    sn60pb40=electrical_resistivity_sn60pb40,
    ss=electrical_resistivity_ss,
)


class StackComponent(StrandComponent):
    """Class that defines StackComponents objects to model HTS stacks of tapes.

    Args:
        StrandComponent (StrandComponent): class that defines general methods for strand and stack objects.
    """

    KIND = "Stack"

    def __init__(
        self,
        simulation,
        sheet,
        icomp: int,
        name: str,
        dict_file_path: dict,
        conductor: object,
    ):
        """Method that makes instance of class StackComponent.

        Args:
            simulation (Simulation): simulation object.
            sheet (Worksheet): worksheet with input data.
            icomp (int): component index.
            name (str): component name.
            dict_file_path (dict): dictionary with paths to load the input files.
            conductor(object): inscance of class Conductor
        """

        self.NAME = name
        # get channels ID consistently with user definition (cdp, 09/2020)
        self.identifier = sheet.cell(row=3, column=4 + icomp).value

        # dictionary declaration
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

        # Call SolidComponent class constructor to deal with StrandMixedComponent time \
        # steps for current, external heating and so on
        SolidComponent(simulation, self)
        if self.operations["IBIFUN"] != -1:
            # Remove key B_field_units.
            del self.operations["B_field_units"]
        # end if

        self.__reorganize_input()
        self.__check_consistecy(conductor)

        # Flag to check if evaluation of homogenized isobaric specific heat can
        # be done or not (depends on homogenized density evaluation).
        self.__stack_density_flag = False

    def __repr__(self):
        return f"{self.__class__.__name__}(Type: {self.NAME}, identifier: {self.identifier})"

    def __str__(self):
        pass

    def __reorganize_input(self):
        """Private method that reorganizes input data stored in dictionary self.inputs to simplify the procedure of properties homogenization."""

        # Create numpy array of string with the identifier of tape material
        self.tape_material = np.array(
            [value.lower() for key, value in self.inputs.items() if key.endswith("material")],
            dtype=str,
        )
        # Get the indexes corresponding to "none" used for consistency check.
        self.__index_material_none = np.nonzero(self.tape_material == "none")[0]
        # Remove "none" items, used to identify not used layer
        self.tape_material = self.tape_material[
            np.nonzero(self.tape_material != "none")[0]
        ]
        # Create numpy array of string with the identifier of tape materials
        # that are not superconducting.
        self.tape_material_not_sc = np.array(
            [
                value.lower()
                for key, value in self.inputs.items()
                if key.endswith("material")
                and key != "HTS_material"
                and value != "none"
            ],
            dtype=str,
        )

        # Create numpy array of float with values of thickness of tape
        # layers in m; order is consistent with values in self.tape_material.
        self.material_thickness = np.array(
            [value for key, value in self.inputs.items() if key.endswith("thickness")],
            dtype=float,
        )
        # Get the indexes corresponding to 0 used for consistency check.
        self.__index_thickness_0 = np.nonzero(self.material_thickness == 0)[0]
        self.material_thickness = self.material_thickness[
            np.nonzero(self.material_thickness)[0]
        ]

        # Total tape thickness in m.
        self.tape_thickness = self.material_thickness.sum()

        # Create numpy array of float with values of thickness of not
        # superconducting tape layers in m; order is consistent with values in self.tape_material_not_sc.
        self.material_thickness_not_sc = np.array(
            [
                value
                for key, value in self.inputs.items()
                if key.endswith("thickness") and key != "HTS_thickness" and value > 0.0
            ],
            dtype=float,
        )
        # Total not superconducting tape thickness in m.
        self.tape_thickness_not_sc = self.material_thickness_not_sc.sum()

        # Create numpy array with density functions according to the tape
        # material; order is consistent with values in self.tape_material.
        self.density_function = np.array(
            [DENSITY_FUNC[key] for key in self.tape_material]
        )

        # Create numpy array with electrical resistivity functions according to
        # the tape material; order is consistent with values in
        # self.tape_material.
        self.electrical_resistivity_function_not_sc = np.array(
            [ELECTRICAL_RESISTIVITY_FUNC[key] for key in self.tape_material_not_sc]
        )

        # Create numpy array with isobaric specific heat functions according to
        # the tape material; order is consistent with values in
        # self.tape_material.
        self.isobaric_specific_heat_function = np.array(
            [ISOBARIC_SPECIFIC_HEAT_FUNC[key] for key in self.tape_material]
        )

        # Create numpy array with thermal conductivity functions according to
        # the tape material; order is consistent with values in
        # self.tape_material.
        self.thermal_conductivity_function = np.array(
            [THERMAL_CONDUCTIVITY_FUNC[key] for key in self.tape_material]
        )

    def __check_consistecy(self, conductor):
        """Private method that checks consistency of stack and or tape user definition.

        Args:
            conductor (Conductor): instance of class Conductor.

        Raises:
            ValueError: if number of tape materials given in input is not consistent with user declared materials.
            ValueError: if number of tape materials given in input is not consistent with not zero user defined material thicknes.
            ValueError: if the indexes of "none" material are not equal to the indexes of thickness equal to 0.
            ValueError: if number of tapes given in input is not consistent with the evaluated one.
            ValueError: if stack cross section given in input is not consistent with the evaluated one.
        """
        # Check that number of tape materials given in input is consistent with
        # user declared materials.
        if self.tape_material.size != self.inputs["N_materal_tape"]:
            # Message to be completed!
            raise ValueError(
                f"{conductor.identifier = } -> {self.identifier = }\nThe number of material constituting the tape ({self.inputs['N_materal_tape'] = }) is inconsistent with the number of defined materials ({self.tape_material.size = }).\nPlease check..."
            )

        # Check that number of tape materials given in input is consistent with
        # not zero user defined material thicknes.
        if self.material_thickness.size != self.inputs["N_materal_tape"]:
            # Message to be completed!
            raise ValueError(
                f"{conductor.identifier = } -> {self.identifier = }\nThe number of material constituting the tape ({self.inputs['N_materal_tape'] = }) is inconsistent with the number of defined thicknesses ({self.material_thickness.size = }).\nPlease check..."
            )

        # Check that the indexes of "none" material are equal to the indexes of
        # thickness equal to 0.
        if any(self.__index_material_none != self.__index_thickness_0):
            # Message to be completed!
            raise ValueError(
                f"{conductor.identifier = } -> {self.identifier = }\nDefined materials and defined thicknesses must be consistent.\nPlease check..."
            )

        # Following quantities are useful to cross check the number of tapes
        # and the stack cross section given in input by the user.

        # Tape cross section in # m^2
        self.tape_cross_section = self.tape_thickness * self.inputs["Stack_width"]
        # Evaluate the number of tapes constituing the stack, used for consistency check.
        self.__tape_number = round(
            self.inputs["Cross_section"] / self.tape_cross_section
        )
        # Evaluate the total cross section of the stack of tapes, used for consistency check.
        self.__cross_section = self.tape_cross_section * self.__tape_number

        # Check that number of tapes given in input is consistent with the
        # evaluated one.
        if self.__tape_number != self.inputs["N_tape"]:
            # Message to be completed!
            raise ValueError(
                f"{conductor.identifier = } -> {self.identifier = }\nInconsistent number of tape: user defines {self.inputs['N_tape'] = } while computed one is {self.__tape_number = }.\nPlease check..."
            )

        # Check that stack cross section given in input is consistent with the
        # evaluated one.
        tol = 1e-3
        if (
            abs(self.__cross_section - self.inputs["Cross_section"])
            / self.inputs["Cross_section"]
            > tol
        ):
            # Message to be completed!
            raise ValueError(
                f"{conductor.identifier = } -> {self.identifier = }\nInconsistent number of tape: user defines {self.inputs['N_tape'] = } while computed one is {self.__tape_number = }.\nPlease check..."
            )

        # Delete no longer useful attributes.
        del (
            self.__index_material_none,
            self.__index_thickness_0,
            self.__tape_number,
            self.__cross_section,
        )

    def stack_density(self, property: dict) -> np.ndarray:
        """Method that evaluates the homogenized denstiy of the stack, which is the same of the tape if the tapes constituting the stack are equals to each other. Homogenization is based on the thickness of tape layers.

        Args:
            property (dict): dictionary with material properties in nodal points or Gauss points according to the value of flag nodal in method eval_sol_comp_properties of class SolidComponent.

        Returns:
            np.ndarray: array with homogenized density of the stack of tapes in kg/m^3.
        """
        # Set fleag to true to allow evaluation of homogenized isobaric
        # specific heat.
        self.__stack_density_flag = True
        density = np.array(
            [func(property["temperature"].size) for func in self.density_function]
        )
        # Evaluate homogenized density of the stack:
        # rho_eq = sum(s_i*rho_i)/s
        self.__density_numerator = np.array(
            list(map(np.multiply, density, self.material_thickness))
        )
        self.__density_numerator_sum = self.__density_numerator.sum(axis=0)
        return self.__density_numerator_sum / self.tape_thickness

    def stack_isobaric_specific_heat(self, property: dict) -> np.ndarray:
        """Method that evaluates homogenized isobaric specific heat of the stack, which is the same of the tape if the tapes constituting the stack are equals to each other. Homogenization is based on the mass of tape layers.

        Args:
            property (dict): dictionary with material properties in nodal points or Gauss points according to the value of flag nodal in method eval_sol_comp_properties of class SolidComponent.

        Returns:
            np.ndarray: array with homogenized isobaric specific heat of the stack of tapes in J/kg/K.
        """
        # Check on homogenized density evaluation before homogenized isobaric
        # specific heat, since some therms are in common and are not evaluated
        # twices.
        if self.__stack_density_flag == False:
            raise ValueError(
                f"Cal method {self.stack_density.__name__} before evaluation of homogenized stack isobaric specific heat.\n"
            )

        # Set flag to false to trigger error in the next homogenized isobaric
        # specific heat evaluation if not done properly.
        self.__stack_density_flag = False
        isobaric_specific_heat = np.array(
            [
                func(property["temperature"])
                for func in self.isobaric_specific_heat_function
            ]
        )
        # Evaluate homogenized isobaric specific heat of the stack:
        # cp_eq = sum(s_i*rho_i*cp_i)/sum(s_i*rho_i)
        return (
            np.array(
                list(map(np.multiply, isobaric_specific_heat, self.__density_numerator))
            ).sum(axis=0)
            / self.__density_numerator_sum
        )

    def stack_thermal_conductivity(self, property: dict) -> np.ndarray:
        """Method that evaluates the homogenized thermal conductivity of the stack, which is the same of the tape if the tapes constituting the stack are equals to each other. Homogenization is based on the thickness of tape layers.

        Args:
            property (dict): dictionary with material properties in nodal points or Gauss points according to the value of flag nodal in method eval_sol_comp_properties of class SolidComponent.

        Returns:
            np.ndarray: array with homogenized thermal conductivity of the stack of tapes in W/m/K.
        """
        thermal_conductivity = np.zeros(
            (self.inputs["N_materal_tape"], property["temperature"].size)
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
        # Evaluate homogenized thermal conductivity of the stack:
        # k_eq = sum(s_i*k_i)/s
        return (
            np.array(
                list(map(np.multiply, thermal_conductivity, self.material_thickness))
            ).sum(axis=0)
            / self.tape_thickness
        )

    def stack_electrical_resistivity_not_sc(self, property: dict) -> np.ndarray:
        """Method that evaluates the homogenized electrical resistivity for not superconducting materials of the stack, which is the same of the tape if the tapes constituting the stack are equals to each other. Homogenization is based on the thickness of tape layers.

        Args:
            property (dict): dictionary with material properties in nodal points or Gauss points according to the value of flag nodal in method eval_sol_comp_properties of class SolidComponent.

        Returns:
            np.ndarray: array with homogenized electrical resistivity of not superconducting materials of the stack of tapes in Ohm*m.
        """
        elestrical_resistivity = np.zeros(
            (self.inputs["N_materal_tape"] - 1, property["temperature"].size)
        )
        for ii, func in enumerate(self.electrical_resistivity_function_not_sc):
            if "cu" in func.__name__:
                elestrical_resistivity[ii, :] = func(
                    property["temperature"],
                    property["B_field"],
                    self.inputs["RRR"],
                )
            else:
                elestrical_resistivity[ii, :] = func(property["temperature"])
        # Evaluate homogenized electrical resistivity of the stack:
        # rho_el_eq = s_not_sc * (sum(s_i/rho_el_i))^-1 for any i not sc
        return (
            np.reciprocal(
                np.array(
                    list(
                        map(
                            np.divide,
                            self.material_thickness_not_sc,
                            elestrical_resistivity,
                        )
                    )
                ).sum(axis=0)
            )
            * self.tape_thickness_not_sc
        )