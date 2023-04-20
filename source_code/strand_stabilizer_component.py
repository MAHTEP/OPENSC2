import numpy as np
import pandas as pd
from solid_component import SolidComponent
from strand_component import StrandComponent

from utility_functions.auxiliary_functions import check_costheta

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

DENSITY_FUNC = dict(
    al=density_al,
    cu=density_cu,
)

THERMAL_CONDUCTIVITY_FUNC = dict(
    al=thermal_conductivity_al,
    cu=thermal_conductivity_cu_nist,
)

ISOBARIC_SPECIFIC_HEAT_FUNC = dict(
    al=isobaric_specific_heat_al,
    cu=isobaric_specific_heat_cu_nist,
)

ELECTRICAL_RESISTIVITY_FUNC = dict(
    al=electrical_resistivity_al,
    cu=electrical_resistivity_cu_nist,
)


class StrandStabilizerComponent(StrandComponent):

    # Class for copper strands objects

    ### INPUT PARAMETERS
    # some are inherited form the parent classes StrandComponent and SolidComponent

    ### THERMOPHYSICAL PROPERTIES
    # inherited from class SolidComponent

    ##### OPERATIONAL PARAMETERS
    # inherited from parents class SolidComponent and StrandComponent

    ### COMPUTED IN INITIALIZATION
    # inherited from class SolidComponent

    ### COMPUTED VECTOR FOR MAGNETIC FIELD
    # inherited from class SolidComponent

    KIND = "StrandStabilizerComponent"

    def __init__(
        self,
        simulation: object,
        sheet,
        icomp: int,
        name: str,
        dict_file_path: dict,
        conductor: object,
    ):
        """Method that makes instance of class StrandMixedComponent.

        Args:
            simulation (object): simulation object.
            sheet (Worksheet): worksheet with input data.
            icomp (int): component index.
            name (str): component name.
            dict_file_path (dict): dictionary with paths to load the input files.
            conductor (object): inscance of class Conductor
        """

        self.name = name
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
        self.time_evol = dict(temperature=dict(), B_field=dict())
        self.time_evol_gauss = dict(
            current_along=dict(),
            voltage_drop_along=dict(),
            linear_power_el_resistance=dict(),
            voltage_drop_along_sum=dict(),
        )
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

        # Check that costheta is in the range (0,1].
        check_costheta(self,dict_file_path["input"],sheet)

        # Call SolidComponent class constructor to deal with StrandStabilizerComponent time \
        # steps for current, external heating and so on (cdp, 11/2020)
        SolidComponent(simulation, self)
        self.inputs["stabilizer_material"] = self.inputs["stabilizer_material"].lower()
        if self.inputs["stabilizer_material"] != "cu":
            # remove key RRR from inputs if stabilizer is not Cu (cdp, 07/2020)
            self.inputs.pop("RRR")
        if self.operations["IBIFUN"] != -1:
            # Remove key B_field_units.
            del self.operations["B_field_units"]

        # Call to method deal_with_flag_IOP_MODE to check and manipulate value
        # of flag self.operations["IOP_MODE"].
        self.deal_with_flag_IOP_MODE()

        # Equivalent radius
        self.radius = np.sqrt(self.inputs["CROSSECTION"] / np.pi)

        # Call method deal_with_fixed_potential to manipulate input about fixed
        # potential values.
        self.deal_with_fixed_potential(conductor.inputs["ZLENGTH"])

    def __repr__(self):
        return f"{self.__class__.__name__}(Type: {self.name}, identifier: {self.identifier})"

    def __str__(self):
        pass

    def strand_density(self, property: dict) -> np.ndarray:
        """Method that evaluates density of the stabilizer.

        Args:
            property (dict): dictionary with material properties in nodal points or Gauss points according to the value of flag nodal in method eval_sol_comp_properties of class SolidComponent.

        Returns:
            np.ndarray: array with density of the stabilizer in kg/m^3.
        """
        return DENSITY_FUNC[self.inputs["stabilizer_material"]](property["temperature"])

    def strand_isobaric_specific_heat(self, property: dict) -> np.ndarray:
        """Method that evaluates isobaric specific heat of the stabilizer.

        Args:
            property (dict): dictionary with material properties in nodal points or Gauss points according to the value of flag nodal in method eval_sol_comp_properties of class SolidComponent.

        Returns:
            np.ndarray: array with isobaric specific heat of the stabilizer in kg/m^3.
        """
        return ISOBARIC_SPECIFIC_HEAT_FUNC[self.inputs["stabilizer_material"]](
            property["temperature"]
        )

    def strand_thermal_conductivity(self, property: dict) -> np.ndarray:
        """Method that evaluates thermal conductivity of the stabilizer.

        Args:
            property (dict): dictionary with material properties in nodal points or Gauss points according to the value of flag nodal in method eval_sol_comp_properties of class SolidComponent.

        Returns:
            np.ndarray: array with thermal conductivity of the stabilizer in W/m/K.
        """
        if self.inputs["stabilizer_material"] == "cu":
            return THERMAL_CONDUCTIVITY_FUNC[self.inputs["stabilizer_material"]](
                property["temperature"],
                property["B_field"],
                self.inputs["RRR"],
            )
        else:
            return THERMAL_CONDUCTIVITY_FUNC[self.inputs["stabilizer_material"]](
                property["temperature"]
            )

    def strand_electrical_resistivity(self, property: dict) -> np.ndarray:
        """Method that evaluates electrical resistivity of the stabilizer.

        Args:
            property (dict): dictionary with material properties in nodal points or Gauss points according to the value of flag nodal in method eval_sol_comp_properties of class SolidComponent.

        Returns:
            np.ndarray: array with electrical resistivity of the stabilizer in Ohm*m.
        """
        if self.inputs["stabilizer_material"] == "cu":
            return ELECTRICAL_RESISTIVITY_FUNC[self.inputs["stabilizer_material"]](
                property["temperature"],
                property["B_field"],
                self.inputs["RRR"],
            )
        else:
            return ELECTRICAL_RESISTIVITY_FUNC[self.inputs["stabilizer_material"]](
                property["temperature"]
            )

    def get_electric_resistance(self, conductor: object) -> np.ndarray:
        f"""Method that evaluate the electric resistance in Gauss node only, used to build the electric_resistance_matrix.

        Args:
            conductor (object): class Conductor object from which node distance is stored to do the calculation.

        Returns:
            np.ndarray: array of electrical resistance in Ohm of length {conductor.grid_input["NELEMS"] = }.
        """
        return (
            self.dict_Gauss_pt["electrical_resistivity_stabilizer"]
            * conductor.node_distance[("StrandComponent", self.identifier)]
            / self.inputs["CROSSECTION"]
        )