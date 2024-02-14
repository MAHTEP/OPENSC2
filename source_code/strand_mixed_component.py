import numpy as np
import pandas as pd
from scipy import optimize
from typing import Tuple, Union

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

# HTS properties
from properties_of_materials.rare_earth_123 import (
    thermal_conductivity_re123,
    isobaric_specific_heat_re123,
    density_re123,
)
# HTS properties
from properties_of_materials.magnesium_diboride import (
    thermal_conductivity_mgb2,
    isobaric_specific_heat_mgb2,
    density_mgb2,
)


DENSITY_FUNC = dict(
    al=density_al,
    cu=density_cu,
    nb3sn=density_nb3sn,
    nbti=density_nbti,
    ybco=density_re123,
    mgb2=density_mgb2,
)

THERMAL_CONDUCTIVITY_FUNC = dict(
    al=thermal_conductivity_al,
    cu=thermal_conductivity_cu_nist,
    nb3sn=thermal_conductivity_nb3sn,
    nbti=thermal_conductivity_nbti,
    ybco=thermal_conductivity_re123,
    mgb2=thermal_conductivity_mgb2,
)

ISOBARIC_SPECIFIC_HEAT_FUNC = dict(
    al=isobaric_specific_heat_al,
    cu=isobaric_specific_heat_cu_nist,
    nb3sn=isobaric_specific_heat_nb3sn,
    nbti=isobaric_specific_heat_nbti,
    ybco=isobaric_specific_heat_re123,
    mgb2=isobaric_specific_heat_mgb2,
)

ELECTRICAL_RESISTIVITY_FUNC = dict(
    al=electrical_resistivity_al,
    cu=electrical_resistivity_cu_nist,
)

INGEGNERISTIC_MODE = 0
PHYSICAL_MODE = 1
MODE_0 = 0
MODE_1 = 1

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
        self.time_evol = dict(
            temperature=dict(),
            B_field=dict(),
            T_cur_sharing=dict(),
            J_critical=dict(),
            )
        self.time_evol_gauss = dict(
            current_along=dict(),
            delta_voltage_along=dict(),
            linear_power_el_resistance=dict(),
            delta_voltage_along_sum=dict(),
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
        self.__compute_cross_section(sheet, dict_file_path["input"])
        self.__get_current_density_cross_section(sheet, dict_file_path["input"])
        # Call SolidComponent class constructor to deal with StrandMixedComponent time \
        # steps for current, external heating and so on (cdp, 11/2020)
        SolidComponent(simulation, self)
        if self.inputs["stabilizer_material"] != "Cu":
            # remove key RRR from inputs if stabilizer is not Cu (cdp, 07/2020)
            self.inputs.pop("RRR")
        if self.operations["IBIFUN"] != -1:
            # Remove key B_field_units.
            del self.operations["B_field_units"]
        # end if (cdp, 07/2020)

        # Call to method deal_with_flag_IOP_MODE to check and manipulate value
        # of flag self.operations["IOP_MODE"].
        self.deal_with_flag_IOP_MODE()

        # Equivalent radius
        self.radius = np.sqrt(self.inputs["CROSSECTION"] / np.pi)

        self.__strand_density_flag = False
        self.__reorganize_input()
        self.__check_consistency(conductor)

        # Call method deal_with_fixed_potential to manipulate input about fixed
        # potential values.
        self.deal_with_fixed_potential(conductor.inputs["ZLENGTH"])

    def __repr__(self):
        return f"{self.__class__.__name__}(Type: {self.name}, identifier: {self.identifier})"

    def __str__(self):
        pass

    def __compute_cross_section(self,sheet, file_path:str):
        """Private method that computes the sloped cross section of the stabilizer and of the superconductor of the StrandMixedComponent.

        Args:
            sheet (Chartsheet | Worksheet | ReadOnlyWorksheet): sheet with input data for StrandMixedComponent definition.
            file_path (str): path to the input file.

        Raises:
            ValueError: if self.inputs["ISTAB_NON_STAB"] =! 0 or self.inputs["ISTAB_NON_STAB"] != 1
        """
        self.cross_section = dict()
        self.cross_section["total"] = self.inputs["CROSSECTION"]
        # Perpendicular superconducting strand cross section.
        self.cross_section["sc_strand"] = self.inputs["N_sc_strand"] * np.pi * self.inputs["d_sc_strand"] ** 2 / 4
        # Sloped superconducting strand cross section.
        self.cross_section["sc_strand_sloped"] = self.cross_section["sc_strand"] / self.inputs["COSTETA"]
        # Perpendicular segregated stabilizer strand cross section.
        self.cross_section["stab_segr"] = self.inputs["N_stab_strand"] * np.pi * self.inputs["d_stab_strand"] ** 2 / 4
        # Sloped segregated stabilizer strand cross section.
        self.cross_section["stab_segr_sloped"] = self.cross_section["stab_segr"] / self.inputs["COSTETA"]

        # Check cross section consistency.
        self.__check_cross_section(sheet, file_path)

        if self.inputs["ISTAB_NON_STAB"] == MODE_0:
            # Superconductor perpendicular cross section in m^2. Remember that 
            # stab_non_stab is referred only to the sc_strand and does not 
            # include the contribution of the segregated stabilizer.

            # Assign the value of STAB_NON_STAB to self.cross_section["alpha_0"]
            self.cross_section["alpha_0"] = self.inputs["STAB_NON_STAB"]
            # Evaluate total superconducting perpendicular cross section.
            self.cross_section["sc"] = self.cross_section["sc_strand"] / (
                1.0 + self.cross_section["alpha_0"]
            )
            # Evaluate the total stab_non_stab ratio, accounting for the 
            # segregated stabilizer strands:
            # alpha_1 = (A_stab_segr + A_stab_matrix) / A_sc
            # A_stab_matrix = A_sc * alpha
            self.cross_section["alpha_1"] = (self.cross_section["stab_segr"] + self.cross_section["alpha_0"] * self.cross_section["sc"]) / self.cross_section["sc"]
        elif self.inputs["ISTAB_NON_STAB"] == MODE_1:
            # Superconductor perpendicular cross section in m^2. Remember that 
            # stab_non_stab includes both the contribution of the stabilizer 
            # constituting the matrix and of segregated stabilizer.

            # Assign the value of STAB_NON_STAB to self.cross_section["alpha_1"]
            self.cross_section["alpha_1"] = self.inputs["STAB_NON_STAB"]
            # Evaluate total superconducting perpendicular cross section.
            self.cross_section["sc"] = self.inputs["CROSSECTION"] / (
                1.0 + self.cross_section["alpha_1"]
            )
            # Evaluate total not segregated stabilizer perpendicular cross 
            # section.
            self.cross_section["stab_no_segr"] = 1. / ((1. + self.cross_section["alpha_1"])) * (self.cross_section["alpha_1"] * self.cross_section["sc_strand"] - self.cross_section["stab_segr"])
            # Evaluate alpha_0 (stabilizer non stabilizer ratio defined wrt the 
            # not segregated stabilizer cross section only).
            self.cross_section["alpha_0"] = self.cross_section["stab_no_segr"] / self.cross_section["sc"]

        else:
            raise ValueError(
                f"Not valid value for flag self.inputs['ISTAB_NON_STAB']. Flag value should be {MODE_0 = } or {MODE_1 = }, current value is {self.inputs['ISTAB_NON_STAB'] = }.\nPlease check in sheet {sheet} of input file {file_path}.\n"
            )

        # Superconductor sloped cross section in m^2. This definition is 
        # independent from flag ISTAB_NON_STAB.
        self.cross_section["sc_sloped"] = self.cross_section["sc"] / (
            self.inputs["COSTETA"]
        )
        # Stabilizer total perpendicular cross section in m^2. It accounts for 
        # the segregated stabilizer strands and for the stabilizer constituting 
        # the matrix in sc strands.
        self.cross_section["stab"] = self.cross_section["alpha_1"] * self.cross_section["sc"]
        # Stabilizer sloped cross section in m^2. It accounts for the
        # contribution of stabilizer strand cross section and for the fraction
        # of stabilizer in sc strands.
        self.cross_section["stab_sloped"] = self.cross_section["alpha_1"] * self.cross_section["sc_sloped"]

    def __check_cross_section(self, sheet, file_path:str):
        """Private method that checks the consistency of the cross section given in input by the user and the one that can be evaluated starting from the diameters and the number of the strands that constitute the StrandMixedComponent.

        Args:
            sheet (Chartsheet | Worksheet | ReadOnlyWorksheet): sheet with input data for StrandMixedComponent definition.
            file_path (str): path to the input file.

        Raises:
            ValueError: if there is not consistency between the total cross section in input and the evaluated one.
        """

        tot_cross_section = self.cross_section["sc_strand"] + self.cross_section["stab_segr"]

        if not np.isclose(self.cross_section["total"],tot_cross_section):
            raise ValueError(f"Total cross section must be consistent with the evaluated total cross sections: total cross section is {self.cross_section['total']} evaluated total cross {tot_cross_section}. Please, check sheet {sheet.title} in input file {file_path}.")

    def __get_current_density_cross_section(self, sheet, file_path:str):
        """Private method that evalutates cross section used to compute the current density and the current sharing temperature according to the definiton of the critical current density scaling parameter c0.

        Args:
            sheet (Chartsheet | Worksheet | ReadOnlyWorksheet): sheet with input data for StrandMixedComponent definition.
            file_path (str): path to the input file.

        Raises:
            ValueError: if self.inputs["C0_MODE"] =! 0 or self.inputs["C0_MODE"] != 1
        """
                
        if self.inputs["C0_MODE"] == INGEGNERISTIC_MODE:
            # Ingegneristic definition for c0 is used, i.e., the value is 
            # normalized with respect to the total perpendicular cross section 
            # of superconducting strand.

            # Convert the ingegneristic definition to the physical definition 
            # (normalized with respect to the total perpendicular cross section 
            # of the superconductor) in order to use the superconducting cross 
            # section in the electric model.
            self.inputs["c0"] = self.inputs["c0"] * (1. + self.cross_section["alpha_0"])
        elif (self.inputs["C0_MODE"] != PHYSICAL_MODE
                and self.inputs["C0_MODE"] != INGEGNERISTIC_MODE
            ):
            raise ValueError(
                f"Not valid value for flag self.inputs['C0_MODE']. Flag value should be {PHYSICAL_MODE = } or {INGEGNERISTIC_MODE = }, current value is {self.inputs['C0_MODE'] = }.\nPlease check in sheet {sheet} of input file {file_path}.\n"
            )

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
                if key.endswith("material") and key != "superconducting_material"
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

    def __check_consistency(self, conductor):
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
                f"{conductor.identifier = } -> {self.identifier = }\nThe number of material constituting the strand mixed ({self.inputs['NUM_MATERIAL_TYPES'] = }) is inconsistent with the number of defined materials ({self.strand_material.size = }).\nPlease check..."
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
            [func(property["temperature"]) for func in self.density_function]
        )
        # Evaluate homogenized density of the strand mixed:
        # rho_eq = (rho_sc + stab_non_stab * rho_stab)/(1 + stab_non_stab)
        self.__density_numerator = density.T * self.homogenization_cefficients
        self.__density_numerator_sum = self.__density_numerator.sum(axis=1)
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
        isobaric_specific_heat = np.zeros(
            (property["temperature"].size, self.inputs["NUM_MATERIAL_TYPES"])
        )
        for ii, func in enumerate(self.isobaric_specific_heat_function):
            if func.__name__ == "isobaric_specific_heat_nb3sn":
                isobaric_specific_heat[:, ii] = func(
                    property["temperature"],
                    property["T_cur_sharing_min"],
                    property["T_critical"],
                    self.inputs["Tc0m"],
                )
            elif func.__name__ == "isobaric_specific_heat_nbti":
                isobaric_specific_heat[:, ii] = func(
                    property["temperature"],
                    property["B_field"],
                    property["T_cur_sharing_min"],
                    property["T_critical"],
                )
            else:
                isobaric_specific_heat[:, ii] = func(property["temperature"])

        # Evaluate homogenized isobaric specific heat of the strand mixed:
        # cp_eq = (cp_sc*rho_sc + stab_non_stab*cp_stab*rho_stab)/(rho_sc + stab_non_stab * rho_stab)
        return (isobaric_specific_heat * self.__density_numerator).sum(
            axis=1
        ) / self.__density_numerator_sum
        # return (isobaric_specific_heat * self.__density_numerator)/self.__density_numerator_sum.reshape(property["temperature"].size,1)

    def strand_thermal_conductivity(self, property: dict) -> np.ndarray:
        """Method that evaluates the homogenized thermal conductivity of the strand mixed, in the case it is made by two materials (stabilizer and superconductor). Homogenization is based on material cross sections.

        Args:
            property (dict): dictionary with material properties in nodal points or Gauss points according to the value of flag nodal in method eval_sol_comp_properties of class SolidComponent.

        Returns:
            np.ndarray: array with homogenized thermal conductivity of the strand mixed in W/m/K.
        """
        thermal_conductivity = np.zeros(
            (property["temperature"].size, self.inputs["NUM_MATERIAL_TYPES"])
        )
        for ii, func in enumerate(self.thermal_conductivity_function):
            if "cu" in func.__name__:
                thermal_conductivity[:, ii] = func(
                    property["temperature"],
                    property["B_field"],
                    self.inputs["RRR"],
                )
            elif "mgb2" in func.__name__:
                thermal_conductivity[:, ii] = func(
                    property["temperature"],
                    property["B_field"],
                )
            else:
                thermal_conductivity[:, ii] = func(property["temperature"])
        # Evaluate homogenized thermal conductivity of the strand mixed:
        # k_eq = (K_sc + stab_non_stab * K_stab)/(1 + stab_non_stab)
        return (thermal_conductivity * self.homogenization_cefficients).sum(
            axis=1
        ) / self.homogenization_cefficients_sum

    def strand_electrical_resistivity_not_sc(self, property: dict) -> np.ndarray:
        """Method that evaluates electrical resisitivity of the stabilizer material of the strand mixed.

        Args:
            property (dict): dictionary with material properties in nodal points or Gauss points according to the value of flag nodal in method eval_sol_comp_properties of class SolidComponent.

        Returns:
            np.ndarray: array with electrical resisitivity of the stabilizer of the strand mixed in Ohm*m.
        """
        electrical_resistivity = np.zeros(
            (property["temperature"].size, self.inputs["NUM_MATERIAL_TYPES"] - 1)
        )
        for ii, func in enumerate(self.electrical_resistivity_function_not_sc):
            if "cu" in func.__name__:
                electrical_resistivity[:, ii] = func(
                    property["temperature"],
                    property["B_field"],
                    self.inputs["RRR"],
                )
            else:
                electrical_resistivity[:, ii] = func(property["temperature"])
        if self.inputs["NUM_MATERIAL_TYPES"] - 1 > 1:
            # Evaluate homogenized electrical resistivity of the strand mixed:
            # rho_el_eq = A_not_sc * (sum(A_i/rho_el_i))^-1 for any i not sc
            return self.cross_section["stab_sloped"] * np.reciprocal((self.cross_section["stab_sloped"] / electrical_resistivity).sum(axis=1))
        if self.inputs["NUM_MATERIAL_TYPES"] - 1 == 1:
            return electrical_resistivity.reshape(property["temperature"].size)

    def superconductor_power_law(
        self,
        current: np.ndarray,
        critical_current: np.ndarray,
        critical_current_density: np.ndarray,
    ) -> np.ndarray:
        """Method that evaluate the electrical resistivity of superconducting material according to the power law.

        Args:
            current (np.ndarray): electric current in superconducting material
            critical_current (np.ndarray): critical current of superconducting material.
            critical_current_density (np.ndarray): critical current density in superconducting material

        Raises:
            ValueError: if arrays current and critical_current does not have the same shape.
            ValueError: if arrays current and critical_current_density does not have the same shape.
            ValueError: if arrays critical_current and critical_current_density does not have the same shape.

        Returns:
            np.ndarray: electrical resistivity of superconducting material in Ohm*m.
        """
        # Check input arrays shape consistency
        if current.shape != critical_current.shape:
            raise ValueError(
                f"Arrays current and critical_current must have the same shape.\n{current.shape = };\n{critical_current.shape = }.\n"
            )
        elif current.shape != critical_current_density.shape:
            raise ValueError(
                f"Arrays current and critical_current_density must have the same shape.\n{current.shape = };\n{critical_current_density.shape = }.\n"
            )
        elif critical_current.shape != critical_current_density.shape:
            raise ValueError(
                f"Arrays critical_current and critical_current_density must have the same shape.\n{critical_current.shape = };\n{critical_current_density.shape = }.\n"
            )

        # Evaluate superconductin electrical resistivity according to the power
        # low scaling:
        # rho_el_sc = E_0 / j_c * (I_sc/I_c)**(n-1) Ohm*m
        return (
            self.inputs["E0"]
            / critical_current_density
            * (current / critical_current) ** (self.inputs["nn"] - 1)
        )

    def solve_current_divider(
        self,
        rho_el_stabilizer: np.ndarray,
        critical_current: np.ndarray,
        current: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Method that solves the not linear system of the current divider between superconduting and stabilizer material in the case of current sharing regime.

        Args:
            rho_el_stabilizer (np.ndarray): array with stabilizer electrical resistivity in Ohm*m.
            critical_current (np.ndarray): array with superconductor critical current in A.
            current (np.ndarray): electric total current array in A.

        Raises:
            ValueError: if arrays rho_el_stabilizer and critical_current does not have the same shape.

        Returns:
            Tuple[np.ndarray, np.ndarray]: superconducting current array in A, stabilizer current array in A.
        """
        # Check array shape
        if rho_el_stabilizer.shape != critical_current.shape:
            raise ValueError(
                f"Arrays rho_el_stabilizer and critical_current must have the same shape.\n {rho_el_stabilizer.shape = };\n{critical_current.shape}.\n"
            )
        if rho_el_stabilizer.shape != current.shape:
            raise ValueError(
                f"Arrays rho_el_stabilizer and current must have the same shape.\n {rho_el_stabilizer.shape = };\n{current.shape}.\n"
            )
        if critical_current.shape != current.shape:
            raise ValueError(
                f"Arrays critical_current and current must have the same shape.\n {critical_current.shape = };\n{current.shape}.\n"
            )

        # Evaluate constant value:
        # psi = rho_el_stab*I_c^n/(E_0*A_stab)
        psi = (
            rho_el_stabilizer
            * critical_current ** self.inputs["nn"]
            / self.inputs["E0"]
            / self.cross_section["stab"]
        )

        # Initialize guess.
        sc_current_guess = np.zeros(current.shape)
        for ii, val in enumerate(current):
            # Evaluate superconducting current guess with bisection method.
            # Set the maximum itaration to 10 and disp to False in order to not
            # rise an error due to not reached convergence.
            sc_current_guess[ii] = optimize.bisect(
                self.__sc_current_residual,
                0.0,
                val,
                args=(psi[ii], val),
                maxiter=10,
                disp=False,
            )
        
        # Tolerance on newton halley increased in case I_critical is very small,
        # to avoid inaccuracies on the divider that could lead to potential 
        # differences between sc and stab
        if min(critical_current)>1e-6:
            # Default tollerance in optimize.newton method
            tollerance = 1.48e-8
        else:
            # Value found trial and error iteration
            tollerance = 1e-12
            # Other possible solution for the correct tollerance
            # tollerance = min(critical_current)/1000

        # Evaluate superconducting with Halley's method
        sc_current = optimize.newton(
            self.__sc_current_residual,
            sc_current_guess,
            args=(psi, current),
            fprime=self.__d_sc_current_residual,
            fprime2=self.__d2_sc_current_residual,
            tol = tollerance,
            maxiter=1000,
        )

        return sc_current, current - sc_current

    def __sc_current_residual(
        self,
        sc_current: Union[float, np.ndarray],
        psi: Union[float, np.ndarray],
        so_current: Union[float, np.ndarray],
    ) -> Union[float, np.ndarray]:
        """Private method that defines the residual function for the evaluation of the superconducting current with bysection and/or Newton-Rampson methods.

        Args:
            sc_current (Union[float, np.ndarray]): superconducting current (guess) in A.
            psi (Union[float, np.ndarray]): costant value in the equation
            so_current (Union[float, np.ndarray]): total current in A.

        Raises:
           Raises:
            ValueError: if arguments sc_current, psi and so_current are not of the same type (float).
            ValueError: if arguments sc_current, psi and so_current are not of the same type (np.ndarray).
            ValueError: if arrays sc_current and psi does not have the same shape.
            ValueError: if arrays sc_current and so_current does not have the same shape.

        Returns:
            Union[float, np.ndarray]: residual value
        """
        # Checks on input arguments.
        if isinstance(sc_current, float):
            if (
                isinstance(psi, float) == False
                or isinstance(so_current, float) == False
            ):
                raise ValueError(
                    f"Arguments sc_current, psi and so_current must be of the same type (float).\n{type(sc_current) = };\n{type(psi) = };\n{type(so_current) = }.\n"
                )
        if isinstance(sc_current, np.ndarray):
            if (
                isinstance(psi, np.ndarray) == False
                or isinstance(so_current, np.ndarray) == False
            ):
                raise ValueError(
                    f"Arguments sc_current, psi and so_current must be of the same type (np.ndarray).\n{type(sc_current) = };\n{type(psi) = };\n{type(so_current) = }.\n"
                )
            if sc_current.shape != psi.shape:
                raise ValueError(
                    f"Arrays sc_current and psi must have the same shape.\n {sc_current.shape = };\n{psi.shape}.\n"
                )
            if sc_current.shape != so_current.shape:
                raise ValueError(
                    f"Arrays sc_current and so_current must have the same shape.\n {sc_current.shape = };\n{so_current.shape}.\n"
                )

        return sc_current ** self.inputs["nn"] + (sc_current - so_current) * psi

    def __d_sc_current_residual(
        self,
        sc_current: Union[float, np.ndarray],
        psi: Union[float, np.ndarray],
        so_current: Union[float, np.ndarray],
    ) -> Union[float, np.ndarray]:
        """Private method that defines the first derivative of residual function wrt sc_current for the evaluation of the superconducting current with Newton-Rampson or Halley's methods.

        Args:
            sc_current (Union[float, np.ndarray]): superconducting current (guess).
            psi (Union[float, np.ndarray]): costant value in the equation.
            so_current (Union[float, np.ndarray]): total current in A, not used but passed by function optimize.newton.

        Returns:
            Union[float, np.ndarray]: residual derivative value
        """

        return self.inputs["nn"] * sc_current ** (self.inputs["nn"] - 1) + psi

    def __d2_sc_current_residual(
        self,
        sc_current: Union[float, np.ndarray],
        psi: Union[float, np.ndarray],
        so_current: Union[float, np.ndarray],
    ) -> Union[float, np.ndarray]:
        """Private method that defines the second derivative of residual function wrt sc_current for the evaluation of the superconducting current with Newton-Rampson or Halley's methods.

        Args:
            sc_current (Union[float, np.ndarray]): superconducting current (guess).
            psi (Union[float, np.ndarray]): costant value in the equation, not needed for this fuction
            so_current (Union[float, np.ndarray]): total current in A, not used but passed by function optimize.newton

        Returns:
            Union[float, np.ndarray]: second derivative of the residual.
        """

        return (
            self.inputs["nn"]
            * (self.inputs["nn"] - 1)
            * sc_current ** (self.inputs["nn"] - 2)
        )

    def get_electric_resistance(self, conductor: object) -> np.ndarray:
        f"""Method that evaluate the electrical resistance in Gauss node only, used to build the electric_resistance_matrix.

        Args:
            conductor (object): class Conductor object from which node distance is stored to do the calculation.

        Returns:
            np.ndarray: array of electrical resistance in Ohm of length {conductor.grid_input["NELEMS"] = }.
        """

        critical_current_gauss = (
            self.cross_section["sc"] * self.dict_Gauss_pt["J_critical"]
        )

        # Make initialization only once for each conductor object.
        if conductor.cond_num_step == 0:
            # Initialize electric resistance arrays in Gauss point; this is the 
            # equivalent electrical resistance, thus it is defined in this way:
            # R_eq = R_sc if superconducting regime
            # R_eq = R_sc * R_stab/(R_sc + R_stab) is normal regime.
            self.dict_Gauss_pt["electric_resistance"] = np.full_like(
                self.dict_Gauss_pt["temperature"], None
            )

            # Initialize array of superconducting electrical resistivit in 
            # Gauss point to None.
            self.dict_Gauss_pt["electrical_resistivity_superconductor"] = np.full_like(
                self.dict_Gauss_pt["temperature"], None
            )

        # Get index for which abs(critical_current_gauss) == 0 (inside normal 
        # zone by definition).
        ind_zero = np.nonzero(abs(critical_current_gauss) == 0)[0]
        # Get index for which abs(critical_current_gauss) > 0 (outside normal 
        # zone by definition).
        ind_not_zero = np.nonzero(abs(critical_current_gauss) > 0)[0]

        # Check if np array ind_zero is not empty: NORMAL REGION BY DEFINITION
        if ind_zero.any():
            # Evaluate electic resistance in normal region (stabilizer only).
            self.dict_Gauss_pt["electric_resistance"][
                ind_zero
            ] = self.electric_resistance(
                conductor, "electrical_resistivity_stabilizer", "stab", ind_zero
            )

        # Check if np array ind_not_zero is not empty: deal with index that 
        # are outside normal zone by definition; however some of them could 
        # still identify a normal region.
        if ind_not_zero.any():

            sc_current_gauss, stab_current_gauss = self.solve_current_divider(
                self.dict_Gauss_pt["electrical_resistivity_stabilizer"][ind_not_zero],
                critical_current_gauss[ind_not_zero],
                self.dict_Gauss_pt["op_current"][ind_not_zero]
                )
            
            self.dict_Gauss_pt["electrical_resistivity_superconductor"][
                ind_not_zero
            ] = self.superconductor_power_law(
                sc_current_gauss,
                critical_current_gauss[ind_not_zero],
                self.dict_Gauss_pt["J_critical"][ind_not_zero],
                )
            # Evaluate the equivalent electric resistance in Ohm.
            self.dict_Gauss_pt["electric_resistance"][
                ind_not_zero] = self.parallel_electric_resistance(
                conductor,
                [
                    "electrical_resistivity_superconductor",
                    "electrical_resistivity_stabilizer",
                ],["sc","stab"],
                ind_not_zero,
            )

            # Compute voltage along stabilizer.
            v_stab = self.dict_Gauss_pt["electrical_resistivity_stabilizer"][
                ind_not_zero
            ] * stab_current_gauss / self.cross_section["stab"]
            # Compute voltage along superconductor
            v_sc = self.inputs["E0"] * (sc_current_gauss / critical_current_gauss[ind_not_zero]) ** self.inputs["nn"]
            # Check that the voltage along stabilizer is equal to the 
            # voltage along superconductor (i.e, check the reliability 
            # of the current divider).
            if all(np.isclose(v_stab,v_sc)) == False:
                raise ValueError(f"Voltage difference along superconductor and stabilizer must be the same.")

        return self.dict_Gauss_pt["electric_resistance"]