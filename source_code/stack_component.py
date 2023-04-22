# Import libraries
import numpy as np
import pandas as pd
from scipy import optimize
from typing import Tuple, Union

# Import classes
from solid_component import SolidComponent
from strand_component import StrandComponent

from utility_functions.auxiliary_functions import check_costheta

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
        simulation: object,
        sheet,
        icomp: int,
        name: str,
        dict_file_path: dict,
        conductor: object,
    ):
        """Method that makes instance of class StackComponent.

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

        # dictionary declaration
        self.inputs = dict()
        self.operations = dict()
        self.dict_node_pt = dict()
        self.dict_Gauss_pt = dict()
        self.dict_num_step = dict()
        self.coordinate = dict()
        # Empty dictionary of list to save variable time evolutions at selected spatial coordinates.
        self.time_evol = dict(temperature=dict(), B_field=dict(), T_cur_sharing=dict())
        self.time_evol_gauss = dict(
            current_along=dict(),
            delta_voltage_along=dict(),
            linear_power_el_resistance=dict(),
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

        # Call SolidComponent class constructor to deal with StrandMixedComponent time \
        # steps for current, external heating and so on
        SolidComponent(simulation, self)
        if self.operations["IBIFUN"] != -1:
            # Remove key B_field_units.
            del self.operations["B_field_units"]
        # end if

        # Call to method deal_with_flag_IOP_MODE to check and manipulate value
        # of flag self.operations["IOP_MODE"].
        self.deal_with_flag_IOP_MODE()

        self.__reorganize_input()
        self.__check_consistency(conductor)

        # Flag to check if evaluation of homogenized isobaric specific heat can
        # be done or not (depends on homogenized density evaluation).
        self.__stack_density_flag = False

        # Call method deal_with_fixed_potential to manipulate input about fixed
        # potential values.
        self.deal_with_fixed_potential(conductor.inputs["ZLENGTH"])

        # Superconductor total cross section in m^2
        self.sc = (
            self.inputs["HTS_thickness"]
            * self.inputs["Stack_width"]
            * self.inputs["Tape_number"]
        )
        # Stabilizer (not sc) total cross section in m^2
        self.stabilizer_cross_section = (
            self.tape_thickness_not_sc
            * self.inputs["Stack_width"]
            * self.inputs["Tape_number"]
        )

    def __repr__(self):
        return f"{self.__class__.__name__}(Type: {self.name}, identifier: {self.identifier})"

    def __str__(self):
        pass

    def __reorganize_input(self):
        """Private method that reorganizes input data stored in dictionary self.inputs to simplify the procedure of properties homogenization."""

        # Create numpy array of string with the identifier of tape material
        self.tape_material = np.array(
            [
                value.lower()
                for key, value in self.inputs.items()
                if key.endswith("material")
            ],
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

    def __check_consistency(self, conductor):
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
        if self.tape_material.size != self.inputs["Material_number"]:
            # Message to be completed!
            raise ValueError(
                f"{conductor.identifier = } -> {self.identifier = }\nThe number of material constituting the tape ({self.inputs['Material_number'] = }) is inconsistent with the number of defined materials ({self.tape_material.size = }).\nPlease check..."
            )

        # Check that number of tape materials given in input is consistent with
        # not zero user defined material thicknes.
        if self.material_thickness.size != self.inputs["Material_number"]:
            # Message to be completed!
            raise ValueError(
                f"{conductor.identifier = } -> {self.identifier = }\nThe number of material constituting the tape ({self.inputs['Material_number'] = }) is inconsistent with the number of defined thicknesses ({self.material_thickness.size = }).\nPlease check..."
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
        self.__tape_number = round(self.inputs["CROSSECTION"] / self.tape_cross_section)
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
            abs(self.__cross_section - self.inputs["CROSSECTION"])
            / self.inputs["CROSSECTION"]
            > tol
        ):
            # Message to be completed!
            raise ValueError(
                f"{conductor.identifier = } -> {self.identifier = }\nInconsistent cross section values: user defines {self.inputs['CROSSECTION'] = } while computed one is {self.__cross_section = }.\nPlease check..."
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
            [func(property["temperature"]) for func in self.density_function]
        )
        # Evaluate homogenized density of the stack:
        # rho_eq = sum(s_i*rho_i)/s
        self.__density_numerator = density.T * self.material_thickness
        self.__density_numerator_sum = self.__density_numerator.sum(axis=1)
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
        return (isobaric_specific_heat.T * self.__density_numerator).sum(
            axis=1
        ) / self.__density_numerator_sum

    def stack_thermal_conductivity(self, property: dict) -> np.ndarray:
        """Method that evaluates the homogenized thermal conductivity of the stack, which is the same of the tape if the tapes constituting the stack are equals to each other. Homogenization is based on the thickness of tape layers.

        Args:
            property (dict): dictionary with material properties in nodal points or Gauss points according to the value of flag nodal in method eval_sol_comp_properties of class SolidComponent.

        Returns:
            np.ndarray: array with homogenized thermal conductivity of the stack of tapes in W/m/K.
        """
        thermal_conductivity = np.zeros(
            (property["temperature"].size, self.inputs["Material_number"])
        )
        for ii, func in enumerate(self.thermal_conductivity_function):
            if "cu" in func.__name__:
                thermal_conductivity[:, ii] = func(
                    property["temperature"],
                    property["B_field"],
                    self.inputs["RRR"],
                )
            else:
                thermal_conductivity[:, ii] = func(property["temperature"])
        # Evaluate homogenized thermal conductivity of the stack:
        # k_eq = sum(s_i*k_i)/s
        return (thermal_conductivity * self.material_thickness).sum(
            axis=1
        ) / self.tape_thickness

    def stack_electrical_resistivity_not_sc(self, property: dict) -> np.ndarray:
        """Method that evaluates the homogenized electrical resistivity for not superconducting materials of the stack, which is the same of the tape if the tapes constituting the stack are equals to each other. Homogenization is based on the thickness of tape layers.

        Args:
            property (dict): dictionary with material properties in nodal points or Gauss points according to the value of flag nodal in method eval_sol_comp_properties of class SolidComponent.

        Returns:
            np.ndarray: array with homogenized electrical resistivity of not superconducting materials of the stack of tapes in Ohm*m.
        """
        electrical_resistivity = np.zeros(
            (property["temperature"].size, self.inputs["Material_number"] - 1)
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
        if self.inputs["Material_number"] - 1 > 1:
            # Evaluate homogenized electrical resistivity of the stack:
            # rho_el_eq = s_not_sc * (sum(s_i/rho_el_i))^-1 for any i not sc
            return self.tape_thickness_not_sc * np.reciprocal(
                (self.material_thickness_not_sc / electrical_resistivity).sum(axis=1)
            )
        elif self.inputs["Material_number"] - 1 == 1:
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
            ind (np.ndarray): array with the index of the location at which electrical resistivity should be evaluated.

        Raises:
            ValueError: if arrays current and critical_current does not have the same shape.
            ValueError: if arrays current and critical_current_density does not have the same shape.
            ValueError: if arrays critical_current and critical_current_density does not have the same shape.

        Returns:
            np.ndarray: electrical resistivity of superconducting material in Ohm*m.
        """
        # Check input arrays shape consistency.
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

        # Evaluate superconducting electrical resistivity according to the power
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

        # Evaluate constant value:
        # psi = rho_el_stab*I_c^n/(E_0*A_stab)
        psi = (
            rho_el_stabilizer
            * critical_current ** self.inputs["nn"]
            / self.inputs["E0"]
            / self.stabilizer_cross_section
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

        # Evaluate superconducting with Halley's method
        sc_current = optimize.newton(
            self.__sc_current_residual,
            sc_current_guess,
            args=(psi, current),
            fprime=self.__d_sc_current_residual,
            fprime2=self.__d2_sc_current_residual,
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
            sc_current (Union(float, np.ndarray)): superconducting current (guess).
            psi (Union(float, np.ndarray)): costant value in the equation
            so_current (Union[float, np.ndarray]): total current in A.

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
            sc_current (Union(float, np.ndarray)): superconducting current (guess).
            psi (Union(float, np.ndarray)): costant value in the equation.
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
            sc_current (Union(float, np.ndarray)): superconducting current (guess).
            psi (Union(float, np.ndarray)): costant value in the equation, not needed for this fuction.
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

        critical_current_node = self.sc_cross_section * self.dict_node_pt["J_critical"]
        critical_current_gauss = (
            self.sc_cross_section * self.dict_Gauss_pt["J_critical"]
        )

        # Get index that correspond to superconducting regime.
        ind_sc_node = np.nonzero(
            self.dict_node_pt["op_current_sc"] / critical_current_node < 0.95
        )[0]
        ind_sc_gauss = np.nonzero(
            self.dict_Gauss_pt["op_current_sc"] / critical_current_gauss < 0.95
        )[0]
        # Get index that correspond to current sharing regime.
        ind_sh_node = np.nonzero(
            self.dict_node_pt["op_current_sc"] / critical_current_node >= 0.95
        )[0]
        ind_sh_gauss = np.nonzero(
            self.dict_Gauss_pt["op_current_sc"] / critical_current_gauss >= 0.95
        )[0]

        # Initialize electric resistance arrays in both nodal and Gauss points;
        # this is the equivalent electrical resistance, thus it is defined in
        # this way:
        # R_eq = R_sc if superconducting regime
        # R_eq = R_sc * R_stab/(R_sc + R_stab) is sharing or normal regime
        self.dict_node_pt["electric_resistance"] = 10.0 * np.ones(
            self.dict_node_pt["temperature"].shape
        )
        self.dict_Gauss_pt["electric_resistance"] = 10.0 * np.ones(
            self.dict_Gauss_pt["temperature"].shape
        )

        ## SUPERCONDUCTING REGIME ##

        # Strand current in superconducting regime is the one carriend by the
        # superconducting material only.
        self.dict_node_pt["op_current"][ind_sc_node] = self.dict_node_pt[
            "op_current_sc"
        ][ind_sc_node]
        self.dict_Gauss_pt["op_current"][ind_sc_gauss] = self.dict_Gauss_pt[
            "op_current_sc"
        ][ind_sc_gauss]

        # Initialize array of superconducting electrical resistivit in nodal and Gauss points to None.
        self.dict_node_pt["electrical_resistivity_superconductor"] = np.full_like(
            self.dict_node_pt["temperature"], None
        )
        self.dict_Gauss_pt["electrical_resistivity_superconductor"] = np.full_like(
            self.dict_Gauss_pt["temperature"], None
        )

        # Compute superconducting electrical resistivity only in index for
        # which the superconducting regime is guaranteed, using the power low.
        self.dict_node_pt["electrical_resistivity_superconductor"][
            ind_sc_node
        ] = self.superconductor_power_law(
            self.dict_node_pt["op_current_sc"][ind_sc_node],
            critical_current_node[ind_sc_node],
            self.dict_node_pt["J_critical"][ind_sc_node],
        )
        self.dict_Gauss_pt["electrical_resistivity_superconductor"][
            ind_sc_gauss
        ] = self.superconductor_power_law(
            self.dict_Gauss_pt["op_current_sc"][ind_sc_gauss],
            critical_current_gauss[ind_sc_gauss],
            self.dict_Gauss_pt["J_critical"][ind_sc_gauss],
        )
        # Evaluate electic resistance in superconducting region (superconductor
        # only).
        self.dict_Gauss_pt["electric_resistance"][
            ind_sc_gauss
        ] = self.electric_resistance(
            conductor, "electrical_resistivity_superconductor", ind_sc_gauss
        )

        ## SHARING OR NORMAL REGIME ##

        # Strand current in sharing regime is the one carried by the both the
        # superconducting and the stabilizer materials.
        # self.dict_node_pt["op_current"][ind_sh_node] = self.dict_node_pt["op_current"][ind_sh_node]
        # self.dict_Gauss_pt["op_current"][ind_sh_node] = self.dict_node_pt["op_current"][ind_sh_gauss]

        # Evaluate how the current is distributed solving the current divider
        # problem in both nodal and Gauss points.
        sc_current_node, stab_current_node = self.solve_current_divider(
            self.dict_node_pt["electrical_resistivity_stabilizer"][ind_sh_node],
            critical_current_node[ind_sh_node],
            self.dict_node_pt["op_current"][ind_sh_node],
        )
        sc_current_gauss, stab_current_gauss = self.solve_current_divider(
            self.dict_Gauss_pt["electrical_resistivity_stabilizer"][ind_sh_gauss],
            critical_current_gauss[ind_sh_gauss],
            self.dict_Gauss_pt["op_current"][ind_sh_gauss],
        )

        # Get index of the normal region, to avoid division by 0 in evaluation
        # of sc electrical resistivity with the power law.
        ind_normal_node = np.nonzero(
            (
                stab_current_node / self.dict_node_pt["op_current"][ind_sh_node]
                > 0.999999
            )
            | (sc_current_node < 1.0)
        )[0]
        ind_normal_gauss = np.nonzero(
            (
                stab_current_gauss / self.dict_Gauss_pt["op_current"][ind_sh_gauss]
                > 0.999999
            )
            | (sc_current_gauss < 1.0)
        )[0]

        ## NORMAL REGIME ONLY ##
        if ind_normal_node.any():
            # Get the index of location of true current sharing region;
            # overwrite ind_sh_node.
            ind_sh_node = np.nonzero(
                (
                    stab_current_node / self.dict_node_pt["op_current"][ind_sh_node]
                    <= 0.999999
                )
                | (sc_current_node >= 1.0)
            )[0]
        if ind_normal_gauss.any():
            # Get the index of location of true current sharing region;
            # overwrite ind_sh_gauss.
            ind_sh_gauss = np.nonzero(
                (
                    stab_current_gauss / self.dict_Gauss_pt["op_current"][ind_sh_gauss]
                    <= 0.999999
                )
                | (sc_current_gauss >= 1.0)
            )[0]
            # Evaluate electic resistance in normal region (stabilizer only).
            self.dict_Gauss_pt["electric_resistance"][
                ind_normal_gauss
            ] = self.electric_resistance(
                conductor, "electrical_resistivity_stabilizer", ind_normal_gauss
            )

        ## SHARING REGIME ONLY ##
        # Evaluate the electrical resistivity of the superconductor according
        # to the power low in both nodal and Gauss points in Ohm*m.
        self.dict_node_pt["electrical_resistivity_superconductor"][
            ind_sh_node
        ] = self.superconductor_power_law(
            sc_current_node[ind_sh_node],
            critical_current_node[ind_sh_node],
            self.dict_node_pt["J_critical"][ind_sh_node],
        )
        self.dict_Gauss_pt["electrical_resistivity_superconductor"][
            ind_sh_gauss
        ] = self.superconductor_power_law(
            sc_current_gauss[ind_sh_gauss],
            critical_current_gauss[ind_sh_gauss],
            self.dict_Gauss_pt["J_critical"][ind_sh_gauss],
        )

        # Evaluate the equivalent electric resistance in Ohm.
        self.dict_Gauss_pt["electric_resistance"][
            ind_sh_gauss
        ] = self.parallel_electric_resistance(
            conductor,
            [
                "electrical_resistivity_superconductor",
                "electrical_resistivity_stabilizer",
            ],
            ind_sh_gauss,
        )

        return self.dict_Gauss_pt["electric_resistance"]
