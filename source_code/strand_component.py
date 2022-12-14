import warnings
from solid_component import SolidComponent
from openpyxl import load_workbook
import numpy as np
import os
from utility_functions.auxiliary_functions import (
    get_from_xlsx,
    load_auxiliary_files,
    build_interpolator,
    do_interpolation,
)

# from utility_functions.InitializationFunctions import Read_input_file
# NbTi properties
from properties_of_materials.niobium_titanium import (
    critical_temperature_nbti,
    critical_current_density_nbti,
    current_sharing_temperature_nbti,
)

# Nb3Sn properties
from properties_of_materials.niobium3_tin import (
    critical_temperature_nb3sn,
    critical_current_density_nb3sn,
    current_sharing_temperature_nb3sn,
)

# RE123 properties
from properties_of_materials.rare_earth_123 import (
    critical_magnetic_field_re123,
    critical_current_density_re123,
    current_sharing_temperature_re123,
)


class StrandComponent(SolidComponent):

    ### INPUT PARAMETERS

    ### OPERATIONAL PARAMETERS
    # inherited from class SolidComponent

    ### COMPUTED IN INITIALIZATION
    # inherited from class SolidComponent

    KIND = "Strand"

    def get_magnetic_field_gradient(self, conductor, nodal=True):

        """
        ############################################################################
        #              Get_alphaB(self, comp)
        ############################################################################
        #
        # Method that initialize magnetic field gradient in python objects of class
        # Strands.
        #
        ############################################################################
        # VARIABLE              I/O    TYPE              DESCRIPTION            UNIT
        # --------------------------------------------------------------------------
        # comp                  I      object            python object of
        #                                                class Strands             -
        # zcoord*               I      np array float    conductor spatial
        #                                                discretization            m
        # I0_OP_MODE*               I      scalar integer    flag to decide how to
        #                                                evaluate current:
        #                                                == 0 -> constant;
        #                                                == 1 -> exponential decay;
        #                                                == -1 -> read from
        #                                                I_file_dummy.xlsx         -
        # IOP_TOT*              I      scalar float      total operation
        #                                                current                   A
        # I0_OP_TOT*             I      scalar float      total operation current
        #                                                @ time = 0                A
        # IALPHAB§              I      scalar integer    flag to define the
        #                                                magnetic field gradient
        #                                                along the strand:
        #                                                == 0 -> no gradient;
        #                                                == -1 -> read from
        #                                                alphab.xlsx               -
        # BASE_PATH*             I      string           path of folder Description
        #                                               of Components              -
        # External_current_path* I     string           name of the external
        #                                               file from which
        #                                               interpolate magnetic
        #                                               field gradient value       -
        # alphaB§               O      np array float   magnetic field
        #                                               gradient                 T/m
        #############################################################################
        # Invoched functions/methods: Get_from_xlsx
        #
        ############################################################################
        # * zcoord, I0_OP_MODE, IOP_TOT, I0_OP_TOT, BASE_PATH and External_current_path
        # are given by conductor.zcoord, conductor.inputs["I0_OP_MODE"], conductor.IOP_TOT, conductor.inputs["I0_OP_TOT"],
        # conductor.BASE_PATH and conductor.file_input["EXTERNAL_CURRENT"].
        # § IALPHAB and alphaB are component attributes: self.IALPHAB, self.dict_node_pt["alpha_B"].
        # N.B. alphaB is a Strands attribute so its value can be assigned
        # directly, it has the same of shape of zcoord and it is a np array.
        ############################################################################
        #
        # Translated and optimized by D. Placido Polito 06/2020
        #
        ############################################################################
        """

        if nodal:
            # compute alpha_B in each node (cdp, 07/2020)
            if self.operations["IALPHAB"] <= -1:  # read from file
                if conductor.cond_time[-1] == 0:
                    # Build file path.
                    file_path = os.path.join(
                        conductor.BASE_PATH, conductor.file_input["EXTERNAL_ALPHAB"]
                    )
                    # Load auxiliary input file.
                    alphab_df, self.flagSpecfield_alpha_b = load_auxiliary_files(
                        file_path, sheetname=self.identifier
                    )
                    # Build interpolator and get the interpolaion flag (space_only,time_only or space_and_time).
                    (
                        self.alphab_interpolator,
                        self.alphab_interp_flag,
                    ) = build_interpolator(
                        alphab_df, self.operations["ALPHAB_INTERPOLATION"]
                    )

                # call load_user_defined_quantity on the component.
                self.dict_node_pt["op_current"] = do_interpolation(
                    self.alphab_interpolator,
                    conductor.grid_features["zcoord"],
                    conductor.cond_time[-1],
                    self.alphab_interp_flag,
                )

                # leggi un file come del campo magnetico
                # controlla se e' per unita' di corrente
                # in caso affermatico moltiplica per IOP_TOT
                if self.flagSpecfield_alpha_b == 2:  # alphaB is per unit of current
                    self.dict_node_pt["alpha_B"] = (
                        self.dict_node_pt["alpha_B"] * conductor.inputs["I0_OP_TOT"]
                    )
                if conductor.inputs["I0_OP_MODE"] < 0:
                    self.dict_node_pt["alpha_B"] = (
                        self.dict_node_pt["alpha_B"]
                        * conductor.inputs["I0_OP_TOT"]
                        / conductor.inputs["I0_OP_TOT"]
                    )
            elif self.operations["IALPHAB"] == 0:
                self.dict_node_pt["alpha_B"] = np.zeros(
                    conductor.grid_features["N_nod"]
                )
        elif nodal == False:
            # compute alpha_B in each Gauss point (cdp, 07/2020)
            self.dict_Gauss_pt["alpha_B"] = (
                np.abs(
                    self.dict_node_pt["alpha_B"][
                        0 : conductor.grid_features["N_nod"] - 1
                    ]
                    + self.dict_node_pt["alpha_B"][
                        1 : conductor.grid_features["N_nod"] + 1
                    ]
                )
                / 2.0
            )

    def get_superconductor_critical_prop(self, conductor, nodal=True):

        """
        Calcola i margini
        usa temperatura strand
        Usa campo magnetico BFIELD
        usa EPSI
        usa alphaB
        :return: Jcritck Tcritchk TCSHRE TCSHREmin
        """

        # Properties evaluation in each nodal point (cdp, 07/2020)
        # Variable where is necessary to correctly evaluate EPSILON calling method \
        # Get_EPS (cdp, 07/2020)
        if nodal:
            self.dict_node_pt = self.eval_critical_properties(self.dict_node_pt)
        # Properties evaluation in each Gauss point (cdp, 07/2020)
        elif nodal == False:
            self.dict_Gauss_pt = self.eval_critical_properties(self.dict_Gauss_pt)

    # End of Method get_superconductor_critical_prop

    def eval_critical_properties(self, dict_dummy):

        if self.inputs["superconducting_material"] == "NbTi":
            dict_dummy["T_critical"] = critical_temperature_nbti(
                dict_dummy["B_field"], self.inputs["Bc20m"], self.inputs["Tc0m"]
            )
            dict_dummy["J_critical"] = critical_current_density_nbti(
                dict_dummy["temperature"],
                dict_dummy["B_field"],
                self.inputs["Bc20m"],
                self.inputs["c0"],
                self.inputs["Tc0m"],
            )
        elif self.inputs["superconducting_material"] == "Nb3Sn":
            dict_dummy["T_critical"] = critical_temperature_nb3sn(
                dict_dummy["B_field"],
                dict_dummy["Epsilon"],
                self.inputs["Tc0m"],
                self.inputs["Bc20m"],
            )
            dict_dummy["J_critical"] = critical_current_density_nb3sn(
                dict_dummy["temperature"],
                dict_dummy["B_field"],
                dict_dummy["Epsilon"],
                self.inputs["Tc0m"],
                self.inputs["Bc20m"],
                self.inputs["c0"],
            )
        elif self.inputs["superconducting_material"] == "YBCO":
            dict_dummy["T_critical"] = self.inputs["Tc0m"] * np.ones(
                dict_dummy["temperature"].shape
            )
            dict_dummy["J_critical"] = critical_current_density_re123(
                dict_dummy["temperature"],
                dict_dummy["B_field"],
                self.inputs["Tc0m"],
                self.inputs["Bc20m"],
                self.inputs["c0"],
            )
        elif self.inputs["superconducting_material"] == "scaling.dat":
            # Get user defined scaling invoking method User_scaling_margin \
            # (cdp, 10/2020)
            self.user_scaling_margin()

        return dict_dummy

    # end method Eval_critical_properties (cdp, 10/2020)

    def get_tcs(self, nodal=True):
        """Method that allows the evaluation of the current sharing temperature in nodal points or in Gauss points, according to flag nodal.

        Args:
            nodal (bool, optional): Flag to evaluate the current sharing temperature in proper location. If True evaluation is on the nodal points, if False evaluation is on the Gauss points. Defaults to True.
        """

        # Current sharing temperature evaluation in each nodal point
        if nodal:
            self.dict_node_pt = self.eval_tcs(self.dict_node_pt)
        # Current sharing temperature evaluation in each Gauss point
        elif nodal == False:
            self.dict_Gauss_pt = self.eval_tcs(self.dict_Gauss_pt)

    # End method get_tcs

    def eval_tcs(self, dict_dummy):

        # Evaluate the current density with a suitable definition of the cross 
        # section according to the definition of the scaling parameter c0. 
        # A_current_density takes already into account of costheta, i.e., is 
        # sloped.
        jop = (
            np.abs(self.dict_node_pt["op_current"])
            / (self.cross_section_current_density)
        )

        bmax = dict_dummy["B_field"] * (1 + dict_dummy["alpha_B"])
        if self.inputs["superconducting_material"] == "NbTi":
            dict_dummy["T_cur_sharing"] = current_sharing_temperature_nbti(
                dict_dummy["B_field"],
                jop,
                self.inputs["Bc20m"],
                self.inputs["c0"],
                self.inputs["Tc0m"],
            )
            dict_dummy["T_cur_sharing_min"] = current_sharing_temperature_nbti(
                bmax,
                jop,
                self.inputs["Bc20m"],
                self.inputs["c0"],
                self.inputs["Tc0m"],
            )
        elif self.inputs["superconducting_material"] == "Nb3Sn":
            dict_dummy["T_cur_sharing"] = current_sharing_temperature_nb3sn(
                dict_dummy["B_field"],
                dict_dummy["Epsilon"],
                jop,
                self.inputs["Tc0m"],
                self.inputs["Bc20m"],
                self.inputs["c0"],
            )
            dict_dummy["T_cur_sharing_min"] = current_sharing_temperature_nb3sn(
                bmax,
                dict_dummy["Epsilon"],
                jop,
                self.inputs["Tc0m"],
                self.inputs["Bc20m"],
                self.inputs["c0"],
            )
        elif self.inputs["superconducting_material"] == "YBCO":
            dict_dummy["T_cur_sharing"] = current_sharing_temperature_re123(
                dict_dummy["B_field"],
                jop,
                self.inputs["Tc0m"],
                self.inputs["Bc20m"],
                self.inputs["c0"],
            )
            dict_dummy["T_cur_sharing_min"] = current_sharing_temperature_re123(
                bmax,
                jop,
                self.inputs["Tc0m"],
                self.inputs["Bc20m"],
                self.inputs["c0"],
            )
        elif self.inputs["superconducting_material"] == "scaling.dat":

            warnings.warn("Still to be understood what to do here!!")

        return dict_dummy

    # End method eval_tcs

    def get_eps(self, conductor, nodal=True):
        # For each strand of type StrandMixedComponent or StackComponent (cdp, 06/2020)
        if nodal:
            # compute Epsilon in each node (cdp, 07/2020)
            if self.operations["IEPS"] < 0:  # strain from file strain.dat

                if conductor.cond_time[-1] == 0:
                    # Build file path.
                    file_path = os.path.join(
                        conductor.BASE_PATH, conductor.file_input["EXTERNAL_STRAIN"]
                    )
                    # Load auxiliary input file.
                    eps_df, self.flagSpecfield_eps = load_auxiliary_files(
                        file_path, sheetname=self.identifier
                    )
                    # Build interpolator and get the interpolaion flag (space_only,time_only or space_and_time).
                    self.eps_interpolator, self.eps_interp_flag = build_interpolator(
                        eps_df, self.operations["IOP_INTERPOLATION"]
                    )

                # call load_user_defined_quantity on the component.
                self.dict_node_pt["Epsilon"] = do_interpolation(
                    self.eps_interpolator,
                    conductor.grid_features["zcoord"],
                    conductor.cond_time[-1],
                    self.eps_interp_flag,
                )

                if self.flagSpecfield_eps == 1:
                   # Add also a logger
                    warnings.warn("Still to be decided what to do here\n")
            elif self.operations["IEPS"] == 0:  # no strain (cdp, 06/2020)
                self.dict_node_pt["Epsilon"] = np.zeros(
                    conductor.grid_features["N_nod"]
                )
            elif self.operations["IEPS"] == 1:
                # constant strain to the value in input file \
                # conductor_i_operation.xlsx (cdp, 06/2020)
                self.dict_node_pt["Epsilon"] = self.operations["EPS"] * np.ones(
                    conductor.grid_features["N_nod"]
                )
        elif nodal == False:
            # compute Epsilon in each Gauss point (cdp, 07/2020)
            self.dict_Gauss_pt["Epsilon"] = (
                self.dict_node_pt["Epsilon"][0 : conductor.grid_features["N_nod"] - 1]
                + self.dict_node_pt["Epsilon"][1 : conductor.grid_features["N_nod"] + 1]
            ) / 2.0

    # end method Get_EPS (cdp, 10/2020)

    def user_scaling_margin(self):

        """
        Method that read file scaling_input.dat to get the user scaling for SuperConductors strands and convert it into an attribute dictionary (cdp, 10/2020)
        """

        # declare dictionary (cdp, 10/2020)
        # construct list of integer values (cdp, 10/2020)
        list_integer = ["emode", "ieavside"]
        # Loop to read file scaling_input.dat by lines and construct a \
        # dictionary (cdp, 10/2020)
        with open("scaling_input.dat", "r") as scaling:
            # Read lines (cdp, 10/2020)
            for line in scaling:
                if line[0] == "#" or line[0] == "\n":
                    # escape comments and void lines (cdp, 10/2020)
                    pass
                else:
                    # split the list into two fields, fields[0] is dictionary key,
                    # fields[1] is the corresponding value. The final \n is ignored \
                    # considering only the first [:-1] characters in string fields[1] \
                    # (cdp, 10/2020)
                    fields = line.split(" = ")
                    if fields[0] in list_integer:
                        # convert to integer (cdp, 10/2020)
                        self.dict_scaling_input[fields[0]] = int(fields[1][:-1])
                    elif fields[0] == "superconducting_material":
                        # flag to
                        self.dict_scaling_input[fields[0]] = str(fields[1][:-1])
                    else:
                        # convert to float (cdp, 10/2020)
                        self.dict_scaling_input[fields[0]] = float(fields[1][:-1])
                    # end if fields[0] (cdp, 10/2020)
                # end if line[0] (cdp, 10/2020)
            # end for line (cdp, 10/2020)
        # end with (cdp, 10/2020)

    # end method User_scaling_margin (cdp, 10/2020)

    def electric_resistance(
        self, conductor: object, electrical_resistivity_key: str, ind: np.ndarray
    ) -> np.ndarray:
        f"""Method that evaluate electric resistance of a single material, such as superconductor, stabilizer, stainless steel and other electric conducting material.

        Args:
            conductor (object): class Conductor object in which distance between consecutive nodes is stored to do the calculation.
            electrical_resistivity_key (str): dictionary key for the electrical resistivity of the material.
            ind (np.ndarray): array with the index of the location in wich electric resistance should be evaluated with this method.

        Returns:
            np.ndarray: array of the electric resistance in Ohm of shape {ind.shape = }. The maximum lenght of the outcome is {conductor.grid_input["NELEMS"] = }.
        """
        return (
            self.dict_Gauss_pt[electrical_resistivity_key][ind]
            * conductor.node_distance[("StrandComponent", self.identifier)][ind]
            / self.inputs["CROSSECTION"]
        )

    def parallel_electric_resistance(
        self, conductor: object, electrical_resistivity_keys: list, ind: np.ndarray
    ) -> np.ndarray:
        f"""Method that evaluate electric resistance in the case of a parallel of two electric conducting materials, as is the case for StackComponent and StrandMixedComponent in current sharing regime.

        Args:
            conductor (object): class Conductor object in which distance between consecutive nodes is stored to do the calculation.
            electrical_resistivity_key (list): list of dictionary key for the electrical resistivity of the materials (typical values for the application of this software are electrical_resistivity_superconductor and electrical_resistivity_stabilizer).
            ind (np.ndarray): array with the index of the location in wich electric resistance should be evaluated with this method.

        Raises:
            ValueError: if list electrical_resistivity_keys does not have exactly two items.
            ValueError: if items in list electrical_resistivity_keys are not of type string.

        Returns:
            np.ndarray: array of the electric resistance in Ohm of shape {ind.shape = }. The maximum lenght of the outcome is {conductor.grid_input["NELEMS"] = }.
        """
        if len(electrical_resistivity_keys) != 2:
            # Check list lenght.
            raise ValueError(
                f"List electrical_resistivity_keys must have 2 items; {len(electrical_resistivity_keys) = }.\n"
            )

        if not all(isinstance(item, str) for item in electrical_resistivity_keys):
            # Check that all items in list are string.
            raise ValueError(
                f"All items in list electrical_resistivity_keys must be of type string. {electrical_resistivity_keys = }.\n"
            )

        # Electric resistance matrix initialization.
        electric_resistances = np.zeros((ind.size, 2))
        # Evaluate electri resistances
        for ii, item in enumerate(electrical_resistivity_keys):
            electric_resistances[:, ii] = self.electric_resistance(conductor, item, ind)

        # Evaluate parallel electric resistance:
        # R_eq = R1*R2/(R1+R2)
        return (
            electric_resistances[:, 0]
            * electric_resistances[:, 1]
            / (electric_resistances.sum(axis=1))
        )


    def __manage_fixed_potental(self, length:float):
        """Method that deals with fixed potentials: converts fixed potential values to array if they are integers or strings and checks the coordinate where fixed potentials are assigned.

        Args:
            length (float): conductor length.
        """
        self.__convert_fixed_potential_to_array()
        self.__checks_fix_potential_coordinate(length)

    def __convert_fixed_potential_to_array(self):
        """Private method that allows to convert keys FIX_POTENTIAL_COORDINATE and FIX_POTENTIAL_VALUE to numpy array.
        """
        # Define dictionary with private methods.
        private_metods = {int: self.__int_or_float_to_array, float: self.__int_or_float_to_array, str: self.__str_to_array}

        # Call private method self.__int_or_float_to_array if 
        # type(self.operations[key]) is int (integer) or float; 
        # self.__str_to_array it type(self.operations[key]) is str (string).
        for key in ["FIX_POTENTIAL_COORDINATE", "FIX_POTENTIAL_VALUE"]:
            private_metods[type(self.operations[key])](key)

    def __int_or_float_to_array(self, key:str):
        """Private method that converts value corresponding to keys FIX_POTENTIAL_COORDINATE and FIX_POTENTIAL_VALUE of dictionary self.operations from integer or float to ndarray of float.

        Args:
            key (str): self.operations key, can be FIX_POTENTIAL_COORDINATE or FIX_POTENTIAL_VALUE.
        """
        # Convert value corresponding to key from int to ndarray of float.
        self.operations[key] = np.array([self.operations[key]], dtype=float)

    def __str_to_array(self, key:str):
        """Private method that converts value corresponding to keys FIX_POTENTIAL_COORDINATE and FIX_POTENTIAL_VALUE of dictionary self.operations from str to ndarray of float.

        Args:
            key (str): self.operations key, can be FIX_POTENTIAL_COORDINATE or FIX_POTENTIAL_VALUE.
        """
        # Convert value corresponding to key from str to ndarray of float.
        self.operations[key] = np.array(self.operations[key].split(","), dtype=float)
    
    def __checks_fix_potential_coordinate(self, length:float):
        """Private method that checks consistency between input values FIX_POTENTIAL_NUMBER, FIX_POTENTIAL_COORDINATE and FIX_POTENTIAL_VALUE after conversion of FIX_POTENTIAL_COORDINATE and FIX_POTENTIAL_VALUE to ndarray.

        Args:
            length (float): length of the conductor.

        Raises:
            ValueError: if length of ndarray FIX_POTENTIAL_COORDINATE or FIX_POTENTIAL_VALUE is not equal to FIX_POTENTIAL_NUMBER.
            ValueError: if maximum value in FIX_POTENTIAL_COORDINATE is larger than the length of the conductor.
        """

        # Build dictionary with error messages.
        message = dict(
            FIX_POTENTIAL_COORDINATE=f"The number of the coordinates of fixed potential must be equal to the number of declared fixed potential:\n{self.operations['FIX_POTENTIAL_COORDINATE']=};\n{self.operations['FIX_POTENTIAL_NUMBER']=}\n",
            FIX_POTENTIAL_VALUE=f"The number of the values of fixed potential must be equal to the number of declared fixed potential:\n{self.operations['FIX_POTENTIAL_VALUE']=};\n{self.operations['FIX_POTENTIAL_NUMBER']=}\n",
        )

        # Check ndarray length consistency.
        for key in message.keys():
            if len(self.operations[key]) != self.operations["FIX_POTENTIAL_NUMBER"]:
                raise ValueError(message[key])

        # Check FIX_POTENTIAL_COORDINATE consistency with conductor length.
        if np.max(self.operations["FIX_POTENTIAL_COORDINATE"]) > length:
            raise ValueError(
                f"Fixed potential coordinate cannot exceed conductor length:\n{np.max(self.operations['FIX_POTENTIAL_COORDINATE'])=}m;\n{length=}m"
            )

    def __delete_fixed_potential_inputs(self, _:float):
        """Method that forces self.operations["FIX_POTENTIAL_NUMBER"] to be zero and deletes no longer useful keys FIX_POTENTIAL_COORDINATE and FIX_POTENTIAL_VALUE from dictionary self.operations.

        Args:
            _ (float): any float varyables, not used.
        """

        # Force self.operations["FIX_POTENTIAL_NUMBER"] to 0.
        self.operations["FIX_POTENTIAL_NUMBER"] = 0
        # Delete no longer useful keys from dictionary self.operations.
        for key in ["FIX_POTENTIAL_COORDINATE", "FIX_POTENTIAL_VALUE"]:
            del self.operations[key]

    def deal_with_fixed_potential(self, length:float):
        """Method that deals with fixed potential input. It converts True value to 1 (odd beaviour) if necessary; defines dictionary with private methods self.__manage_fixed_potental and self.__delete_fixed_potential_inputs, and calls them according to input value FIX_POTENTIAL_FLAG.

        Args:
            length (float): conductor length.
        """

        # Temporary solution to mangage input file loading, strange behavior: 1
        # are converted to True but 0 not converted to False.
        for key in [
            "FIX_POTENTIAL_NUMBER",
            "FIX_POTENTIAL_COORDINATE",
            "FIX_POTENTIAL_VALUE",
        ]:
            if self.operations[key] == True:
                self.operations[key] = 1
        
        # Define dictionary with methods inherited from class Strand used to 
        # deal with fixed potential.
        methods = {
            True: self.__manage_fixed_potental,
            False: self.__delete_fixed_potential_inputs,
        }

        # Calls method self.manage_equipotential_surfaces_index if
        # FIX_POTENTIAL_FLAG is True; self.delete_fixed_potential_inputs if 
        # FIX_POTENTIAL_FLAG is False.
        methods[self.operations["FIX_POTENTIAL_FLAG"]](length)