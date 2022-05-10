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
        # xcoord*               I      np array float    conductor spatial
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
        # * xcoord, I0_OP_MODE, IOP_TOT, I0_OP_TOT, BASE_PATH and External_current_path
        # are given by conductor.xcoord, conductor.inputs["I0_OP_MODE"], conductor.IOP_TOT, conductor.inputs["I0_OP_TOT"],
        # conductor.BASE_PATH and conductor.file_input["EXTERNAL_CURRENT"].
        # § IALPHAB and alphaB are component attributes: self.IALPHAB, self.dict_node_pt["alpha_B"].
        # N.B. alphaB is a Strands attribute so its value can be assigned
        # directly, it has the same of shape of xcoord and it is a np array.
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
                    alphab_df, flagSpecfield = load_auxiliary_files(
                        file_path, sheetname=self.ID
                    )
                    # Build interpolator and get the interpolaion flag (space_only,time_only or space_and_time).
                    (
                        self.alphab_interpolator,
                        self.alphab_interp_flag,
                    ) = build_interpolator(
                        alphab_df, self.operations["ALPHAB_INTERPOLATION"]
                    )

                # call load_user_defined_quantity on the component.
                self.dict_node_pt["IOP"] = do_interpolation(
                    self.alphab_interpolator,
                    conductor.dict_discretization["xcoord"],
                    conductor.cond_time[-1],
                    self.alphab_interp_flag,
                )

                # leggi un file come del campo magnetico
                # controlla se e' per unita' di corrente
                # in caso affermatico moltiplica per IOP_TOT
                if flagSpecfield == 2:  # alphaB is per unit of current
                    self.dict_node_pt["alpha_B"] = (
                        self.dict_node_pt["alpha_B"] * conductor.IOP_TOT
                    )
                if conductor.inputs["I0_OP_MODE"] < 0:
                    self.dict_node_pt["alpha_B"] = (
                        self.dict_node_pt["alpha_B"]
                        * conductor.IOP_TOT
                        / conductor.inputs["I0_OP_TOT"]
                    )
            elif self.operations["IALPHAB"] == 0:
                self.dict_node_pt["alpha_B"] = np.zeros(
                    conductor.dict_discretization["N_nod"]
                )
        elif nodal == False:
            # compute alpha_B in each Gauss point (cdp, 07/2020)
            self.dict_Gauss_pt["alpha_B"] = (
                np.abs(
                    self.dict_node_pt["alpha_B"][
                        0 : conductor.dict_discretization["N_nod"] - 1
                    ]
                    + self.dict_node_pt["alpha_B"][
                        1 : conductor.dict_discretization["N_nod"] + 1
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

        if self.inputs["ISUPERCONDUCTOR"] == "NbTi":
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
        elif self.inputs["ISUPERCONDUCTOR"] == "Nb3Sn":
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
        elif self.inputs["ISUPERCONDUCTOR"] == "HTS":
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
        elif self.inputs["ISUPERCONDUCTOR"] == "scaling.dat":
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

        jop = (
            np.abs(self.dict_node_pt["IOP"][0])
            / (self.ASC / self.inputs["COSTETA"])
            * np.ones(dict_dummy["B_field"].shape)
        )

        bmax = dict_dummy["B_field"] * (1 + dict_dummy["alpha_B"])
        if self.inputs["ISUPERCONDUCTOR"] == "NbTi":
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
        elif self.inputs["ISUPERCONDUCTOR"] == "Nb3Sn":
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
        elif self.inputs["ISUPERCONDUCTOR"] == "HTS":
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
        elif self.inputs["ISUPERCONDUCTOR"] == "scaling.dat":

            warnings.warn("Still to be understood what to do here!!")

        return dict_dummy

    # End method eval_tcs

    def get_eps(self, conductor, nodal=True):
        # For each strand of type StrandMixedComponent or StrandSuperconductorComponent (cdp, 06/2020)
        if nodal:
            # compute Epsilon in each node (cdp, 07/2020)
            if self.operations["IEPS"] < 0:  # strain from file strain.dat

                if conductor.cond_time[-1] == 0:
                    # Build file path.
                    file_path = os.path.join(
                        conductor.BASE_PATH, conductor.file_input["EXTERNAL_STRAIN"]
                    )
                    # Load auxiliary input file.
                    eps_df, flagSpecfield = load_auxiliary_files(
                        file_path, sheetname=self.ID
                    )
                    # Build interpolator and get the interpolaion flag (space_only,time_only or space_and_time).
                    self.eps_interpolator, self.eps_interp_flag = build_interpolator(
                        eps_df, self.operations["IOP_INTERPOLATION"]
                    )

                # call load_user_defined_quantity on the component.
                self.dict_node_pt["Epsilon"] = do_interpolation(
                    self.eps_interpolator,
                    conductor.dict_discretization["xcoord"],
                    conductor.cond_time[-1],
                    self.eps_interp_flag,
                )

                if flagSpecfield == 1:
                    print("still to be decided what to do here\n")
            elif self.operations["IEPS"] == 0:  # no strain (cdp, 06/2020)
                self.dict_node_pt["Epsilon"] = np.zeros(
                    conductor.dict_discretization["N_nod"]
                )
            elif self.operations["IEPS"] == 1:
                # constant strain to the value in input file \
                # conductor_i_operation.xlsx (cdp, 06/2020)
                self.dict_node_pt["Epsilon"] = self.operations["EPS"] * np.ones(
                    conductor.dict_discretization["N_nod"]
                )
        elif nodal == False:
            # compute Epsilon in each Gauss point (cdp, 07/2020)
            self.dict_Gauss_pt["Epsilon"] = (
                self.dict_node_pt["Epsilon"][
                    0 : conductor.dict_discretization["N_nod"] - 1
                ]
                + self.dict_node_pt["Epsilon"][
                    1 : conductor.dict_discretization["N_nod"] + 1
                ]
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
                    elif fields[0] == "ISUPERCONDUCTOR":
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
