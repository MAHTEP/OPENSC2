import numpy as np
import os
import warnings

from utility_functions.auxiliary_functions import (
    load_auxiliary_files,
    build_interpolator,
    do_interpolation,
)

from conductor_flags import (
    IOP_NOT_DEFINED,
    IOP_CONSTANT,
    IOP_FROM_FILE,
    IOP_FROM_EXT_FUNCTION,
)

class SolidComponent:
    def __init__(self, simulation, s_comp):

        """
        Constructor method of class SolidComponent (cdp, 11/2020)
        """

        # Questa è una bozza, quando e se si dovranno considerare altri flag come \
        # IBFUN o IALPHAB, valutare se è il caso di sfruttare un metodo per \
        # evitare di scrivere lo stesso codice più volte (cdp, 11/2020)
        if s_comp.operations["IQFUN"] > 0:
            # External heating parameters given in file operation.xlsx \
            # (cdp, 11/2020)
            # The heating will be on at some times (cdp, 01/2021)
            s_comp.flag_heating = "On"
            if simulation.transient_input["IADAPTIME"] == 0:
                # Time adaptivity off and (cdp, 11/2020)
                s_comp.dict_num_step["IQFUN"] = dict(
                    ON=int(
                        s_comp.operations["TQBEG"]
                        / simulation.transient_input["STPMIN"]
                    ),
                    OFF=int(
                        s_comp.operations["TQEND"]
                        / simulation.transient_input["STPMIN"]
                    ),
                )
            else:
                # Time adaptivity on
                print("Still to be decided what to do there (cdp, 11/2020)\n")
        # End s_comp.operations["IQFUN"].

    # end method __init__ (cdp, 11/2020)

    def eval_sol_comp_properties(self, inventory, nodal=True):

        """
        Method that evaluate total_density, specific_heat and thermal conductivity of SolidComponent class objects in both nodal points and Gauss points according to **options input parameter (cdp, 07/2020)
        """

        # Properties evaluation in each nodal point (cdp, 07/2020)
        if nodal:
            dict_dummy = self.dict_node_pt
            self.dict_node_pt = self.eval_properties(dict_dummy, inventory)
        # Properties evaluation in each Gauss point (cdp, 07/2020)
        elif nodal == False:
            dict_dummy = self.dict_Gauss_pt
            self.dict_Gauss_pt = self.eval_properties(dict_dummy, inventory)

    def eval_properties(self, dict_dummy: dict, inventory: dict) -> dict:
        """Method that actually evaluate total_density, specific_heat and thermal conductivity of SolidComponent class objects regardless of the location (nodal or Gauss points).

        Args:
            dict_dummy (dict): dictionary with material properties in nodal points or Gauss points according to the value of flag nodal in method eval_sol_comp_properties of class SolidComponent.

        Returns:
            dict: dictionary with updated material properties in nodal points or Gauss points according to the value of flag nodal in method eval_sol_comp_properties of class SolidComponent.
        """

        if (
            self.name == inventory["StrandMixedComponent"].name
            or self.name == inventory["StrandStabilizerComponent"].name
        ):
            dict_dummy.update(total_density=self.strand_density(dict_dummy))
            dict_dummy.update(
                total_isobaric_specific_heat=self.strand_isobaric_specific_heat(
                    dict_dummy
                )
            )
            dict_dummy.update(
                total_thermal_conductivity=self.strand_thermal_conductivity(dict_dummy)
            )
            if self.name == inventory["StrandMixedComponent"].name:
                dict_dummy.update(
                    electrical_resistivity_stabilizer=self.strand_electrical_resistivity_not_sc(
                        dict_dummy
                    )
                )
            elif self.name == inventory["StrandStabilizerComponent"].name:
                dict_dummy.update(
                    electrical_resistivity_stabilizer=self.strand_electrical_resistivity(
                        dict_dummy
                    )
                )
        elif self.name == inventory["StackComponent"].name:
            dict_dummy.update(total_density=self.stack_density(dict_dummy))
            dict_dummy.update(
                total_isobaric_specific_heat=self.stack_isobaric_specific_heat(
                    dict_dummy
                )
            )
            dict_dummy.update(
                total_thermal_conductivity=self.stack_thermal_conductivity(dict_dummy)
            )
            dict_dummy.update(
                electrical_resistivity_stabilizer=self.stack_electrical_resistivity_not_sc(
                    dict_dummy
                )
            )
        elif inventory["JacketComponent"].name:
            dict_dummy.update(total_density=self.jacket_density(dict_dummy))
            dict_dummy.update(
                total_isobaric_specific_heat=self.jacket_isobaric_specific_heat(
                    dict_dummy
                )
            )
            dict_dummy.update(
                total_thermal_conductivity=self.jacket_thermal_conductivity(dict_dummy)
            )
            dict_dummy.update(
                total_electrical_resistivity=self.jacket_electrical_resistivity(
                    dict_dummy
                )
            )

        return dict_dummy

    def get_current_fractions(
        self, total_sc_area: float, total_so_area: float, inventory: dict
    ):
        """Method that evaluates: 1) fraction of the total current that flows in superconductor cross section of each strand or stack object if in superconducting regime (op_current_fraction_sc); 2) fraction of the total current that flows in the total cross section (superconductor and not superconductor or stabilizer materials) of each strand or stack object if in current sharing regime (op_current_fraction_sh). Both are methods of the generic object.

        Note: for the time being both fractions are set to 0.0 for objects of kind JacketComponent.

        Args:
            total_sc_area (float): total superconductor cross section of the conductor.
            total_so_area (float): total cross section of strands and or stacks of the conductor.
            inventory (dict): dictionary with all the conductor components.
        """

        # Fraction of the total current that goes in the superconductor cross
        # section in superconducting regime.
        if (
            self.name == inventory["StrandMixedComponent"].name
            or self.name == inventory["StackComponent"].name
        ):
            self.op_current_fraction_sc = self.cross_section["sc"] / total_sc_area
        elif (
            self.name == inventory["JacketComponent"].name
            or self.name == inventory["StrandStabilizerComponent"].name
        ):
            self.op_current_fraction_sc = 0.0

        # Fraction of the total current that goes in the strand or thape cross
        # section in current sharing regime.
        if (
            self.name == inventory["StackComponent"].name
            or self.name == inventory["StrandMixedComponent"].name
            or self.name == inventory["StrandStabilizerComponent"].name
        ):
            self.op_current_fraction_sh = self.inputs["CROSSECTION"] / total_so_area
        else:  # jacket object.
            self.op_current_fraction_sh = 0.0

    def __check_current_mode(self, conductor: object):
        """Private method that checks consistency between flags conductor.inputs['I0_OP_MODE'] and self.operations['IOP_MODE'] to deal with current definition.

        Args:
            conductor (object): ConductorComponent object with all informations to make the check.

        Raises:
            ValueError: self.operations["IOP_MODE"] != None and conductor.inputs["I0_OP_MODE"] == -1 and self.operations["IOP_MODE"] != -1.
            ValueError: self.operations["IOP_MODE"] != None and conductor.inputs["I0_OP_MODE"] == 0 and self.operations["IOP_MODE"] != 0.
        """

        # Initialize dictionary with error message to be printed.
        message_switch = {
            -1: f"{conductor.inputs['I0_OP_MODE']=} implies that current carried by object {self.identifier = } should be read from file. Flag self.operations['IOP_MODE'] should be set to -1; current value is {self.operations['IOP_MODE']=}. Please check sheet {self.identifier} of input file conductor_operation.xlsx.\n",
            0: f"{conductor.inputs['I0_OP_MODE']=} implies that current carried by object {self.identifier = } is evaluated from the code since the total current carried by the conductor is assigned. Flag self.operations['IOP_MODE'] should be set to 0; current value is {self.operations['IOP_MODE']=}. Please check sheet {self.identifier} of input file conductor_operation.xlsx.\n",
        }

        # Check consistency between flags conductor.inputs['I0_OP_MODE'] and
        # self.operations['IOP_MODE'].
        if self.operations["IOP_MODE"] != None:
            if (
                conductor.inputs["I0_OP_MODE"] == -1
                and self.operations["IOP_MODE"] != -1
            ):
                raise ValueError(message_switch[conductor.inputs["I0_OP_MODE"]])
            elif (
                conductor.inputs["I0_OP_MODE"] == 0 and self.operations["IOP_MODE"] != 0
            ):
                raise ValueError(message_switch[conductor.inputs["I0_OP_MODE"]])

    def get_current(self, conductor: object):
        """Method that evaluates the current source therm at each thermal time step according to flag I0_OP_MODE set as input by the user.

        Args:
            conductor (object): ConductorComponent object with all informations to make the calculation.

        Raises:
            ValueError: if I0_OP_MODE is -2 (still to be implemented).
            ValueError: if a not valid value is given to flag I0_OP_MODE.
        """

        # Check consistency between flags conductor.inputs['I0_OP_MODE']
        # and self.operations['IOP_MODE'] only the first time.
        if conductor.cond_time[-1] == 0:
            self.__check_current_mode(conductor)

        # Get current.
        if self.operations["IOP_MODE"] != None:
            # The object carryes a current and its value is defied as below.
            if conductor.inputs["I0_OP_MODE"] == IOP_FROM_FILE:

                if conductor.cond_time[-1] == 0:
                    # Build file path.
                    file_path = os.path.join(
                        conductor.BASE_PATH, conductor.file_input["EXTERNAL_CURRENT"]
                    )
                    # Load auxiliary input file.
                    current_df, self.flagSpecfield_current = load_auxiliary_files(
                        file_path, sheetname=self.identifier
                    )
                    # Build interpolator and get the interpolaion flag (space_only,time_only or space_and_time).
                    (
                        self.current_interpolator,
                        self.current_interp_flag,
                    ) = build_interpolator(
                        current_df, self.operations["IOP_INTERPOLATION"]
                    )

                # Evaluate current of generic solid component object by
                # interpolation.
                self.dict_node_pt["op_current"] = do_interpolation(
                    self.current_interpolator,
                    conductor.grid_features["zcoord"],
                    conductor.electric_time,
                    self.current_interp_flag,
                )

                if self.current_interp_flag == "time_only":
                    # Convert to array
                    self.dict_node_pt["op_current"] = self.dict_node_pt[
                        "op_current"
                    ] * np.ones(conductor.grid_features["N_nod"])
                # Evaluate current in the Gauss nodal points.
                self.dict_Gauss_pt["op_current"] = (
                    self.dict_node_pt["op_current"][:-1]
                    + self.dict_node_pt["op_current"][1:]
                ) / 2.0
                # This is exploited in the electric resistance evaluation.
                if (
                    self.name == conductor.inventory["StackComponent"].name
                    or self.name == conductor.inventory["StrandMixedComponent"].name
                ):
                    # Build an alias for convenience when dealing with electric
                    # resistance evaluation.
                    self.dict_node_pt["op_current_sc"] = self.dict_node_pt["op_current"]
                    self.dict_Gauss_pt["op_current_sc"] = self.dict_Gauss_pt[
                        "op_current"
                    ]

                if self.flagSpecfield_current == 2:
                    # Add also a logger
                    warnings.warn("Still to be decided what to do here\n")
            elif conductor.inputs["I0_OP_MODE"] == IOP_CONSTANT:
                # Evaluate both attributes self.dict_node_pt["op_current"] and
                # self.dict_node_pt["op_current_sc"] for convenience in the evaluation of
                # electrical resistivity.
                self.dict_node_pt["op_current"] = (
                    conductor.inputs["I0_OP_TOT"]
                    * self.op_current_fraction_sh
                    * np.ones(conductor.grid_features["N_nod"])
                )
                self.dict_Gauss_pt["op_current"] = (
                    self.dict_node_pt["op_current"][:-1]
                    + self.dict_node_pt["op_current"][1:]
                ) / 2.0
                if (
                    self.name == conductor.inventory["StackComponent"].name
                    or self.name == conductor.inventory["StrandMixedComponent"].name
                ):
                    self.dict_node_pt["op_current_sc"] = (
                        conductor.inputs["I0_OP_TOT"]
                        * self.op_current_fraction_sc
                        * np.ones(conductor.grid_features["N_nod"])
                    )
                    self.dict_Gauss_pt["op_current_sc"] = (
                        self.dict_node_pt["op_current_sc"][:-1]
                        + self.dict_node_pt["op_current_sc"][1:]
                    ) / 2.0

            elif conductor.inputs["I0_OP_MODE"] == IOP_NOT_DEFINED:
                # User does not specify a current: set current carrient 
                # operating current to zero only the first time in both nodal 
                # and gauss points.
                if conductor.cond_num_step == 0:
                    self.dict_node_pt["op_current"] = np.zeros(conductor.grid_features["N_nod"])
                    self.dict_Gauss_pt["op_current"] = np.zeros(conductor.grid_input["NELEMS"])
            
            elif conductor.inputs["I0_OP_MODE"] == IOP_FROM_EXT_FUNCTION:
                raise ValueError(
                    f"Current from external function to be implemented!"
                )
            else:
                raise ValueError(
                    f"Not defined value for flag I0_OP_MODE: {conductor.inputs['I0_OP_MODE']=}.\n"
                )
        else:
            # The object does not carry a current; arrays are initialized to 0.
            # Initialize array op_current to 0 in dictionary dict_node_pt to
            # avoid error.
            self.dict_node_pt["op_current"] = np.zeros(conductor.grid_features["N_nod"])
            # Initialize array op_current to 0 in dictionary dict_Gauss_pt to
            # avoid error.
            self.dict_Gauss_pt["op_current"] = np.zeros(conductor.grid_input["NELEMS"])
            # This is exploited in the electric resistance evaluation.
            if (
                self.name == conductor.inventory["StackComponent"].name
                or self.name == conductor.inventory["StrandMixedComponent"].name
            ):
                # Build an alias for convenience when dealing with electric
                # resistance evaluation.
                self.dict_node_pt["op_current_sc"] = self.dict_node_pt["op_current"]
                self.dict_Gauss_pt["op_current_sc"] = self.dict_Gauss_pt["op_current"]

    # end Get_I

    def get_magnetic_field(self, conductor, nodal=True):
        if nodal:
            # compute B_field in each node (cdp, 07/2020)
            if self.operations["IBIFUN"] < 0:
                # cza to enable other negative \
                # (read from file) flags -> ibifun.eq.-3, see below (August 29, 2018)

                if conductor.cond_time[-1] == 0:
                    # Build file path.
                    file_path = os.path.join(
                        conductor.BASE_PATH, conductor.file_input["EXTERNAL_BFIELD"]
                    )
                    # Load auxiliary input file.
                    bfield_df, _ = load_auxiliary_files(
                        file_path, sheetname=self.identifier
                    )
                    # Build interpolator and get the interpolaion flag (space_only,time_only or space_and_time).
                    (
                        self.bfield_interpolator,
                        self.bfield_interp_flag,
                    ) = build_interpolator(
                        bfield_df, self.operations["B_INTERPOLATION"]
                    )

                # call load_user_defined_quantity on the component.
                self.dict_node_pt["B_field"] = do_interpolation(
                    self.bfield_interpolator,
                    conductor.grid_features["zcoord"],
                    conductor.electric_time,
                    self.bfield_interp_flag,
                )
                if self.operations["B_field_units"] == "T/A":
                    # BFIELD is per unit of current
                    self.dict_node_pt["B_field"] = (
                        self.dict_node_pt["B_field"] * conductor.inputs["I0_OP_TOT"]
                    )
                if (
                    conductor.inputs["I0_OP_MODE"] != 0
                    and conductor.inputs["I0_OP_TOT"] > 0
                ):
                    #### bfield e' un self e' un vettore
                    self.dict_node_pt["B_field"] = (
                        self.dict_node_pt["B_field"]
                        * conductor.inputs["I0_OP_TOT"]
                        / conductor.inputs["I0_OP_TOT"]
                    )
            elif self.operations["IBIFUN"] == 0:
                self.dict_node_pt["B_field"] = np.linspace(
                    self.operations["BISS"],
                    self.operations["BOSS"],
                    conductor.grid_features["N_nod"],
                )
            elif self.operations["IBIFUN"] == 1:
                self.dict_node_pt["B_field"] = np.linspace(
                    self.operations["BISS"],
                    self.operations["BOSS"],
                    conductor.grid_features["N_nod"],
                ) + conductor.inputs["I0_OP_TOT"] / conductor.inputs[
                    "I0_OP_TOT"
                ] * np.linspace(
                    self.operations["BITR"],
                    self.operations["BOTR"],
                    conductor.grid_features["N_nod"],
                )
        elif nodal == False:
            # compute B_field in each Gauss point (cdp, 07/2020)
            self.dict_Gauss_pt["B_field"] = (
                np.abs(
                    self.dict_node_pt["B_field"][:-1] + self.dict_node_pt["B_field"][1:]
                )
                / 2.0
            )

    # end Get_B_field

    # HERE STARTS THE DEFINITION OF MODULES USEFUL TO INITIALIZE THE DRIVERS FOR \
    # THE EXTERNAL HEATING. D. Placido (06/2020)

    def get_heat(self, conductor):

        """
        Method that evaluates the external heating according to the value of flag IQUFN, thaing unto account the chosen solution method (cdp, 11/2020)
        """

        # START INITIALIZATION (cdp, 10/2020)
        if conductor.cond_num_step == 0:
            if self.operations["IQFUN"] == 0:
                # Initialization is done always in the same way ragardless of the \
                # solution method: a column vector to exploit the array smart notation \
                # in Conductor class method Eval_Gauss_point. It is the only times at \
                # which this method is invoked (cdp, 11/2020)
                self.dict_node_pt["EXTFLX"] = np.zeros(
                    (conductor.grid_features["N_nod"], 1)
                )
            else:
                # Initialization is done always in the same way ragardless of the \
                # value of the flag IQFUN and coherently with the chosen solution \
                # algorithm (cdp, 10/2020)
                if (
                    conductor.inputs["METHOD"] == "BE"
                    or conductor.inputs["METHOD"] == "CN"
                ):
                    # Backward Euler or Crank-Nicolson (cdp, 10/2020)
                    self.dict_node_pt["EXTFLX"] = np.zeros(
                        (conductor.grid_features["N_nod"], 2)
                    )
                elif conductor.inputs["METHOD"] == "AM4":
                    # Adams-Moulton 4 (cdp, 10/2020)
                    self.dict_node_pt["EXTFLX"] = np.zeros(
                        (conductor.grid_features["N_nod"], 4)
                    )
            # end if self.operations["IQFUN"] (cdp, 11/2020)
        # end if conductor.cond_num_step (cdp, 10/2020)
        # END INITIALIZATION (cdp, 10/2020)
        # criteria to decide how to evaluate the external heating (cdp, 10/2020)
        if self.operations["IQFUN"] > 0:
            # invoke method Q0_where to evaluate the external heating according to \
            # the function corresponding to the value of flag IQFUN. It is not \
            # necessary to distinguish according to "Method" options at this point \
            # (cdp, 10/2020)
            if (
                conductor.cond_time[-1] >= self.operations["TQBEG"]
                and conductor.cond_time[-1] <= self.operations["TQEND"]
            ):
                self.heat0_where(conductor)
            elif (
                conductor.cond_time[-1] > self.operations["TQEND"]
                and self.flag_heating == "On"
            ):
                self.dict_node_pt["EXTFLX"][:, 0] = 0.0
                self.flag_heating = "Off"
            # end if (cdp, 10/2020)
        elif self.operations["IQFUN"] == -1:
            if conductor.cond_time[-1] == 0:
                # Build file path.
                file_path = os.path.join(
                    conductor.BASE_PATH, conductor.file_input["EXTERNAL_HEAT"]
                )
                # Load auxiliary input file.
                heat_df, _ = load_auxiliary_files(file_path, sheetname=self.identifier)
                # Build interpolator and get the interpolaion flag (space_only,time_only or space_and_time).
                self.heat_interpolator, self.heat_interp_flag = build_interpolator(
                    heat_df, self.operations["Q_INTERPOLATION"]
                )

                # compute external heating at conductor initialization calling function do_interpolation.

                self.dict_node_pt["EXTFLX"][:, 0] = do_interpolation(
                    self.heat_interpolator,
                    conductor.grid_features["zcoord"],
                    conductor.cond_time[-1],
                    self.heat_interp_flag,
                )
            elif conductor.cond_num_step > 0:
                if conductor.cond_num_step == 1:
                    # Store the old values only immediately after the initializzation, since after that the whole SYSLOD array is saved and there is no need to compute twice the same values.
                    self.dict_node_pt["EXTFLX"][:, 1] = self.dict_node_pt["EXTFLX"][
                        :, 0
                    ].copy()
                # end if conductor.cond_num_step (cdp, 10/2020)
                # call method load_user_defined_quantity to compute heat and overwrite the previous values.
                self.dict_node_pt["EXTFLX"][:, 0] = do_interpolation(
                    self.heat_interpolator,
                    conductor.grid_features["zcoord"],
                    conductor.cond_time[-1],
                    self.heat_interp_flag,
                )
            # end if conductor.cond_num_step (cdp, 10/2020)
        elif self.operations["IQFUN"] == -2:
            # AM2 part to be implemented
            if conductor.cond_time[-1] == 0:
                self.user_heat_function(conductor)
            elif conductor.cond_num_step > 0:
                if conductor.cond_num_step == 1:
                    # Store the old values only immediately after the initializzation, since after that the whole SYSLOD array is saved and there is no need to compute twice the same values.
                    self.dict_node_pt["EXTFLX"][:, 1] = self.dict_node_pt["EXTFLX"][
                        :, 0
                    ].copy()
                # end if conductor.cond_num_step (cdp, 10/2020)
                # call method load_user_defined_quantity to compute heat and overwrite the previous values.
                self.user_heat_function(conductor)
            # end if conductor.cond_num_step (cdp, 10/2020)
        # end self.operations["IQFUN"] (cdp, 10/2020)

    # end Get_Q

    def heat0_where(self, conductor):

        # Compute at each time step since the mesh can change
        lower_bound = np.min(
            np.nonzero(conductor.grid_features["zcoord"] >= self.operations["XQBEG"])
        )
        upper_bound = np.max(
            np.nonzero(conductor.grid_features["zcoord"] <= self.operations["XQEND"])
        )
        if self.operations["IQFUN"] == 1:
            # Square wave in time and space (cdp, 11/2020)
            if conductor.inputs["METHOD"] == "BE" or conductor.inputs["METHOD"] == "CN":
                # Backward Euler or Crank-Nicolson (cdp, 10/2020)
                if conductor.cond_num_step == 0:
                    # Initialization to Q0 value: this occurs when TQBEG = 0.0 s, i.e. \
                    # heating starts at the beginning of the simulation (cdp, 10/2020)
                    self.dict_node_pt["EXTFLX"][
                        lower_bound : upper_bound + 1, 0
                    ] = self.operations["Q0"]
                elif conductor.cond_num_step > 0:
                    if conductor.cond_num_step == 1:
                        # Store the old values only immediately after the initializzation, \
                        # since after that the whole SYSLOD array is saved and there is no \
                        # need to compute twice the same values (cdp, 10/2020)
                        self.dict_node_pt["EXTFLX"][:, 1] = self.dict_node_pt["EXTFLX"][
                            :, 0
                        ].copy()
                    self.dict_node_pt["EXTFLX"][
                        lower_bound : upper_bound + 1, 0
                    ] = self.operations["Q0"]
                # end if (cdp, 10/2020)
            elif conductor.inputs["METHOD"] == "AM4":
                # Adams-Moulton 4 (cdp, 10/2020)
                if conductor.cond_num_step == 0:
                    self.dict_node_pt["EXTFLX"][
                        lower_bound : upper_bound + 1, 0
                    ] = self.operations["Q0"]
                    for cc in range(1, 4):
                        # initialize the other columns to the same array: dummy steady \
                        # state (cdp, 10/2020)
                        self.dict_node_pt["EXTFLX"][:, cc] = self.dict_node_pt[
                            "EXTFLX"
                        ][:, 0].copy()
                    # end for cc (cdp, 10/2020)
                elif conductor.cond_num_step > 0:
                    self.dict_node_pt["EXTFLX"][:, 1:4] = self.dict_node_pt["EXTFLX"][
                        :, 0:3
                    ].copy()
                    self.dict_node_pt["EXTFLX"][
                        lower_bound : upper_bound + 1, 0
                    ] = self.operations["Q0"]
                # end if (cdp, 10/2020)
            # end if conductor.inputs["METHOD"] (cdp, 10/2020)
        # end if self.operations["IQFUN"] (cdp, 11/2020)

    # end Q0_where

    def user_heat_function(self, arg):
        # Method that allows user to define an arbitrary function for heat
        # load.
        # To be edited.
        pass

    def jhtflx_new_0(self, conductor):  # tesded: ok (cdp, 06/2020)

        """
        ############################################################################
        #              JHTFLX_new_0(self, conductor)
        ############################################################################
        #
        # Method that initialize to zero the Joule heating flux in python objects
        # of class SolidComponent.
        #
        ############################################################################
        # VARIABLE    I/O    TYPE              DESCRIPTION                      UNIT
        # --------------------------------------------------------------------------
        # self        I      object            python object of
        #                                      class SolidComponent           -
        # zcoord*     I      np array float    conductor spatial
        #                                      discretization                  m
        # JHTFLX      O      np array float    Joule heating flux
        #                                      vector                          W/m^2
        #############################################################################
        # Invoched functions/methods: none
        #
        ############################################################################
        # * zcoord is given by conductor.zcoord
        # N.B. JHTFLX is a SolidComponent attribute so its value can be assigned
        # directly, it has the same of shape of zcoord and it is a np array.
        ############################################################################
        #
        # Author D. Placido Polito 06/2020
        #
        ############################################################################
        """

        # Method JHTFLX_new_0 starts here. (cdp, 06/2020)
        if conductor.inputs["METHOD"] == "BE" or conductor.inputs["METHOD"] == "CN":
            # Backward Euler or Crank-Nicolson (cdp, 10/2020)
            if conductor.cond_time[-1] == 0:
                # Initialization (cdp, 10/2020)
                self.dict_node_pt["JHTFLX"] = np.zeros(
                    (conductor.grid_features["N_nod"], 2)
                )
            elif conductor.cond_time[-1] > 0:
                if conductor.cond_num_step == 1:
                    # Store the old values only immediately after the initializzation, \
                    # since after that the whole SYSLOD array is saved and there is no \
                    # need to compute twice the same values (cdp, 10/2020)
                    self.dict_node_pt["JHTFLX"][:, 1] = self.dict_node_pt["JHTFLX"][
                        :, 0
                    ].copy()
                # Update value at the current time step (cdp, 10/2020)
                self.dict_node_pt["JHTFLX"][:, 0] = 0.0
            # end if conductor.cond_time[-1] (cdp, 10/2020)
        elif conductor.inputs["METHOD"] == "AM4":
            # Adams-Moulton 4 (cdp, 10/2020)
            if conductor.cond_time[-1] == 0:
                # Initialization (cdp, 10/2020)
                self.dict_node_pt["JHTFLX"] = np.zeros(
                    (conductor.grid_features["N_nod"], 4)
                )
            elif conductor.cond_time[-1] > 0:
                self.dict_node_pt["JHTFLX"][:, 1:4] = self.dict_node_pt["JHTFLX"][
                    :, 0:3
                ].copy()
                # Update value at the current time step (cdp, 10/2020)
                self.dict_node_pt["JHTFLX"][:, 0] = 0.0
                # end if conductor.cond_time[-1] (cdp, 10/2020)
        # end if conductor.inputs["METHOD"] (cdp, 10/2020)

    # end JHTFLX_new_0

    def initialize_electric_quantities(self, conductor):
        """Method that initializes to zero some arrays that are an outcome of the electric method for each SolidComponent object:

        * self.dict_Gauss_pt["current_along"]
        * self.dict_Gauss_pt["delta_voltage_along"]
        * self.dict_Gauss_pt["delta_voltage_along_sum"]
        * self.dict_Gauss_pt["integral_power_el_res"]
        * self.dict_Gauss_pt["linear_power_el_resistance"]
        * self.dict_node_pt["total_linear_power_el_cond"]
        * self.dict_node_pt["total_power_el_cond"]
        * self.dict_node_pt["integral_power_el_cond"]
        """

        # Aliases
        n_nod = conductor.grid_features["N_nod"]
        n_elem = conductor.grid_input["NELEMS"]
        method = conductor.inputs["METHOD"]

        self.dict_Gauss_pt["current_along"] = np.zeros(n_elem)
        self.dict_Gauss_pt["delta_voltage_along"] = np.zeros(n_elem)
        self.dict_Gauss_pt["delta_voltage_along_sum"] = np.zeros(n_elem)

        # Introduced to evaluate the integral value Joule power due to electric 
        # resistance along the cable in the electric module as Delta V * I.
        self.dict_Gauss_pt["integral_power_el_res"] = np.zeros(n_elem)

        self.dict_node_pt["total_power_el_cond"] = np.zeros(n_nod)
        # Introduced to evaluate the integral value Joule power due to electric 
        # conductance across the cable in the electric module.
        self.dict_node_pt["integral_power_el_cond"] = np.zeros(n_nod)

        if method == "BE" or method == "CN":
            # Backward Euler or Crank-Nicolson.
            self.dict_node_pt["total_linear_power_el_cond"] = np.zeros(
                (n_nod,2)
            )
            self.dict_Gauss_pt["linear_power_el_resistance"] = np.zeros(
                    (n_elem,2)
                )
        elif method == "AM4":
            # Adams-Moulton 4.
            self.dict_node_pt["total_linear_power_el_cond"] = np.zeros(
                (n_nod,4)
            )
            self.dict_Gauss_pt["linear_power_el_resistance"] = np.zeros(
                    (n_elem,4)
                )

    def set_power_array_to_zeros(self, conductor):
        """Method that sets to zeros arrays used in the evaluation of the Joule power afther the initialization carried out with method initialize_electric_quantities.

        Args:
            conductor (Conductor): object with all the information useful for the computation.
        """

        # Aliases
        n_nod = conductor.grid_features["N_nod"]
        n_elem = conductor.grid_input["NELEMS"]

        # Set the following arrays to zeros
        self.dict_Gauss_pt["integral_power_el_res"] = np.zeros(n_elem)
        self.dict_node_pt["integral_power_el_cond"] = np.zeros(n_nod)

    def get_joule_power_along(self, conductor: object):
        """Method that evaluate the contribution to the integral power in the element of Joule power (in W/m) due to the voltage difference (due to electric resistance and inductance) along the SolidComponent objects. This approach was discussed with prof. Zach Hartwig and Dr. Nicolò Riva and is more general and coservative that the evaluation that takes into account only the contribution of the electric resistance.
        This method should be called in the electric method, when the transient solution is used and only for current carriers.

        Args:
            conductor (object): ConductorComponent object with all informations to make the calculation.
        """

        # Alias
        method = conductor.inputs["METHOD"]

        # Finalize the evaluation of the integral value of the Joule power due 
        # to electric resistance along the current carrier started in method 
        # conductor.eval_integral_joule_power. Array integral_power_el_res 
        # stores the numerator (energy in J) that here is divided by the 
        # thermal hydraulic to get again a power (W), which is further divided 
        # by the length of the discreticazion element corrected by cos(theta) 
        # to get a linear power density (W/m):
        # P_Joule = sum_1^N_em P_{Joule,i} / (Delta_t_TH * Delta_z * costheta)
        integral_j_pow_along = (
            self.dict_Gauss_pt["integral_power_el_res"]
            / (conductor.grid_features["delta_z"] * self.inputs["COSTETA"]
            * conductor.time_step)
        )

        if method == "BE" or method == "CN":
            if conductor.cond_num_step == 1:
                # Store the old values only immediately after the
                # initializzation, since after that the whole SYSLOD array
                # is saved and there is no need to compute twice the same
                # values.
                self.dict_Gauss_pt["linear_power_el_resistance"][
                    :,1
                ] = self.dict_Gauss_pt["linear_power_el_resistance"][
                        :,0
                    ].copy()

            # Assign the integral value of the linear power density of the 
            # Joule power along current carrier.
            self.dict_Gauss_pt["linear_power_el_resistance"][:, 0] = (
                integral_j_pow_along
            )
        elif method == "AM4":
            # Adams-Moulton 4.
            self.dict_Gauss_pt["linear_power_el_resistance"][
                :, 1:4
            ] = self.dict_Gauss_pt["linear_power_el_resistance"][:, 0:3].copy()
            # Assign the integral value of the linear power density of the 
            # Joule power along current carrier.
            self.dict_Gauss_pt["linear_power_el_resistance"][:, 0] = (
                integral_j_pow_along
            )

    def get_joule_power_across(self, conductor: object):
        """Method that evaluates the contribution to the total power in the nodes of Joule power (in W/m) due to the electic conductance across the SolidComponent objects.
        This method should be called in the electric method, when the transient solution is used and only for current carriers.

        Args:
            conductor (object): ConductorComponent object with all informations to make the calculation.
        """

        # Alias
        method = conductor.inputs["METHOD"]

        # Finalize the evaluation of the integral value of the Joule power due 
        # to electric conductance between the current carrier started in method 
        # conductor.eval_integral_joule_power. Array integral_power_el_cond 
        # stores the numerator (energy in J) that here is divided by the 
        # thermal hydraulic to get again a power (W), which is further divided 
        # by a suitable length to get a linear power density (W/m):
        # P_Joule = sum_1^N_em P_{Joule,i} / (Delta_t_TH * Delta_z_tilde)
        integral_j_pow_across = (
            self.dict_node_pt["integral_power_el_cond"]
            / (conductor.grid_features["delta_z_tilde"] * conductor.time_step)
        )

        if method == "BE" or method == "CN":
            # Backward Euler or Crank-Nicolson.
            if conductor.cond_num_step == 1:
                # Store the old values only immediately after the
                # initializzation, since after that the whole SYSLOD array
                # is saved and there is no need to compute twice the same
                # values.
                self.dict_node_pt["total_linear_power_el_cond"][
                    :, 1
                ] = self.dict_node_pt["total_linear_power_el_cond"][:, 0].copy()
            # Assign the integral value of the linear power density of the 
            # Joule power across current carrier.
            self.dict_node_pt["total_linear_power_el_cond"][:, 0] = (
                integral_j_pow_across
            )
        elif method == "AM4":
            # Adams-Moulton 4.
            self.dict_node_pt["total_linear_power_el_cond"][
                :, 1:4
            ] = self.dict_node_pt["total_linear_power_el_cond"][:, 0:3].copy()
            # Assign the integral value of the linear power density of the 
            # Joule power across current carrier.
            self.dict_node_pt["total_linear_power_el_cond"][:, 0] = (
                integral_j_pow_across
            )

    def get_joule_power_along_steady(self, conductor: object):
        """Method that evaluate the contribution to the total power in the element of Joule power (in W/m) due to the voltage difference (due to electric resistance and inductance) along the SolidComponent objects.
        This approach was discussed with prof. Zach Hartwig and Dr. Nicolò Riva and is more general and coservative that the evaluation that takes into account only the contribution of the electric resistance.
        This method should be called in the electric method, when the transient solution is used and for the solid component that actually carries a current. It works when the steady state solution for the electric module is computed.

        Args:
            conductor (object): ConductorComponent object with all informations to make the calculation.
        """

        # Alias
        current = self.dict_Gauss_pt["current_along"]
        voltage = self.dict_Gauss_pt["delta_voltage_along"]
        d_z_tilde = conductor.grid_features["delta_z"] * self.inputs["COSTETA"]

        # Evaluate Joule linear power along the strand in W, due
        # to electric resistances only for current carriers:
        # P_along = Delta_Phi_along * I_along
        # N.B. this evaluation accounts aslo for the voltage due to the 
        # inductance and is a conservative an more general approach. Discussed 
        # with prof. Zach Hartwig and Dr. Nicolò Riva.
        self.dict_Gauss_pt["integral_power_el_res"] = voltage * current

        # Convert W in W/m keping into account the cos(theta).
        # This is independent of the time integration method since at the 
        # initialization (when this function is used) BE, CN and AM4 should 
        # fill only the 0 index column.
        self.dict_Gauss_pt["linear_power_el_resistance"][:, 0] = (
            self.dict_Gauss_pt["integral_power_el_res"] / d_z_tilde
        )

    def get_joule_power_across_steady(self, conductor: object):
        """Method that evaluates the contribution to the total power in the nodes of Joule power (in W/m) due to the electic conductance across the SolidComponent objects.
        This method should be called in the electric method, when the transient solution is used. It works when the steady state solution for the electric module is computed.

        Args:
            conductor (object): ConductorComponent object with all informations to make the calculation.
        """

        # Evaluate total Joule linear power across the strand in W/m, due to 
        # electric conductance only for current carriers:
        # P_l_t = P_t / Delta_z_tilde
        self.dict_node_pt["total_linear_power_el_cond"][:, 0] = (
            self.dict_node_pt["total_power_el_cond"]
            / conductor.grid_features["delta_z_tilde"]
        )

    def set_energy_counters(self, conductor):
        # tesded: ok (cdp, 06/2020)

        """
        ############################################################################
        #              Set_energy_counters(self, conductor)
        ############################################################################
        #
        # Method that initialize to zero the external energy and Joule heating in
        # python objects of class SolidComponent.
        #
        ############################################################################
        # VARIABLE    I/O    TYPE              DESCRIPTION                      UNIT
        # --------------------------------------------------------------------------
        # self        I      object            python object of
        #                                      class SolidComponent           -
        # zcoord*     I      np array float    conductor spatial
        #                                      discretization                  m
        # EEXT        O      np array float    external heating vector         MJ
        # EJHT        O      np array float    external Joule heating
        #                                      vector                          MJ
        #############################################################################
        # Invoched functions/methods: none
        #
        ############################################################################
        # * zcoord is given by conductor.zcoord
        # N.B. EEXT and EJHT are SolidComponent attributes so therir value can be
        # assigned directly they have the same of shape of zcoord and they are np
        # arrays.
        ############################################################################
        #
        # Author D. Placido Polito 06/2020
        #
        ############################################################################
        """

        # Method Set_energy_counters starts here. (cdp, 06/2020)

        if conductor.inputs["METHOD"] == "BE" or conductor.inputs["METHOD"] == "CN":
            # Backward Euler or Crank-Nicolson (cdp, 10/2020)
            if conductor.cond_time[-1] == 0:
                # Initialization (cdp, 10/2020)
                self.dict_node_pt["EEXT"] = np.zeros(
                    (conductor.grid_features["N_nod"], 2)
                )
                self.dict_node_pt["EJHT"] = np.zeros(
                    (conductor.grid_features["N_nod"], 2)
                )
            elif conductor.cond_time[-1] > 0:
                if conductor.cond_num_step == 1:
                    # Store the old values only immediately after the initializzation, \
                    # since after that the whole SYSLOD array is saved and there is no \
                    # need to compute twice the same values (cdp, 10/2020)
                    self.dict_node_pt["EEXT"][:, 1] = self.dict_node_pt["EEXT"][
                        :, 0
                    ].copy()
                self.dict_node_pt["EJHT"][:, 1] = self.dict_node_pt["EJHT"][:, 0].copy()
                # Update value at the current time step (cdp, 10/2020)
                self.dict_node_pt["EEXT"][:, 0] = 0.0
                self.dict_node_pt["EJHT"][:, 0] = 0.0
            # end if conductor.cond_time[-1] (cdp, 10/2020)
        elif conductor.inputs["METHOD"] == "AM4":
            # Adams-Moulton 4 (cdp, 10/2020)
            if conductor.cond_time[-1] == 0:
                # Initialization (cdp, 10/2020)
                self.dict_node_pt["EEXT"] = np.zeros(
                    (conductor.grid_features["N_nod"], 4)
                )
                self.dict_node_pt["EJHT"] = np.zeros(
                    (conductor.grid_features["N_nod"], 4)
                )
            elif conductor.cond_time[-1] > 0:
                self.dict_node_pt["EEXT"][:, 1:4] = self.dict_node_pt["EEXT"][
                    :, 0:3
                ].copy()
                self.dict_node_pt["EJHT"][:, 1:4] = self.dict_node_pt["EJHT"][
                    :, 0:3
                ].copy()
                # Update value at the current time step (cdp, 10/2020)
                self.dict_node_pt["EEXT"][:, 0] = 0.0
                self.dict_node_pt["EJHT"][:, 0] = 0.0
            # end if conductor.cond_time[-1] (cdp, 10/2020)
        # end if conductor.inputs["METHOD"] (cdp, 10/2020)

    # end Set_energy_counters

    def deal_with_flag_IOP_MODE(self):
        """Method that checks and converts values assigned to flag IOP_MODE.

        Raises:
            ValueError: if self.operations['IOP_MODE'] is a string different from 'none'.
        """
        # Check if self.operations["IOP_MODE"] is a string
        if type(self.operations["IOP_MODE"]) == str:
            # Convert string to lower case.
            self.operations["IOP_MODE"] = self.operations["IOP_MODE"].lower()
            # Check if string value is "none".
            if self.operations["IOP_MODE"] == "none":
                # Convert string "none" to None.
                self.operations["IOP_MODE"] = None
            else:
                raise ValueError(
                    f"Not valid value to flag self.operations['IOP_MODE']. Possible values are -1, 0 or 'none', current value is {self.operations['IOP_MODE']=}. Check sheet {self.identifier} of input file conuctor_operation.xlsx.\n"
                )
        else:
            # Temporary solution to mangage input file loading, strange
            # behavior: 1 are converted to True but 0 not converted to False.
            if self.operations["IOP_MODE"] == True:
                self.operations["IOP_MODE"] = 1
            elif self.operations["IOP_MODE"] == False:
                self.operations["IOP_MODE"] = 0