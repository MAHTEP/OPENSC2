"""The FlidComponents is a class SHOULD BE DONE AFTER A MEETING
properties:
-number: Number of fluid channels (M)
-type: Type of fluid	To be chosen among He or N2 – thermophysical properties to be evaluated, according to the choice
-crossSection: Cross section of any of the M channels	Typically in [mm^2]
-costheta: Cos(teta) of any of the M channels	Angle of inclination wrt the jacket axis
-contactPerimeterSC: Contact perimeter (equivalent to area per unit length) with any of the Z SC strands and K stabilizer strands, or N SC + stabilizer strands	Typically in [mm]
-htc: Heat transfer coefficient (equivalent to area per unit length) with any of the Z SC strands and K stabilizer strands, or N SC + stabilizer strands	Typically in [W/m2K] NON SOLO NUMERI QUI, MA ANCHE FUNZIONI
-contactPerimeterJ: Contact perimeter (equivalent to area per unit length) with the jacket	Typically in [mm]
-htcJ: Heat transfer coefficient (equivalent to area per unit length) with the jacket	Typically in [W/m2K] NON SOLO NUMERI QUI, MA ANCHE FUNZIONI
-perimeterTransfer: Perimeter (equivalent to area per unit length) available for the mass transfer with the (M-1) fluid components	Typically in [mm]
-perimeterConductiveFluid: Perimeter (equivalent to area per unit length) available for the pure conductive heat transfer with the (M-1) fluid components	Typically in [mm]
-htc: Heat transfer coefficient (equivalent to area per unit length) with the (M-1) fluid components	Typically in [W/m2K] NON SOLO NUMERI QUI, MA ANCHE FUNZIONI
"""

import pandas as pd
from collections import namedtuple
from typing import Union, NamedTuple
from utility_functions.auxiliary_functions import get_from_xlsx

class FluidComponentInput:
    """Interface class used to get the input data for fluid components objects. Attributes inputs and operations are inherited from this class by Coolant and Channel class. Being an interface, this class has only the constructror (__init__) method."""

    def __init__(self, sheet, sheetOpar, dict_file_path, identifier):
        """[summary]

        Args:
            sheet ([type]): [description]
            sheetOpar ([type]): [description]
            cond ([type]): [description]
            identifier ([type]): [description]
        """
        # Dictionary declaration (cdp, 11/2020)
        self.inputs = dict()
        self.operations = dict()
        # Dictionary initialization: inputs.
        self.inputs = pd.read_excel(
            dict_file_path["input"],
            sheet_name=sheet.title,
            skiprows=2,
            header=0,
            index_col=0,
            usecols=["Variable name", identifier],
        )[identifier].to_dict()
        # Dictionary initialization: operations.
        self.operations = pd.read_excel(
            dict_file_path["operation"],
            sheet_name=sheetOpar.title,
            skiprows=2,
            header=0,
            index_col=0,
            usecols=["Variable name", identifier],
        )[identifier].to_dict()
        # Tuning input and operational parameters according to input flags value
        if self.inputs["ISRECTANGULAR"] == False:
            # Remove keys SIDE1 and SIDE2 from inputs
            del self.inputs["SIDE1"]
            del self.inputs["SIDE2"]

    # End  method __init__.


# End class FluidComponentInput.

# Import classes Channel and Coolant: done here to avoid circular import error
from channel import Channel
from coolant import Coolant


class FluidComponent:

    # class variable shared by all instances
    KIND = "Fluid_component"

    def __init__(self, sheet, sheetOpar, icomp, dict_file_path):
        """[summary]

        Args:
            sheet ([type]): [description]
            sheetOpar ([type]): [description]
            icomp ([type]): [description]
            cond ([type]): [description]
        """
        # Get channels ID consistently with user definition (cdp, 09/2020)
        self.identifier = sheet.cell(row=3, column=4 + icomp).value
        # Instance of class Channel (build a coolant object)
        self.channel = Channel(sheet, sheetOpar, dict_file_path, self.identifier)
        # Instance of class Coolant (build a coolant object)
        self.coolant = Coolant(sheet, sheetOpar, dict_file_path, self.identifier)
        self.coordinate = dict()

    # End method __init__.

    def __repr__(self):
        """[summary]

        Returns:
            [type]: [description]
        """
        return f"{self.__class__.__name__}(Type: {self.name}, identifier: {self.identifier})"

    # End method __repr__.

    def __str__(self):
        pass

    # End method __str__.

    def build_th_bc_index(
        self,
        conductor:object,
    ):
        """Method that builds and initialize the data structure to store the values of the index used to assign the boundary conditions at inlet and outlet for the thermal hydraulic problem. The data structure is a namedtuple with three keys (velocity pressure and temperature), each key has a namedtuple with two keys (forward and backward) to properly assing the boundary conditions in case of coolants in counter current flow.

        Args:
            conductor (Conductor): object with all the information of the conductor.
        """

        # ALIAS
        ndf = conductor.dict_N_equation["NODOFS"]
        # eq_idx (NamedTuple): collection of fluid equation index (velocity, 
        # pressure and temperaure equations).
        eq_idx = conductor.equation_index[self.identifier]

        # Namedtuple constuctor
        BC_idx = namedtuple("BC_idx",("velocity","pressure","temperature"))
        Flow_dir = namedtuple("Flow_dir",("forward","backward"))
        
        # Build namedtuple with the index used to assign the inlet BC.
        self.inl_idx = BC_idx(
            # Inlet velocity index (forward and backward flow).
            velocity=Flow_dir(
                forward=eq_idx.velocity, # first node
                backward=- ndf + eq_idx.velocity, # last node
            ),
            # Inlet pressure index (forward and backward flow).
            pressure=Flow_dir(
                forward=eq_idx.pressure, # first node
                backward=- ndf + eq_idx.pressure, # last node
            ),
            # Inlet temperature index (forward and backward flow).
            temperature=Flow_dir(
                forward=eq_idx.temperature, # first node
                backward=- ndf + eq_idx.temperature, # last node
            ),
        )

        # Build namedtuple with the index used to assign the outlet BC.
        self.out_idx = BC_idx(
            # Outlet velocity index (forward and backward flow).
            velocity=Flow_dir(
                forward=- ndf + eq_idx.velocity, # last node
                backward=eq_idx.velocity, # first node
            ),
            # Outlet pressure index (forward and backward flow).
            pressure=Flow_dir(
                forward=- ndf + eq_idx.pressure, # last node
                backward=eq_idx.pressure, # first node
            ),
            # Outlet temperature index (forward and backward flow).
            temperature=Flow_dir(
                forward=- ndf + eq_idx.temperature, # last node
                backward=eq_idx.temperature, # first node
            ),
        )


    def impose_pressure_drop(
        self,
        ndarrays:tuple,
        conductor:object,
        path:str,
        )->tuple:
        
        """Method that imposes a pressure drop as boundary conditions for the thermal hydraulic problem. Imposed quantities are:
            * inlet pressure
            * inlet temperature (or outlet temperature if back flow)
            * outlet pressure
        The method suitably accounts for forward (left to right)/backward (right to left) flow as well as for back flow (part of the fluid moves towards the inlet due to a pressure build up).

        Args:
            ndarrays (tuple): collection of known term vector (Known) and system matrix (SYSMAT) np.ndarray on wich impose the boundary conditions.
            conductor (Conductor): object with all the information of the conductor.
            path (str): path of the auxiliary input file with the values for the boundary conditions (if INTIAL = -1).

        Returns:
            tuple: collection of known term vector (Known) and system matrix (SYSMAT) np.ndarray with imposed boundary conditions.
        """
        
        known,sysmat = ndarrays
        # Get boundary condition values.
        if self.coolant.operations["INTIAL"] == 1:
            # inlet pressure.
            p_inl = self.coolant.operations["PREINL"]
            # inlet temperature.
            T_inl = self.coolant.operations["TEMINL"]
            # outlet pressure.
            p_out = self.coolant.operations["PREOUT"]
            # outlet temperature: to be assigned if outlet velocity is negative.
            T_out = self.coolant.operations["TEMOUT"]
        else:
            # get from file with interpolation in time.
            [flow_par, flagSpecfield] = get_from_xlsx(
                conductor,
                path,
                self,
                "INTIAL",
                self.coolant.operations["INTIAL"]
            )
            print(
                f"""flagSpecfield == {flagSpecfield}: still to be decided if 
            it is useful and if yes still to be defined\n"""
            )
            # inlet pressure.
            p_inl = flow_par[2]
            # inlet temperature.
            T_inl = flow_par[0]
            # outlet pressure.
            p_out = flow_par[3]
            # outlet temperature: to be assigned if outlet velocity is negative.
            T_out = flow_par[1]
        
        # ALIAS
        _, inl_p_idx, inl_t_idx = self.inl_idx # tuple unpack
        _, out_p_idx, out_t_idx = self.out_idx # tuple unpack
        main_d_idx = conductor.dict_band["Main_diag"]
        flow_dir = self.channel.flow_dir[0]
        velocity = self.coolant.dict_node_pt["velocity"]
        
        # Assign BC
        if flow_dir == "forward":
            # p_inl
            sysmat[:,inl_p_idx.forward] = 0.0
            # main diagonal.
            sysmat[main_d_idx,inl_p_idx.forward] = 1.0
            known[inl_p_idx.forward] = p_inl
            # p_out
            sysmat[:,out_p_idx.forward] = 0.0
            # main diagonal.
            sysmat[main_d_idx,out_p_idx.forward] = 1.0
            known[out_p_idx.forward] = p_out
            # T_inl
            if velocity[0] > 0:
                sysmat[:,inl_t_idx.forward] = 0.0
                # main diagonal.
                sysmat[main_d_idx,inl_t_idx.forward] = 1.0
                known[inl_t_idx.forward] = T_inl
            # T_out
            if velocity[-1] < 0:
                sysmat[:,out_t_idx.forward] = 0.0
                # main diagonal.
                sysmat[main_d_idx,out_t_idx.forward] = 1.0
                known[out_t_idx.forward] = T_out
        elif flow_dir == "backward":
            # p_inl
            sysmat[:,inl_p_idx.backward] = 0.0
            # main diagonal.
            sysmat[main_d_idx,inl_p_idx.backward] = 1.0
            known[inl_p_idx.backward] = p_inl
            # p_out
            sysmat[:,out_p_idx.backward] = 0.0
            # main diagonal.
            sysmat[main_d_idx,out_p_idx.backward] = 1.0
            known[out_p_idx.backward] = p_out
            # T_inl
            if velocity[-1] < 0:
                sysmat[:,inl_t_idx.backward] = 0.0
                # main diagonal.
                sysmat[main_d_idx,inl_t_idx.backward] = 1.0
                known[inl_t_idx.backward] = T_inl
            # T_out
            if velocity[0] > 0:
                sysmat[:,out_t_idx.backward] = 0.0
                # main diagonal.
                sysmat[main_d_idx,out_t_idx.backward] = 1.0
                known[out_t_idx.backward] = T_out

        return known,sysmat

    def impose_inl_p_out_v(
        self,
        ndarrays:tuple,
        conductor:object,
        path:str,
        )->tuple:
        
        """Method that imposes a inlet pressure and outlet velocity as boundary conditions for the thermal hydraulic problem. Imposed quantities are:
            * inlet pressure
            * inlet temperature (or outlet temperature if backflow)
            * outlet velocity
        The method suitably accounts for forward (left to right)/backward (right to left) flow as well as for back flow (part of the fluid moves towards the inlet due to a pressure build up).

        Args:
            ndarrays (tuple): collection of known term vector (Known) and system matrix (SYSMAT) np.ndarray on wich impose the boundary conditions.
            conductor (Conductor): object with all the information of the conductor.
            path (str): path of the auxiliary input file with the values for the boundary conditions (if INTIAL = -2).

        Returns:
            tuple: collection of known term vector (Known) and system matrix (SYSMAT) np.ndarray with imposed boundary conditions.
        """

        # INLET AND OUTLET RESERVOIRS, INLET CONDITIONS AND FLOW SPECIFIED
        known,sysmat = ndarrays
        # Get boundary condition values.
        if self.coolant.operations["INTIAL"] == 2:
            # outlet mass flow rate.
            mfr_out = self.coolant.operations["MDTOUT"]
            # inlet pressure.
            p_inl = self.coolant.operations["PREINL"]
            # inlet temperature.
            T_inl = self.coolant.operations["TEMINL"]
            # outlet temperature: to be assigned if outlet velocity is negative.
            T_out = self.coolant.operations["TEMOUT"]
        else:  # N.B va aggiustato per renderlo conforme caso positivo!
            # all values from flow_dummy.xlsx: call get_from_xlsx.
            [flow_par, flagSpecfield] = get_from_xlsx(
                conductor,
                path,
                self,
                "INTIAL",
                self.coolant.operations["INTIAL"],
            )
            print(
                f"""flagSpecfield == {flagSpecfield}: still to be decided if it 
            useful and if yes still to be defined\n"""
            )
            mfr_out = flow_par[3]  # inlet mass flow rate.
            p_inl = flow_par[2]  # inlet pressure.
            T_inl = flow_par[0]  # inlet temperature.
            # outlet temperature: to be assigned if outlet velocity is negative.
            T_out = flow_par[1]

        # ALIAS
        _, inl_p_idx, inl_t_idx = self.inl_idx # tuple unpack
        out_v_idx, _, out_t_idx = self.out_idx # tuple unpack
        main_d_idx = conductor.dict_band["Main_diag"]
        flow_dir = self.channel.flow_dir[0]
        velocity = self.coolant.dict_node_pt["velocity"]
        density = self.coolant.dict_node_pt["total_density"]
        cross_section = self.channel.inputs["CROSSECTION"]
        
        # Assign BC
        if flow_dir == "forward":
            # Flow direction from x = 0 to x = L.
            # v_out
            sysmat[:,out_v_idx.forward] = 0.0
            # main diagonal.
            sysmat[main_d_idx,out_v_idx.forward] = 1.0
            known[out_v_idx.forward] = (mfr_out / density[-1] / cross_section)
            # p_inl
            sysmat[0:, inl_p_idx.forward] = 0.0
            # main diagonal.
            sysmat[main_d_idx, inl_p_idx.forward] = 1.0
            known[inl_p_idx.forward] = p_inl
            # T_inl
            if velocity[0] > 0:
                sysmat[:,inl_t_idx.forward] = 0.0
                # main diagonal.
                sysmat[main_d_idx,inl_t_idx.forward] = 1.0
                known[inl_t_idx.forward] = T_inl
            # T_out (T_inl if mfr_out < 0)
            if velocity[-1] < 0:
                sysmat[:,out_t_idx.forward] = 0.0
                # main diagonal.
                sysmat[main_d_idx,out_t_idx.forward] = 1.0
                known[out_t_idx.forward] = T_out
        elif flow_dir == "backward":
            # Flow direction from x = L to x = 0.
            # v_out
            sysmat[:,out_v_idx.backward] = 0.0
            # main diagonal.
            sysmat[main_d_idx,out_v_idx.backward] = 1.0
            known[out_v_idx.backward] = (mfr_out / density[0] / cross_section)
            # p_inl
            sysmat[0:, inl_p_idx.backward] = 0.0
            # main diagonal.
            sysmat[main_d_idx, inl_p_idx.backward] = 1.0
            known[inl_p_idx.backward] = p_inl
            # T_inl
            if velocity[-1] < 0:
                sysmat[:,inl_t_idx.backward] = 0.0
                # main diagonal.
                sysmat[main_d_idx,inl_t_idx.backward] = 1.0
                known[inl_t_idx.backward] = T_inl
            # T_out
            if velocity[0] > 0:
                sysmat[:,out_t_idx.backward] = 0.0
                # main diagonal.
                sysmat[main_d_idx,out_t_idx.backward] = 1.0
                known[out_t_idx.backward] = T_out
        
        return known,sysmat

    def impose_inl_v_out_p(
        self,
        ndarrays:tuple,
        conductor:object,
        path:str,
        )->tuple:
        
        """Method that imposes a inlet velocity and outlet pressure as boundary conditions for the thermal hydraulic problem. Imposed quantities are:
            * inlet velocity
            * inlet temperature (or outlet temperature if backflow)
            * outlet pressure
        The method suitably accounts for forward (left to right)/backward (right to left) flow as well as for back flow (part of the fluid moves towards the inlet due to a pressure build up).

        Args:
            ndarrays (tuple): collection of known term vector (Known) and system matrix (SYSMAT) np.ndarray on wich impose the boundary conditions.
            conductor (Conductor): object with all the information of the conductor.
            path (str): path of the auxiliary input file with the values for the boundary conditions (if INTIAL = -3).

        Returns:
            tuple: collection of known term vector (Known) and system matrix (SYSMAT) np.ndarray with imposed boundary conditions.
        """

        known,sysmat = ndarrays
        # Get boundary condition values.
        if self.coolant.operations["INTIAL"] == 3:
            # outlet mass flow rate.
            mfr_inl = self.coolant.operations["MDTIN"]
            # outlet pressure.
            p_out = self.coolant.operations["PREOUT"]
            # inlet temperature.
            T_inl = self.coolant.operations["TEMINL"]
            # outlet temperature: to be assigned if outlet velocity is negative.
            T_out = self.coolant.operations["TEMOUT"]
        else:  # N.B va aggiustato per renderlo conforme caso positivo!
            # all values from flow_dummy.xlsx: call get_from_xlsx.
            [flow_par, flagSpecfield] = get_from_xlsx(
                conductor,
                path,
                self,
                "INTIAL",
                self.coolant.operations["INTIAL"],
            )
            print(
                f"""flagSpecfield == {flagSpecfield}: still to be decided if it 
            useful and if yes still to be defined\n"""
            )
            mfr_inl = flow_par[3]  # inlet mass flow rate.
            p_out = flow_par[2]  # outlet pressure.
            T_inl = flow_par[0]  # inlet temperature.
            # outlet temperature: to be assigned if outlet velocity is negative.
            T_out = flow_par[1]

        # ALIAS
        inl_v_idx, _, inl_t_idx = self.inl_idx # tuple unpack
        _, out_p_idx, out_t_idx = self.out_idx # tuple unpack
        main_d_idx = conductor.dict_band["Main_diag"]
        flow_dir = self.channel.flow_dir[0]
        velocity = self.coolant.dict_node_pt["velocity"]
        density = self.coolant.dict_node_pt["total_density"]
        cross_section = self.channel.inputs["CROSSECTION"]
        
        # Assign BC
        if flow_dir == "forward":
            # Flow direction from x = 0 to x = L.
            # v_inl
            sysmat[:,inl_v_idx.forward] = 0.0
            # main diagonal.
            sysmat[main_d_idx,inl_v_idx.forward] = 1.0
            known[inl_v_idx.forward] = (mfr_inl / density[0] / cross_section)
            # p_out
            sysmat[0:, out_p_idx.forward] = 0.0
            # main diagonal.
            sysmat[main_d_idx, out_p_idx.forward] = 1.0
            known[out_p_idx.forward] = p_out
            # T_inl
            if velocity[0] > 0:
                sysmat[:,inl_t_idx.forward] = 0.0
                # main diagonal.
                sysmat[main_d_idx,inl_t_idx.forward] = 1.0
                known[inl_t_idx.forward] = T_inl
            # T_out (T_inl if mfr_inl < 0)
            if velocity[-1] < 0:
                sysmat[:,out_t_idx.forward] = 0.0
                # main diagonal.
                sysmat[main_d_idx,out_t_idx.forward] = 1.0
                known[out_t_idx.forward] = T_out
        elif flow_dir == "backward":
            # Flow direction from x = L to x = 0.
            # v_inl
            sysmat[:,inl_v_idx.backward] = 0.0
            # main diagonal.
            sysmat[main_d_idx,inl_v_idx.backward] = 1.0
            known[inl_v_idx.backward] = (mfr_inl / density[-1] / cross_section)
            # p_out
            sysmat[0:, out_p_idx.backward] = 0.0
            # main diagonal.
            sysmat[main_d_idx, out_p_idx.backward] = 1.0
            known[out_p_idx.backward] = p_out
            # T_inl
            if velocity[-1] < 0:
                sysmat[:,inl_t_idx.backward] = 0.0
                # main diagonal.
                sysmat[main_d_idx,inl_t_idx.backward] = 1.0
                known[inl_t_idx.backward] = T_inl
            # T_out
            if velocity[0] > 0:
                sysmat[:,out_t_idx.backward] = 0.0
                # main diagonal.
                sysmat[main_d_idx,out_t_idx.backward] = 1.0
                known[out_t_idx.backward] = T_out
        
        return known,sysmat

# End class FluidComponent.
