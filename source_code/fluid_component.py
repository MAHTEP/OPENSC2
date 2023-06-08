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

# End class FluidComponent.
