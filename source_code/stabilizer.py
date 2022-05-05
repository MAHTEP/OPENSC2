from solid_components import SolidComponents
from strands import Strands
import pandas as pd


class Stabilizer(Strands):

    # Class for copper strands objects

    ### INPUT PARAMETERS
    # some are inherited form the parent classes Strands and SolidComponents

    ### THERMOPHYSICAL PROPERTIES
    # inherited from class SolidComponents

    ##### OPERATIONAL PARAMETERS
    # inherited from parents class SolidComponents and Strands

    ### COMPUTED IN INITIALIZATION
    # inherited from class SolidComponents

    ### COMPUTED VECTOR FOR MAGNETIC FIELD
    # inherited from class SolidComponents

    KIND = "Stabilizer"

    def __init__(self, simulation, sheet, icomp, name, dict_file_path):

        self.NAME = name
        # get channels ID consistently with user definition (cdp, 09/2020)
        self.ID = sheet.cell(row=3, column=4 + icomp).value

        # dictionary declaration (cdp, 11/2020)
        self.inputs = dict()
        self.dict_operation = dict()
        self.dict_node_pt = dict()
        self.dict_Gauss_pt = dict()
        self.dict_num_step = dict()
        # Empty dictionary of list to save variable time evolutions at selected spatial coordinates.
        self.time_evol = dict(temperature=dict(), B_field=dict())
        self.dict_scaling_input = dict()
        # Dictionary initialization: inputs.
        self.inputs = pd.read_excel(
            dict_file_path["input"],
            sheet_name=sheet.title,
            skiprows=2,
            header=0,
            index_col=0,
            usecols=["Variable name", self.ID],
        )[self.ID].to_dict()
        # Dictionary initialization: dict_operation.
        self.dict_operation = pd.read_excel(
            dict_file_path["operation"],
            sheet_name=sheet.title,
            skiprows=2,
            header=0,
            index_col=0,
            usecols=["Variable name", self.ID],
        )[self.ID].to_dict()

        # Call SolidComponents class constructor to deal with Stabilizer time \
        # steps for current, external heating and so on (cdp, 11/2020)
        SolidComponents(simulation, self)
        if self.inputs["ISTABILIZER"] != "Cu":
            # remove key RRR from inputs if stabilizer is not Cu (cdp, 07/2020)
            self.inputs.pop("RRR")
        if self.dict_operation["IBIFUN"] != -1:
            # Remove key B_field_units.
            del self.dict_operation["B_field_units"]

    def __repr__(self):
        return f"{self.__class__.__name__}(Type: {self.NAME}, ID: {self.ID})"

    def __str__(self):
        pass

    # agiungere metodi per calcolo proprietà nel gauss point
