from solid_components import SolidComponents
from strands import Strands
import pandas as pd

# Modified by D.Placido PoliTo
# data 05/05/2020
class SuperConductor(Strands):

    # Class for superconductor strands objects

    ### INPUT PARAMETERS
    # some are inherited form the parent classes Strands and SolidComponents

    ### THERMOPHYSICAL PROPERTIES
    # inherited from class SolidComponents

    ### OPERATIONAL PARAMETERS
    # inherited from parent classes Strands and SolidComponents

    ### COMPUTED IN INITIALIZATION
    # inherited from class SolidComponents

    ### COMPUTED VECTOR FOR MAGNETIC FIELD
    # inherited from class SolidComponents

    KIND = "Super_conductor"

    def __init__(self, simulation, sheet, icomp, name, dict_file_path):

        self.NAME = name
        # get channels ID consistently with user definition (cdp, 09/2020)
        self.ID = sheet.cell(row=3, column=4 + icomp).value

        # dictionary declaration (cdp, 11/2020)
        self.dict_input = dict()
        self.dict_operation = dict()
        self.dict_node_pt = dict()
        self.dict_Gauss_pt = dict()
        self.dict_num_step = dict()
        # Empty dictionary of list to save variable time evolutions at selected spatial coordinates.
        self.time_evol = dict(temperature=dict(), B_field=dict(), T_cur_sharing=dict())
        self.dict_scaling_input = dict()
        # Dictionary initialization: dict_input.
        self.dict_input = pd.read_excel(
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
        self.ASC = self.dict_input["CROSSECTION"]
        if self.dict_operation["IBIFUN"] != -1:
            # Remove key B_field_units.
            del self.dict_operation["B_field_units"]
        # Call SolidComponents class constructor to deal with SuperConductor time \
        # steps for current, external heating and so on (cdp, 11/2020)
        SolidComponents(simulation, self)

    def __repr__(self):
        return f"{self.__class__.__name__}(Type: {self.NAME}, ID: {self.ID})"

    def __str__(self):
        pass
