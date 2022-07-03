import numpy as np
import pandas as pd

from solid_component import SolidComponent
from strand_component import StrandComponent

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

DENSITY_FUNC = dict(
    al=density_al,
    cu=density_cu,
    nb3sn=density_nb3sn,
    nbti=density_nbti,
)

THERMAL_CONDUCTIVITY_FUNC = dict(
    al=thermal_conductivity_al,
    cu=thermal_conductivity_cu_nist,
    nb3sn=thermal_conductivity_nb3sn,
    nbti=thermal_conductivity_nbti,
)

ISOBARIC_SPECIFIC_HEAT_FUNC = dict(
    al=isobaric_specific_heat_al,
    cu=isobaric_specific_heat_cu_nist,
    nb3sn=isobaric_specific_heat_nb3sn,
    nbti=isobaric_specific_heat_nbti,
)

ELECTRICAL_RESISTIVITY_FUNC = dict(
    al=electrical_resistivity_al,
    cu=electrical_resistivity_cu_nist,
)

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

    def __init__(self, simulation, sheet, icomp, name, dict_file_path):

        self.NAME = name
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
        self.time_evol = dict(temperature=dict(), B_field=dict(), T_cur_sharing=dict())
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

        self.ASC = self.inputs["CROSSECTION"] / (1.0 + self.inputs["STAB_NON_STAB"])
        # Call SolidComponent class constructor to deal with StrandMixedComponent time \
        # steps for current, external heating and so on (cdp, 11/2020)
        SolidComponent(simulation, self)
        if self.inputs["ISTABILIZER"] != "Cu":
            # remove key RRR from inputs if stabilizer is not Cu (cdp, 07/2020)
            self.inputs.pop("RRR")
        if self.operations["IBIFUN"] != -1:
            # Remove key B_field_units.
            del self.operations["B_field_units"]
        # end if (cdp, 07/2020)

    def __repr__(self):
        return f"{self.__class__.__name__}(Type: {self.NAME}, identifier: {self.identifier})"

    def __str__(self):
        pass
