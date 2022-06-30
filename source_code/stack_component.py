# Import libraries
import numpy as np
import pandas as pd

# Import classes
from solid_component import SolidComponent
from strand_component import StrandComponent

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
    cu=thermal_conductivity_cu,
    ge=thermal_conductivity_ge,
    hc276=thermal_conductivity_hc276,
    re123=thermal_conductivity_re123,
    sn60pb40=thermal_conductivity_sn60pb40,
    ss=thermal_conductivity_ss,
)

ISOBARIC_SPECIFIC_HEAT_FUNC = dict(
    ag=isobaric_specific_heat_ag,
    al=isobaric_specific_heat_al,
    cu=isobaric_specific_heat_cu,
    ge=isobaric_specific_heat_ge,
    hc276=isobaric_specific_heat_hc276,
    re123=isobaric_specific_heat_re123,
    sn60pb40=isobaric_specific_heat_sn60pb40,
    ss=isobaric_specific_heat_ss,
)

ELECTRICAL_RESISTIVITY_FUNC = dict(
    ag=electrical_resistivity_ag,
    al=electrical_resistivity_al,
    cu=electrical_resistivity_cu,
    # ge=electrical_resistivity_ge, not defined
    hc276=electrical_resistivity_hc276,
    re123=electrical_resistivity_re123,
    sn60pb40=electrical_resistivity_sn60pb40,
    ss=electrical_resistivity_ss,
)

class StackComponent(StrandComponent:StrandComponent):
    """Class that defines StackComponents objects to model HTS stacks of tapes.

    Args:
        StrandComponent (StrandComponent): class that defines general methods for strand and stack objects.
    """
    KIND = "Stack"
