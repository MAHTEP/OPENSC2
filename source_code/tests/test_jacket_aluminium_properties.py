import numpy as np

from jacket_component import JacketComponent
from properties_of_materials.aluminium import (
    density_al,
    electrical_resistivity_al,
    isobaric_specific_heat_al,
    thermal_conductivity_al,
)
from properties_of_materials.glass_epoxy import (
    density_ge,
    electrical_resistivity_ge,
    isobaric_specific_heat_ge,
    thermal_conductivity_ge,
)


def _build_jacket_for_property_test(
    *,
    jacket_material,
    insulation_material="none",
    jacket_cross_section=2.0,
    insulation_cross_section=0.0,
):
    jacket = object.__new__(JacketComponent)
    jacket.inputs = {
        "jacket_material": jacket_material,
        "insulation_material": insulation_material,
        "jacket_cross_section": jacket_cross_section,
        "insulation_cross_section": insulation_cross_section,
        "NUM_MATERIAL_TYPES": (
            1 if insulation_material == "none" else 2
        ),
        "CROSSECTION": jacket_cross_section + insulation_cross_section,
    }
    jacket._JacketComponent__reorganize_input()
    jacket._JacketComponent__jacket_density_flag = False
    return jacket


def test_pure_aluminium_jacket_uses_aluminium_property_functions():
    jacket = _build_jacket_for_property_test(jacket_material="al")
    temperature = np.array([4.5, 20.0, 80.0, 120.0])
    properties = {"temperature": temperature}

    np.testing.assert_allclose(
        jacket.jacket_density(properties),
        density_al(temperature),
    )
    np.testing.assert_allclose(
        jacket.jacket_isobaric_specific_heat(properties),
        isobaric_specific_heat_al(temperature),
    )
    np.testing.assert_allclose(
        jacket.jacket_thermal_conductivity(properties),
        thermal_conductivity_al(temperature),
    )
    np.testing.assert_allclose(
        jacket.jacket_electrical_resistivity(properties),
        electrical_resistivity_al(temperature),
    )


def test_aluminium_jacket_thermal_conductivity_follows_local_temperature():
    jacket = _build_jacket_for_property_test(jacket_material="al")
    temperature = np.array([4.5, 20.0, 40.0, 80.0, 100.0])

    conductivity = jacket.jacket_thermal_conductivity(
        {"temperature": temperature}
    )

    np.testing.assert_allclose(
        conductivity,
        thermal_conductivity_al(temperature),
    )
    assert np.ptp(conductivity) > 100.0


def test_aluminium_and_glass_epoxy_jacket_homogenization_remains_consistent():
    jacket = _build_jacket_for_property_test(
        jacket_material="al",
        insulation_material="ge",
        jacket_cross_section=2.0,
        insulation_cross_section=1.0,
    )
    temperature = np.array([4.5, 20.0, 80.0])
    properties = {"temperature": temperature}

    rho_al = density_al(temperature)
    rho_ge = density_ge(temperature)
    cp_al = isobaric_specific_heat_al(temperature)
    cp_ge = isobaric_specific_heat_ge(temperature)
    k_al = thermal_conductivity_al(temperature)
    k_ge = thermal_conductivity_ge(temperature)
    rho_el_al = electrical_resistivity_al(temperature)
    rho_el_ge = electrical_resistivity_ge(temperature)

    expected_density = (2.0 * rho_al + rho_ge) / 3.0
    np.testing.assert_allclose(
        jacket.jacket_density(properties),
        expected_density,
    )

    expected_cp = (
        2.0 * rho_al * cp_al + rho_ge * cp_ge
    ) / (2.0 * rho_al + rho_ge)
    np.testing.assert_allclose(
        jacket.jacket_isobaric_specific_heat(properties),
        expected_cp,
    )

    expected_k = (2.0 * k_al + k_ge) / 3.0
    np.testing.assert_allclose(
        jacket.jacket_thermal_conductivity(properties),
        expected_k,
    )

    expected_rho_el = 3.0 / (
        2.0 / rho_el_al + 1.0 / rho_el_ge
    )
    np.testing.assert_allclose(
        jacket.jacket_electrical_resistivity(properties),
        expected_rho_el,
    )
