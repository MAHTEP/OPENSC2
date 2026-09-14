import numpy as np

from properties_of_materials.silver import (
    density_ag,
    density_ag_cryosoft,
    electrical_resistivity_ag,
    electrical_resistivity_ag_cryosoft,
    isobaric_specific_heat_ag_cryosoft,
    thermal_conductivity_ag,
    thermal_conductivity_ag_cryosoft,
)
from properties_of_materials.copper import (
    thermal_conductivity_cu_nist,
    electrical_resistivity_cu_nist,
)
from stack_component import StackComponent
from properties_of_materials.rare_earth_123 import thermal_conductivity_re123


def test_cryosoft_silver_reference_values_against_fortran():
    temperature = np.array([4.5, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0, 150.0])
    magnetic_field = 17.0
    rrr = 12.8

    expected_k = np.array([
        6.0481350e01,
        1.3290617e02,
        2.4195525e02,
        3.2108044e02,
        3.6890332e02,
        3.8662924e02,
        3.9936749e02,
        4.1337726e02,
    ])
    expected_rho = np.array([
        1.8160712e-09,
        1.8172006e-09,
        1.8474121e-09,
        2.3825661e-09,
        3.5386392e-09,
        4.7578692e-09,
        5.9375243e-09,
        8.9337151e-09,
    ])
    expected_cp = np.array([
        1.7233650e-01,
        1.7977964e00,
        1.5675276e01,
        7.8556503e01,
        1.3181342e02,
        1.6567111e02,
        1.8712546e02,
        2.1387729e02,
    ])

    np.testing.assert_allclose(
        thermal_conductivity_ag_cryosoft(
            temperature, magnetic_field, rrr
        ),
        expected_k,
        rtol=2e-6,
        atol=0.0,
    )
    np.testing.assert_allclose(
        electrical_resistivity_ag_cryosoft(
            temperature, magnetic_field, rrr
        ),
        expected_rho,
        rtol=2e-6,
        atol=0.0,
    )
    np.testing.assert_allclose(
        isobaric_specific_heat_ag_cryosoft(temperature),
        expected_cp,
        rtol=2e-6,
        atol=0.0,
    )


def test_cryosoft_silver_density_reference_value():
    temperature = np.array([4.5, 20.0, 100.0])
    np.testing.assert_allclose(
        density_ag_cryosoft(temperature),
        10490.0 * np.ones(temperature.shape),
    )
    assert not np.allclose(
        density_ag_cryosoft(temperature),
        density_ag(temperature),
    )


def test_legacy_silver_functions_remain_separate():
    temperature = np.array([4.5, 20.0, 60.0, 100.0])
    legacy_k = thermal_conductivity_ag(temperature)
    cryosoft_k = thermal_conductivity_ag_cryosoft(
        temperature, 17.0, 12.8
    )
    legacy_rho = electrical_resistivity_ag(temperature)
    cryosoft_rho = electrical_resistivity_ag_cryosoft(
        temperature, 17.0, 12.8
    )

    assert np.max(np.abs(legacy_k - cryosoft_k)) > 1.0e3
    assert np.max(cryosoft_rho / legacy_rho) > 100.0


def test_stack_dispatch_uses_independent_cu_and_ag_rrr():
    # Material_number includes one superconducting layer; the electrical
    # resistivity array contains only the non-superconducting layers.
    component = object.__new__(StackComponent)
    component.inputs = {
        "Material_number": 3,
        "RRR": 50.0,
        "RRR_Ag": 12.8,
    }

    component.tape_material = np.array(["ybco", "cu", "ag_cryosoft"])
    component.tape_material_not_sc = np.array(["cu", "ag_cryosoft"])

    component.material_thickness = np.array([0.10, 0.45, 0.45])
    component.material_thickness_not_sc = np.array([0.45, 0.45])
    component.tape_thickness = 1.0
    component.tape_thickness_not_sc = 0.90

    component.thermal_conductivity_function = np.array(
        [
            thermal_conductivity_re123,
            thermal_conductivity_cu_nist,
            thermal_conductivity_ag_cryosoft,
        ],
        dtype=object,
    )
    component.electrical_resistivity_function_not_sc = np.array(
        [
            electrical_resistivity_cu_nist,
            electrical_resistivity_ag_cryosoft,
        ],
        dtype=object,
    )

    prop = {
        "temperature": np.array([4.5, 20.0, 100.0]),
        "B_field": np.array([17.0, 17.0, 17.0]),
    }

    expected_k = (
        0.10 * thermal_conductivity_re123(prop["temperature"])
        + 0.45
        * thermal_conductivity_cu_nist(
            prop["temperature"], prop["B_field"], 50.0
        )
        + 0.45
        * thermal_conductivity_ag_cryosoft(
            prop["temperature"], prop["B_field"], 12.8
        )
    )
    np.testing.assert_allclose(
        component.stack_thermal_conductivity(prop), expected_k
    )

    rho_cu = electrical_resistivity_cu_nist(
        prop["temperature"], prop["B_field"], 50.0
    )
    rho_ag = electrical_resistivity_ag_cryosoft(
        prop["temperature"], prop["B_field"], 12.8
    )
    expected_rho = 0.90 / (0.45 / rho_cu + 0.45 / rho_ag)

    np.testing.assert_allclose(
        component.stack_electrical_resistivity_not_sc(prop), expected_rho
    )
