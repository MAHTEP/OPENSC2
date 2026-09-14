import numpy as np

from properties_of_materials.aluminium import (
    density_al,
    density_al_cryosoft,
    electrical_resistivity_al,
    electrical_resistivity_al_cryosoft,
    isobaric_specific_heat_al_cryosoft,
    thermal_conductivity_al,
    thermal_conductivity_al_cryosoft,
)
from strand_stabilizer_component import StrandStabilizerComponent
from jacket_component import JacketComponent


def test_cryosoft_aluminium_reference_values_at_benchmark_conditions():
    temperature = np.array([4.5, 20.0, 60.0, 100.0, 150.0])
    expected_k = np.array([
        667.0182837569719,
        1780.199606429915,
        595.8186633286475,
        296.72505532775267,
        243.83465961984385,
    ])
    expected_rho = np.array([
        1.6421707970129425e-10,
        1.9937689701419716e-10,
        1.152006785645088e-09,
        4.623103502834715e-09,
        1.0207345042599918e-08,
    ])
    np.testing.assert_allclose(
        thermal_conductivity_al_cryosoft(temperature, 17.0, 300.0),
        expected_k, rtol=2e-12, atol=0.0,
    )
    np.testing.assert_allclose(
        electrical_resistivity_al_cryosoft(temperature, 17.0, 300.0),
        expected_rho, rtol=2e-12, atol=0.0,
    )


def test_cryosoft_aluminium_cp_and_density_reference_values():
    temperature = np.array([4.5, 20.0, 60.0, 100.0])
    expected_cp = np.array([
        0.310155529875,
        8.853871065000003,
        213.01710611220915,
        483.49094768768555,
    ])
    np.testing.assert_allclose(
        isobaric_specific_heat_al_cryosoft(temperature),
        expected_cp, rtol=2e-12, atol=0.0,
    )
    np.testing.assert_allclose(
        density_al_cryosoft(temperature),
        density_al(temperature),
        rtol=0.0, atol=0.0,
    )


def test_legacy_aluminium_functions_remain_separate():
    temperature = np.array([4.5, 20.0, 60.0, 100.0])
    assert not np.allclose(
        thermal_conductivity_al(temperature),
        thermal_conductivity_al_cryosoft(temperature, 17.0, 300.0),
    )
    assert np.max(np.abs(
        electrical_resistivity_al(temperature)
        - electrical_resistivity_al_cryosoft(temperature, 17.0, 300.0)
    )) > 1e-9


def test_strand_stabilizer_dispatches_cryosoft_aluminium_with_b_and_rrr():
    component = object.__new__(StrandStabilizerComponent)
    component.inputs = {"stabilizer_material": "al_cryosoft", "RRR": 300.0}
    prop = {
        "temperature": np.array([4.5, 20.0, 100.0]),
        "B_field": np.array([17.0, 17.0, 17.0]),
    }
    np.testing.assert_allclose(
        component.strand_thermal_conductivity(prop),
        thermal_conductivity_al_cryosoft(
            prop["temperature"], prop["B_field"], 300.0
        ),
    )
    np.testing.assert_allclose(
        component.strand_electrical_resistivity(prop),
        electrical_resistivity_al_cryosoft(
            prop["temperature"], prop["B_field"], 300.0
        ),
    )


def test_jacket_dispatches_cryosoft_aluminium_with_b_and_rrr():
    component = object.__new__(JacketComponent)
    component.inputs = {"NUM_MATERIAL_TYPES": 1, "CROSSECTION": 1.0, "RRR": 300.0}
    component.materials = np.array(["al_cryosoft"])
    component.cross_sections = np.array([1.0])
    component.thermal_conductivity_function = np.array(
        [thermal_conductivity_al_cryosoft]
    )
    component.electrical_resistivity_function = np.array(
        [electrical_resistivity_al_cryosoft]
    )
    prop = {
        "temperature": np.array([4.5, 20.0, 100.0]),
        "B_field": np.array([17.0, 17.0, 17.0]),
    }
    np.testing.assert_allclose(
        component.jacket_thermal_conductivity(prop),
        thermal_conductivity_al_cryosoft(
            prop["temperature"], prop["B_field"], 300.0
        ),
    )
    np.testing.assert_allclose(
        component.jacket_electrical_resistivity(prop),
        electrical_resistivity_al_cryosoft(
            prop["temperature"], prop["B_field"], 300.0
        ),
    )
