from types import SimpleNamespace

import numpy as np
import pytest
from scipy.optimize import brentq

from stack_component import StackComponent
from strand_mixed_component import StrandMixedComponent


E0 = 1.0e-4
N_VALUE = 15
STABILIZER_AREA = 1.784e-6
STABILIZER_RESISTIVITY = 2.11e-9
TOTAL_CURRENT = 809.708


def build_component(component_class=StackComponent):
    component = component_class.__new__(component_class)
    component.inputs = {"E0": E0, "nn": N_VALUE}
    component.cross_section = {"stab": STABILIZER_AREA}
    component.identifier = component_class.__name__
    return component


def solve(critical_current):
    critical_current = np.asarray(critical_current, dtype=float)
    rho = np.full_like(critical_current, STABILIZER_RESISTIVITY)
    current = np.full_like(critical_current, TOTAL_CURRENT)
    return build_component().solve_current_divider(
        rho,
        critical_current,
        current,
    )


def test_dimensionless_solver_matches_independent_well_conditioned_roots():
    critical_current = np.array([100.0, 500.0, 1000.0, 2000.0])
    sc_current, stab_current = solve(critical_current)

    reference = np.array(
        [
            brentq(
                lambda candidate: (
                    E0 * (candidate / critical) ** N_VALUE
                    - STABILIZER_RESISTIVITY
                    * (TOTAL_CURRENT - candidate)
                    / STABILIZER_AREA
                ),
                0.0,
                TOTAL_CURRENT,
            )
            for critical in critical_current
        ]
    )

    np.testing.assert_allclose(sc_current, reference, rtol=1.0e-11)
    np.testing.assert_allclose(sc_current + stab_current, TOTAL_CURRENT)


def test_dimensionless_solver_survives_critical_current_power_underflow():
    critical_current = np.array([1.0e-21])
    dimensional_psi = (
        STABILIZER_RESISTIVITY
        * critical_current**N_VALUE
        / E0
        / STABILIZER_AREA
    )
    assert dimensional_psi[0] == 0.0

    sc_current, stab_current = solve(critical_current)
    v_sc, v_stab = build_component()._evaluate_current_divider_voltage(
        np.array([STABILIZER_RESISTIVITY]),
        critical_current,
        sc_current,
        stab_current,
    )

    assert np.all(np.isfinite(sc_current))
    assert sc_current[0] > 0.0
    np.testing.assert_allclose(sc_current + stab_current, TOTAL_CURRENT)
    np.testing.assert_allclose(v_sc, v_stab, rtol=1.0e-12, atol=0.0)


def test_solution_converges_continuously_to_stabilizer_only_normal_state():
    positive_critical_current = np.logspace(-2, -30, 29)
    sc_current, stab_current = solve(
        np.concatenate((positive_critical_current, [0.0]))
    )

    assert np.all(np.diff(sc_current) <= 0.0)
    assert np.all(np.diff(stab_current) >= 0.0)
    assert sc_current[-1] == 0.0
    assert stab_current[-1] == TOTAL_CURRENT
    np.testing.assert_allclose(sc_current + stab_current, TOTAL_CURRENT)


def test_solution_is_symmetric_for_reversed_current():
    component = build_component()
    critical_current = np.array([1000.0])
    rho = np.array([STABILIZER_RESISTIVITY])
    positive_current = np.array([TOTAL_CURRENT])
    positive_sc, positive_stab = component.solve_current_divider(
        rho,
        critical_current,
        positive_current,
    )
    negative_sc, negative_stab = component.solve_current_divider(
        rho,
        critical_current,
        -positive_current,
    )

    np.testing.assert_allclose(negative_sc, -positive_sc)
    np.testing.assert_allclose(negative_stab, -positive_stab)
    negative_v_sc, negative_v_stab = (
        component._evaluate_current_divider_voltage(
            rho,
            critical_current,
            negative_sc,
            negative_stab,
        )
    )
    np.testing.assert_allclose(negative_v_sc, negative_v_stab)


@pytest.mark.parametrize(
    ("component_class", "extra_arguments"),
    (
        (StackComponent, (1,)),
        (StrandMixedComponent, ()),
    ),
)
def test_both_component_types_use_stable_current_divider(
    component_class,
    extra_arguments,
):
    component = build_component(component_class)
    rho = np.array([STABILIZER_RESISTIVITY])
    critical_current = np.array([1.0e-21])
    current = np.array([TOTAL_CURRENT])

    sc_current, stab_current = component.solve_current_divider(
        rho,
        critical_current,
        current,
        *extra_arguments,
    )

    np.testing.assert_allclose(sc_current + stab_current, current)
    v_sc, v_stab = component._evaluate_current_divider_voltage(
        rho,
        critical_current,
        sc_current,
        stab_current,
    )
    np.testing.assert_allclose(v_sc, v_stab, rtol=1.0e-12, atol=0.0)


@pytest.mark.parametrize("component_class", (StackComponent, StrandMixedComponent))
def test_normal_region_at_gauss_index_zero_is_not_skipped(component_class):
    component = component_class.__new__(component_class)
    component.cross_section = {"sc": 1.0, "stab": STABILIZER_AREA}
    component.dict_Gauss_pt = {
        "J_critical": np.array([0.0]),
        "temperature": np.array([70.0]),
        "electrical_resistivity_stabilizer": np.array(
            [STABILIZER_RESISTIVITY]
        ),
    }
    component.electric_resistance = (
        lambda conductor, resistivity_key, cross_section_key, indices: np.array(
            [42.0]
        )
    )
    conductor = SimpleNamespace(cond_num_step=0, grid_input={"NELEMS": 1})

    resistance = component.get_electric_resistance(conductor)

    np.testing.assert_allclose(resistance.astype(float), [42.0])


def test_voltage_check_reports_actionable_worst_point_diagnostics():
    component = build_component()
    component.identifier = "STACK_1"
    with pytest.raises(ValueError) as error:
        component._assert_current_divider_voltage_balance(
            np.array([STABILIZER_RESISTIVITY]),
            np.array([1.0e-21]),
            np.array([TOTAL_CURRENT]),
            np.array([0.0]),
            np.array([TOTAL_CURRENT]),
            gauss_indices=np.array([734]),
            temperature=np.array([66.251964861917]),
            thermal_hydraulic_time=611.25,
            electric_time=611.24775,
        )

    message = str(error.value)
    assert "component = 'STACK_1'" in message
    assert "gauss_index = 734" in message
    assert "thermal_hydraulic_time = 611.25 s" in message
    assert "critical_current" in message
    assert "dimensionless_residual" in message


@pytest.mark.parametrize(
    "invalid_critical_current",
    (np.array([-1.0]), np.array([np.nan])),
)
def test_solver_rejects_invalid_critical_current(invalid_critical_current):
    with pytest.raises(ValueError):
        build_component().solve_current_divider(
            np.array([STABILIZER_RESISTIVITY]),
            invalid_critical_current,
            np.array([TOTAL_CURRENT]),
        )
