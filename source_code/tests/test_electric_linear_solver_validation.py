import numpy as np
import pytest
from scipy import sparse

import utility_functions.electric_auxiliary_functions as electric_aux


def _identity_system():
    matrix = sparse.eye(2, format="csr")
    right_hand_side = np.array([1.0, -2.0])
    return matrix, right_hand_side


def test_valid_pardiso_solution_is_returned_without_reset(monkeypatch):
    matrix, right_hand_side = _identity_system()
    expected = right_hand_side.copy()
    solve_calls = []
    reset_calls = []

    def fake_solve(received_matrix, received_rhs):
        solve_calls.append((received_matrix, received_rhs))
        return expected.copy()

    monkeypatch.setattr(electric_aux, "pardiso_spsolve", fake_solve)
    monkeypatch.setattr(
        electric_aux.pypardiso_solver,
        "free_memory",
        lambda **kwargs: reset_calls.append(kwargs),
    )

    result = electric_aux._solve_electric_linear_system(
        matrix,
        right_hand_side,
        context="unit test",
    )

    np.testing.assert_array_equal(result, expected)
    assert len(solve_calls) == 1
    assert solve_calls[0][0] is matrix
    assert solve_calls[0][1] is right_hand_side
    assert reset_calls == []


def test_bad_solution_is_retried_on_identical_system_after_reset(
    monkeypatch,
):
    matrix, right_hand_side = _identity_system()
    candidates = iter(
        (
            np.array([1.0e12, -1.0e12]),
            right_hand_side.copy(),
        )
    )
    solve_calls = []
    reset_calls = []
    warning_messages = []

    def fake_solve(received_matrix, received_rhs):
        solve_calls.append((received_matrix, received_rhs))
        return next(candidates)

    monkeypatch.setattr(electric_aux, "pardiso_spsolve", fake_solve)
    monkeypatch.setattr(
        electric_aux.pypardiso_solver,
        "free_memory",
        lambda **kwargs: reset_calls.append(kwargs),
    )
    monkeypatch.setattr(
        electric_aux.LOGGER,
        "warning",
        lambda message, *args: warning_messages.append(message % args),
    )

    result = electric_aux._solve_electric_linear_system(
        matrix,
        right_hand_side,
        context="TH step=1001, electric step=4",
    )

    np.testing.assert_array_equal(result, right_hand_side)
    assert len(solve_calls) == 2
    assert all(call[0] is matrix for call in solve_calls)
    assert all(call[1] is right_hand_side for call in solve_calls)
    assert reset_calls == [{"everything": True}]
    assert len(warning_messages) == 2
    assert "unacceptable electric solution" in warning_messages[0]
    assert "recovered after one full" in warning_messages[1]


def test_nonfinite_solution_is_retried(monkeypatch):
    matrix, right_hand_side = _identity_system()
    candidates = iter(
        (
            np.array([np.nan, 0.0]),
            right_hand_side.copy(),
        )
    )
    reset_calls = []

    monkeypatch.setattr(
        electric_aux,
        "pardiso_spsolve",
        lambda matrix, rhs: next(candidates),
    )
    monkeypatch.setattr(
        electric_aux.pypardiso_solver,
        "free_memory",
        lambda **kwargs: reset_calls.append(kwargs),
    )

    result = electric_aux._solve_electric_linear_system(
        matrix,
        right_hand_side,
        context="unit test",
    )

    np.testing.assert_array_equal(result, right_hand_side)
    assert reset_calls == [{"everything": True}]


def test_all_bad_solutions_fail_before_returning_invalid_state(monkeypatch):
    matrix, right_hand_side = _identity_system()
    bad_solution = np.array([1.0e12, -1.0e12])
    reset_calls = []

    monkeypatch.setattr(
        electric_aux,
        "pardiso_spsolve",
        lambda matrix, rhs: bad_solution.copy(),
    )
    monkeypatch.setattr(
        electric_aux.pypardiso_solver,
        "free_memory",
        lambda **kwargs: reset_calls.append(kwargs),
    )
    monkeypatch.setattr(
        electric_aux.sparse.linalg,
        "spsolve",
        lambda matrix, rhs: bad_solution.copy(),
    )

    with pytest.raises(RuntimeError) as error:
        electric_aux._solve_electric_linear_system(
            matrix,
            right_hand_side,
            context="unit-test",
        )

    message = str(error.value)
    assert "SciPy spsolve fallback" in message
    assert "The invalid electric solution was not applied" in message
    assert reset_calls == [
        {"everything": True},
        {"everything": True},
    ]


def test_zero_rhs_exact_solution_passes_validation():
    matrix = sparse.eye(2, format="csr")
    right_hand_side = np.zeros(2)
    solution = np.zeros(2)

    diagnostics = electric_aux._electric_linear_solution_diagnostics(
        matrix,
        right_hand_side,
        solution,
    )

    assert diagnostics["acceptable"] is True
    assert diagnostics["rhs_relative_residual"] == 0.0
