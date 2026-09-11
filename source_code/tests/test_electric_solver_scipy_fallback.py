import numpy as np
import pytest
from scipy import sparse

import utility_functions.electric_auxiliary_functions as eaf


def _simple_system():
    matrix = sparse.csr_matrix(
        np.array(
            [
                [2.0, 0.0],
                [0.0, 4.0],
            ]
        )
    )
    rhs = np.array([2.0, 8.0])
    exact = np.array([1.0, 2.0])
    return matrix, rhs, exact


def test_scipy_fallback_recovers_after_two_invalid_pardiso_solutions(monkeypatch):
    matrix, rhs, exact = _simple_system()

    pardiso_calls = []
    free_memory_calls = []
    scipy_calls = []

    def fake_pardiso(matrix_arg, rhs_arg):
        pardiso_calls.append(1)
        return np.zeros_like(rhs_arg, dtype=float)

    def fake_free_memory(*, everything):
        free_memory_calls.append(everything)

    def fake_scipy_spsolve(matrix_arg, rhs_arg):
        scipy_calls.append(1)
        return exact.copy()

    monkeypatch.setattr(eaf, "pardiso_spsolve", fake_pardiso)
    monkeypatch.setattr(
        eaf.pypardiso_solver,
        "free_memory",
        fake_free_memory,
    )
    monkeypatch.setattr(
        eaf.sparse.linalg,
        "spsolve",
        fake_scipy_spsolve,
    )

    solution = eaf._solve_electric_linear_system(
        matrix,
        rhs,
        context="unit-test",
    )

    np.testing.assert_allclose(solution, exact)
    assert len(pardiso_calls) == 2
    assert free_memory_calls == [True, True]
    assert len(scipy_calls) == 1


def test_scipy_fallback_is_rejected_if_residual_is_also_invalid(monkeypatch):
    matrix, rhs, _ = _simple_system()

    free_memory_calls = []

    def fake_pardiso(matrix_arg, rhs_arg):
        return np.zeros_like(rhs_arg, dtype=float)

    def fake_free_memory(*, everything):
        free_memory_calls.append(everything)

    def fake_scipy_spsolve(matrix_arg, rhs_arg):
        return np.zeros_like(rhs_arg, dtype=float)

    monkeypatch.setattr(eaf, "pardiso_spsolve", fake_pardiso)
    monkeypatch.setattr(
        eaf.pypardiso_solver,
        "free_memory",
        fake_free_memory,
    )
    monkeypatch.setattr(
        eaf.sparse.linalg,
        "spsolve",
        fake_scipy_spsolve,
    )

    with pytest.raises(
        RuntimeError,
        match="SciPy spsolve fallback",
    ):
        eaf._solve_electric_linear_system(
            matrix,
            rhs,
            context="unit-test",
        )

    assert free_memory_calls == [True, True]
