from types import SimpleNamespace
from unittest.mock import Mock, patch

import numpy as np

from simulation import _save_spatial_distribution_if_due


def _conductor(target_time, current_time, time_step=0.1):
    return SimpleNamespace(
        Space_save=np.asarray([0.0, target_time]),
        i_save=1,
        i_save_max=2,
        cond_time=[0.0, current_time],
        time_step=time_step,
        store_spatial_distributions=Mock(),
        store_interp_spatial_distributions=Mock(),
    )


def test_exact_requested_time_is_stored_without_interpolation():
    conductor = _conductor(
        target_time=0.3,
        current_time=0.1 + 0.2,
    )

    with patch("simulation.save_simulation_space") as save_space:
        _save_spatial_distribution_if_due(
            conductor,
            "spatial-output",
        )

    conductor.store_spatial_distributions.assert_called_once_with("t_save")
    conductor.store_interp_spatial_distributions.assert_not_called()
    save_space.assert_called_once_with(conductor, "spatial-output")
    assert conductor.i_save == 2
    assert not hasattr(conductor, "t_save_left")
    assert not hasattr(conductor, "t_save_right")


def test_requested_time_approached_from_left_stores_left_endpoint():
    conductor = _conductor(
        target_time=1.0,
        current_time=0.95,
    )

    with patch("simulation.save_simulation_space") as save_space:
        _save_spatial_distribution_if_due(
            conductor,
            "spatial-output",
        )

    assert conductor.t_save_left == 0.95
    conductor.store_spatial_distributions.assert_called_once_with()
    conductor.store_interp_spatial_distributions.assert_not_called()
    save_space.assert_not_called()
    assert conductor.i_save == 1


def test_crossed_requested_time_uses_interpolation():
    conductor = _conductor(
        target_time=1.0,
        current_time=1.05,
    )
    conductor.t_save_left = 0.95

    with patch("simulation.save_simulation_space") as save_space:
        _save_spatial_distribution_if_due(
            conductor,
            "spatial-output",
        )

    assert conductor.t_save_right == 1.05
    conductor.store_spatial_distributions.assert_not_called()
    conductor.store_interp_spatial_distributions.assert_called_once_with()
    save_space.assert_called_once_with(conductor, "spatial-output")
    assert conductor.i_save == 2


def test_completed_schedule_is_ignored_without_indexing_space_save():
    conductor = SimpleNamespace(
        Space_save=np.asarray([0.0]),
        i_save=1,
        i_save_max=1,
    )

    with patch("simulation.save_simulation_space") as save_space:
        _save_spatial_distribution_if_due(
            conductor,
            "spatial-output",
        )

    save_space.assert_not_called()
