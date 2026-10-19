from types import SimpleNamespace

import numpy as np
import pytest

from conductor import Conductor


def _simulation_fixed_dt(dt=0.005):
    return SimpleNamespace(
        transient_input={
            "IADAPTIME": 0,
            "TIME_STEP": dt,
            "STPMIN": dt,
        }
    )


def _bare_conductor():
    # The two methods under test only use their explicit arguments, so a full
    # Conductor initialization is unnecessary and would turn these unit tests
    # into integration tests.
    return object.__new__(Conductor)


def test_aux_event_intervals_exactly_ten_steps_accept_float_roundoff():
    conductor = _bare_conductor()
    simulation = _simulation_fixed_dt(0.005)

    # Nominal spacing is exactly 0.05 s = 10 * dt. With absolute times near
    # 610 s, binary subtraction can produce values infinitesimally below 0.05.
    event_times = np.array(
        [610.65, 610.70, 610.75, 610.80],
        dtype=float,
    )

    conductor._Conductor__check_event_time_aux_input(
        event_times,
        simulation,
        "conductor_current_dump.xlsx",
    )


def test_aux_event_interval_genuinely_shorter_than_ten_steps_is_rejected():
    conductor = _bare_conductor()
    simulation = _simulation_fixed_dt(0.005)

    event_times = np.array(
        [610.65, 610.699],
        dtype=float,
    )

    with pytest.raises(ValueError, match="at least 10"):
        conductor._Conductor__check_event_time_aux_input(
            event_times,
            simulation,
            "conductor_current_dump.xlsx",
        )


def test_main_input_heating_interval_exactly_ten_steps_accept_float_roundoff():
    conductor = _bare_conductor()
    simulation = _simulation_fixed_dt(0.005)

    # Layout expected by __check_event_time_main_input:
    # first half = TQBEG values, second half = TQEND values.
    event_times = [610.70, 610.75]

    conductor._Conductor__check_event_time_main_input(
        event_times,
        simulation,
        "conductor_operation.xlsx",
    )


def test_main_input_heating_interval_genuinely_shorter_is_rejected():
    conductor = _bare_conductor()
    simulation = _simulation_fixed_dt(0.005)

    event_times = [610.65, 610.699]

    with pytest.raises(ValueError, match="at least 10"):
        conductor._Conductor__check_event_time_main_input(
            event_times,
            simulation,
            "conductor_operation.xlsx",
        )
