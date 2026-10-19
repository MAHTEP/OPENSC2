from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import utility_functions.checkpoint as checkpoint_module
from utility_functions.checkpoint_schedule import (
    CheckpointBoundary,
    CheckpointSchedule,
)


_SCHEDULE_NOT_SET = object()


def make_runtime(
    *,
    step,
    interval,
    current_time,
    schedule=_SCHEDULE_NOT_SET,
):
    runtime = SimpleNamespace(
        transient_input={
            "CHECKPOINT_EVERY_N_STEPS": interval,
        },
        num_step=step,
        simulation_time=[current_time],
        epsilon=1.0e-6,
        dict_path={"Checkpoint_dir": "checkpoint-dir"},
    )
    if schedule is not _SCHEDULE_NOT_SET:
        runtime.checkpoint_schedule = schedule
    return runtime


def schedule_at(time, trigger):
    return CheckpointSchedule(
        user_enabled=True,
        boundaries=(CheckpointBoundary(time=time, trigger=trigger),),
    )


class CheckpointTriggerArbitrationTests(unittest.TestCase):
    def assert_checkpoint_written(self, runtime, expected_trigger):
        expected_path = Path("written-checkpoint.h5")
        with patch.object(
            checkpoint_module,
            "write_checkpoint",
            return_value=expected_path,
        ) as writer:
            result = checkpoint_module.write_checkpoint_if_due(runtime)

        self.assertEqual(result, expected_path)
        writer.assert_called_once_with(
            runtime,
            "checkpoint-dir",
            trigger=expected_trigger,
        )

    def test_requested_boundary_overrides_coincident_periodic_checkpoint(self):
        runtime = make_runtime(
            step=10,
            interval=5,
            current_time=0.25,
            schedule=schedule_at(0.25, "requested"),
        )

        self.assert_checkpoint_written(runtime, "requested")

    def test_final_boundary_overrides_coincident_periodic_checkpoint(self):
        runtime = make_runtime(
            step=10,
            interval=5,
            current_time=1.0,
            schedule=schedule_at(1.0, "final"),
        )

        self.assert_checkpoint_written(runtime, "final")

    def test_periodic_checkpoint_remains_active_before_final_boundary(self):
        runtime = make_runtime(
            step=10,
            interval=5,
            current_time=0.25,
            schedule=schedule_at(1.0, "final"),
        )

        self.assert_checkpoint_written(runtime, "periodic")

    def test_requested_checkpoint_works_when_periodic_mode_is_disabled(self):
        runtime = make_runtime(
            step=11,
            interval=0,
            current_time=0.25,
            schedule=schedule_at(0.25, "requested"),
        )

        self.assert_checkpoint_written(runtime, "requested")

    def test_no_due_trigger_returns_none_without_writing(self):
        runtime = make_runtime(
            step=11,
            interval=5,
            current_time=0.25,
            schedule=schedule_at(1.0, "final"),
        )

        with patch.object(checkpoint_module, "write_checkpoint") as writer:
            result = checkpoint_module.write_checkpoint_if_due(runtime)

        self.assertIsNone(result)
        writer.assert_not_called()

    def test_legacy_runtime_without_schedule_keeps_periodic_behavior(self):
        runtime = make_runtime(
            step=10,
            interval=5,
            current_time=0.25,
        )

        self.assert_checkpoint_written(runtime, "periodic")


if __name__ == "__main__":
    unittest.main()
