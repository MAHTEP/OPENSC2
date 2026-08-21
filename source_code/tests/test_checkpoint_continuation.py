from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

from test_checkpoint import make_restore_target, make_simulation
from utility_functions.checkpoint import (
    CheckpointValidationError,
    apply_checkpoint_to_runtime,
    read_checkpoint,
    write_checkpoint,
)


class CheckpointContinuationRestoreTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        self.input_dir = self.root / "input"
        self.input_dir.mkdir()
        (self.input_dir / "transitory_input.xlsx").write_bytes(
            b"input-content"
        )

    def _checkpoint_with_current_inputs(self):
        simulation = make_simulation(self.input_dir)
        checkpoint_path = write_checkpoint(
            simulation,
            self.root / "checkpoints",
            trigger="requested",
        )
        return read_checkpoint(checkpoint_path)

    def test_continuation_keeps_new_time_policy_and_event_timeline(self):
        checkpoint = self._checkpoint_with_current_inputs()
        target = make_restore_target(self.input_dir)
        target.transient_input = {
            "IADAPTIME": 0,
            "TIME_STEP": 0.025,
            "STPMIN": 0.0025,
            "TEND": 0.4,
        }

        conductor = target.list_of_Conductors[0]
        conductor.time_step = 0.025
        conductor.events_time = np.array([0.05, 0.1, 0.15, 0.3])
        conductor.i_event = 0
        conductor.Space_save = np.array([0.0, 0.1, 0.15, 0.3, 0.4])
        conductor.num_step_save = np.zeros(5, dtype=int)
        conductor.i_save_max = 4

        compatibility = SimpleNamespace(
            is_compatible=True,
            blocking_reasons=(),
            warnings=(),
        )

        with patch(
            "utility_functions.checkpoint.evaluate_restart_compatibility",
            return_value=compatibility,
        ):
            apply_checkpoint_to_runtime(
                checkpoint,
                target,
                mode="continuation",
            )

        self.assertEqual(target.simulation_time, [0.0, 0.1])
        self.assertEqual(target.num_step, 1)
        self.assertTrue(target.restored_from_checkpoint)
        self.assertEqual(conductor.time_step, 0.025)
        np.testing.assert_allclose(
            conductor.events_time,
            [0.05, 0.1, 0.15, 0.3],
        )
        self.assertEqual(conductor.i_event, 2)
        np.testing.assert_allclose(conductor.Space_save, [0.15, 0.3, 0.4])
        np.testing.assert_array_equal(conductor.num_step_save, [0, 0, 0])
        self.assertEqual(conductor.i_save, 0)
        self.assertEqual(conductor.i_save_max, 2)

    def test_continuation_accepts_adaptive_time_policy(self):
        checkpoint = self._checkpoint_with_current_inputs()
        target = make_restore_target(self.input_dir)
        target.transient_input = {
            "IADAPTIME": 1,
            "TIME_STEP": 0.025,
            "STPMIN": 0.0025,
            "STPMAX": 0.05,
            "MLT_INCREASE": 1.2,
            "MLT_DECREASE": 0.5,
            "TEND": 0.4,
        }
        conductor = target.list_of_Conductors[0]
        conductor.time_step = 0.025
        conductor.i_event = 0
        conductor.Space_save = np.array([0.0, 0.1, 0.2, 0.4])
        conductor.num_step_save = np.zeros(4, dtype=int)
        conductor.i_save_max = 3

        compatibility = SimpleNamespace(
            is_compatible=True,
            blocking_reasons=(),
            warnings=(),
        )

        with patch(
            "utility_functions.checkpoint.evaluate_restart_compatibility",
            return_value=compatibility,
        ):
            apply_checkpoint_to_runtime(
                checkpoint,
                target,
                mode="continuation",
            )

        self.assertEqual(target.simulation_time, [0.0, 0.1])
        self.assertEqual(target.num_step, 1)
        self.assertTrue(target.restored_from_checkpoint)
        self.assertEqual(conductor.time_step, 0.025)
        self.assertEqual(conductor.cond_time, [0.0, 0.1])
        self.assertEqual(conductor.cond_num_step, 1)
        np.testing.assert_allclose(conductor.Space_save, [0.2, 0.4])
        np.testing.assert_array_equal(conductor.num_step_save, [0, 0])

    def test_continuation_accepts_a_different_spatial_output_schedule(self):
        checkpoint = self._checkpoint_with_current_inputs()
        target = make_restore_target(self.input_dir)
        target.transient_input = {
            "IADAPTIME": 0,
            "TIME_STEP": 0.0025,
            "STPMIN": 0.0025,
            "TEND": 0.4,
        }
        conductor = target.list_of_Conductors[0]
        conductor.time_step = 0.0025
        conductor.i_event = 0
        conductor.Space_save = np.array(
            [0.0, 0.1, 0.15, 0.2, 0.25, 0.4]
        )
        conductor.num_step_save = np.zeros(6, dtype=np.int64)
        conductor.i_save_max = 5

        compatibility = SimpleNamespace(
            is_compatible=True,
            blocking_reasons=(),
            warnings=(),
        )

        with patch(
            "utility_functions.checkpoint.evaluate_restart_compatibility",
            return_value=compatibility,
        ):
            apply_checkpoint_to_runtime(
                checkpoint,
                target,
                mode="continuation",
            )

        np.testing.assert_allclose(
            conductor.Space_save,
            [0.15, 0.2, 0.25, 0.4],
        )
        np.testing.assert_array_equal(
            conductor.num_step_save,
            np.zeros(4, dtype=np.int64),
        )
        self.assertEqual(conductor.i_save, 0)
        self.assertEqual(conductor.i_save_max, 3)
        self.assertTrue(target.restored_from_checkpoint)

    def test_continuation_rejects_output_schedule_without_future_times(self):
        checkpoint = self._checkpoint_with_current_inputs()
        target = make_restore_target(self.input_dir)
        target.transient_input = {
            "IADAPTIME": 0,
            "TIME_STEP": 0.025,
            "STPMIN": 0.0025,
            "TEND": 0.4,
        }
        conductor = target.list_of_Conductors[0]
        conductor.i_event = 0
        conductor.Space_save = np.array([0.0, 0.1])

        compatibility = SimpleNamespace(
            is_compatible=True,
            blocking_reasons=(),
            warnings=(),
        )

        with patch(
            "utility_functions.checkpoint.evaluate_restart_compatibility",
            return_value=compatibility,
        ):
            with self.assertRaisesRegex(
                CheckpointValidationError,
                "Space_save contains no time after the checkpoint",
            ):
                apply_checkpoint_to_runtime(
                    checkpoint,
                    target,
                    mode="continuation",
                )

        self.assertFalse(target.restored_from_checkpoint)

    def test_continuation_rejects_end_time_at_checkpoint_before_mutation(self):
        checkpoint = self._checkpoint_with_current_inputs()
        target = make_restore_target(self.input_dir)
        target.transient_input = {
            "IADAPTIME": 0,
            "TIME_STEP": 0.025,
            "STPMIN": 0.0025,
            "TEND": 0.1,
        }
        conductor = target.list_of_Conductors[0]
        conductor.i_event = 0

        compatibility = SimpleNamespace(
            is_compatible=True,
            blocking_reasons=(),
            warnings=(),
        )

        with patch(
            "utility_functions.checkpoint.evaluate_restart_compatibility",
            return_value=compatibility,
        ):
            with self.assertRaisesRegex(
                CheckpointValidationError,
                "TEND.*checkpoint",
            ):
                apply_checkpoint_to_runtime(
                    checkpoint,
                    target,
                    mode="continuation",
                )

        self.assertEqual(target.simulation_time, [0.0])
        self.assertEqual(target.num_step, 0)
        self.assertFalse(target.restored_from_checkpoint)
        self.assertEqual(conductor.cond_time, [0.0])
        self.assertEqual(conductor.cond_num_step, 0)

    def test_recovery_still_restores_saved_clock_policy(self):
        checkpoint = self._checkpoint_with_current_inputs()
        target = make_restore_target(self.input_dir)
        conductor = target.list_of_Conductors[0]
        conductor.time_step = 0.025

        apply_checkpoint_to_runtime(checkpoint, target)

        self.assertEqual(conductor.time_step, 0.1)
        self.assertEqual(conductor.i_event, 1)
        np.testing.assert_allclose(conductor.events_time, [0.01, 1.0])


if __name__ == "__main__":
    unittest.main()
