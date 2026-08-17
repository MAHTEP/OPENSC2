import copy
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from test_checkpoint import make_restore_target, make_simulation
from utility_functions.checkpoint import (
    apply_checkpoint_to_runtime,
    evaluate_restart_compatibility,
    read_checkpoint,
    write_checkpoint,
)


class CheckpointContinuationCompatibilityIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        self.input_dir = self.root / "input"
        self.input_dir.mkdir()
        (self.input_dir / "transitory_input.xlsx").write_bytes(
            b"input-content"
        )

    @staticmethod
    def _transient_input():
        return {
            "IADAPTIME": 0,
            "TIME_STEP": 0.1,
            "STPMIN": 0.01,
            "TEND": 1.0,
        }

    def _checkpoint(self):
        simulation = make_simulation(self.input_dir)
        simulation.transient_input = self._transient_input()
        checkpoint_path = write_checkpoint(
            simulation,
            self.root / "checkpoints",
            trigger="requested",
        )
        return read_checkpoint(checkpoint_path)

    def test_identical_profile_allows_explicit_continuation(self):
        checkpoint = self._checkpoint()
        runtime_profile = copy.deepcopy(checkpoint.continuation_profile)

        report = evaluate_restart_compatibility(
            checkpoint,
            self.input_dir,
            mode="continuation",
            runtime_profile=runtime_profile,
        )

        self.assertTrue(report.is_compatible)
        self.assertEqual(report.mode, "continuation")
        self.assertTrue(report.manifest_comparison.is_match)
        self.assertEqual(report.blocking_reasons, ())
        self.assertEqual(report.warnings, ())
        self.assertIsNotNone(report.continuation_comparison)
        self.assertTrue(report.continuation_comparison.is_compatible)
        self.assertEqual(
            report.continuation_comparison.immutable_differences,
            (),
        )
        self.assertEqual(
            report.continuation_comparison.time_policy_differences,
            (),
        )
        self.assertEqual(
            report.continuation_comparison.driver_differences,
            (),
        )

    def test_allowed_changes_are_reported_without_blocking(self):
        checkpoint = self._checkpoint()
        runtime_profile = copy.deepcopy(checkpoint.continuation_profile)
        runtime_profile.time_policy["TIME_STEP"] = 0.025
        runtime_profile.time_policy["TEND"] = 2.0
        runtime_profile.drivers["COND_1"] = {
            "current": {
                "source": {"kind": "canonical_input"},
                "parameters": {
                    "I0_OP_MODE": 0,
                    "I0_OP_TOT": 15000.0,
                },
                "components": {},
            },
        }

        report = evaluate_restart_compatibility(
            checkpoint,
            self.input_dir,
            mode="continuation",
            runtime_profile=runtime_profile,
        )

        self.assertTrue(report.is_compatible)
        self.assertEqual(report.blocking_reasons, ())
        self.assertEqual(
            report.continuation_comparison.time_policy_differences,
            (
                "time_policy.TEND",
                "time_policy.TIME_STEP",
            ),
        )
        self.assertTrue(report.continuation_comparison.driver_differences)
        warning_text = "\n".join(report.warnings)
        self.assertIn("time_policy.TIME_STEP", warning_text)
        self.assertIn("drivers.COND_1.current", warning_text)

    def test_immutable_change_blocks_continuation(self):
        checkpoint = self._checkpoint()
        runtime_profile = copy.deepcopy(checkpoint.continuation_profile)
        runtime_profile.immutable["IADAPTIME"] = 1

        report = evaluate_restart_compatibility(
            checkpoint,
            self.input_dir,
            mode="continuation",
            runtime_profile=runtime_profile,
        )

        self.assertFalse(report.is_compatible)
        self.assertFalse(report.continuation_comparison.is_compatible)
        self.assertEqual(
            report.continuation_comparison.immutable_differences,
            ("immutable.IADAPTIME",),
        )
        self.assertIn(
            "immutable.IADAPTIME",
            "\n".join(report.blocking_reasons),
        )

    def test_missing_runtime_profile_blocks_without_fallback(self):
        checkpoint = self._checkpoint()

        report = evaluate_restart_compatibility(
            checkpoint,
            self.input_dir,
            mode="continuation",
        )

        self.assertFalse(report.is_compatible)
        self.assertIsNone(report.continuation_comparison)
        self.assertIn(
            "current continuation profile",
            "\n".join(report.blocking_reasons).lower(),
        )

    def test_apply_uses_continuation_semantics_without_recovery_fallback(self):
        checkpoint = self._checkpoint()
        target = make_restore_target(self.input_dir)
        target.transient_input = self._transient_input()
        conductor = target.list_of_Conductors[0]
        conductor.events_time = np.array([0.05, 0.1, 0.2, 1.0])
        conductor.i_event = 0

        apply_checkpoint_to_runtime(
            checkpoint,
            target,
            mode="continuation",
        )

        self.assertTrue(target.restored_from_checkpoint)
        self.assertEqual(target.simulation_time, [0.0, 0.1])
        self.assertEqual(target.num_step, 1)
        self.assertEqual(conductor.time_step, 0.1)
        np.testing.assert_allclose(
            conductor.events_time,
            [0.05, 0.1, 0.2, 1.0],
        )
        self.assertEqual(conductor.i_event, 2)


if __name__ == "__main__":
    unittest.main()
