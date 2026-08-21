from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from test_checkpoint import make_restore_target, make_simulation
from utility_functions.checkpoint import (
    CheckpointValidationError,
    apply_checkpoint_to_runtime,
    build_continuation_profile,
    evaluate_restart_compatibility,
    read_checkpoint,
    write_checkpoint,
)


class CheckpointContinuationAdaptivePolicyTests(unittest.TestCase):
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
    def _time_policy(iadaptime):
        return {
            "IADAPTIME": iadaptime,
            "TIME_STEP": 0.025,
            "STPMIN": 0.0025,
            "STPMAX": 0.05,
            "MLT_INCREASE": 1.2,
            "MLT_DECREASE": 0.5,
            "TIMEREF": 0.2,
            "TAUREF": 0.05,
            "TEND": 0.4,
        }

    def _checkpoint(self, iadaptime=0):
        simulation = make_simulation(self.input_dir)
        simulation.transient_input = self._time_policy(iadaptime)
        checkpoint_path = write_checkpoint(
            simulation,
            self.root / "checkpoints",
            trigger="requested",
        )
        return read_checkpoint(checkpoint_path)

    def test_profile_classifies_iadaptime_as_time_policy(self):
        simulation = make_simulation(self.input_dir)
        simulation.transient_input = self._time_policy(0)

        profile = build_continuation_profile(simulation)

        self.assertNotIn("IADAPTIME", profile.immutable)
        self.assertEqual(profile.time_policy["IADAPTIME"], 0)

    def test_fixed_to_adaptive_change_is_allowed_and_reported(self):
        checkpoint = self._checkpoint(iadaptime=0)
        target = make_restore_target(self.input_dir)
        target.transient_input = self._time_policy(1)
        runtime_profile = build_continuation_profile(target)

        report = evaluate_restart_compatibility(
            checkpoint,
            self.input_dir,
            mode="continuation",
            runtime_profile=runtime_profile,
        )

        self.assertTrue(report.is_compatible)
        self.assertEqual(report.blocking_reasons, ())
        self.assertEqual(
            report.continuation_comparison.immutable_differences,
            (),
        )
        self.assertEqual(
            report.continuation_comparison.time_policy_differences,
            ("time_policy.IADAPTIME",),
        )
        self.assertEqual(
            report.warnings,
            (
                "Continuation time-policy input differs: "
                "time_policy.IADAPTIME.",
            ),
        )

    def test_restore_keeps_fresh_adaptive_policy_and_restores_state(self):
        checkpoint = self._checkpoint(iadaptime=0)
        target = make_restore_target(self.input_dir)
        target.transient_input = self._time_policy(1)
        conductor = target.list_of_Conductors[0]
        conductor.time_step = target.transient_input["TIME_STEP"]
        conductor.events_time = np.array([0.05, 0.1, 0.15, 0.4])
        conductor.i_event = 0
        conductor.Space_save = np.array([0.0, 0.1, 0.2, 0.4])
        conductor.num_step_save = np.zeros(4, dtype=int)
        conductor.i_save_max = 3
        conductor.EQTEIG[:] = -1.0
        conductor.dict_Step["SYSVAR"][:] = -1.0

        try:
            apply_checkpoint_to_runtime(
                checkpoint,
                target,
                mode="continuation",
            )
        except CheckpointValidationError as exc:
            self.fail(
                "A valid fixed-to-adaptive continuation was rejected: "
                f"{exc}"
            )

        self.assertTrue(target.restored_from_checkpoint)
        self.assertEqual(target.transient_input["IADAPTIME"], 1)
        self.assertEqual(target.simulation_time, [0.0, 0.1])
        self.assertEqual(target.num_step, 1)
        self.assertEqual(conductor.time_step, 0.025)
        np.testing.assert_allclose(
            conductor.events_time,
            [0.05, 0.1, 0.15, 0.4],
        )
        self.assertEqual(conductor.i_event, 2)
        np.testing.assert_allclose(conductor.Space_save, [0.2, 0.4])
        np.testing.assert_array_equal(conductor.num_step_save, [0, 0])
        np.testing.assert_allclose(conductor.EQTEIG, [1.0, 2.0])
        np.testing.assert_allclose(
            conductor.dict_Step["SYSVAR"],
            np.arange(16.0).reshape(8, 2),
        )

    def test_adaptive_to_fixed_change_is_allowed_and_reported(self):
        checkpoint = self._checkpoint(iadaptime=1)
        target = make_restore_target(self.input_dir)
        target.transient_input = self._time_policy(0)

        report = evaluate_restart_compatibility(
            checkpoint,
            self.input_dir,
            mode="continuation",
            runtime_profile=build_continuation_profile(target),
        )

        self.assertTrue(report.is_compatible)
        self.assertEqual(report.blocking_reasons, ())
        self.assertEqual(
            report.continuation_comparison.time_policy_differences,
            ("time_policy.IADAPTIME",),
        )

    def test_adaptive_mode_change_is_allowed_and_reported(self):
        checkpoint = self._checkpoint(iadaptime=1)
        target = make_restore_target(self.input_dir)
        target.transient_input = self._time_policy(2)

        report = evaluate_restart_compatibility(
            checkpoint,
            self.input_dir,
            mode="continuation",
            runtime_profile=build_continuation_profile(target),
        )

        self.assertTrue(report.is_compatible)
        self.assertEqual(report.blocking_reasons, ())
        self.assertEqual(
            report.continuation_comparison.time_policy_differences,
            ("time_policy.IADAPTIME",),
        )

    def _assert_policy_rejected_before_mutation(
        self,
        transient_input,
        message,
    ):
        checkpoint = self._checkpoint(iadaptime=0)
        target = make_restore_target(self.input_dir)
        target.transient_input = transient_input
        conductor = target.list_of_Conductors[0]

        original_time = list(target.simulation_time)
        original_step = target.num_step
        original_time_step = conductor.time_step
        original_events = conductor.events_time.copy()
        original_sysvar = conductor.dict_Step["SYSVAR"].copy()

        with self.assertRaisesRegex(CheckpointValidationError, message):
            apply_checkpoint_to_runtime(
                checkpoint,
                target,
                mode="continuation",
            )

        self.assertFalse(target.restored_from_checkpoint)
        self.assertEqual(target.simulation_time, original_time)
        self.assertEqual(target.num_step, original_step)
        self.assertEqual(conductor.time_step, original_time_step)
        np.testing.assert_array_equal(
            conductor.events_time,
            original_events,
        )
        np.testing.assert_array_equal(
            conductor.dict_Step["SYSVAR"],
            original_sysvar,
        )

    def test_invalid_iadaptime_is_rejected_before_mutation(self):
        policy = self._time_policy(3)

        self._assert_policy_rejected_before_mutation(
            policy,
            r"Continuation IADAPTIME must be one of",
        )

    def test_unimplemented_iadaptime_is_rejected_before_mutation(self):
        policy = self._time_policy(-1)

        self._assert_policy_rejected_before_mutation(
            policy,
            r"Continuation IADAPTIME=-1 is not implemented",
        )

    def test_adaptive_policy_requires_stpmax_before_mutation(self):
        policy = self._time_policy(1)
        del policy["STPMAX"]

        self._assert_policy_rejected_before_mutation(
            policy,
            r"Continuation STPMAX must be a finite numeric scalar",
        )

    def test_adaptive_policy_rejects_inverted_bounds_before_mutation(self):
        policy = self._time_policy(1)
        policy["STPMIN"] = 0.06

        self._assert_policy_rejected_before_mutation(
            policy,
            r"Continuation STPMIN must not exceed STPMAX",
        )

    def test_adaptive_policy_rejects_initial_step_outside_bounds(self):
        policy = self._time_policy(2)
        policy["TIME_STEP"] = 0.06

        self._assert_policy_rejected_before_mutation(
            policy,
            r"Continuation TIME_STEP must lie between STPMIN and STPMAX",
        )

    def test_adaptive_policy_rejects_nonfinite_multiplier(self):
        policy = self._time_policy(1)
        policy["MLT_INCREASE"] = np.inf

        self._assert_policy_rejected_before_mutation(
            policy,
            r"Continuation MLT_INCREASE must be a finite numeric scalar",
        )


if __name__ == "__main__":
    unittest.main()
