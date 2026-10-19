from dataclasses import FrozenInstanceError
import unittest

import numpy as np

from utility_functions.checkpoint_schedule import (
    CheckpointBoundary,
    build_checkpoint_schedule,
    next_checkpoint_boundary,
    parse_user_checkpoints_flag,
)


class CheckpointScheduleTests(unittest.TestCase):
    def build(
        self,
        transient_input,
        times=(),
        *,
        sheet_present=True,
        epsilon=1.0e-6,
    ):
        return build_checkpoint_schedule(
            transient_input,
            times,
            sheet_present=sheet_present,
            epsilon=epsilon,
        )

    def test_missing_flag_preserves_periodic_only_backward_compatibility(self):
        schedule = self.build({}, sheet_present=False)

        self.assertFalse(schedule.user_enabled)
        self.assertEqual(schedule.boundaries, ())

    def test_flag_accepts_only_explicit_true_and_false_values(self):
        accepted = (
            (True, True),
            (False, False),
            (np.bool_(True), True),
            (np.bool_(False), False),
            ("TRUE", True),
            (" true ", True),
            ("FALSE", False),
            (" false ", False),
        )
        for raw_value, expected in accepted:
            with self.subTest(raw_value=raw_value):
                self.assertIs(
                    parse_user_checkpoints_flag(
                        {"USER_CHECKPOINTS": raw_value}
                    ),
                    expected,
                )

        rejected = (0, 1, "ON", "OFF", "yes", "no", None, np.nan)
        for raw_value in rejected:
            with self.subTest(raw_value=raw_value):
                with self.assertRaisesRegex(ValueError, "USER_CHECKPOINTS"):
                    parse_user_checkpoints_flag(
                        {"USER_CHECKPOINTS": raw_value}
                    )

    def test_disabled_mode_ignores_nonempty_sheet_with_warning(self):
        with self.assertWarnsRegex(UserWarning, "CHECKPOINTS.*ignored"):
            schedule = self.build(
                {"USER_CHECKPOINTS": False, "TEND": 0.65},
                (0.2, 0.4),
            )

        self.assertFalse(schedule.user_enabled)
        self.assertEqual(schedule.boundaries, ())

    def test_enabled_mode_requires_checkpoints_sheet(self):
        with self.assertRaisesRegex(ValueError, "CHECKPOINTS.*required"):
            self.build(
                {"USER_CHECKPOINTS": True, "TEND": 0.65},
                sheet_present=False,
            )

    def test_enabled_empty_sheet_schedules_only_final_checkpoint(self):
        schedule = self.build(
            {"USER_CHECKPOINTS": True, "TEND": 0.65},
            (),
        )

        self.assertTrue(schedule.user_enabled)
        self.assertEqual(
            schedule.boundaries,
            (CheckpointBoundary(time=0.65, trigger="final"),),
        )

    def test_requested_times_are_sorted_and_tend_is_appended(self):
        with self.assertWarnsRegex(UserWarning, "sorted"):
            schedule = self.build(
                {"USER_CHECKPOINTS": "TRUE", "TEND": 0.65},
                (0.4, 0.2),
            )

        self.assertEqual(
            schedule.boundaries,
            (
                CheckpointBoundary(time=0.2, trigger="requested"),
                CheckpointBoundary(time=0.4, trigger="requested"),
                CheckpointBoundary(time=0.65, trigger="final"),
            ),
        )

    def test_duplicate_times_within_tolerance_are_collapsed(self):
        with self.assertWarnsRegex(UserWarning, "duplicate"):
            schedule = self.build(
                {"USER_CHECKPOINTS": True, "TEND": 0.65},
                (0.2, 0.2000004, 0.4),
            )

        self.assertEqual(
            schedule.boundaries,
            (
                CheckpointBoundary(time=0.2, trigger="requested"),
                CheckpointBoundary(time=0.4, trigger="requested"),
                CheckpointBoundary(time=0.65, trigger="final"),
            ),
        )

    def test_explicit_time_within_tolerance_of_tend_has_requested_priority(self):
        schedule = self.build(
            {"USER_CHECKPOINTS": True, "TEND": 0.65},
            (0.6499995,),
        )

        self.assertEqual(
            schedule.boundaries,
            (CheckpointBoundary(time=0.65, trigger="requested"),),
        )

    def test_invalid_requested_times_are_rejected(self):
        invalid_cases = (
            (("not-a-number",), "finite numeric"),
            ((None,), "finite numeric"),
            ((np.nan,), "finite numeric"),
            ((np.inf,), "finite numeric"),
            ((0.0,), "greater than zero"),
            ((-0.1,), "greater than zero"),
            ((0.6501,), "cannot exceed TEND"),
        )
        for values, message in invalid_cases:
            with self.subTest(values=values):
                with self.assertRaisesRegex(ValueError, message):
                    self.build(
                        {"USER_CHECKPOINTS": True, "TEND": 0.65},
                        values,
                    )

    def test_restart_selects_first_boundary_strictly_after_current_time(self):
        schedule = self.build(
            {"USER_CHECKPOINTS": True, "TEND": 0.65},
            (0.2, 0.4),
        )

        boundary = next_checkpoint_boundary(
            schedule,
            current_time=0.2000005,
            epsilon=1.0e-6,
        )

        self.assertEqual(
            boundary,
            CheckpointBoundary(time=0.4, trigger="requested"),
        )

    def test_schedule_and_boundaries_are_immutable(self):
        schedule = self.build(
            {"USER_CHECKPOINTS": True, "TEND": 0.65},
            (0.2,),
        )

        with self.assertRaises(FrozenInstanceError):
            schedule.user_enabled = False
        with self.assertRaises(FrozenInstanceError):
            schedule.boundaries[0].time = 0.3


if __name__ == "__main__":
    unittest.main()
