import unittest
from copy import deepcopy
from contextlib import redirect_stdout
from io import StringIO
from types import SimpleNamespace

from utility_functions.time_step_planning import (
    apply_time_step_plan,
    plan_next_time_step,
)


class DummyConductor(SimpleNamespace):
    def append_time(self, time):
        self.cond_time.append(time)
        self.cond_num_step += 1

    def move_to_next_event(self):
        if self.i_event < self.i_event_max:
            self.i_event += 1


def make_conductor(*, current_time=1.0, time_step=0.2, events=(10.0,)):
    return DummyConductor(
        cond_time=[current_time],
        cond_num_step=4,
        time_step=time_step,
        events_time=list(events),
        i_event=0,
        i_event_max=len(events) - 1,
        force_next_tstep_flag=False,
        next_time_step=0.0,
        identifier="dummy",
    )


class PlanNextTimeStepTests(unittest.TestCase):
    def apply_plan(self, conductor, plan):
        with redirect_stdout(StringIO()):
            apply_time_step_plan(conductor, plan)

    def test_normal_step_is_planned_without_mutating_conductor(self):
        conductor = make_conductor()
        state_before = deepcopy(vars(conductor))

        plan = plan_next_time_step(conductor, epsilon=1e-6, t_step_min=0.05)

        self.assertEqual(vars(conductor), state_before)
        self.assertAlmostEqual(plan.target_time, 1.2)
        self.assertEqual(plan.time_step, 0.2)
        self.assertFalse(plan.event_reached)
        self.assertFalse(plan.time_step_forced)
        self.assertIsNone(plan.following_time_step)

    def test_event_within_epsilon_preserves_nominal_time_step(self):
        conductor = make_conductor(events=(1.2000005, 2.0))

        plan = plan_next_time_step(conductor, epsilon=1e-6, t_step_min=0.05)
        self.apply_plan(conductor, plan)

        self.assertEqual(conductor.cond_time, [1.0, 1.2000005])
        self.assertEqual(conductor.cond_num_step, 5)
        self.assertEqual(conductor.time_step, 0.2)
        self.assertEqual(conductor.i_event, 1)
        self.assertTrue(conductor.force_next_tstep_flag)
        self.assertEqual(conductor.next_time_step, 0.05)

    def test_event_within_epsilon_is_symmetric(self):
        conductor = make_conductor(events=(1.1999995, 2.0))

        plan = plan_next_time_step(conductor, epsilon=1e-6, t_step_min=0.05)
        self.apply_plan(conductor, plan)

        self.assertEqual(conductor.cond_time, [1.0, 1.1999995])
        self.assertEqual(conductor.time_step, 0.2)
        self.assertEqual(conductor.i_event, 1)
        self.assertEqual(conductor.next_time_step, 0.05)

    def test_crossed_event_forces_current_step_and_minimum_following_step(self):
        conductor = make_conductor(time_step=0.3, events=(1.2, 2.0))

        plan = plan_next_time_step(conductor, epsilon=1e-6, t_step_min=0.05)
        self.apply_plan(conductor, plan)

        self.assertEqual(conductor.cond_time, [1.0, 1.2])
        self.assertAlmostEqual(conductor.time_step, 0.2)
        self.assertTrue(plan.time_step_forced)
        self.assertEqual(conductor.next_time_step, 0.05)
        self.assertEqual(conductor.i_event, 1)

    def test_forced_step_below_minimum_is_reused_after_event(self):
        conductor = make_conductor(time_step=0.2, events=(1.02, 2.0))

        plan = plan_next_time_step(conductor, epsilon=1e-6, t_step_min=0.05)
        self.apply_plan(conductor, plan)

        self.assertAlmostEqual(conductor.time_step, 0.02)
        self.assertAlmostEqual(conductor.next_time_step, 0.02)

    def test_last_event_index_does_not_advance_past_timeline(self):
        conductor = make_conductor(events=(1.2,))

        plan = plan_next_time_step(conductor, epsilon=1e-6, t_step_min=0.05)
        self.apply_plan(conductor, plan)

        self.assertEqual(conductor.i_event, 0)

    def test_consecutive_events_are_planned_one_iteration_at_a_time(self):
        conductor = make_conductor(events=(1.2, 1.25))

        first_plan = plan_next_time_step(
            conductor, epsilon=1e-6, t_step_min=0.05
        )
        self.apply_plan(conductor, first_plan)

        # This reproduces get_time_step consuming the one-shot override after
        # the first event step has been solved.
        conductor.time_step = conductor.next_time_step
        conductor.force_next_tstep_flag = False

        second_plan = plan_next_time_step(
            conductor, epsilon=1e-6, t_step_min=0.05
        )
        self.apply_plan(conductor, second_plan)

        self.assertEqual(conductor.cond_time, [1.0, 1.2, 1.25])
        self.assertEqual(conductor.cond_num_step, 6)
        self.assertEqual(conductor.i_event, 1)


if __name__ == "__main__":
    unittest.main()
