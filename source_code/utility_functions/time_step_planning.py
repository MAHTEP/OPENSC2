"""Plan and commit the next conductor time without mutating completed history."""

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class NextTimeStepPlan:
    """Decision required to begin one thermal-hydraulic time step.

    ``time_step`` is the numerical increment used by the solver.  It can
    differ slightly from ``target_time - current_time`` when an event is
    synchronized within the accepted time tolerance; this preserves the
    historical OPENSC2 behaviour.
    """

    target_time: float
    time_step: float
    event_reached: bool = False
    time_step_forced: bool = False
    following_time_step: Optional[float] = None


def plan_next_time_step(conductor, epsilon: float, t_step_min: float) -> NextTimeStepPlan:
    """Return the next time-step decision without changing ``conductor``."""

    current_time = conductor.cond_time[-1]
    proposed_time_step = conductor.time_step
    proposed_target_time = current_time + proposed_time_step
    event_time = conductor.events_time[conductor.i_event]

    if proposed_target_time - epsilon <= event_time <= proposed_target_time + epsilon:
        return NextTimeStepPlan(
            target_time=event_time,
            time_step=proposed_time_step,
            event_reached=True,
            following_time_step=t_step_min,
        )

    if current_time < event_time < proposed_target_time:
        forced_time_step = event_time - current_time
        following_time_step = (
            forced_time_step if forced_time_step < t_step_min else t_step_min
        )
        return NextTimeStepPlan(
            target_time=event_time,
            time_step=forced_time_step,
            event_reached=True,
            time_step_forced=True,
            following_time_step=following_time_step,
        )

    return NextTimeStepPlan(
        target_time=proposed_target_time,
        time_step=proposed_time_step,
    )


def apply_time_step_plan(conductor, plan: NextTimeStepPlan):
    """Atomically install ``plan`` at the beginning of the next iteration."""

    if plan.time_step_forced:
        print(
            f"Forced {conductor.identifier} time step: {plan.time_step} s\n"
        )

    conductor.append_time(plan.target_time)
    conductor.time_step = plan.time_step

    if plan.following_time_step is not None:
        conductor.force_next_tstep_flag = True
        conductor.next_time_step = plan.following_time_step

    if plan.event_reached:
        conductor.move_to_next_event()

    return conductor
