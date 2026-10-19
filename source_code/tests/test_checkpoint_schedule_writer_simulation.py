from types import SimpleNamespace
import unittest
from unittest.mock import patch

from simulation import Simulation
from utility_functions.checkpoint_schedule import (
    CheckpointBoundary,
    CheckpointSchedule,
)


class CheckpointWriterSimulationTests(unittest.TestCase):
    def test_solution_delegates_checkpoint_selection_to_single_arbiter(self):
        simulation = Simulation.__new__(Simulation)
        simulation.checkpoint_schedule = CheckpointSchedule(
            user_enabled=True,
            boundaries=(
                CheckpointBoundary(time=0.1, trigger="requested"),
            ),
        )
        simulation.epsilon = 1.0e-6
        simulation.simulation_time = [0.0]
        simulation.simulation_time_step = 0.0
        simulation.transient_input = {
            "SIMULATION": "schedule-writer-test",
            "TEND": 0.1,
            "STPMIN": 0.1,
        }
        simulation.num_step = 0
        simulation.numObj = 1
        simulation.n_digit_time = 6
        simulation.restored_from_checkpoint = True
        simulation.environment = SimpleNamespace()
        simulation.starter_file_path = "transitory_input.xlsx"
        simulation.dict_qsource = {"COND_1": None}
        simulation.dict_path = {}

        conductor = SimpleNamespace(
            identifier="COND_1",
            cond_time=[0.0],
            time_step=0.1,
            inputs={"I0_OP_MODE": 0},
            inventory={
                "FluidComponent": SimpleNamespace(collection=[]),
            },
            Space_save=[10.0],
            i_save=0,
            i_save_max=0,
            compute_radiative_heat_exhange_jk=lambda: None,
            compute_heat_exchange_jk_env=lambda environment: None,
            electric_method=lambda: None,
            operating_conditions_em=lambda: None,
            operating_conditions_th=lambda runtime: None,
            build_heat_source=lambda runtime: None,
        )
        simulation.list_of_Conductors = [conductor]

        def apply_plan(runtime_conductor, plan):
            runtime_conductor.time_step = 0.1
            runtime_conductor.cond_time.append(0.1)

        with (
            patch(
                "simulation.IOP_NOT_DEFINED",
                0,
            ),
            patch(
                "simulation.plan_next_time_step",
                return_value=SimpleNamespace(),
            ) as planner,
            patch(
                "simulation.apply_time_step_plan",
                side_effect=apply_plan,
            ),
            patch("simulation.step"),
            patch("simulation.update_real_time_plots"),
            patch("simulation.save_simulation_time"),
            patch(
                "simulation.get_time_step",
                return_value=0.1,
            ),
            patch(
                "simulation.write_checkpoint_if_due",
                create=True,
                return_value="checkpoint-result",
            ) as arbiter,
            patch(
                "simulation.write_periodic_checkpoint_if_due",
                create=True,
                return_value=None,
            ) as legacy_writer,
            patch("builtins.print") as output,
        ):
            simulation.conductor_solution(gui=None)

        planner.assert_called_once_with(
            conductor,
            simulation.epsilon,
            simulation.transient_input["STPMIN"],
            scheduled_boundary_time=0.1,
        )
        arbiter.assert_called_once_with(simulation)
        legacy_writer.assert_not_called()
        output.assert_any_call(
            "Checkpoint saved: checkpoint-result"
        )
        self.assertEqual(simulation.num_step, 1)
        self.assertEqual(simulation.simulation_time, [0.0, 0.1])


if __name__ == "__main__":
    unittest.main()
