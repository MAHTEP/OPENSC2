import importlib.util
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace
import unittest
from unittest.mock import Mock, call, patch


def _module(name, **attributes):
    module = ModuleType(name)
    module.__dict__.update(attributes)
    return module


def _load_simulation_module():
    """Load simulation.py with lightweight stand-ins for unrelated imports."""

    noop = lambda *args, **kwargs: None
    utility_package = _module("utility_functions")
    utility_package.__path__ = []
    stubs = {
        "line_profiler": _module("line_profiler", LineProfiler=object),
        "conductor": _module("conductor", Conductor=object),
        "conductor_flags": _module(
            "conductor_flags", IOP_NOT_DEFINED=None, SHEET_NAME=""
        ),
        "environment": _module("environment", Environment=object),
        "utility_functions": utility_package,
        "utility_functions.auxiliary_functions": _module(
            "utility_functions.auxiliary_functions",
            check_repeated_headings=noop,
            check_object_number=noop,
            check_flag_value=noop,
            check_sheet_names=noop,
            with_read_csv=noop,
            with_read_excel=noop,
        ),
        "utility_functions.transient_solution_functions": _module(
            "utility_functions.transient_solution_functions",
            get_time_step=noop,
            step=noop,
        ),
        "utility_functions.time_step_planning": _module(
            "utility_functions.time_step_planning",
            apply_time_step_plan=noop,
            plan_next_time_step=noop,
        ),
        "utility_functions.checkpoint": _module(
            "utility_functions.checkpoint",
            write_periodic_checkpoint_if_due=noop,
            write_checkpoint_if_due=noop,
        ),
        "utility_functions.checkpoint_schedule": _module(
            "utility_functions.checkpoint_schedule",
            load_checkpoint_schedule=noop,
            next_checkpoint_boundary=noop,
        ),
        "utility_functions.output": _module(
            "utility_functions.output",
            save_simulation_space=Mock(),
            reorganize_spatial_distribution=noop,
            save_simulation_time=noop,
            save_properties=noop,
        ),
        "utility_functions.plots": _module(
            "utility_functions.plots",
            plot_properties=noop,
            make_plots=noop,
            create_real_time_plots=noop,
            update_real_time_plots=noop,
        ),
        "simulation_global_info": _module(
            "simulation_global_info", MLT_DEFAULT_VALUE=0
        ),
        "utility_functions.utils_global_info": _module(
            "utility_functions.utils_global_info", VALID_FLAG_VALUES=()
        ),
    }

    module_path = Path(__file__).parents[1] / "simulation.py"
    spec = importlib.util.spec_from_file_location(
        "simulation_restart_preamble_test_target", module_path
    )
    module = importlib.util.module_from_spec(spec)

    with patch.dict(sys.modules, stubs):
        spec.loader.exec_module(module)

    return module


class RestartPreambleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.simulation_module = _load_simulation_module()

    def make_runtime(self, *, restored):
        conductor = SimpleNamespace(
            identifier="COND_1",
            i_save=7,
            compute_radiative_heat_exhange_jk=Mock(),
            compute_heat_exchange_jk_env=Mock(),
            store_spatial_distributions_t0=Mock(),
            store_spatial_distributions=Mock(),
            post_processing=Mock(),
        )
        simulation = self.simulation_module.Simulation.__new__(
            self.simulation_module.Simulation
        )
        simulation.list_of_Conductors = [conductor]
        simulation.environment = object()
        simulation.dict_path = {
            "Output_Spatial_distribution_COND_1_dir": "spatial",
            "Output_Solution_COND_1_dir": "solution",
        }
        simulation.n_digit_time = 6
        simulation.simulation_time = [0.5]
        simulation.transient_input = {
            "TEND": 0.5,
            "STPMIN": 0.01,
            "SIMULATION": "restart-preamble-test",
        }
        simulation.restored_from_checkpoint = restored
        return simulation, conductor

    def test_fresh_run_preserves_t0_output_preamble(self):
        simulation, conductor = self.make_runtime(restored=False)
        save_space = self.simulation_module.save_simulation_space
        save_space.reset_mock()

        simulation.conductor_solution(None)

        conductor.store_spatial_distributions_t0.assert_called_once_with(
            "t_save"
        )
        save_space.assert_called_once()
        self.assertEqual(conductor.i_save, 8)
        conductor.compute_radiative_heat_exhange_jk.assert_called_once_with()
        conductor.compute_heat_exchange_jk_env.assert_called_once_with(
            simulation.environment
        )

    def test_post_processing_refreshes_tend_buffer_before_writing(self):
        simulation, conductor = self.make_runtime(restored=True)
        save_space = self.simulation_module.save_simulation_space
        save_space.reset_mock()

        ordered_calls = Mock()
        ordered_calls.attach_mock(
            conductor.store_spatial_distributions,
            "refresh",
        )
        ordered_calls.attach_mock(save_space, "write")

        simulation.conductor_post_processing()

        conductor.post_processing.assert_called_once_with(simulation)
        conductor.store_spatial_distributions.assert_called_once_with(
            t_save_key="t_save"
        )
        ordered_calls.assert_has_calls(
            [
                call.refresh(t_save_key="t_save"),
                call.write(conductor, "spatial"),
            ]
        )
        conductor.store_spatial_distributions_t0.assert_not_called()

    def test_restart_skips_only_t0_output_preamble(self):
        simulation, conductor = self.make_runtime(restored=True)
        save_space = self.simulation_module.save_simulation_space
        save_space.reset_mock()

        simulation.conductor_solution(None)

        conductor.store_spatial_distributions_t0.assert_not_called()
        save_space.assert_not_called()
        self.assertEqual(conductor.i_save, 7)
        conductor.compute_radiative_heat_exhange_jk.assert_called_once_with()
        conductor.compute_heat_exchange_jk_env.assert_called_once_with(
            simulation.environment
        )


if __name__ == "__main__":
    unittest.main()
