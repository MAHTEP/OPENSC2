from types import SimpleNamespace
from unittest import TestCase
from unittest.mock import patch

import numpy as np
import pandas as pd

from conductor import Conductor
from solid_component import SolidComponent
from strand_component import StrandComponent
from utility_functions.step_matrix_construction import (
    eval_transport_coefficients,
)
from utility_functions.utils_global_info import IADAPTIME_VALUES


class InputControlSemanticsTests(TestCase):
    def test_eps_interpolation_is_independent_from_current_interpolation(self):
        strand = StrandComponent.__new__(StrandComponent)
        strand.identifier = "STACK_1"
        strand.operations = {
            "IEPS": -1,
            "IOP_INTERPOLATION": "linear",
            "EPS_INTERPOLATION": "cubic",
        }
        strand.dict_node_pt = {}
        conductor = SimpleNamespace(
            cond_time=[0.0],
            electric_time=0.0,
            BASE_PATH=".",
            file_input={"EXTERNAL_STRAIN": "strain.xlsx"},
            grid_features={"zcoord": np.array([0.0, 1.0]), "N_nod": 2},
        )

        with (
            patch(
                "strand_component.load_auxiliary_files",
                return_value=(object(), 0),
            ),
            patch(
                "strand_component.build_interpolator",
                return_value=(object(), "space_only"),
            ) as build,
            patch(
                "strand_component.do_interpolation",
                return_value=np.array([0.0, 0.0]),
            ),
        ):
            strand.get_eps(conductor)

        self.assertEqual(build.call_args.args[1], "cubic")

    def test_default_electric_subcycling_uses_ten_steps(self):
        conductor = Conductor.__new__(Conductor)
        conductor.time_step = 0.025
        conductor.inputs = {"ELECTRIC_TIME_STEP": None}
        conductor.operations = {"MAXIMUM_ITERATION_NUMBER": 1000}

        conductor._Conductor__get_electric_time_step()

        self.assertEqual(conductor.electric_time_step_number, 10)
        self.assertAlmostEqual(conductor.electric_time_step, 0.0025)

    def test_explicit_electric_step_exactly_partitions_th_interval(self):
        conductor = Conductor.__new__(Conductor)
        conductor.time_step = 1.0
        conductor.inputs = {"ELECTRIC_TIME_STEP": 0.3}
        conductor.operations = {"MAXIMUM_ITERATION_NUMBER": 4}

        conductor._Conductor__get_electric_time_step()

        self.assertEqual(conductor.electric_time_step_number, 4)
        self.assertAlmostEqual(conductor.electric_time_step, 0.25)

    def test_exact_decimal_ratio_does_not_create_an_extra_substep(self):
        conductor = Conductor.__new__(Conductor)
        conductor.time_step = 0.0025
        conductor.inputs = {"ELECTRIC_TIME_STEP": 0.00025}
        conductor.operations = {"MAXIMUM_ITERATION_NUMBER": 1000}

        conductor._Conductor__get_electric_time_step()

        self.assertEqual(conductor.electric_time_step_number, 10)

    def test_maximum_electric_substep_count_is_enforced(self):
        conductor = Conductor.__new__(Conductor)
        conductor.time_step = 1.0
        conductor.inputs = {"ELECTRIC_TIME_STEP": 0.2}
        conductor.operations = {"MAXIMUM_ITERATION_NUMBER": 4}

        with self.assertRaisesRegex(ValueError, "exceeds"):
            conductor._Conductor__get_electric_time_step()

    def test_transverse_transport_multiplier_scales_all_coefficients(self):
        comp_1 = SimpleNamespace(
            identifier="CHAN_1",
            coolant=SimpleNamespace(
                dict_Gauss_pt={
                    "velocity": np.array([2.0]),
                    "total_density": np.array([8.0]),
                    "total_enthalpy": np.array([10.0]),
                }
            ),
        )
        comp_2 = SimpleNamespace(identifier="CHAN_2")
        interface = SimpleNamespace(
            interf_name="CHAN_1_CHAN_2",
            comp_1=comp_1,
            comp_2=comp_2,
        )
        multiplier = pd.DataFrame(
            [[0.0, 0.5], [0.5, 0.0]],
            index=["CHAN_1", "CHAN_2"],
            columns=["CHAN_1", "CHAN_2"],
        )
        conductor = SimpleNamespace(
            dict_interf_peri={
                "ch_ch": {
                    "Open": {
                        "Gauss": {"CHAN_1_CHAN_2": np.array([3.0])}
                    }
                }
            },
            dict_df_coupling={"trans_transp_multiplier": multiplier},
            interface=SimpleNamespace(fluid_fluid=[interface]),
            dict_Gauss_pt={
                "K1": {"CHAN_1_CHAN_2": np.zeros(1)},
                "K2": {"CHAN_1_CHAN_2": np.zeros(1)},
                "K3": {"CHAN_1_CHAN_2": np.zeros(1)},
            },
            k_loc=1.0,
            lambda_v=1.0,
        )

        eval_transport_coefficients(
            conductor,
            "CHAN_1_CHAN_2",
            comp_1,
            np.array([0]),
            np.array([4.0]),
        )

        expected_k1 = 0.5 * 3.0 * np.sqrt(2.0 * 8.0 / 4.0)
        self.assertAlmostEqual(
            conductor.dict_Gauss_pt["K1"]["CHAN_1_CHAN_2"][0],
            expected_k1,
        )
        self.assertAlmostEqual(
            conductor.dict_Gauss_pt["K2"]["CHAN_1_CHAN_2"][0],
            expected_k1 * 2.0,
        )
        self.assertAlmostEqual(
            conductor.dict_Gauss_pt["K3"]["CHAN_1_CHAN_2"][0],
            expected_k1 * 12.0,
        )

    def test_legacy_bitr_botr_mode_is_explicitly_archived(self):
        solid = SolidComponent.__new__(SolidComponent)
        solid.operations = {"IBIFUN": 1}

        with self.assertRaisesRegex(NotImplementedError, "BITR/BOTR"):
            solid.get_magnetic_field(SimpleNamespace(), nodal=True)

    def test_user_function_adaptive_mode_is_not_exposed(self):
        self.assertNotIn(-2, IADAPTIME_VALUES)


if __name__ == "__main__":
    import unittest

    unittest.main()
