import hashlib
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from test_checkpoint import make_simulation
from utility_functions.checkpoint import build_continuation_profile


class CheckpointContinuationDriverProfileTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        self.input_dir = self.root / "input"
        self.input_dir.mkdir()
        (self.input_dir / "transitory_input.xlsx").write_bytes(
            b"input-content"
        )

    def _simulation(self):
        simulation = make_simulation(self.input_dir)
        simulation.transient_input = {
            "IADAPTIME": 0,
            "TIME_STEP": 0.025,
            "STPMIN": 0.0025,
            "STPMAX": 0.05,
            "MLT_INCREASE": 1.2,
            "MLT_DECREASE": 0.5,
            "TIMEREF": 0.2,
            "TAUREF": 0.05,
            "TEND": 0.4,
            "CHECKPOINT_EVERY_N_STEPS": 50,
            "USER_CHECKPOINTS": True,
        }
        conductor = simulation.list_of_Conductors[0]
        conductor.inputs.update(
            I0_OP_MODE=0,
            I0_OP_TOT=12000.0,
            ELECTRIC_TIME_STEP=0.001,
        )
        conductor.file_input = {
            "EXTERNAL_CURRENT": "none",
            "EXTERNAL_BFIELD": "none",
            "EXTERNAL_ALPHAB": "none",
            "EXTERNAL_HEAT": "none",
        }
        solid = conductor.inventory["SolidComponent"].collection[0]
        solid.operations = {
            "IOP_MODE": 0,
            "IOP_INTERPOLATION": "linear",
            "IBIFUN": 1,
            "BISS": 2.0,
            "BOSS": 3.0,
            "BITR": 0.5,
            "BOTR": 0.7,
            "B_INTERPOLATION": "linear",
            "IALPHAB": 0,
            "ALPHAB_INTERPOLATION": "linear",
            "IQFUN": 1,
            "Q_INTERPOLATION": "linear",
            "XQBEG": 1.0,
            "XQEND": 3.0,
            "Q0": 250.0,
            "TQBEG": 0.1,
            "TQEND": 0.2,
        }
        return simulation, conductor, solid

    def _write_auxiliary(self, relative_path, content):
        path = self.input_dir / relative_path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
        return {
            "kind": "auxiliary_file",
            "path": relative_path.as_posix(),
            "sha256": hashlib.sha256(content).hexdigest(),
            "size": len(content),
        }

    def test_canonical_driver_profile_is_identifier_keyed(self):
        simulation, _, _ = self._simulation()

        profile = build_continuation_profile(simulation)

        self.assertEqual(
            profile.time_policy["conductors"],
            {"COND_1": {"ELECTRIC_TIME_STEP": 0.001}},
        )
        self.assertEqual(
            profile.drivers,
            {
                "COND_1": {
                    "current": {
                        "source": {"kind": "canonical_input"},
                        "parameters": {
                            "I0_OP_MODE": 0,
                            "I0_OP_TOT": 12000.0,
                        },
                        "components": {
                            "STACK_1": {"IOP_MODE": 0},
                        },
                    },
                    "magnetic_field": {
                        "components": {
                            "STACK_1": {
                                "source": {"kind": "canonical_input"},
                                "parameters": {
                                    "IBIFUN": 1,
                                    "BISS": 2.0,
                                    "BOSS": 3.0,
                                    "BITR": 0.5,
                                    "BOTR": 0.7,
                                },
                            },
                        },
                    },
                    "magnetic_field_gradient": {
                        "components": {
                            "STACK_1": {
                                "source": {"kind": "disabled"},
                                "parameters": {"IALPHAB": 0},
                            },
                        },
                    },
                    "external_heat": {
                        "components": {
                            "STACK_1": {
                                "source": {"kind": "canonical_input"},
                                "parameters": {
                                    "IQFUN": 1,
                                    "XQBEG": 1.0,
                                    "XQEND": 3.0,
                                    "Q0": 250.0,
                                    "TQBEG": 0.1,
                                    "TQEND": 0.2,
                                },
                            },
                        },
                    },
                },
            },
        )

    def test_auxiliary_driver_profile_records_content_identity(self):
        simulation, conductor, solid = self._simulation()
        expected_current = self._write_auxiliary(
            Path("drivers/current.tsv"),
            b"current-data",
        )
        expected_field = self._write_auxiliary(
            Path("drivers/field.tsv"),
            b"field-data",
        )
        expected_gradient = self._write_auxiliary(
            Path("drivers/gradient.tsv"),
            b"gradient-data",
        )
        expected_heat = self._write_auxiliary(
            Path("drivers/heat.tsv"),
            b"heat-data",
        )
        conductor.inputs["I0_OP_MODE"] = -1
        conductor.file_input.update(
            EXTERNAL_CURRENT=expected_current["path"],
            EXTERNAL_BFIELD=expected_field["path"],
            EXTERNAL_ALPHAB=expected_gradient["path"],
            EXTERNAL_HEAT=expected_heat["path"],
        )
        solid.operations.update(
            IOP_MODE=-1,
            IOP_INTERPOLATION="cubic",
            IBIFUN=-1,
            B_INTERPOLATION="cubic",
            B_field_units="T/A",
            IALPHAB=-1,
            ALPHAB_INTERPOLATION="cubic",
            IQFUN=-1,
            Q_INTERPOLATION="cubic",
        )

        profile = build_continuation_profile(simulation)
        drivers = profile.drivers["COND_1"]

        self.assertEqual(drivers["current"]["source"], expected_current)
        self.assertEqual(
            drivers["current"]["components"]["STACK_1"],
            {"IOP_MODE": -1, "IOP_INTERPOLATION": "cubic"},
        )
        self.assertEqual(
            drivers["magnetic_field"]["components"]["STACK_1"],
            {
                "source": expected_field,
                "parameters": {
                    "IBIFUN": -1,
                    "B_INTERPOLATION": "cubic",
                    "B_field_units": "T/A",
                },
            },
        )
        self.assertEqual(
            drivers["magnetic_field_gradient"]["components"][
                "STACK_1"
            ],
            {
                "source": expected_gradient,
                "parameters": {
                    "IALPHAB": -1,
                    "ALPHAB_INTERPOLATION": "cubic",
                },
            },
        )
        self.assertEqual(
            drivers["external_heat"]["components"]["STACK_1"],
            {
                "source": expected_heat,
                "parameters": {
                    "IQFUN": -1,
                    "Q_INTERPOLATION": "cubic",
                    "TQBEG": 0.1,
                    "TQEND": 0.2,
                },
            },
        )

    def test_function_driver_profile_reserves_callable_identity(self):
        simulation, conductor, solid = self._simulation()
        conductor.inputs["I0_OP_MODE"] = -2
        solid.operations["IQFUN"] = -2

        profile = build_continuation_profile(simulation)
        drivers = profile.drivers["COND_1"]

        self.assertEqual(
            drivers["current"]["source"],
            {
                "kind": "python_function",
                "callable": (
                    "utility_functions.electric_auxiliary_functions."
                    "custom_current_function"
                ),
            },
        )
        self.assertEqual(
            drivers["external_heat"]["components"]["STACK_1"][
                "source"
            ],
            {
                "kind": "python_function",
                "callable": (
                    "solid_component.SolidComponent.user_heat_function"
                ),
            },
        )


if __name__ == "__main__":
    unittest.main()
