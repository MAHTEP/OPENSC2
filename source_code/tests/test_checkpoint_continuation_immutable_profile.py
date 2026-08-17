import copy
import hashlib
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

from test_checkpoint import make_simulation
from utility_functions.checkpoint import (
    build_continuation_profile,
    compare_continuation_profiles,
)


class CheckpointContinuationImmutableProfileTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        self.input_dir = self.root / "input"
        self.input_dir.mkdir()
        (self.input_dir / "transitory_input.xlsx").write_bytes(
            b"input-content"
        )

    def _write_input(self, name, content):
        path = self.input_dir / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
        return {
            "path": name,
            "sha256": hashlib.sha256(content).hexdigest(),
            "size": len(content),
        }

    def _simulation(self):
        grid = self._write_input("grid.dat", b"grid-content")
        coupling = self._write_input(
            "coupling.dat",
            b"coupling-content",
        )
        structure = self._write_input(
            "structure.dat",
            b"structure-content",
        )
        (self.input_dir / "operation.xlsx").write_bytes(
            b"mixed-driver-and-static-content"
        )

        simulation = make_simulation(self.input_dir)
        simulation.transient_input = {
            "IADAPTIME": 0,
            "TIME_STEP": 0.025,
            "STPMIN": 0.0025,
            "TEND": 0.4,
            "MAGNET": "magnet.xlsx",
            "ENVIRONMENT": "environment.xlsx",
            "SIMULATION": "continuation_case",
        }
        simulation.environment = SimpleNamespace(
            inputs={
                "Medium": "helium",
                "Temperature": 4.5,
            },
        )

        conductor = simulation.list_of_Conductors[0]
        conductor.inputs.update(
            ZLENGTH=10.0,
            NELEMS=100,
            I0_OP_MODE=0,
            I0_OP_TOT=12000.0,
            ELECTRIC_TIME_STEP=0.001,
        )
        conductor.operations = {
            "ELECTRIC_SOLVER": 0,
            "SELF_INDUCTANCE": 1.0e-6,
        }
        conductor.file_input = {
            "GRID_DEFINITION": grid["path"],
            "STRUCTURE_COUPLING": coupling["path"],
            "STRUCTURE_ELEMENTS": structure["path"],
            "OPERATION": "operation.xlsx",
            "EXTERNAL_CURRENT": "drivers/current.tsv",
            "EXTERNAL_BFIELD": "drivers/field.tsv",
            "EXTERNAL_ALPHAB": "drivers/gradient.tsv",
            "EXTERNAL_HEAT": "drivers/heat.tsv",
        }

        fluid = conductor.inventory["FluidComponent"].collection[0]
        fluid.inputs = {
            "HYDIAMETER": 0.01,
            "CROSSECTION": 1.0e-4,
        }
        fluid.operations = {
            "INTIAL": 1,
            "TEMINL": 4.5,
        }

        solid = conductor.inventory["SolidComponent"].collection[0]
        solid.inputs = {
            "CROSSECTION": 2.0e-4,
            "stabilizer_material": "cu",
            "RRR": 100.0,
        }
        solid.operations = {
            "TCS_EVALUATION": True,
            "FIX_POTENTIAL_FLAG": 0,
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

        return simulation, conductor, fluid, solid, {
            "GRID_DEFINITION": grid,
            "STRUCTURE_COUPLING": coupling,
            "STRUCTURE_ELEMENTS": structure,
        }

    def test_profile_captures_non_driver_runtime_input_semantics(self):
        simulation, _, _, _, _ = self._simulation()

        profile = build_continuation_profile(simulation)
        immutable = profile.immutable

        self.assertNotIn("IADAPTIME", immutable)
        self.assertEqual(profile.time_policy["IADAPTIME"], 0)
        self.assertEqual(
            immutable["transient_input"],
            {
                "ENVIRONMENT": "environment.xlsx",
                "MAGNET": "magnet.xlsx",
                "SIMULATION": "continuation_case",
            },
        )
        self.assertEqual(
            immutable["environment"],
            {
                "inputs": {
                    "Medium": "helium",
                    "Temperature": 4.5,
                },
            },
        )

        conductor = immutable["conductors"]["COND_1"]
        self.assertEqual(
            conductor["inputs"],
            {
                "ELECTRIC_METHOD": "CN",
                "METHOD": "AM4",
                "NELEMS": 100,
                "ZLENGTH": 10.0,
            },
        )
        self.assertEqual(
            conductor["operations"],
            {
                "ELECTRIC_SOLVER": 0,
                "SELF_INDUCTANCE": 1.0e-6,
            },
        )
        self.assertEqual(
            conductor["components"]["CHAN_1"],
            {
                "kind": "fluid",
                "inputs": {
                    "CROSSECTION": 1.0e-4,
                    "HYDIAMETER": 0.01,
                },
                "operations": {
                    "INTIAL": 1,
                    "TEMINL": 4.5,
                },
            },
        )
        self.assertEqual(
            conductor["components"]["STACK_1"],
            {
                "kind": "solid",
                "inputs": {
                    "CROSSECTION": 2.0e-4,
                    "RRR": 100.0,
                    "stabilizer_material": "cu",
                },
                "operations": {
                    "FIX_POTENTIAL_FLAG": 0,
                    "TCS_EVALUATION": True,
                },
            },
        )

    def test_profile_records_only_static_file_content_identity(self):
        simulation, _, _, _, expected_files = self._simulation()

        profile = build_continuation_profile(simulation)
        conductor = profile.immutable["conductors"]["COND_1"]

        self.assertEqual(conductor["static_files"], expected_files)
        static_text = repr(conductor["static_files"])
        self.assertNotIn("EXTERNAL_CURRENT", static_text)
        self.assertNotIn("EXTERNAL_BFIELD", static_text)
        self.assertNotIn("EXTERNAL_ALPHAB", static_text)
        self.assertNotIn("EXTERNAL_HEAT", static_text)
        self.assertNotIn("OPERATION", static_text)

    def test_fixed_parameter_changes_are_blocking(self):
        first, _, _, _, _ = self._simulation()
        second = copy.deepcopy(first)
        second_conductor = second.list_of_Conductors[0]
        second_conductor.inputs["ZLENGTH"] = 12.0
        second_solid = second_conductor.inventory[
            "SolidComponent"
        ].collection[0]
        second_solid.inputs["CROSSECTION"] = 3.0e-4

        comparison = compare_continuation_profiles(
            build_continuation_profile(first),
            build_continuation_profile(second),
        )

        self.assertFalse(comparison.is_compatible)
        self.assertEqual(
            comparison.immutable_differences,
            (
                (
                    "immutable.conductors.COND_1.components."
                    "STACK_1.inputs.CROSSECTION"
                ),
                "immutable.conductors.COND_1.inputs.ZLENGTH",
            ),
        )

    def test_time_and_driver_changes_do_not_enter_immutable_profile(self):
        first, _, _, _, _ = self._simulation()
        second = copy.deepcopy(first)
        second.transient_input["IADAPTIME"] = 1
        second.transient_input["TIME_STEP"] = 0.005
        second.transient_input["TEND"] = 1.0
        second_conductor = second.list_of_Conductors[0]
        second_conductor.inputs["I0_OP_TOT"] = 15000.0
        second_solid = second_conductor.inventory[
            "SolidComponent"
        ].collection[0]
        second_solid.operations["Q0"] = 500.0

        first_profile = build_continuation_profile(first)
        second_profile = build_continuation_profile(second)
        comparison = compare_continuation_profiles(
            first_profile,
            second_profile,
        )

        self.assertEqual(first_profile.immutable, second_profile.immutable)
        self.assertTrue(comparison.is_compatible)
        self.assertIn(
            "time_policy.IADAPTIME",
            comparison.time_policy_differences,
        )
        self.assertTrue(comparison.driver_differences)

    def test_profile_build_is_deterministic_and_non_mutating(self):
        simulation, _, _, _, _ = self._simulation()
        transient_before = copy.deepcopy(simulation.transient_input)
        environment_before = copy.deepcopy(simulation.environment.inputs)
        conductor = simulation.list_of_Conductors[0]
        conductor_inputs_before = copy.deepcopy(conductor.inputs)
        conductor_operations_before = copy.deepcopy(conductor.operations)

        first = build_continuation_profile(simulation)
        second = build_continuation_profile(simulation)

        self.assertEqual(first, second)
        self.assertEqual(simulation.transient_input, transient_before)
        self.assertEqual(simulation.environment.inputs, environment_before)
        self.assertEqual(conductor.inputs, conductor_inputs_before)
        self.assertEqual(conductor.operations, conductor_operations_before)


if __name__ == "__main__":
    unittest.main()
