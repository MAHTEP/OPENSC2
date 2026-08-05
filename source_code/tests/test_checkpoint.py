import hashlib
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import h5py
import numpy as np

from utility_functions.checkpoint import (
    CheckpointValidationError,
    checkpoint_interval,
    SCHEMA_VERSION,
    build_input_manifest,
    validate_checkpoint_state,
    write_checkpoint,
    write_periodic_checkpoint_if_due,
)


def _collection(*items):
    return SimpleNamespace(collection=list(items))


def make_simulation(base_path, *, force_next_tstep_flag=False):
    solid = SimpleNamespace(
        identifier="STACK_1",
        dict_node_pt={
            "EEXT": np.arange(6.0).reshape(3, 2),
            "EJHT": np.arange(6.0, 12.0).reshape(3, 2),
        },
        store_sd_node={"temperature": {"t_save_left": np.array([4.2, 4.3])}},
        time_evol={"temperature": [4.2, 4.3]},
        time_evol_gauss={},
        time_evol_max_temperature={"time (s)": [0.0, 0.1]},
    )
    coolant = SimpleNamespace(
        time_evol={"pressure": [1.0e5, 1.01e5]},
        time_evol_io={"time (s)": [0.0, 0.1]},
        time_evol_max_temperature={},
    )
    fluid = SimpleNamespace(identifier="CHAN_1", coolant=coolant)

    conductor = SimpleNamespace(
        identifier="COND_1",
        cond_time=[0.0, 0.1],
        cond_num_step=1,
        time_step=0.1,
        EQTEIG=np.array([1.0, 2.0]),
        events_time=np.array([0.01, 1.0]),
        i_event=1,
        force_next_tstep_flag=force_next_tstep_flag,
        dict_Step={
            "SYSVAR": np.arange(6.0).reshape(3, 2),
            "SYSLOD": np.arange(12.0).reshape(3, 4),
            "AM4_AA": np.arange(24.0).reshape(4, 2, 3),
        },
        inputs={
            "METHOD": "AM4",
            "ELECTRIC_METHOD": "CN",
            "I0_OP_MODE": 0,
        },
        electric_solution=np.array([10.0, 20.0]),
        electric_solution_steady=np.array([9.0, 19.0]),
        electric_time=0.1,
        cond_el_num_step=10,
        enthalpy_balance=1.0,
        enthalpy_inl=2.0,
        enthalpy_out=3.0,
        E_sol_ini=4.0,
        E_str_ini=5.0,
        E_jk_ini=6.0,
        i_save=1,
        num_step_save=np.array([0, 1]),
        t_save_left=0.1,
        store_sd_node={"zcoord": {"t_save_left": np.array([0.0, 1.0])}},
        store_sd_gauss={},
        inventory={
            "SolidComponent": _collection(solid),
            "all_component": _collection(fluid, solid),
        },
    )
    return SimpleNamespace(
        basePath=str(base_path),
        simulation_time=[0.0, 0.1],
        num_step=1,
        list_of_Conductors=[conductor],
    )


class CheckpointTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        self.input_dir = self.root / "input"
        self.input_dir.mkdir()
        (self.input_dir / "transitory_input.xlsx").write_bytes(b"input-content")

    def test_writer_persists_schema_and_restart_state(self):
        simulation = make_simulation(self.input_dir)

        checkpoint = write_checkpoint(
            simulation, self.root / "checkpoints", trigger="periodic"
        )

        self.assertEqual(checkpoint.name, "checkpoint_step_000001.h5")
        self.assertTrue(checkpoint.is_file())
        self.assertFalse(Path(f"{checkpoint}.tmp").exists())

        with h5py.File(checkpoint, "r") as h5file:
            metadata = h5file["metadata"]
            self.assertEqual(metadata.attrs["schema_version"], SCHEMA_VERSION)
            self.assertEqual(metadata.attrs["trigger"], "periodic")
            self.assertTrue(bool(metadata.attrs["complete"]))

            np.testing.assert_allclose(
                h5file["simulation/simulation_time"][:], [0.0, 0.1]
            )
            self.assertEqual(h5file["simulation/num_step"][()], 1)

            conductor = h5file["conductors/COND_1"]
            np.testing.assert_allclose(
                conductor["th_history/AM4_AA"][:],
                simulation.list_of_Conductors[0].dict_Step["AM4_AA"],
            )
            np.testing.assert_allclose(
                conductor["electric/electric_solution"][:], [10.0, 20.0]
            )
            np.testing.assert_allclose(
                conductor[
                    "components/STACK_1/energy_history/EEXT"
                ][:],
                np.arange(6.0).reshape(3, 2),
            )
            np.testing.assert_allclose(
                conductor[
                    "output_state/buffers/components/CHAN_1/coolant/"
                    "time_evol/pressure"
                ][:],
                [1.0e5, 1.01e5],
            )

    def test_manifest_is_deterministic_and_uses_relative_paths(self):
        nested = self.input_dir / "auxiliary"
        nested.mkdir()
        data_file = nested / "current.tsv"
        data_file.write_bytes(b"0\t1\n")
        simulation = make_simulation(self.input_dir)

        manifest = build_input_manifest(simulation)

        self.assertEqual(
            [entry["path"] for entry in manifest],
            ["auxiliary/current.tsv", "transitory_input.xlsx"],
        )
        self.assertEqual(
            manifest[0]["sha256"], hashlib.sha256(b"0\t1\n").hexdigest()
        )

    def test_invalid_boundary_does_not_create_final_checkpoint(self):
        simulation = make_simulation(
            self.input_dir, force_next_tstep_flag=True
        )
        checkpoint_dir = self.root / "checkpoints"

        with self.assertRaisesRegex(
            CheckpointValidationError, "force_next_tstep_flag must be False"
        ):
            write_checkpoint(simulation, checkpoint_dir, trigger="periodic")

        self.assertFalse((checkpoint_dir / "checkpoint_step_000001.h5").exists())

    def test_validation_rejects_inconsistent_conductor_history(self):
        simulation = make_simulation(self.input_dir)
        simulation.list_of_Conductors[0].cond_num_step = 2

        with self.assertRaisesRegex(
            CheckpointValidationError, "must equal len\\(cond_time\\) - 1"
        ):
            validate_checkpoint_state(simulation)

    def test_writer_accepts_conductor_without_active_electric_model(self):
        simulation = make_simulation(self.input_dir)
        conductor = simulation.list_of_Conductors[0]
        conductor.inputs["I0_OP_MODE"] = None
        del conductor.electric_solution
        del conductor.electric_solution_steady

        checkpoint = write_checkpoint(
            simulation, self.root / "checkpoints", trigger="final"
        )

        with h5py.File(checkpoint, "r") as h5file:
            self.assertEqual(h5file["metadata"].attrs["trigger"], "final")
            self.assertEqual(len(h5file["conductors/COND_1/electric"]), 0)

    def test_validation_rejects_non_integer_global_step(self):
        simulation = make_simulation(self.input_dir)
        simulation.num_step = 1.5

        with self.assertRaisesRegex(
            CheckpointValidationError, "num_step must be an integer"
        ):
            validate_checkpoint_state(simulation)

    def test_invalid_trigger_is_rejected(self):
        simulation = make_simulation(self.input_dir)

        with self.assertRaisesRegex(ValueError, "Invalid checkpoint trigger"):
            write_checkpoint(
                simulation, self.root / "checkpoints", trigger="manual"
            )

    def test_periodic_checkpoint_uses_default_interval(self):
        simulation = make_simulation(self.input_dir)
        simulation.num_step = 100
        simulation.transient_input = {}
        simulation.dict_path = {"Checkpoint_dir": self.root / "checkpoints"}

        with patch(
            "utility_functions.checkpoint.write_checkpoint",
            return_value=self.root / "checkpoints" / "checkpoint_step_000100.h5",
        ) as writer:
            checkpoint = write_periodic_checkpoint_if_due(simulation)

        self.assertIsNotNone(checkpoint)
        writer.assert_called_once_with(
            simulation,
            self.root / "checkpoints",
            trigger="periodic",
        )

    def test_periodic_checkpoint_skips_non_boundary_and_disabled_mode(self):
        simulation = make_simulation(self.input_dir)
        simulation.num_step = 99
        simulation.transient_input = {"CHECKPOINT_EVERY_N_STEPS": 100}
        simulation.dict_path = {"Checkpoint_dir": self.root / "checkpoints"}

        with patch("utility_functions.checkpoint.write_checkpoint") as writer:
            self.assertIsNone(write_periodic_checkpoint_if_due(simulation))
            simulation.num_step = 100
            simulation.transient_input["CHECKPOINT_EVERY_N_STEPS"] = 0
            self.assertIsNone(write_periodic_checkpoint_if_due(simulation))

        writer.assert_not_called()

    def test_checkpoint_interval_accepts_integral_excel_number(self):
        self.assertEqual(
            checkpoint_interval({"CHECKPOINT_EVERY_N_STEPS": 25.0}),
            25,
        )

    def test_checkpoint_interval_rejects_invalid_value(self):
        for value in (-1, 2.5, "not-a-number", True):
            with self.subTest(value=value):
                with self.assertRaisesRegex(
                    ValueError, "must be a non-negative integer"
                ):
                    checkpoint_interval({"CHECKPOINT_EVERY_N_STEPS": value})


if __name__ == "__main__":
    unittest.main()
