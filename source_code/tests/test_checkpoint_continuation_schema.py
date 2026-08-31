from dataclasses import FrozenInstanceError
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import h5py

from test_checkpoint import make_simulation
from utility_functions.checkpoint import (
    CheckpointReadError,
    SCHEMA_VERSION,
    SUPPORTED_SCHEMA_VERSIONS,
    read_checkpoint,
    write_checkpoint,
)


class CheckpointContinuationSchemaTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        self.input_dir = self.root / "input"
        self.input_dir.mkdir()
        (self.input_dir / "transitory_input.xlsx").write_bytes(
            b"input-content"
        )

    def _simulation_with_continuation_policy(self):
        simulation = make_simulation(self.input_dir)
        simulation.transient_input = {
            "IADAPTIME": 0,
            "TIME_STEP": 0.025,
            "STPMIN": 0.0025,
            "TEND": 0.4,
            "CHECKPOINT_EVERY_N_STEPS": 50,
            "USER_CHECKPOINTS": True,
        }
        return simulation

    def _write_checkpoint_with_schema(
        self,
        schema_version,
        *,
        remove_continuation_profile=False,
        remove_spatial_schedule=False,
    ):
        checkpoint_path = write_checkpoint(
            self._simulation_with_continuation_policy(),
            self.root / "checkpoints",
            trigger="requested",
        )
        with h5py.File(checkpoint_path, "r+") as h5file:
            metadata = h5file["metadata"]
            metadata.attrs["schema_version"] = schema_version
            if (
                remove_continuation_profile
                and "continuation_profile" in metadata
            ):
                del metadata["continuation_profile"]
            if remove_spatial_schedule:
                for output_state in h5file["conductors"].values():
                    del output_state["output_state/Space_save"]
                    del output_state["output_state/i_save_max"]
        return checkpoint_path

    def test_schema_1_2_persists_required_continuation_profile(self):
        self.assertEqual(SCHEMA_VERSION, "1.2")
        self.assertEqual(
            SUPPORTED_SCHEMA_VERSIONS,
            frozenset(("1.1", "1.2")),
        )

        checkpoint_path = write_checkpoint(
            self._simulation_with_continuation_policy(),
            self.root / "checkpoints",
            trigger="requested",
        )

        with h5py.File(checkpoint_path, "r") as h5file:
            profile = h5file["metadata/continuation_profile"]
            self.assertEqual(profile.attrs["profile_version"], "1.1")
            self.assertEqual(profile["time_policy/IADAPTIME"][()], 0)
            self.assertEqual(
                profile["time_policy/TIME_STEP"][()],
                0.025,
            )
            self.assertEqual(
                profile["time_policy/CHECKPOINT_EVERY_N_STEPS"][()],
                50,
            )
            self.assertTrue(
                bool(profile["time_policy/USER_CHECKPOINTS"][()])
            )
            self.assertEqual(len(profile["drivers"]), 0)

        checkpoint = read_checkpoint(checkpoint_path)
        self.assertEqual(checkpoint.schema_version, "1.2")
        self.assertEqual(
            checkpoint.continuation_profile.immutable,
            {
                "conductors": {
                    "COND_1": {
                        "components": {
                            "CHAN_1": {"kind": "fluid"},
                            "STACK_1": {"kind": "solid"},
                        },
                        "inputs": {
                            "ELECTRIC_METHOD": "CN",
                            "METHOD": "AM4",
                        },
                    },
                },
            },
        )
        self.assertEqual(
            checkpoint.continuation_profile.time_policy,
            {
                "CHECKPOINT_EVERY_N_STEPS": 50,
                "IADAPTIME": 0,
                "STPMIN": 0.0025,
                "TEND": 0.4,
                "TIME_STEP": 0.025,
                "USER_CHECKPOINTS": True,
            },
        )
        self.assertEqual(checkpoint.continuation_profile.drivers, {})
        with self.assertRaises(FrozenInstanceError):
            checkpoint.continuation_profile.time_policy = {}

    def test_reader_rejects_development_schema_1_1_without_profile(self):
        checkpoint_path = self._write_checkpoint_with_schema(
            "1.1",
            remove_continuation_profile=True,
        )

        with self.assertRaisesRegex(
            CheckpointReadError,
            r"continuation_profile",
        ):
            read_checkpoint(checkpoint_path)

    def test_reader_accepts_schema_1_1_with_continuation_profile(self):
        checkpoint_path = self._write_checkpoint_with_schema(
            "1.1",
            remove_spatial_schedule=True,
        )

        checkpoint = read_checkpoint(checkpoint_path)

        self.assertEqual(checkpoint.schema_version, "1.1")

    def test_reader_rejects_schema_1_2_without_spatial_schedule(self):
        checkpoint_path = self._write_checkpoint_with_schema(
            "1.2",
            remove_spatial_schedule=True,
        )

        with self.assertRaisesRegex(
            CheckpointReadError,
            r"output state is missing 'Space_save'",
        ):
            read_checkpoint(checkpoint_path)

    def test_reader_rejects_unknown_schema_1_3(self):
        checkpoint_path = self._write_checkpoint_with_schema("1.3")

        with self.assertRaisesRegex(
            CheckpointReadError,
            r"Unsupported checkpoint schema '1\.3'",
        ):
            read_checkpoint(checkpoint_path)


if __name__ == "__main__":
    unittest.main()
