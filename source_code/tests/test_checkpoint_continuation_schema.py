from dataclasses import FrozenInstanceError
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import h5py

from test_checkpoint import make_simulation
from utility_functions.checkpoint import (
    SCHEMA_VERSION,
    evaluate_restart_compatibility,
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

    def _write_legacy_schema_checkpoint(self):
        checkpoint_path = write_checkpoint(
            self._simulation_with_continuation_policy(),
            self.root / "checkpoints",
            trigger="requested",
        )
        with h5py.File(checkpoint_path, "r+") as h5file:
            metadata = h5file["metadata"]
            metadata.attrs["schema_version"] = "1.1"
            if "continuation_profile" in metadata:
                del metadata["continuation_profile"]
        return checkpoint_path

    def test_schema_1_2_persists_detached_continuation_profile(self):
        self.assertEqual(SCHEMA_VERSION, "1.2")

        checkpoint_path = write_checkpoint(
            self._simulation_with_continuation_policy(),
            self.root / "checkpoints",
            trigger="requested",
        )

        with h5py.File(checkpoint_path, "r") as h5file:
            profile = h5file["metadata/continuation_profile"]
            self.assertEqual(profile.attrs["profile_version"], "1.0")
            self.assertEqual(profile["immutable/IADAPTIME"][()], 0)
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
                "IADAPTIME": 0,
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
                "STPMIN": 0.0025,
                "TEND": 0.4,
                "TIME_STEP": 0.025,
                "USER_CHECKPOINTS": True,
            },
        )
        self.assertEqual(checkpoint.continuation_profile.drivers, {})
        with self.assertRaises(FrozenInstanceError):
            checkpoint.continuation_profile.time_policy = {}

    def test_reader_accepts_legacy_schema_1_1_for_recovery(self):
        checkpoint = read_checkpoint(
            self._write_legacy_schema_checkpoint()
        )

        self.assertEqual(checkpoint.schema_version, "1.1")
        self.assertIsNone(checkpoint.continuation_profile)

        report = evaluate_restart_compatibility(
            checkpoint,
            self.input_dir,
            mode="recovery",
        )
        self.assertTrue(report.is_compatible)
        self.assertEqual(report.blocking_reasons, ())

    def test_legacy_schema_1_1_cannot_start_continuation(self):
        checkpoint = read_checkpoint(
            self._write_legacy_schema_checkpoint()
        )

        report = evaluate_restart_compatibility(
            checkpoint,
            self.input_dir,
            mode="continuation",
        )

        self.assertFalse(report.is_compatible)
        self.assertEqual(len(report.blocking_reasons), 1)
        self.assertRegex(
            report.blocking_reasons[0],
            r"schema 1\.1.*continuation profile",
        )


if __name__ == "__main__":
    unittest.main()
