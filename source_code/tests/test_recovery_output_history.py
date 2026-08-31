from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

import numpy as np
import pandas as pd

from utility_functions.recovery_output import seed_recovery_output_history


class RecoveryOutputHistoryTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        self.source_run = self.root / "source" / "BE" / "case"
        self.checkpoint_directory = self.source_run / "Checkpoints"
        self.checkpoint_directory.mkdir(parents=True)
        self.checkpoint_path = self.checkpoint_directory / "checkpoint.h5"
        self.checkpoint_path.write_bytes(b"checkpoint")

        self.source_spatial = (
            self.source_run
            / "Output"
            / "Spatial_distribution"
            / "COND_1"
        )
        self.source_time = (
            self.source_run / "Output" / "Time_evolution" / "COND_1"
        )
        self.source_spatial.mkdir(parents=True)
        self.source_time.mkdir(parents=True)

        self.target_spatial = self.root / "target" / "spatial"
        self.target_time = self.root / "target" / "time"
        self.target_spatial.mkdir(parents=True)
        self.target_time.mkdir(parents=True)

    @staticmethod
    def _write_time_file(path, times):
        pd.DataFrame(
            {
                "time (s)": times,
                "zcoord = 0.0 (m)": np.asarray(times) + 1.0,
            }
        ).to_csv(path, sep="\t", index=False)

    def _checkpoint(self, *, buffers):
        saved = SimpleNamespace(
            clock={"cond_time": [0.0, 10.0]},
            output_state={
                "i_save": 3,
                "num_step_save": np.asarray(
                    [710, 810, 910, 0, 0],
                    dtype=int,
                ),
                "buffers": buffers,
            },
        )
        return SimpleNamespace(
            path=self.checkpoint_path,
            conductors={"COND_1": saved},
        )

    def _simulation(self):
        return SimpleNamespace(
            list_of_Conductors=[SimpleNamespace(identifier="COND_1")],
            dict_path={
                "Output_Spatial_distribution_COND_1_dir": str(
                    self.target_spatial
                ),
                "Output_Time_evolution_COND_1_dir": str(self.target_time),
            },
        )

    def test_seeds_only_spatial_files_recorded_before_checkpoint(self):
        for step in (710, 810, 910, 1010):
            (self.source_spatial / f"CHAN_1_({step})_sd.tsv").write_text(
                f"step {step}",
                encoding="utf-8",
            )
            (
                self.source_spatial
                / f"STACK_1_({step})_gauss_sd.tsv"
            ).write_text(
                f"gauss step {step}",
                encoding="utf-8",
            )

        self._write_time_file(
            self.source_time / "CHAN_1_temperature_te.tsv",
            [0.0, 8.0, 9.0, 10.0, 11.0],
        )
        self._write_time_file(
            self.target_time / "CHAN_1_temperature_te.tsv",
            [],
        )

        checkpoint = self._checkpoint(
            buffers={
                "components": {
                    "CHAN_1": {
                        "coolant": {
                            "time_evol": {
                                "temperature": {
                                    "time (s)": [9.0, 10.0],
                                }
                            }
                        }
                    }
                }
            }
        )
        reports = seed_recovery_output_history(
            checkpoint,
            self._simulation(),
        )

        for step in (710, 810, 910):
            self.assertTrue(
                (self.target_spatial / f"CHAN_1_({step})_sd.tsv").is_file()
            )
            self.assertTrue(
                (
                    self.target_spatial
                    / f"STACK_1_({step})_gauss_sd.tsv"
                ).is_file()
            )
        self.assertFalse(
            (self.target_spatial / "CHAN_1_(1010)_sd.tsv").exists()
        )
        self.assertEqual(reports[0]["spatial_steps"], (710, 810, 910))
        self.assertEqual(reports[0]["spatial_files"], 6)

    def test_time_prefix_ends_before_earliest_restored_buffer_value(self):
        for step in (710, 810, 910):
            (self.source_spatial / f"CHAN_1_({step})_sd.tsv").write_text(
                f"step {step}",
                encoding="utf-8",
            )
        self._write_time_file(
            self.source_time / "CHAN_1_temperature_te.tsv",
            [0.0, 8.0, 9.0, 10.0, 11.0],
        )
        self._write_time_file(
            self.target_time / "CHAN_1_temperature_te.tsv",
            [],
        )

        checkpoint = self._checkpoint(
            buffers={
                "components": {
                    "CHAN_1": {
                        "coolant": {
                            "time_evol": {
                                "temperature": {
                                    "time (s)": [9.0, 10.0],
                                }
                            }
                        }
                    }
                }
            }
        )
        seed_recovery_output_history(checkpoint, self._simulation())

        restored = pd.read_csv(
            self.target_time / "CHAN_1_temperature_te.tsv",
            delimiter="\t",
        )
        np.testing.assert_allclose(restored["time (s)"], [0.0, 8.0])

        buffered = pd.DataFrame(
            {
                "time (s)": [9.0, 10.0],
                "zcoord = 0.0 (m)": [10.0, 11.0],
            }
        )
        buffered.to_csv(
            self.target_time / "CHAN_1_temperature_te.tsv",
            sep="\t",
            mode="a",
            header=False,
            index=False,
        )
        complete = pd.read_csv(
            self.target_time / "CHAN_1_temperature_te.tsv",
            delimiter="\t",
        )
        np.testing.assert_allclose(
            complete["time (s)"],
            [0.0, 8.0, 9.0, 10.0],
        )

    def test_empty_buffers_keep_history_through_checkpoint(self):
        for step in (710, 810, 910):
            (self.source_spatial / f"CHAN_1_({step})_sd.tsv").write_text(
                f"step {step}",
                encoding="utf-8",
            )
        self._write_time_file(
            self.source_time / "CHAN_1_temperature_te.tsv",
            [0.0, 9.0, 10.0, 11.0],
        )
        self._write_time_file(
            self.target_time / "CHAN_1_temperature_te.tsv",
            [],
        )

        checkpoint = self._checkpoint(buffers={})
        seed_recovery_output_history(checkpoint, self._simulation())

        restored = pd.read_csv(
            self.target_time / "CHAN_1_temperature_te.tsv",
            delimiter="\t",
        )
        np.testing.assert_allclose(restored["time (s)"], [0.0, 9.0, 10.0])

    def test_rejects_a_reused_recovery_time_directory(self):
        for step in (710, 810, 910):
            (self.source_spatial / f"CHAN_1_({step})_sd.tsv").write_text(
                f"step {step}",
                encoding="utf-8",
            )
        self._write_time_file(
            self.source_time / "CHAN_1_temperature_te.tsv",
            [0.0, 9.0, 10.0],
        )
        self._write_time_file(
            self.target_time / "CHAN_1_temperature_te.tsv",
            [0.0],
        )

        with self.assertRaisesRegex(
            FileExistsError,
            "Select a new output directory",
        ):
            seed_recovery_output_history(
                self._checkpoint(buffers={}),
                self._simulation(),
            )


if __name__ == "__main__":
    unittest.main()
