from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

import numpy as np
import pandas as pd

from utility_functions.output import (
    reorganize_heat_sd,
    reorganize_spatial_distribution,
)


class SpatialOutputReorganizationTests(unittest.TestCase):
    @staticmethod
    def write_fluid_output(output_directory, step, offset):
        columns = [
            "zcoord (m)",
            "velocity (m/s)",
            "pressure (Pa)",
            "temperature (K)",
            "total_density (kg/m^3)",
            "friction_factor (~)",
        ]
        values = np.asarray(
            [
                [0.0, 1.0, 2.0, offset + 3.0, 4.0, 5.0],
                [1.0, 6.0, 7.0, offset + 8.0, 9.0, 10.0],
            ]
        )
        pd.DataFrame(values, columns=columns).to_csv(
            output_directory / f"CHAN_1_({step})_sd.tsv",
            sep="\t",
            index=False,
        )

    @staticmethod
    def make_fluid_conductor():
        return SimpleNamespace(
            Space_save=np.asarray([0.0, 0.9, 1.0]),
            num_step_save=np.asarray([0, 9, 10], dtype=int),
            inventory={
                "FluidComponent": SimpleNamespace(
                    collection=[SimpleNamespace(identifier="CHAN_1")]
                ),
                "SolidComponent": SimpleNamespace(collection=[]),
            },
        )

    def test_reorganization_uses_only_files_available_after_restart(self):
        conductor = self.make_fluid_conductor()

        with TemporaryDirectory() as temporary_directory:
            output_directory = Path(temporary_directory)
            for step, offset in ((9, 90.0), (10, 100.0)):
                self.write_fluid_output(output_directory, step, offset)

            self.assertFalse(
                (output_directory / "CHAN_1_(0)_sd.tsv").exists()
            )

            reorganize_spatial_distribution(
                conductor,
                str(output_directory),
                n_digit_time=6,
            )

            zcoord = pd.read_csv(
                output_directory / "zcoord.tsv",
                sep="\t",
            )
            temperature = pd.read_csv(
                output_directory / "CHAN_1_temperature_sd.tsv",
                sep="\t",
            )
            expected_columns = ["time = 0.9 (s)", "time = 1.0 (s)"]
            self.assertEqual(list(zcoord.columns), expected_columns)
            self.assertEqual(list(temperature.columns), expected_columns)
            np.testing.assert_allclose(
                temperature["time = 0.9 (s)"],
                [93.0, 98.0],
            )
            np.testing.assert_allclose(
                temperature["time = 1.0 (s)"],
                [103.0, 108.0],
            )

    def test_reorganization_preserves_complete_continuous_history(self):
        conductor = self.make_fluid_conductor()
        with TemporaryDirectory() as temporary_directory:
            output_directory = Path(temporary_directory)
            for step, offset in ((0, 0.0), (9, 90.0), (10, 100.0)):
                self.write_fluid_output(output_directory, step, offset)

            reorganize_spatial_distribution(
                conductor,
                str(output_directory),
                n_digit_time=6,
            )

            temperature = pd.read_csv(
                output_directory / "CHAN_1_temperature_sd.tsv",
                sep="\t",
            )
            self.assertEqual(
                list(temperature.columns),
                [
                    "time = 0.0 (s)",
                    "time = 0.9 (s)",
                    "time = 1.0 (s)",
                ],
            )

    def test_heat_reorganization_can_start_after_checkpoint(self):
        conductor = self.make_fluid_conductor()
        with TemporaryDirectory() as temporary_directory:
            output_directory = Path(temporary_directory)
            for step, values in (
                (9, [9.0, 9.5]),
                (10, [10.0, 10.5]),
            ):
                pd.DataFrame({"JACKET_1": values}).to_csv(
                    output_directory / f"Heat_rad_inner_({step})_sd.tsv",
                    sep="\t",
                    index=False,
                )

            reorganize_heat_sd(
                conductor,
                str(output_directory),
                "Heat_rad_inner",
                "Heat_rad",
                n_digit_time=6,
            )

            heat = pd.read_csv(
                output_directory / "Heat_rad_JACKET_1_sd.tsv",
                sep="\t",
            )
            self.assertEqual(
                list(heat.columns),
                ["time = 0.9 (s)", "time = 1.0 (s)"],
            )
            np.testing.assert_allclose(
                heat["time = 0.9 (s)"],
                [9.0, 9.5],
            )
            np.testing.assert_allclose(
                heat["time = 1.0 (s)"],
                [10.0, 10.5],
            )


if __name__ == "__main__":
    unittest.main()
