from datetime import datetime, timezone
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from openpyxl import Workbook

from test_checkpoint import make_simulation
from utility_functions.checkpoint import (
    build_continuation_profile,
    compare_continuation_profiles,
)


class CheckpointContinuationXlsxIdentityTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)

    @staticmethod
    def _write_workbook(path, value, modified):
        path.parent.mkdir(parents=True, exist_ok=True)
        workbook = Workbook()
        worksheet = workbook.active
        worksheet.title = "Spatial_distribution"
        worksheet["A1"] = "Variable"
        worksheet["B2"] = value
        workbook.properties.creator = modified.isoformat()
        workbook.properties.modified = modified
        workbook.save(path)
        workbook.close()

    @staticmethod
    def _write_binary(path, content):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)

    @staticmethod
    def _simulation(input_directory, static_name):
        simulation = make_simulation(input_directory)
        simulation.transient_input = {
            "IADAPTIME": 0,
            "TIME_STEP": 0.025,
            "STPMIN": 0.0025,
            "TEND": 0.4,
        }
        conductor = simulation.list_of_Conductors[0]
        conductor.file_input = {"OUTPUT": static_name}
        return simulation

    def test_semantically_equal_xlsx_rewrites_are_compatible(self):
        first_dir = self.root / "first"
        second_dir = self.root / "second"
        first_path = first_dir / "diagnostic.xlsx"
        second_path = second_dir / "diagnostic.xlsx"
        self._write_workbook(
            first_path,
            value=4.5,
            modified=datetime(2026, 8, 16, 10, 0, tzinfo=timezone.utc),
        )
        self._write_workbook(
            second_path,
            value=4.5,
            modified=datetime(2026, 8, 16, 20, 0, tzinfo=timezone.utc),
        )
        self.assertNotEqual(first_path.read_bytes(), second_path.read_bytes())

        comparison = compare_continuation_profiles(
            build_continuation_profile(
                self._simulation(first_dir, "diagnostic.xlsx")
            ),
            build_continuation_profile(
                self._simulation(second_dir, "diagnostic.xlsx")
            ),
        )

        self.assertTrue(comparison.is_compatible)
        self.assertEqual(comparison.immutable_differences, ())

    def test_xlsx_cell_change_remains_blocking(self):
        first_dir = self.root / "first_changed"
        second_dir = self.root / "second_changed"
        self._write_workbook(
            first_dir / "diagnostic.xlsx",
            value=4.5,
            modified=datetime(2026, 8, 16, 10, 0, tzinfo=timezone.utc),
        )
        self._write_workbook(
            second_dir / "diagnostic.xlsx",
            value=4.6,
            modified=datetime(2026, 8, 16, 20, 0, tzinfo=timezone.utc),
        )

        comparison = compare_continuation_profiles(
            build_continuation_profile(
                self._simulation(first_dir, "diagnostic.xlsx")
            ),
            build_continuation_profile(
                self._simulation(second_dir, "diagnostic.xlsx")
            ),
        )

        self.assertFalse(comparison.is_compatible)
        self.assertIn(
            "immutable.conductors.COND_1.static_files.OUTPUT.sha256",
            comparison.immutable_differences,
        )

    def test_non_xlsx_static_files_remain_byte_sensitive(self):
        first_dir = self.root / "first_binary"
        second_dir = self.root / "second_binary"
        self._write_binary(first_dir / "geometry.dat", b"geometry-v1")
        self._write_binary(second_dir / "geometry.dat", b"geometry-v2")

        comparison = compare_continuation_profiles(
            build_continuation_profile(
                self._simulation(first_dir, "geometry.dat")
            ),
            build_continuation_profile(
                self._simulation(second_dir, "geometry.dat")
            ),
        )

        self.assertFalse(comparison.is_compatible)
        self.assertIn(
            "immutable.conductors.COND_1.static_files.OUTPUT.sha256",
            comparison.immutable_differences,
        )


if __name__ == "__main__":
    unittest.main()
