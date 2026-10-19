from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from openpyxl import Workbook

from utility_functions.checkpoint_schedule import (
    CheckpointBoundary,
    load_checkpoint_schedule,
)


class CheckpointScheduleExcelTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        self.workbook_path = self.root / "transitory_input.xlsx"

    def write_workbook(
        self,
        *,
        sheet_name="CHECKPOINTS",
        header="Time (s)",
        times=(),
        extra_cells=(),
    ):
        workbook = Workbook()
        worksheet = workbook.active
        worksheet.title = sheet_name
        if header is not None:
            worksheet["A1"] = header
        for row_index, value in enumerate(times, start=2):
            worksheet.cell(row=row_index, column=1, value=value)
        for coordinate, value in extra_cells:
            worksheet[coordinate] = value
        workbook.save(self.workbook_path)
        workbook.close()
        return self.workbook_path

    def load(self, transient_input, *, epsilon=1.0e-6):
        return load_checkpoint_schedule(
            self.workbook_path,
            transient_input,
            epsilon=epsilon,
        )

    def test_enabled_mode_reads_checkpoint_times_from_workbook(self):
        self.write_workbook(times=(0.2, 0.4))

        schedule = self.load(
            {"USER_CHECKPOINTS": True, "TEND": 0.65}
        )

        self.assertEqual(
            schedule.boundaries,
            (
                CheckpointBoundary(time=0.2, trigger="requested"),
                CheckpointBoundary(time=0.4, trigger="requested"),
                CheckpointBoundary(time=0.65, trigger="final"),
            ),
        )

    def test_enabled_header_only_sheet_schedules_final_checkpoint(self):
        self.write_workbook()

        schedule = self.load(
            {"USER_CHECKPOINTS": "TRUE", "TEND": 0.65}
        )

        self.assertEqual(
            schedule.boundaries,
            (CheckpointBoundary(time=0.65, trigger="final"),),
        )

    def test_enabled_mode_rejects_missing_checkpoints_sheet(self):
        self.write_workbook(sheet_name="OTHER")

        with self.assertRaisesRegex(ValueError, "CHECKPOINTS.*required"):
            self.load({"USER_CHECKPOINTS": True, "TEND": 0.65})

    def test_disabled_mode_allows_legacy_workbook_without_sheet(self):
        self.write_workbook(sheet_name="TRANSIENT")

        schedule = self.load(
            {"USER_CHECKPOINTS": False, "TEND": 0.65}
        )

        self.assertFalse(schedule.user_enabled)
        self.assertEqual(schedule.boundaries, ())

    def test_disabled_mode_warns_when_populated_sheet_is_ignored(self):
        self.write_workbook(times=(0.2, 0.4))

        with self.assertWarnsRegex(UserWarning, "CHECKPOINTS.*ignored"):
            schedule = self.load(
                {"USER_CHECKPOINTS": False, "TEND": 0.65}
            )

        self.assertFalse(schedule.user_enabled)
        self.assertEqual(schedule.boundaries, ())

    def test_reader_rejects_incorrect_column_header(self):
        self.write_workbook(header="Checkpoint time", times=(0.2,))

        with self.assertRaisesRegex(ValueError, "Time \\(s\\)"):
            self.load({"USER_CHECKPOINTS": True, "TEND": 0.65})

    def test_reader_ignores_blank_rows_inside_time_column(self):
        self.write_workbook(times=(0.2, None, 0.4))

        schedule = self.load(
            {"USER_CHECKPOINTS": True, "TEND": 0.65}
        )

        self.assertEqual(
            schedule.boundaries,
            (
                CheckpointBoundary(time=0.2, trigger="requested"),
                CheckpointBoundary(time=0.4, trigger="requested"),
                CheckpointBoundary(time=0.65, trigger="final"),
            ),
        )

    def test_reader_rejects_content_outside_time_column(self):
        self.write_workbook(
            times=(0.2,),
            extra_cells=(("B2", "unexpected"),),
        )

        with self.assertRaisesRegex(
            ValueError,
            "CHECKPOINTS.*column A",
        ):
            self.load({"USER_CHECKPOINTS": True, "TEND": 0.65})

    def test_reader_rejects_missing_workbook(self):
        missing_path = self.root / "missing_transitory_input.xlsx"

        with self.assertRaisesRegex(
            FileNotFoundError,
            "checkpoint schedule workbook does not exist",
        ):
            load_checkpoint_schedule(
                missing_path,
                {"USER_CHECKPOINTS": True, "TEND": 0.65},
                epsilon=1.0e-6,
            )


if __name__ == "__main__":
    unittest.main()
