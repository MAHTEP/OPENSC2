import hashlib
import os
from dataclasses import FrozenInstanceError, replace
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import h5py
import numpy as np

from utility_functions.checkpoint import (
    CheckpointReadError,
    CheckpointValidationError,
    apply_checkpoint_to_runtime,
    checkpoint_interval,
    SCHEMA_VERSION,
    build_input_manifest,
    compare_input_manifest,
    evaluate_restart_compatibility,
    read_checkpoint,
    validate_runtime_restore_target,
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
        Space_save=np.array([0.0, 0.1]),
        t_save_left=0.1,
        store_sd_node={"zcoord": {"t_save_left": np.array([0.0, 1.0])}},
        store_sd_gauss={},
        inventory={
            "FluidComponent": _collection(fluid),
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


def make_restore_target(base_path):
    simulation = make_simulation(base_path)
    simulation.simulation_time = [0.0]
    simulation.num_step = 0
    conductor = simulation.list_of_Conductors[0]
    conductor.cond_time = [0.0]
    conductor.cond_num_step = 0
    conductor.electric_time = 0.0
    conductor.cond_el_num_step = 0
    conductor.i_save = 0
    conductor.num_step_save[:] = 0
    return simulation


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
            self.assertEqual(
                conductor["components/CHAN_1"].attrs["kind"], "fluid"
            )
            self.assertEqual(
                conductor["components/STACK_1"].attrs["kind"], "solid"
            )
            self.assertEqual(
                set(conductor["components"]), {"CHAN_1", "STACK_1"}
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

    def test_manifest_comparison_accepts_identical_inputs(self):
        simulation = make_simulation(self.input_dir)
        checkpoint_path = write_checkpoint(
            simulation, self.root / "checkpoints", trigger="periodic"
        )
        checkpoint = read_checkpoint(checkpoint_path)

        comparison = compare_input_manifest(checkpoint, self.input_dir)

        self.assertTrue(comparison.is_match)
        self.assertEqual(comparison.checkpoint_entries, checkpoint.input_manifest)
        self.assertEqual(comparison.current_entries, checkpoint.input_manifest)
        self.assertEqual(comparison.missing, ())
        self.assertEqual(comparison.added, ())
        self.assertEqual(comparison.modified, ())

    def test_manifest_comparison_reports_all_difference_categories(self):
        removed = self.input_dir / "removed.tsv"
        removed.write_bytes(b"removed")
        modified = self.input_dir / "modified.tsv"
        modified.write_bytes(b"before")
        simulation = make_simulation(self.input_dir)
        checkpoint_path = write_checkpoint(
            simulation, self.root / "checkpoints", trigger="periodic"
        )
        checkpoint = read_checkpoint(checkpoint_path)

        removed.unlink()
        modified.write_bytes(b"after-and-longer")
        (self.input_dir / "added.tsv").write_bytes(b"added")

        comparison = compare_input_manifest(checkpoint, self.input_dir)

        self.assertFalse(comparison.is_match)
        self.assertEqual(
            [entry.path for entry in comparison.missing], ["removed.tsv"]
        )
        self.assertEqual([entry.path for entry in comparison.added], ["added.tsv"])
        self.assertEqual(
            [change.checkpoint.path for change in comparison.modified],
            ["modified.tsv"],
        )
        change = comparison.modified[0]
        self.assertEqual(
            change.checkpoint.sha256,
            hashlib.sha256(b"before").hexdigest(),
        )
        self.assertEqual(
            change.current.sha256,
            hashlib.sha256(b"after-and-longer").hexdigest(),
        )
        self.assertEqual(change.checkpoint.size, len(b"before"))
        self.assertEqual(change.current.size, len(b"after-and-longer"))

    def test_manifest_comparison_is_deterministic(self):
        simulation = make_simulation(self.input_dir)
        checkpoint_path = write_checkpoint(
            simulation, self.root / "checkpoints", trigger="periodic"
        )
        checkpoint = read_checkpoint(checkpoint_path)
        (self.input_dir / "z_added.tsv").write_bytes(b"z")
        (self.input_dir / "a_added.tsv").write_bytes(b"a")

        comparison = compare_input_manifest(checkpoint, self.input_dir)

        self.assertEqual(
            [entry.path for entry in comparison.added],
            ["a_added.tsv", "z_added.tsv"],
        )

    def test_manifest_comparison_detects_same_size_content_change(self):
        target = self.input_dir / "same_size.tsv"
        target.write_bytes(b"before")
        simulation = make_simulation(self.input_dir)
        checkpoint_path = write_checkpoint(
            simulation, self.root / "checkpoints", trigger="periodic"
        )
        checkpoint = read_checkpoint(checkpoint_path)

        target.write_bytes(b"after!")
        comparison = compare_input_manifest(checkpoint, self.input_dir)

        self.assertEqual(len(comparison.modified), 1)
        self.assertEqual(comparison.modified[0].checkpoint.size, 6)
        self.assertEqual(comparison.modified[0].current.size, 6)
        self.assertNotEqual(
            comparison.modified[0].checkpoint.sha256,
            comparison.modified[0].current.sha256,
        )

    def test_manifest_comparison_rejects_missing_input_directory(self):
        simulation = make_simulation(self.input_dir)
        checkpoint_path = write_checkpoint(
            simulation, self.root / "checkpoints", trigger="periodic"
        )
        checkpoint = read_checkpoint(checkpoint_path)

        with self.assertRaisesRegex(
            CheckpointValidationError, "input directory does not exist"
        ):
            compare_input_manifest(checkpoint, self.root / "missing")

    def test_manifest_comparison_requires_detached_checkpoint_data(self):
        with self.assertRaisesRegex(TypeError, "CheckpointData"):
            compare_input_manifest(object(), self.input_dir)

    def test_recovery_compatibility_accepts_identical_inputs(self):
        simulation = make_simulation(self.input_dir)
        checkpoint_path = write_checkpoint(
            simulation, self.root / "checkpoints", trigger="periodic"
        )
        checkpoint = read_checkpoint(checkpoint_path)

        report = evaluate_restart_compatibility(
            checkpoint, self.input_dir, mode="recovery"
        )

        self.assertEqual(report.mode, "recovery")
        self.assertTrue(report.is_compatible)
        self.assertTrue(report.manifest_comparison.is_match)
        self.assertEqual(report.blocking_reasons, ())
        self.assertEqual(report.warnings, ())
        with self.assertRaises(FrozenInstanceError):
            report.is_compatible = False

    def test_recovery_compatibility_rejects_missing_input(self):
        checkpoint = self._checkpoint_with_current_inputs()
        (self.input_dir / "transitory_input.xlsx").unlink()

        report = evaluate_restart_compatibility(checkpoint, self.input_dir)

        self.assertFalse(report.is_compatible)
        self.assertEqual(
            report.blocking_reasons,
            ("Input file is missing: transitory_input.xlsx.",),
        )

    def test_recovery_compatibility_rejects_added_input(self):
        checkpoint = self._checkpoint_with_current_inputs()
        (self.input_dir / "added.tsv").write_bytes(b"added")

        report = evaluate_restart_compatibility(checkpoint, self.input_dir)

        self.assertFalse(report.is_compatible)
        self.assertEqual(
            report.blocking_reasons,
            ("Unexpected input file was added: added.tsv.",),
        )

    def test_recovery_compatibility_rejects_modified_input(self):
        checkpoint = self._checkpoint_with_current_inputs()
        (self.input_dir / "transitory_input.xlsx").write_bytes(b"changed")

        report = evaluate_restart_compatibility(checkpoint, self.input_dir)

        self.assertFalse(report.is_compatible)
        self.assertEqual(
            report.blocking_reasons,
            ("Input file was modified: transitory_input.xlsx.",),
        )

    def test_recovery_compatibility_accumulates_all_incompatibilities(self):
        (self.input_dir / "missing.tsv").write_bytes(b"missing")
        (self.input_dir / "modified.tsv").write_bytes(b"before")
        checkpoint = self._checkpoint_with_current_inputs()
        (self.input_dir / "missing.tsv").unlink()
        (self.input_dir / "modified.tsv").write_bytes(b"after")
        (self.input_dir / "added.tsv").write_bytes(b"added")

        report = evaluate_restart_compatibility(checkpoint, self.input_dir)

        self.assertFalse(report.is_compatible)
        self.assertEqual(
            report.blocking_reasons,
            (
                "Input file is missing: missing.tsv.",
                "Unexpected input file was added: added.tsv.",
                "Input file was modified: modified.tsv.",
            ),
        )

    def test_continuation_is_explicitly_blocked(self):
        checkpoint = self._checkpoint_with_current_inputs()

        report = evaluate_restart_compatibility(
            checkpoint, self.input_dir, mode="continuation"
        )

        self.assertFalse(report.is_compatible)
        self.assertTrue(report.manifest_comparison.is_match)
        self.assertEqual(
            report.blocking_reasons,
            ("Restart mode 'continuation' is not supported yet.",),
        )

    def test_unknown_restart_mode_is_rejected(self):
        checkpoint = self._checkpoint_with_current_inputs()

        with self.assertRaisesRegex(ValueError, "Invalid restart mode"):
            evaluate_restart_compatibility(
                checkpoint, self.input_dir, mode="resume"
            )

    def test_compatibility_evaluation_does_not_mutate_checkpoint(self):
        checkpoint = self._checkpoint_with_current_inputs()
        original_manifest = checkpoint.input_manifest
        original_time = checkpoint.simulation_time.copy()

        evaluate_restart_compatibility(checkpoint, self.input_dir)

        self.assertIs(checkpoint.input_manifest, original_manifest)
        np.testing.assert_array_equal(checkpoint.simulation_time, original_time)

    def test_compatibility_evaluation_is_deterministic(self):
        checkpoint = self._checkpoint_with_current_inputs()
        (self.input_dir / "z_added.tsv").write_bytes(b"z")
        (self.input_dir / "a_added.tsv").write_bytes(b"a")

        first = evaluate_restart_compatibility(checkpoint, self.input_dir)
        second = evaluate_restart_compatibility(checkpoint, self.input_dir)

        self.assertEqual(first, second)

    def test_runtime_restore_target_accepts_matching_fresh_runtime(self):
        checkpoint = self._checkpoint_with_current_inputs()
        target = make_restore_target(self.input_dir)

        report = validate_runtime_restore_target(checkpoint, target)

        self.assertTrue(report.is_valid)
        self.assertEqual(report.blocking_reasons, ())
        self.assertEqual(report.warnings, ())
        with self.assertRaises(FrozenInstanceError):
            report.is_valid = False

    def test_runtime_restore_target_requires_detached_checkpoint(self):
        with self.assertRaisesRegex(TypeError, "CheckpointData"):
            validate_runtime_restore_target(
                object(), make_restore_target(self.input_dir)
            )

    def test_runtime_restore_target_reports_inconsistent_checkpoint_clock(self):
        checkpoint = self._checkpoint_with_current_inputs()
        checkpoint = replace(
            checkpoint,
            simulation_time=np.array([0.0, 0.2]),
            num_step=2,
        )

        report = validate_runtime_restore_target(
            checkpoint, make_restore_target(self.input_dir)
        )

        self.assertIn(
            "Checkpoint global num_step does not equal "
            "len(simulation_time) - 1.",
            report.blocking_reasons,
        )
        self.assertIn(
            "Conductor 'COND_1' checkpoint step does not match the global "
            "checkpoint step.",
            report.blocking_reasons,
        )
        self.assertIn(
            "Checkpoint global time does not equal the latest conductor time.",
            report.blocking_reasons,
        )

    def test_runtime_restore_target_rejects_non_fresh_runtime(self):
        checkpoint = self._checkpoint_with_current_inputs()
        target = make_simulation(self.input_dir)

        report = validate_runtime_restore_target(checkpoint, target)

        self.assertFalse(report.is_valid)
        self.assertIn(
            "Runtime simulation num_step must be zero before restore.",
            report.blocking_reasons,
        )
        self.assertIn(
            "Conductor 'COND_1' runtime cond_num_step must be zero before restore.",
            report.blocking_reasons,
        )

    def test_runtime_restore_target_reports_conductor_inventory_mismatch(self):
        checkpoint = self._checkpoint_with_current_inputs()
        target = make_restore_target(self.input_dir)
        target.list_of_Conductors[0].identifier = "OTHER"

        report = validate_runtime_restore_target(checkpoint, target)

        self.assertEqual(
            report.blocking_reasons,
            (
                "Runtime is missing conductor 'COND_1'.",
                "Runtime contains unexpected conductor 'OTHER'.",
            ),
        )

    def test_runtime_restore_target_reports_method_history_and_event_mismatch(self):
        checkpoint = self._checkpoint_with_current_inputs()
        target = make_restore_target(self.input_dir)
        conductor = target.list_of_Conductors[0]
        conductor.inputs["METHOD"] = "BE"
        conductor.inputs["ELECTRIC_METHOD"] = "BE"
        conductor.dict_Step["SYSVAR"] = np.zeros((7, 2))
        conductor.events_time = np.array([0.02, 1.0])

        report = validate_runtime_restore_target(checkpoint, target)

        self.assertFalse(report.is_valid)
        self.assertIn("Conductor 'COND_1' TH method differs", report.blocking_reasons[0])
        self.assertTrue(
            any("SYSVAR" in reason and "shape differs" in reason
                for reason in report.blocking_reasons)
        )
        self.assertIn(
            "Conductor 'COND_1' event timeline differs from the checkpoint.",
            report.blocking_reasons,
        )

    def test_runtime_restore_target_reports_component_kind_mismatch(self):
        checkpoint = self._checkpoint_with_current_inputs()
        target = make_restore_target(self.input_dir)
        conductor = target.list_of_Conductors[0]
        fluid = conductor.inventory["FluidComponent"].collection.pop()
        conductor.inventory["SolidComponent"].collection.insert(0, fluid)

        report = validate_runtime_restore_target(checkpoint, target)

        self.assertIn(
            "Conductor 'COND_1' component 'CHAN_1' kind differs: "
            "checkpoint='fluid', runtime='solid'.",
            report.blocking_reasons,
        )

    def test_runtime_restore_target_reports_electric_activation_mismatch(self):
        checkpoint = self._checkpoint_with_current_inputs()
        target = make_restore_target(self.input_dir)
        target.list_of_Conductors[0].inputs["I0_OP_MODE"] = None

        report = validate_runtime_restore_target(checkpoint, target)

        self.assertIn(
            "Conductor 'COND_1' electric activation differs: "
            "checkpoint=True, runtime=False.",
            report.blocking_reasons,
        )

    def test_runtime_restore_target_reports_output_incompatibilities(self):
        checkpoint = self._checkpoint_with_current_inputs()
        target = make_restore_target(self.input_dir)
        conductor = target.list_of_Conductors[0]
        conductor.Space_save = np.array([0.0])
        del conductor.inventory["SolidComponent"].collection[0].time_evol

        report = validate_runtime_restore_target(checkpoint, target)

        self.assertIn(
            "Conductor 'COND_1' saved i_save is outside runtime Space_save.",
            report.blocking_reasons,
        )
        self.assertTrue(
            any("missing buffer 'time_evol'" in reason
                for reason in report.blocking_reasons)
        )

    def test_runtime_restore_validation_is_deterministic_and_non_mutating(self):
        checkpoint = self._checkpoint_with_current_inputs()
        target = make_restore_target(self.input_dir)
        original_sysvar = target.list_of_Conductors[0].dict_Step["SYSVAR"].copy()
        original_time = list(target.simulation_time)

        first = validate_runtime_restore_target(checkpoint, target)
        second = validate_runtime_restore_target(checkpoint, target)

        self.assertEqual(first, second)
        self.assertEqual(target.simulation_time, original_time)
        np.testing.assert_array_equal(
            target.list_of_Conductors[0].dict_Step["SYSVAR"], original_sysvar
        )

    def test_checkpoint_application_restores_persisted_runtime_state(self):
        checkpoint = self._checkpoint_with_current_inputs()
        target = make_restore_target(self.input_dir)
        conductor = target.list_of_Conductors[0]
        solid = conductor.inventory["SolidComponent"].collection[0]
        coolant = conductor.inventory["FluidComponent"].collection[0].coolant

        conductor.dict_Step["SYSVAR"][:] = -1.0
        conductor.electric_solution[:] = -2.0
        conductor.enthalpy_balance = -3.0
        solid.dict_node_pt["EEXT"][:] = -4.0
        coolant.time_evol["pressure"] = [0.0]
        del conductor.t_save_left

        result = apply_checkpoint_to_runtime(checkpoint, target)

        self.assertIsNone(result)
        self.assertEqual(target.simulation_time, [0.0, 0.1])
        self.assertEqual(target.num_step, 1)
        self.assertEqual(conductor.cond_time, [0.0, 0.1])
        self.assertEqual(conductor.cond_num_step, 1)
        self.assertEqual(conductor.time_step, 0.1)
        np.testing.assert_allclose(conductor.EQTEIG, [1.0, 2.0])
        self.assertEqual(conductor.i_event, 1)
        self.assertEqual(conductor.electric_time, 0.1)
        self.assertEqual(conductor.cond_el_num_step, 10)
        np.testing.assert_allclose(
            conductor.dict_Step["SYSVAR"],
            checkpoint.conductors["COND_1"].th_history["SYSVAR"],
        )
        np.testing.assert_allclose(conductor.electric_solution, [10.0, 20.0])
        self.assertEqual(conductor.enthalpy_balance, 1.0)
        np.testing.assert_allclose(
            solid.dict_node_pt["EEXT"], np.arange(6.0).reshape(3, 2)
        )
        self.assertEqual(conductor.i_save, 1)
        np.testing.assert_array_equal(conductor.num_step_save, [0, 1])
        self.assertEqual(conductor.t_save_left, 0.1)
        self.assertEqual(coolant.time_evol["pressure"], [1.0e5, 1.01e5])

    def test_checkpoint_application_preserves_appendable_history_types(self):
        checkpoint = self._checkpoint_with_current_inputs()
        target = make_restore_target(self.input_dir)
        conductor = target.list_of_Conductors[0]
        coolant = conductor.inventory["FluidComponent"].collection[0].coolant
        original_events = conductor.events_time
        original_space_save = conductor.Space_save

        apply_checkpoint_to_runtime(checkpoint, target)

        self.assertIsInstance(target.simulation_time, list)
        self.assertIsInstance(conductor.cond_time, list)
        self.assertIsInstance(coolant.time_evol["pressure"], list)
        self.assertIsInstance(coolant.time_evol_io["time (s)"], list)
        self.assertIs(conductor.events_time, original_events)
        self.assertIs(conductor.Space_save, original_space_save)
        target.simulation_time.append(0.2)
        conductor.cond_time.append(0.2)
        coolant.time_evol_io["time (s)"].append(0.2)

    def test_checkpoint_application_uses_independent_deep_copies(self):
        checkpoint = self._checkpoint_with_current_inputs()
        target = make_restore_target(self.input_dir)

        apply_checkpoint_to_runtime(checkpoint, target)

        conductor = target.list_of_Conductors[0]
        saved = checkpoint.conductors["COND_1"]
        conductor.dict_Step["SYSVAR"][0, 0] = -99.0
        conductor.inventory["SolidComponent"].collection[0].dict_node_pt[
            "EEXT"
        ][0, 0] = -98.0
        conductor.inventory["FluidComponent"].collection[0].coolant.time_evol[
            "pressure"
        ][0] = -97.0

        self.assertEqual(saved.th_history["SYSVAR"][0, 0], 0.0)
        self.assertEqual(
            saved.components["STACK_1"]["energy_history"]["EEXT"][0, 0],
            0.0,
        )
        self.assertEqual(
            saved.output_state["buffers"]["components"]["CHAN_1"]
            ["coolant"]["time_evol"]["pressure"][0],
            1.0e5,
        )

    def test_checkpoint_application_rejects_invalid_target_before_mutation(self):
        checkpoint = self._checkpoint_with_current_inputs()
        target = make_restore_target(self.input_dir)
        target.num_step = 7
        original_time = list(target.simulation_time)
        original_sysvar = target.list_of_Conductors[0].dict_Step[
            "SYSVAR"
        ].copy()

        with self.assertRaisesRegex(
            CheckpointValidationError, "Runtime restore target is not valid"
        ):
            apply_checkpoint_to_runtime(checkpoint, target)

        self.assertEqual(target.simulation_time, original_time)
        self.assertEqual(target.num_step, 7)
        np.testing.assert_array_equal(
            target.list_of_Conductors[0].dict_Step["SYSVAR"], original_sysvar
        )

    def test_checkpoint_application_rechecks_manifest_before_mutation(self):
        checkpoint = self._checkpoint_with_current_inputs()
        target = make_restore_target(self.input_dir)
        (self.input_dir / "transitory_input.xlsx").write_bytes(b"changed")

        with self.assertRaisesRegex(
            CheckpointValidationError,
            "Checkpoint inputs are not compatible with recovery",
        ):
            apply_checkpoint_to_runtime(checkpoint, target)

        self.assertEqual(target.simulation_time, [0.0])
        self.assertEqual(target.num_step, 0)

    def test_runtime_validation_covers_mutating_restore_destinations(self):
        checkpoint = self._checkpoint_with_current_inputs()
        target = make_restore_target(self.input_dir)
        conductor = target.list_of_Conductors[0]
        solid = conductor.inventory["SolidComponent"].collection[0]
        coolant = conductor.inventory["FluidComponent"].collection[0].coolant
        del conductor.enthalpy_balance
        solid.dict_node_pt["EEXT"] = np.zeros((1, 1))
        coolant.time_evol["runtime_only"] = [0.0]

        report = validate_runtime_restore_target(checkpoint, target)

        self.assertIn(
            "Conductor 'COND_1' is missing balance 'enthalpy_balance'.",
            report.blocking_reasons,
        )
        self.assertTrue(
            any("STACK_1" in reason and "EEXT" in reason
                and "shape differs" in reason
                for reason in report.blocking_reasons)
        )
        self.assertTrue(
            any("runtime_only" in reason and "unexpected runtime key" in reason
                for reason in report.blocking_reasons)
        )

    def test_checkpoint_application_restores_lazily_created_output_keys(self):
        checkpoint = self._checkpoint_with_current_inputs()
        target = make_restore_target(self.input_dir)
        conductor = target.list_of_Conductors[0]
        saved_store = checkpoint.conductors["COND_1"].output_state[
            "buffers"
        ]["conductor"]["store_sd_node"]
        saved_store["htc_ch_sol"] = {
            "t_save": {"CHAN_1_STACK_1": np.array([10.0, 11.0])},
            "t_save_left": {"CHAN_1_STACK_1": np.array([8.0, 9.0])},
        }
        conductor.store_sd_node["htc_ch_sol"] = {
            "t_save": {},
            "t_save_left": {},
        }

        report = validate_runtime_restore_target(checkpoint, target)

        self.assertTrue(report.is_valid, report.blocking_reasons)
        apply_checkpoint_to_runtime(checkpoint, target)
        restored = conductor.store_sd_node["htc_ch_sol"]["t_save"][
            "CHAN_1_STACK_1"
        ]
        np.testing.assert_array_equal(restored, [10.0, 11.0])
        self.assertFalse(
            np.shares_memory(
                restored,
                saved_store["htc_ch_sol"]["t_save"]["CHAN_1_STACK_1"],
            )
        )

    def test_checkpoint_application_accepts_inactive_electric_model(self):
        source = make_simulation(self.input_dir)
        conductor = source.list_of_Conductors[0]
        conductor.inputs["I0_OP_MODE"] = None
        del conductor.electric_solution
        del conductor.electric_solution_steady
        checkpoint_path = write_checkpoint(
            source, self.root / "checkpoints", trigger="periodic"
        )
        checkpoint = read_checkpoint(checkpoint_path)
        source.simulation_time = [0.0]
        source.num_step = 0
        conductor.cond_time = [0.0]
        conductor.cond_num_step = 0
        conductor.i_save = 0
        conductor.num_step_save[:] = 0

        apply_checkpoint_to_runtime(checkpoint, source)

        self.assertEqual(source.simulation_time, [0.0, 0.1])
        self.assertEqual(conductor.cond_time, [0.0, 0.1])
        self.assertFalse(hasattr(conductor, "electric_solution"))
        self.assertFalse(hasattr(conductor, "electric_solution_steady"))

    def _checkpoint_with_current_inputs(self):
        simulation = make_simulation(self.input_dir)
        checkpoint_path = write_checkpoint(
            simulation, self.root / "checkpoints", trigger="periodic"
        )
        return read_checkpoint(checkpoint_path)

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

    def test_reader_returns_detached_intermediate_state(self):
        simulation = make_simulation(self.input_dir)
        checkpoint = write_checkpoint(
            simulation, self.root / "checkpoints", trigger="periodic"
        )

        loaded = read_checkpoint(checkpoint)

        self.assertEqual(loaded.path, checkpoint.resolve())
        self.assertEqual(loaded.schema_version, SCHEMA_VERSION)
        self.assertEqual(loaded.trigger, "periodic")
        self.assertEqual(loaded.num_step, 1)
        self.assertEqual(loaded.input_manifest[0].path, "transitory_input.xlsx")
        self.assertEqual(loaded.input_manifest[0].size, len(b"input-content"))
        np.testing.assert_allclose(loaded.simulation_time, [0.0, 0.1])

        conductor = loaded.conductors["COND_1"]
        self.assertEqual(conductor.identifier, "COND_1")
        self.assertEqual(conductor.th_method, "AM4")
        np.testing.assert_allclose(
            conductor.th_history["AM4_AA"],
            simulation.list_of_Conductors[0].dict_Step["AM4_AA"],
        )
        np.testing.assert_allclose(
            conductor.electric["electric_solution"], [10.0, 20.0]
        )
        np.testing.assert_allclose(
            conductor.components["STACK_1"]["energy_history"]["EJHT"],
            np.arange(6.0, 12.0).reshape(3, 2),
        )
        self.assertEqual(conductor.components["CHAN_1"], {"kind": "fluid"})
        self.assertEqual(conductor.components["STACK_1"]["kind"], "solid")
        np.testing.assert_allclose(
            conductor.output_state["buffers"]["components"]["CHAN_1"]
            ["coolant"]["time_evol"]["pressure"],
            [1.0e5, 1.01e5],
        )

        # The file is closed on return: its contents can be replaced and the
        # detached data remain usable.
        checkpoint.unlink()
        np.testing.assert_allclose(conductor.clock["cond_time"], [0.0, 0.1])

    def test_reader_does_not_mutate_source_simulation(self):
        simulation = make_simulation(self.input_dir)
        checkpoint = write_checkpoint(
            simulation, self.root / "checkpoints", trigger="periodic"
        )
        original_time = list(simulation.simulation_time)
        original_sysvar = simulation.list_of_Conductors[0].dict_Step[
            "SYSVAR"
        ].copy()

        read_checkpoint(checkpoint)

        self.assertEqual(simulation.simulation_time, original_time)
        np.testing.assert_array_equal(
            simulation.list_of_Conductors[0].dict_Step["SYSVAR"],
            original_sysvar,
        )

    def test_reader_accepts_checkpoint_without_active_electric_model(self):
        simulation = make_simulation(self.input_dir)
        conductor = simulation.list_of_Conductors[0]
        conductor.inputs["I0_OP_MODE"] = None
        del conductor.electric_solution
        del conductor.electric_solution_steady
        checkpoint = write_checkpoint(
            simulation, self.root / "checkpoints", trigger="final"
        )

        loaded = read_checkpoint(checkpoint)

        self.assertEqual(loaded.conductors["COND_1"].electric, {})

    def test_reader_rejects_incomplete_checkpoint(self):
        simulation = make_simulation(self.input_dir)
        checkpoint = write_checkpoint(
            simulation, self.root / "checkpoints", trigger="periodic"
        )
        with h5py.File(checkpoint, "r+") as h5file:
            h5file["metadata"].attrs["complete"] = False

        with self.assertRaisesRegex(CheckpointReadError, "incomplete"):
            read_checkpoint(checkpoint)

    def test_reader_rejects_incompatible_schema(self):
        simulation = make_simulation(self.input_dir)
        checkpoint = write_checkpoint(
            simulation, self.root / "checkpoints", trigger="periodic"
        )
        with h5py.File(checkpoint, "r+") as h5file:
            h5file["metadata"].attrs["schema_version"] = "999.0"

        with self.assertRaisesRegex(
            CheckpointReadError, "Unsupported checkpoint schema"
        ):
            read_checkpoint(checkpoint)

    def test_reader_rejects_missing_required_dataset(self):
        simulation = make_simulation(self.input_dir)
        checkpoint = write_checkpoint(
            simulation, self.root / "checkpoints", trigger="periodic"
        )
        with h5py.File(checkpoint, "r+") as h5file:
            del h5file["simulation/num_step"]

        with self.assertRaisesRegex(
            CheckpointReadError, "missing required entries: num_step"
        ):
            read_checkpoint(checkpoint)

    def test_reader_rejects_corrupt_non_hdf5_file(self):
        checkpoint = self.root / "corrupt.h5"
        checkpoint.write_bytes(b"not an HDF5 checkpoint")

        with self.assertRaisesRegex(CheckpointReadError, "Cannot read checkpoint"):
            read_checkpoint(checkpoint)

    def test_reader_rejects_invalid_manifest_hash(self):
        simulation = make_simulation(self.input_dir)
        checkpoint = write_checkpoint(
            simulation, self.root / "checkpoints", trigger="periodic"
        )
        with h5py.File(checkpoint, "r+") as h5file:
            h5file["metadata/input_manifest/sha256"][0] = "invalid"

        with self.assertRaisesRegex(CheckpointReadError, "Invalid SHA-256"):
            read_checkpoint(checkpoint)

    def test_reader_rejects_empty_conductor_group(self):
        simulation = make_simulation(self.input_dir)
        checkpoint = write_checkpoint(
            simulation, self.root / "checkpoints", trigger="periodic"
        )
        with h5py.File(checkpoint, "r+") as h5file:
            del h5file["conductors/COND_1"]

        with self.assertRaisesRegex(CheckpointReadError, "contains no conductors"):
            read_checkpoint(checkpoint)

    def test_reader_rejects_inconsistent_saved_conductor_clock(self):
        simulation = make_simulation(self.input_dir)
        checkpoint = write_checkpoint(
            simulation, self.root / "checkpoints", trigger="periodic"
        )
        with h5py.File(checkpoint, "r+") as h5file:
            h5file["conductors/COND_1/clock/cond_num_step"][...] = 2

        with self.assertRaisesRegex(
            CheckpointReadError, "must equal len\\(cond_time\\) - 1"
        ):
            read_checkpoint(checkpoint)

    def test_validation_rejects_incomplete_component_inventory(self):
        simulation = make_simulation(self.input_dir)
        conductor = simulation.list_of_Conductors[0]
        conductor.inventory["all_component"].collection.pop()

        with self.assertRaisesRegex(
            CheckpointValidationError,
            "all_component must contain exactly the fluid and solid components",
        ):
            validate_checkpoint_state(simulation)

    def test_reader_rejects_invalid_component_kind(self):
        simulation = make_simulation(self.input_dir)
        checkpoint = write_checkpoint(
            simulation, self.root / "checkpoints", trigger="periodic"
        )
        with h5py.File(checkpoint, "r+") as h5file:
            h5file["conductors/COND_1/components/CHAN_1"].attrs["kind"] = "gas"

        with self.assertRaisesRegex(CheckpointReadError, "invalid kind"):
            read_checkpoint(checkpoint)

    def test_reader_rejects_component_buffer_mismatch(self):
        simulation = make_simulation(self.input_dir)
        checkpoint = write_checkpoint(
            simulation, self.root / "checkpoints", trigger="periodic"
        )
        with h5py.File(checkpoint, "r+") as h5file:
            del h5file[
                "conductors/COND_1/output_state/buffers/components/CHAN_1"
            ]

        with self.assertRaisesRegex(
            CheckpointReadError, "component inventory does not match"
        ):
            read_checkpoint(checkpoint)


if __name__ == "__main__":
    unittest.main()
