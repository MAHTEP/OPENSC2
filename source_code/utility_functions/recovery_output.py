"""Carry valid on-disk output history into a recovery run.

Recovery resumes the same logical trajectory represented by its checkpoint.
When the target output directory is new, already-flushed history must be
copied from the run containing the checkpoint before the solver appends the
restored in-memory buffers and newly computed values.
"""

from collections.abc import Mapping
from pathlib import Path
import re
from shutil import copy2

import numpy as np
import pandas as pd


_RAW_SPATIAL_OUTPUT = re.compile(
    r"_\((?P<step>\d+)\)_(?:gauss_)?sd\.tsv$"
)
_TIME_COLUMNS = frozenset(("time", "time (s)"))


def seed_recovery_output_history(checkpoint, simulation):
    """Seed a fresh recovery directory with history valid at the checkpoint.

    The source run is inferred from the checkpoint location:
    ``<run>/Checkpoints/<checkpoint>.h5``. Spatial files are selected using
    the saved ``i_save`` cursor and ``num_step_save`` values, so files created
    by a failed source attempt *after* the checkpoint are not imported. Time
    histories are truncated before the earliest restored buffer entry; the
    buffer will later be flushed by the resumed solver without duplication.
    """

    checkpoint_path = Path(checkpoint.path)
    source_run_directory = checkpoint_path.parent.parent
    if checkpoint_path.parent.name != "Checkpoints":
        raise ValueError(
            "Cannot infer the recovery output source from checkpoint path "
            f"{checkpoint_path!s}; expected a Checkpoints directory."
        )

    runtime_by_identifier = {
        conductor.identifier: conductor
        for conductor in simulation.list_of_Conductors
    }
    checkpoint_identifiers = set(checkpoint.conductors)
    runtime_identifiers = set(runtime_by_identifier)
    if checkpoint_identifiers != runtime_identifiers:
        raise ValueError(
            "Checkpoint and runtime conductor inventories differ while "
            "preparing recovery output history."
        )

    reports = []
    for identifier in sorted(checkpoint_identifiers):
        saved = checkpoint.conductors[identifier]
        runtime = runtime_by_identifier[identifier]
        checkpoint_time = float(np.asarray(saved.clock["cond_time"])[-1])

        saved_i_save = int(saved.output_state["i_save"])
        saved_steps = np.asarray(
            saved.output_state["num_step_save"],
            dtype=int,
        )
        valid_spatial_steps = tuple(
            int(step) for step in saved_steps[:saved_i_save]
        )

        source_spatial = (
            source_run_directory
            / "Output"
            / "Spatial_distribution"
            / identifier
        )
        target_spatial = Path(
            simulation.dict_path[
                f"Output_Spatial_distribution_{identifier}_dir"
            ]
        )
        source_time = (
            source_run_directory
            / "Output"
            / "Time_evolution"
            / identifier
        )
        target_time = Path(
            simulation.dict_path[f"Output_Time_evolution_{identifier}_dir"]
        )

        _require_distinct_directories(source_spatial, target_spatial)
        _require_distinct_directories(source_time, target_time)

        spatial_count = _copy_valid_spatial_prefix(
            source_spatial,
            target_spatial,
            valid_spatial_steps,
        )
        buffered_time_start = _earliest_buffered_time(
            saved.output_state["buffers"]
        )
        time_count = _copy_valid_time_prefix(
            source_time,
            target_time,
            checkpoint_time,
            buffered_time_start,
        )
        reports.append(
            {
                "identifier": runtime.identifier,
                "checkpoint_time": checkpoint_time,
                "spatial_steps": valid_spatial_steps,
                "spatial_files": spatial_count,
                "time_files": time_count,
                "buffered_time_start": buffered_time_start,
            }
        )

    return tuple(reports)


def _require_distinct_directories(source, target):
    if source.resolve() == target.resolve():
        raise ValueError(
            "Recovery output history cannot be seeded in place. Select a "
            "new output directory so the interrupted run remains unchanged."
        )


def _copy_valid_spatial_prefix(source, target, valid_steps):
    if not source.is_dir():
        raise FileNotFoundError(
            f"Source spatial-output directory not found: {source!s}"
        )
    target.mkdir(parents=True, exist_ok=True)

    valid_step_set = set(valid_steps)
    selected_files = []
    copied_steps = set()
    for source_path in sorted(source.iterdir()):
        if not source_path.is_file():
            continue
        match = _RAW_SPATIAL_OUTPUT.search(source_path.name)
        if match is None:
            continue
        step = int(match.group("step"))
        if step not in valid_step_set:
            continue

        target_path = target / source_path.name
        if target_path.exists():
            raise FileExistsError(
                "Recovery target already contains spatial output "
                f"{target_path!s}."
            )
        selected_files.append((source_path, target_path))
        copied_steps.add(step)

    missing_steps = valid_step_set - copied_steps
    if missing_steps:
        raise FileNotFoundError(
            "Source run is missing raw spatial outputs for checkpoint steps "
            f"{sorted(missing_steps)!r}."
        )

    for source_path, target_path in selected_files:
        copy2(source_path, target_path)
    return len(selected_files)


def _copy_valid_time_prefix(
    source,
    target,
    checkpoint_time,
    buffered_time_start,
):
    if not source.is_dir():
        raise FileNotFoundError(
            f"Source time-output directory not found: {source!s}"
        )
    target.mkdir(parents=True, exist_ok=True)

    source_files = {
        path.name: path
        for path in source.glob("*.tsv")
        if path.is_file()
    }
    target_files = {
        path.name: path
        for path in target.glob("*.tsv")
        if path.is_file()
    }
    missing_sources = set(target_files) - set(source_files)
    if missing_sources:
        raise FileNotFoundError(
            "Interrupted run is missing time-output files expected by the "
            f"recovery runtime: {sorted(missing_sources)!r}."
        )

    prepared_files = []
    for file_name, source_path in sorted(source_files.items()):
        source_values = pd.read_csv(source_path, delimiter="\t")
        if source_values.columns.empty:
            raise ValueError(
                f"Time-output file has no headings: {source_path!s}"
            )
        time_column = source_values.columns[0]
        if time_column not in _TIME_COLUMNS:
            raise ValueError(
                f"Unexpected time heading {time_column!r} in "
                f"{source_path!s}."
            )

        target_path = target / file_name
        if target_path.exists():
            target_values = pd.read_csv(
                target_path,
                delimiter="\t",
            )
            if target_values.columns.tolist() != source_values.columns.tolist():
                raise ValueError(
                    "Recovery and source time-output headings differ for "
                    f"{file_name!r}."
                )
            if not target_values.empty:
                raise FileExistsError(
                    "Recovery target time-output file is not empty: "
                    f"{target_path!s}. Select a new output directory."
                )

        times = pd.to_numeric(source_values[time_column], errors="raise")
        if buffered_time_start is None:
            keep = (times < checkpoint_time) | np.isclose(
                times,
                checkpoint_time,
                rtol=1.0e-12,
                atol=1.0e-12,
            )
        else:
            keep = (times < buffered_time_start) & ~np.isclose(
                times,
                buffered_time_start,
                rtol=1.0e-12,
                atol=1.0e-12,
            )
        keep &= (times < checkpoint_time) | np.isclose(
            times,
            checkpoint_time,
            rtol=1.0e-12,
            atol=1.0e-12,
        )

        prepared_files.append((target_path, source_values.loc[keep]))

    for target_path, values in prepared_files:
        values.to_csv(target_path, sep="\t", index=False)
    return len(prepared_files)


def _earliest_buffered_time(buffers):
    times = []

    def collect(value):
        if not isinstance(value, Mapping):
            return
        for key, item in value.items():
            if key in _TIME_COLUMNS and not isinstance(item, Mapping):
                try:
                    values = np.asarray(item, dtype=float).reshape(-1)
                except (TypeError, ValueError):
                    continue
                finite = values[np.isfinite(values)]
                if finite.size:
                    times.append(float(np.min(finite)))
            else:
                collect(item)

    collect(buffers)
    return min(times) if times else None
