"""HDF5 checkpoint writing utilities.

This module intentionally implements only the persistence side of restart.
Restoring a checkpoint into live OPENSC2 objects is handled by a later
increment, after the initialization path has been made restart-aware.
"""

from collections.abc import Mapping
from datetime import datetime, timezone
import hashlib
import os
from pathlib import Path
import subprocess
from urllib.parse import quote

import h5py
import numpy as np


SCHEMA_VERSION = "1.0"
VALID_TRIGGERS = frozenset(("periodic", "requested", "final"))
DEFAULT_CHECKPOINT_EVERY_N_STEPS = 100
CHECKPOINT_INTERVAL_INPUT = "CHECKPOINT_EVERY_N_STEPS"

_BALANCE_ATTRIBUTES = (
    "enthalpy_balance",
    "enthalpy_inl",
    "enthalpy_out",
    "E_sol_ini",
    "E_str_ini",
    "E_jk_ini",
)

_OUTPUT_STATE_ATTRIBUTES = (
    "i_save",
    "num_step_save",
    "t_save_left",
)

_OUTPUT_BUFFER_ATTRIBUTES = (
    "store_sd_node",
    "store_sd_gauss",
    "time_evol",
    "time_evol_gauss",
    "time_evol_io",
    "time_evol_max_temperature",
)


class CheckpointValidationError(ValueError):
    """Raised when a simulation is not at a safe checkpoint boundary."""


def checkpoint_interval(transient_input):
    """Return the configured periodic checkpoint interval.

    A missing input keeps the historical default of one checkpoint every 100
    global TH iterations.  Zero disables periodic checkpoints.
    """

    raw_value = transient_input.get(
        CHECKPOINT_INTERVAL_INPUT,
        DEFAULT_CHECKPOINT_EVERY_N_STEPS,
    )
    if isinstance(raw_value, (bool, np.bool_)):
        raise ValueError(
            f"{CHECKPOINT_INTERVAL_INPUT} must be a non-negative integer."
        )

    try:
        numeric_value = float(raw_value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"{CHECKPOINT_INTERVAL_INPUT} must be a non-negative integer."
        ) from exc

    if (
        not np.isfinite(numeric_value)
        or numeric_value < 0.0
        or not numeric_value.is_integer()
    ):
        raise ValueError(
            f"{CHECKPOINT_INTERVAL_INPUT} must be a non-negative integer."
        )
    return int(numeric_value)


def write_periodic_checkpoint_if_due(simulation):
    """Write a periodic checkpoint when the global step reaches its interval.

    Return the written path, or ``None`` when periodic checkpointing is
    disabled or the current global iteration is not a checkpoint boundary.
    """

    interval = checkpoint_interval(simulation.transient_input)
    if interval == 0 or simulation.num_step % interval != 0:
        return None

    try:
        checkpoint_dir = simulation.dict_path["Checkpoint_dir"]
    except KeyError as exc:
        raise CheckpointValidationError(
            "Missing simulation.dict_path['Checkpoint_dir']."
        ) from exc

    return write_checkpoint(simulation, checkpoint_dir, trigger="periodic")


def write_checkpoint(simulation, checkpoint_dir, trigger):
    """Write an atomic HDF5 checkpoint and return its final path.

    The caller must invoke this function only after every conductor has
    completed its current TH iteration, including output updates and the
    local call to ``get_time_step``.
    """

    validate_checkpoint_state(simulation)
    if trigger not in VALID_TRIGGERS:
        raise ValueError(
            f"Invalid checkpoint trigger {trigger!r}; expected one of "
            f"{sorted(VALID_TRIGGERS)}."
        )

    checkpoint_dir = Path(checkpoint_dir)
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    filename = f"checkpoint_step_{int(simulation.num_step):06d}.h5"
    final_path = checkpoint_dir / filename
    temporary_path = checkpoint_dir / f"{filename}.tmp"

    manifest = build_input_manifest(simulation, excluded_dir=checkpoint_dir)

    with h5py.File(temporary_path, "w") as h5file:
        _write_metadata(h5file, simulation, trigger, manifest)
        _write_simulation_state(h5file, simulation)
        h5file.flush()

    # Mark completion only after the whole payload has been flushed and the
    # first file handle has been closed.
    with h5py.File(temporary_path, "r+") as h5file:
        h5file["metadata"].attrs["complete"] = True
        h5file.flush()

    os.replace(temporary_path, final_path)
    return final_path


def validate_checkpoint_state(simulation):
    """Validate invariants required by the selected checkpoint boundary."""

    simulation_time = _finite_time_history(
        getattr(simulation, "simulation_time", None), "simulation_time"
    )
    if not hasattr(simulation, "num_step"):
        raise CheckpointValidationError("Missing simulation.num_step.")
    if not isinstance(simulation.num_step, (int, np.integer)):
        raise CheckpointValidationError("simulation.num_step must be an integer.")
    if simulation.num_step < 0:
        raise CheckpointValidationError("simulation.num_step cannot be negative.")

    conductors = list(getattr(simulation, "list_of_Conductors", ()))
    if not conductors:
        raise CheckpointValidationError("The simulation has no conductors.")

    identifiers = [getattr(item, "identifier", None) for item in conductors]
    if any(not identifier for identifier in identifiers):
        raise CheckpointValidationError("Every conductor needs an identifier.")
    if len(set(identifiers)) != len(identifiers):
        raise CheckpointValidationError("Conductor identifiers must be unique.")

    for conductor in conductors:
        _validate_conductor(conductor)

    return simulation_time


def build_input_manifest(simulation, excluded_dir=None):
    """Return deterministic SHA-256 entries for files under ``basePath``."""

    base_path = Path(getattr(simulation, "basePath", ""))
    if not base_path.is_dir():
        raise CheckpointValidationError(
            f"Simulation input directory does not exist: {base_path!s}."
        )

    base_path = base_path.resolve()
    excluded_dir = Path(excluded_dir).resolve() if excluded_dir else None
    manifest = []
    for path in sorted(base_path.rglob("*")):
        if not path.is_file():
            continue
        resolved_path = path.resolve()
        if excluded_dir is not None and _is_relative_to(resolved_path, excluded_dir):
            continue
        manifest.append(
            {
                "path": path.relative_to(base_path).as_posix(),
                "sha256": _sha256(path),
                "size": path.stat().st_size,
            }
        )
    return manifest


def _validate_conductor(conductor):
    label = f"conductor {conductor.identifier!r}"
    cond_time = _finite_time_history(
        getattr(conductor, "cond_time", None), f"{label}.cond_time"
    )

    cond_num_step = getattr(conductor, "cond_num_step", None)
    if cond_num_step != len(cond_time) - 1:
        raise CheckpointValidationError(
            f"{label}: cond_num_step ({cond_num_step}) must equal "
            f"len(cond_time) - 1 ({len(cond_time) - 1})."
        )

    time_step = getattr(conductor, "time_step", None)
    if time_step is None or not np.isfinite(time_step) or time_step < 0.0:
        raise CheckpointValidationError(
            f"{label}: time_step must be a finite non-negative value."
        )

    if getattr(conductor, "force_next_tstep_flag", None) is not False:
        raise CheckpointValidationError(
            f"{label}: force_next_tstep_flag must be False after get_time_step()."
        )

    events_time = np.asarray(getattr(conductor, "events_time", ()), dtype=float)
    if not np.isfinite(events_time).all():
        raise CheckpointValidationError(
            f"{label}: events_time contains non-finite values."
        )
    i_event = getattr(conductor, "i_event", None)
    if events_time.size:
        if i_event is None or not 0 <= int(i_event) < events_time.size:
            raise CheckpointValidationError(
                f"{label}: i_event is outside the reconstructed event timeline."
            )
    elif i_event not in (None, 0):
        raise CheckpointValidationError(
            f"{label}: i_event must be 0 or None when there are no events."
        )

    for attribute in ("EQTEIG", "dict_Step"):
        if not hasattr(conductor, attribute):
            raise CheckpointValidationError(f"{label}: missing {attribute}.")
    for key in ("SYSVAR", "SYSLOD"):
        if key not in conductor.dict_Step:
            raise CheckpointValidationError(
                f"{label}: dict_Step is missing required key {key!r}."
            )

    electric_active = getattr(conductor, "inputs", {}).get("I0_OP_MODE") is not None
    if electric_active:
        for attribute in (
            "electric_solution",
            "electric_solution_steady",
            "electric_time",
            "cond_el_num_step",
        ):
            if not hasattr(conductor, attribute):
                raise CheckpointValidationError(
                    f"{label}: active electric model is missing {attribute}."
                )
        if cond_num_step > 0 and not np.isclose(
            conductor.electric_time,
            cond_time[-1],
            rtol=1.0e-12,
            atol=1.0e-12,
        ):
            raise CheckpointValidationError(
                f"{label}: electric_time must coincide with cond_time[-1]."
            )

    for attribute in _BALANCE_ATTRIBUTES:
        if not hasattr(conductor, attribute):
            raise CheckpointValidationError(f"{label}: missing {attribute}.")

    for attribute in ("i_save", "num_step_save"):
        if not hasattr(conductor, attribute):
            raise CheckpointValidationError(f"{label}: missing {attribute}.")

    solid_components = _collection(conductor, "SolidComponent")
    for component in solid_components:
        component_label = getattr(component, "identifier", repr(component))
        node_state = getattr(component, "dict_node_pt", {})
        for key in ("EEXT", "EJHT"):
            if key not in node_state:
                raise CheckpointValidationError(
                    f"{label}, solid component {component_label!r}: missing {key}."
                )


def _write_metadata(h5file, simulation, trigger, manifest):
    metadata = h5file.create_group("metadata")
    metadata.attrs["schema_version"] = SCHEMA_VERSION
    metadata.attrs["created_utc"] = datetime.now(timezone.utc).isoformat()
    metadata.attrs["trigger"] = trigger
    metadata.attrs["git_commit"] = _git_commit()
    metadata.attrs["complete"] = False

    input_manifest = metadata.create_group("input_manifest")
    string_type = h5py.string_dtype(encoding="utf-8")
    input_manifest.create_dataset(
        "paths",
        data=np.asarray([entry["path"] for entry in manifest], dtype=object),
        dtype=string_type,
    )
    input_manifest.create_dataset(
        "sha256",
        data=np.asarray([entry["sha256"] for entry in manifest], dtype=object),
        dtype=string_type,
    )
    input_manifest.create_dataset(
        "sizes", data=np.asarray([entry["size"] for entry in manifest], dtype=np.int64)
    )


def _write_simulation_state(h5file, simulation):
    simulation_group = h5file.create_group("simulation")
    simulation_group.create_dataset(
        "simulation_time", data=np.asarray(simulation.simulation_time, dtype=float)
    )
    simulation_group.create_dataset("num_step", data=int(simulation.num_step))

    conductors_group = h5file.create_group("conductors")
    for conductor in simulation.list_of_Conductors:
        conductor_group = conductors_group.create_group(
            _safe_name(conductor.identifier)
        )
        conductor_group.attrs["identifier"] = conductor.identifier
        _write_conductor(conductor_group, conductor)


def _write_conductor(group, conductor):
    metadata = group.create_group("metadata")
    metadata.attrs["th_method"] = getattr(conductor, "inputs", {}).get(
        "METHOD", "unknown"
    )
    metadata.attrs["electric_method"] = getattr(conductor, "inputs", {}).get(
        "ELECTRIC_METHOD", "unknown"
    )

    clock = group.create_group("clock")
    for attribute in (
        "cond_time",
        "cond_num_step",
        "time_step",
        "EQTEIG",
        "i_event",
        "events_time",
        "cond_el_num_step",
        "electric_time",
    ):
        if hasattr(conductor, attribute):
            _write_value(clock, attribute, getattr(conductor, attribute))

    th_history = group.create_group("th_history")
    for key, value in conductor.dict_Step.items():
        _write_value(th_history, key, value)

    electric = group.create_group("electric")
    for attribute in ("electric_solution", "electric_solution_steady"):
        if hasattr(conductor, attribute):
            _write_value(electric, attribute, getattr(conductor, attribute))

    balances = group.create_group("balances")
    for attribute in _BALANCE_ATTRIBUTES:
        _write_value(balances, attribute, getattr(conductor, attribute))

    components = group.create_group("components")
    for component in _collection(conductor, "SolidComponent"):
        identifier = getattr(component, "identifier", component.__class__.__name__)
        component_group = components.create_group(_safe_name(identifier))
        component_group.attrs["identifier"] = identifier
        energy_history = component_group.create_group("energy_history")
        _write_value(energy_history, "EEXT", component.dict_node_pt["EEXT"])
        _write_value(energy_history, "EJHT", component.dict_node_pt["EJHT"])

    output_state = group.create_group("output_state")
    for attribute in _OUTPUT_STATE_ATTRIBUTES:
        if hasattr(conductor, attribute):
            _write_value(output_state, attribute, getattr(conductor, attribute))
    output_buffers = output_state.create_group("buffers")
    for owner_path, owner in _iter_output_owners(conductor):
        owner_group = _require_path(output_buffers, owner_path)
        for attribute in _OUTPUT_BUFFER_ATTRIBUTES:
            if hasattr(owner, attribute):
                _write_value(owner_group, attribute, getattr(owner, attribute))


def _iter_output_owners(conductor):
    yield "conductor", conductor
    for component in _collection(conductor, "all_component"):
        identifier = getattr(component, "identifier", component.__class__.__name__)
        base_path = f"components/{_safe_name(identifier)}"
        yield f"{base_path}/object", component
        if hasattr(component, "coolant"):
            yield f"{base_path}/coolant", component.coolant
        if hasattr(component, "channel"):
            yield f"{base_path}/channel", component.channel


def _write_value(parent, name, value):
    safe_name = _safe_name(name)
    if isinstance(value, Mapping):
        child = parent.create_group(safe_name)
        child.attrs["python_type"] = "mapping"
        child.attrs["original_name"] = str(name)
        for key, item in value.items():
            _write_value(child, key, item)
        return

    if value is None:
        child = parent.create_group(safe_name)
        child.attrs["python_type"] = "none"
        child.attrs["original_name"] = str(name)
        return

    if isinstance(value, (list, tuple)):
        array = np.asarray(value)
        if array.dtype == object:
            child = parent.create_group(safe_name)
            child.attrs["python_type"] = "tuple" if isinstance(value, tuple) else "list"
            child.attrs["original_name"] = str(name)
            for index, item in enumerate(value):
                _write_value(child, str(index), item)
            return
        value = array

    if isinstance(value, str):
        dataset = parent.create_dataset(
            safe_name, data=value, dtype=h5py.string_dtype(encoding="utf-8")
        )
    else:
        array = np.asarray(value)
        if array.dtype == object:
            raise TypeError(
                f"Cannot serialize object-valued checkpoint entry {name!r}."
            )
        if array.dtype.kind in ("U", "S"):
            dataset = parent.create_dataset(
                safe_name,
                data=array.astype(object),
                dtype=h5py.string_dtype(encoding="utf-8"),
            )
        else:
            dataset = parent.create_dataset(safe_name, data=value)
    dataset.attrs["original_name"] = str(name)


def _collection(conductor, key):
    inventory = getattr(conductor, "inventory", {})
    item = inventory.get(key) if isinstance(inventory, Mapping) else None
    return list(getattr(item, "collection", ()))


def _require_path(parent, path):
    group = parent
    for part in path.split("/"):
        group = group.require_group(part)
    return group


def _finite_time_history(values, label):
    if values is None:
        raise CheckpointValidationError(f"Missing {label}.")
    array = np.asarray(values, dtype=float)
    if array.ndim != 1 or not array.size:
        raise CheckpointValidationError(f"{label} must be a non-empty 1-D history.")
    if not np.isfinite(array).all():
        raise CheckpointValidationError(f"{label} contains non-finite values.")
    if np.any(np.diff(array) < 0.0):
        raise CheckpointValidationError(f"{label} must be non-decreasing.")
    return array


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _safe_name(value):
    return quote(str(value), safe="_-.")


def _is_relative_to(path, parent):
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True


def _git_commit():
    source_root = Path(__file__).resolve().parents[1]
    result = subprocess.run(
        ("git", "rev-parse", "HEAD"),
        cwd=source_root,
        capture_output=True,
        check=False,
        text=True,
    )
    return result.stdout.strip() if result.returncode == 0 else "unknown"
