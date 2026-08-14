"""HDF5 checkpoint persistence and recovery utilities.

The reader returns a detached intermediate representation.  Runtime validation
and state application are separate operations so every incompatibility can be
reported before a live OPENSC2 object is mutated.
"""

import copy
from collections.abc import Mapping
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import os
from pathlib import Path
import re
import subprocess
from urllib.parse import quote

import h5py
import numpy as np


SCHEMA_VERSION = "1.1"
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


class CheckpointReadError(ValueError):
    """Raised when a checkpoint cannot be safely read or validated."""


@dataclass(frozen=True)
class InputManifestEntry:
    """One input file recorded in a checkpoint manifest."""

    path: str
    sha256: str
    size: int


@dataclass(frozen=True)
class InputManifestChange:
    """One input file whose saved and current contents differ."""

    checkpoint: InputManifestEntry
    current: InputManifestEntry


@dataclass(frozen=True)
class InputManifestComparison:
    """Detached report comparing saved and current simulation inputs."""

    checkpoint_entries: tuple[InputManifestEntry, ...]
    current_entries: tuple[InputManifestEntry, ...]
    missing: tuple[InputManifestEntry, ...]
    added: tuple[InputManifestEntry, ...]
    modified: tuple[InputManifestChange, ...]

    @property
    def is_match(self):
        """Return ``True`` only when both manifests are identical."""

        return not (self.missing or self.added or self.modified)


@dataclass(frozen=True)
class RestartCompatibilityReport:
    """Non-destructive decision about whether a restart may proceed."""

    mode: str
    is_compatible: bool
    manifest_comparison: InputManifestComparison
    blocking_reasons: tuple[str, ...]
    warnings: tuple[str, ...]


@dataclass(frozen=True)
class RuntimeRestoreValidationReport:
    """Non-destructive validation of a runtime restore destination."""

    is_valid: bool
    blocking_reasons: tuple[str, ...]
    warnings: tuple[str, ...]


@dataclass(frozen=True)
class ConductorCheckpointData:
    """Detached persistent state of one conductor."""

    identifier: str
    th_method: str
    electric_method: str
    clock: dict
    th_history: dict
    electric: dict
    balances: dict
    components: dict
    output_state: dict


@dataclass(frozen=True)
class CheckpointData:
    """Validated checkpoint contents, independent of an open HDF5 file."""

    path: Path
    schema_version: str
    created_utc: str
    trigger: str
    git_commit: str
    input_manifest: tuple[InputManifestEntry, ...]
    simulation_time: np.ndarray
    num_step: int
    conductors: dict[str, ConductorCheckpointData]


def read_checkpoint(checkpoint_path):
    """Read and validate a checkpoint without mutating a simulation.

    All HDF5 datasets are copied into ordinary Python values or NumPy arrays
    before the file is closed.  The returned object therefore owns no live
    HDF5 handles.
    """

    checkpoint_path = Path(checkpoint_path)
    if not checkpoint_path.is_file():
        raise CheckpointReadError(
            f"Checkpoint file does not exist: {checkpoint_path!s}."
        )

    try:
        with h5py.File(checkpoint_path, "r") as h5file:
            return _read_checkpoint_file(h5file, checkpoint_path.resolve())
    except CheckpointReadError:
        raise
    except (KeyError, OSError, TypeError, ValueError) as exc:
        raise CheckpointReadError(
            f"Cannot read checkpoint {checkpoint_path!s}: {exc}"
        ) from exc


def compare_input_manifest(checkpoint, input_directory):
    """Compare a checkpoint manifest with the current input directory.

    The function is deliberately diagnostic: it reports every difference but
    does not decide whether a recovery or continuation may proceed.  It also
    does not mutate the detached checkpoint representation.
    """

    if not isinstance(checkpoint, CheckpointData):
        raise TypeError("checkpoint must be a CheckpointData instance.")

    checkpoint_entries = tuple(checkpoint.input_manifest)
    current_entries = _build_input_manifest_entries(input_directory)
    checkpoint_by_path = {entry.path: entry for entry in checkpoint_entries}
    current_by_path = {entry.path: entry for entry in current_entries}

    missing = tuple(
        checkpoint_by_path[path]
        for path in sorted(checkpoint_by_path.keys() - current_by_path.keys())
    )
    added = tuple(
        current_by_path[path]
        for path in sorted(current_by_path.keys() - checkpoint_by_path.keys())
    )
    modified = tuple(
        InputManifestChange(
            checkpoint=checkpoint_by_path[path],
            current=current_by_path[path],
        )
        for path in sorted(checkpoint_by_path.keys() & current_by_path.keys())
        if (
            checkpoint_by_path[path].sha256 != current_by_path[path].sha256
            or checkpoint_by_path[path].size != current_by_path[path].size
        )
    )

    return InputManifestComparison(
        checkpoint_entries=checkpoint_entries,
        current_entries=current_entries,
        missing=missing,
        added=added,
        modified=modified,
    )


def evaluate_restart_compatibility(
    checkpoint,
    input_directory,
    mode="recovery",
):
    """Evaluate restart compatibility without mutating runtime state.

    ``recovery`` is accepted only when the current input manifest is exactly
    equal to the one persisted in the checkpoint.  ``continuation`` is kept as
    an explicit interface value, but is blocked until semantic input
    compatibility has been implemented.  Unknown modes are programming or
    user-interface errors and are rejected immediately.
    """

    valid_modes = ("recovery", "continuation")
    if mode not in valid_modes:
        raise ValueError(
            f"Invalid restart mode {mode!r}; expected one of {valid_modes!r}."
        )

    comparison = compare_input_manifest(checkpoint, input_directory)
    blocking_reasons = []

    if mode == "continuation":
        blocking_reasons.append(
            "Restart mode 'continuation' is not supported yet."
        )
    else:
        blocking_reasons.extend(
            f"Input file is missing: {entry.path}."
            for entry in comparison.missing
        )
        blocking_reasons.extend(
            f"Unexpected input file was added: {entry.path}."
            for entry in comparison.added
        )
        blocking_reasons.extend(
            f"Input file was modified: {change.checkpoint.path}."
            for change in comparison.modified
        )

    return RestartCompatibilityReport(
        mode=mode,
        is_compatible=not blocking_reasons,
        manifest_comparison=comparison,
        blocking_reasons=tuple(blocking_reasons),
        warnings=(),
    )


def validate_runtime_restore_target(checkpoint, simulation):
    """Validate a freshly initialized runtime before applying a checkpoint.

    The function checks identities, methods, component kinds, numerical-state
    shapes, event timelines, and output ownership.  It deliberately performs
    no assignments: every incompatibility is accumulated before the immutable
    report is returned.
    """

    if not isinstance(checkpoint, CheckpointData):
        raise TypeError("checkpoint must be a CheckpointData instance.")

    reasons = []
    _validate_fresh_runtime_clock(simulation, reasons)
    runtime_conductors = list(
        getattr(simulation, "list_of_Conductors", ())
    )
    runtime_by_identifier = {}
    for conductor in runtime_conductors:
        identifier = getattr(conductor, "identifier", None)
        if not identifier:
            reasons.append("Runtime conductor is missing an identifier.")
            continue
        if identifier in runtime_by_identifier:
            reasons.append(
                f"Runtime conductor identifier is duplicated: {identifier!r}."
            )
            continue
        runtime_by_identifier[identifier] = conductor

    checkpoint_ids = set(checkpoint.conductors)
    runtime_ids = set(runtime_by_identifier)
    reasons.extend(
        f"Runtime is missing conductor {identifier!r}."
        for identifier in sorted(checkpoint_ids - runtime_ids)
    )
    reasons.extend(
        f"Runtime contains unexpected conductor {identifier!r}."
        for identifier in sorted(runtime_ids - checkpoint_ids)
    )

    _validate_checkpoint_global_clock(checkpoint, reasons)
    for identifier in sorted(checkpoint_ids & runtime_ids):
        _validate_runtime_conductor_target(
            checkpoint.conductors[identifier],
            runtime_by_identifier[identifier],
            reasons,
        )

    return RuntimeRestoreValidationReport(
        is_valid=not reasons,
        blocking_reasons=tuple(reasons),
        warnings=(),
    )


def apply_checkpoint_to_runtime(checkpoint, simulation):
    """Apply validated checkpoint state to a freshly initialized runtime.

    This function restores state persisted in schema 1.1 and synchronizes the
    fluid and solid primary nodal variables from ``SYSVAR``.  It does not
    rebuild derived properties, touch output files, or enter the transient
    loop.

    All replacement values are prepared before the first assignment.  Runtime
    container types are preserved where they are operationally significant;
    in particular, time histories remain lists so the solver can append the
    next completed time.
    """

    compatibility = evaluate_restart_compatibility(
        checkpoint,
        getattr(simulation, "basePath", None),
        mode="recovery",
    )
    if not compatibility.is_compatible:
        details = "\n".join(
            f"- {reason}" for reason in compatibility.blocking_reasons
        )
        raise CheckpointValidationError(
            "Checkpoint inputs are not compatible with recovery:\n" + details
        )

    report = validate_runtime_restore_target(checkpoint, simulation)
    if not report.is_valid:
        details = "\n".join(f"- {reason}" for reason in report.blocking_reasons)
        raise CheckpointValidationError(
            "Runtime restore target is not valid:\n" + details
        )

    runtime_by_identifier = {
        conductor.identifier: conductor
        for conductor in simulation.list_of_Conductors
    }

    prepared_simulation_time = _copy_for_runtime(
        checkpoint.simulation_time,
        simulation.simulation_time,
    )
    prepared_conductors = []
    for identifier in sorted(checkpoint.conductors):
        saved = checkpoint.conductors[identifier]
        runtime = runtime_by_identifier[identifier]
        runtime_components = {
            component.identifier: component
            for component, _kind in _component_inventory(
                runtime, f"conductor {identifier!r}"
            )
        }

        clock = {
            attribute: _copy_for_runtime(
                saved.clock[attribute], getattr(runtime, attribute)
            )
            for attribute in (
                "cond_time",
                "cond_num_step",
                "time_step",
                "EQTEIG",
                "i_event",
            )
        }
        for attribute in ("electric_time", "cond_el_num_step"):
            if attribute in saved.clock:
                clock[attribute] = _copy_for_runtime(
                    saved.clock[attribute], getattr(runtime, attribute)
                )

        th_history = _copy_for_runtime(saved.th_history, runtime.dict_Step)
        component_primary_state = _prepare_component_primary_state(
            runtime, th_history
        )
        electric = {
            attribute: _copy_for_runtime(value, getattr(runtime, attribute))
            for attribute, value in saved.electric.items()
        }
        balances = {
            attribute: _copy_for_runtime(value, getattr(runtime, attribute))
            for attribute, value in saved.balances.items()
        }

        component_energy = {}
        for component_identifier, component_state in saved.components.items():
            if component_state["kind"] != "solid":
                continue
            runtime_component = runtime_components[component_identifier]
            component_energy[component_identifier] = {
                key: _copy_for_runtime(
                    value, runtime_component.dict_node_pt[key]
                )
                for key, value in component_state["energy_history"].items()
            }

        output_attributes = {
            attribute: _copy_for_runtime(
                value, getattr(runtime, attribute, None)
            )
            for attribute, value in saved.output_state.items()
            if attribute != "buffers"
        }
        output_buffers = []
        for owner_path, owner in _runtime_output_owners(runtime):
            saved_owner = _mapping_path(
                saved.output_state["buffers"], owner_path
            )
            for attribute, value in saved_owner.items():
                output_buffers.append(
                    (
                        owner,
                        attribute,
                        _copy_for_runtime(value, getattr(owner, attribute)),
                    )
                )

        prepared_conductors.append(
            (
                runtime,
                clock,
                th_history,
                electric,
                balances,
                runtime_components,
                component_primary_state,
                component_energy,
                output_attributes,
                output_buffers,
            )
        )

    simulation.simulation_time = prepared_simulation_time
    simulation.num_step = int(checkpoint.num_step)
    for (
        runtime,
        clock,
        th_history,
        electric,
        balances,
        runtime_components,
        component_primary_state,
        component_energy,
        output_attributes,
        output_buffers,
    ) in prepared_conductors:
        for attribute, value in clock.items():
            setattr(runtime, attribute, value)
        runtime.dict_Step = th_history
        for owner, values in component_primary_state:
            for key, value in values.items():
                owner.dict_node_pt[key] = value
        for attribute, value in electric.items():
            setattr(runtime, attribute, value)
        for attribute, value in balances.items():
            setattr(runtime, attribute, value)
        for component_identifier, energy_history in component_energy.items():
            node_state = runtime_components[component_identifier].dict_node_pt
            for key, value in energy_history.items():
                node_state[key] = value
        for attribute, value in output_attributes.items():
            setattr(runtime, attribute, value)
        for owner, attribute, value in output_buffers:
            setattr(owner, attribute, value)

    # This flag is deliberately the last mutation. Invalid checkpoints and
    # failures while preparing the replacement state must leave a fresh
    # runtime on the normal t=0 path.
    simulation.restored_from_checkpoint = True


def _prepare_component_primary_state(conductor, th_history):
    """Prepare detached primary nodal state reconstructed from ``SYSVAR``."""

    sysvar = np.asarray(th_history["SYSVAR"])
    if sysvar.ndim != 2 or sysvar.shape[1] < 1:
        raise CheckpointValidationError(
            f"Conductor {conductor.identifier!r} SYSVAR must be a two-"
            "dimensional array with at least one history column."
        )

    try:
        ndf = int(conductor.dict_N_equation["NODOFS"])
        equation_index = conductor.equation_index
    except (AttributeError, KeyError, TypeError, ValueError) as error:
        raise CheckpointValidationError(
            f"Conductor {conductor.identifier!r} lacks a valid runtime "
            "equation mapping."
        ) from error

    if ndf <= 0:
        raise CheckpointValidationError(
            f"Conductor {conductor.identifier!r} NODOFS must be positive."
        )

    prepared = []
    for component in conductor.inventory["FluidComponent"].collection:
        try:
            indices = equation_index[component.identifier]
            values = {
                name: sysvar[getattr(indices, name)::ndf, 0].copy()
                for name in ("velocity", "pressure", "temperature")
            }
            owner = component.coolant
        except (AttributeError, KeyError, TypeError) as error:
            raise CheckpointValidationError(
                f"Fluid component {component.identifier!r} lacks a valid "
                "runtime equation mapping or nodal state."
            ) from error
        _validate_prepared_primary_shapes(component.identifier, owner, values)
        prepared.append((owner, values))

    for component in conductor.inventory["SolidComponent"].collection:
        try:
            values = {
                "temperature": sysvar[
                    equation_index[component.identifier]::ndf, 0
                ].copy()
            }
        except (KeyError, TypeError) as error:
            raise CheckpointValidationError(
                f"Solid component {component.identifier!r} lacks a valid "
                "runtime equation mapping or nodal state."
            ) from error
        _validate_prepared_primary_shapes(
            component.identifier, component, values
        )
        prepared.append((component, values))

    return prepared


def _validate_prepared_primary_shapes(identifier, owner, values):
    """Reject invalid SYSVAR slices before any runtime state is mutated."""

    node_state = getattr(owner, "dict_node_pt", None)
    if not isinstance(node_state, Mapping):
        raise CheckpointValidationError(
            f"Component {identifier!r} lacks a nodal state mapping."
        )

    for key, value in values.items():
        if key not in node_state:
            raise CheckpointValidationError(
                f"Component {identifier!r} nodal state is missing {key!r}."
            )
        expected_shape = np.asarray(node_state[key]).shape
        if value.shape != expected_shape:
            raise CheckpointValidationError(
                f"Component {identifier!r} {key!r} reconstructed from "
                f"SYSVAR has shape {value.shape}, expected {expected_shape}."
            )


def _copy_for_runtime(saved, runtime_template):
    """Deep-copy saved data while retaining mutable runtime container types."""

    if isinstance(runtime_template, Mapping):
        # Legacy schema-1.1 HDF5 groups may expose mapping keys in
        # lexicographic order.  Rebuild shared entries in the order of the
        # already initialized runtime mapping, then append saved-only keys
        # that are legitimately created lazily during the transient.
        restored = {
            key: _copy_for_runtime(saved[key], value)
            for key, value in runtime_template.items()
            if key in saved
        }
        restored.update(
            {
                key: copy.deepcopy(value)
                for key, value in saved.items()
                if key not in runtime_template
            }
        )
        return restored
    if isinstance(runtime_template, list):
        values = saved.tolist() if isinstance(saved, np.ndarray) else list(saved)
        return copy.deepcopy(values)
    if isinstance(runtime_template, tuple):
        values = saved.tolist() if isinstance(saved, np.ndarray) else list(saved)
        return tuple(copy.deepcopy(values))
    if isinstance(runtime_template, np.ndarray):
        return np.array(saved, dtype=runtime_template.dtype, copy=True)
    if isinstance(runtime_template, np.generic):
        return np.asarray(saved, dtype=runtime_template.dtype).item()
    return copy.deepcopy(saved)


def _runtime_output_owners(conductor):
    """Yield runtime output owners using decoded mapping path components."""

    yield ("conductor",), conductor
    for component, _kind in _component_inventory(
        conductor, f"conductor {conductor.identifier!r}"
    ):
        base_path = ("components", component.identifier)
        yield base_path + ("object",), component
        if hasattr(component, "coolant"):
            yield base_path + ("coolant",), component.coolant
        if hasattr(component, "channel"):
            yield base_path + ("channel",), component.channel


def _mapping_path(mapping, path):
    value = mapping
    for key in path:
        value = value[key]
    return value


def _validate_fresh_runtime_clock(simulation, reasons):
    simulation_time = np.asarray(
        getattr(simulation, "simulation_time", ()), dtype=float
    )
    if (
        simulation_time.shape != (1,)
        or not np.isclose(simulation_time[0], 0.0)
    ):
        reasons.append(
            "Runtime simulation must be freshly initialized at time zero."
        )
    if getattr(simulation, "num_step", None) != 0:
        reasons.append("Runtime simulation num_step must be zero before restore.")


def _validate_checkpoint_global_clock(checkpoint, reasons):
    simulation_time = np.asarray(checkpoint.simulation_time)
    if checkpoint.num_step != simulation_time.size - 1:
        reasons.append(
            "Checkpoint global num_step does not equal "
            "len(simulation_time) - 1."
        )

    conductor_times = []
    for conductor in checkpoint.conductors.values():
        cond_time = np.asarray(conductor.clock.get("cond_time", ()))
        if cond_time.ndim == 1 and cond_time.size:
            conductor_times.append(float(cond_time[-1]))
        if conductor.clock.get("cond_num_step") != checkpoint.num_step:
            reasons.append(
                f"Conductor {conductor.identifier!r} checkpoint step does not "
                "match the global checkpoint step."
            )

    if conductor_times and not np.isclose(
        float(simulation_time[-1]),
        max(conductor_times),
        rtol=1.0e-12,
        atol=1.0e-12,
    ):
        reasons.append(
            "Checkpoint global time does not equal the latest conductor time."
        )


def _validate_runtime_conductor_target(saved, runtime, reasons):
    label = f"Conductor {saved.identifier!r}"
    runtime_cond_time = np.asarray(
        getattr(runtime, "cond_time", ()), dtype=float
    )
    if (
        runtime_cond_time.shape != (1,)
        or not np.isclose(runtime_cond_time[0], 0.0)
    ):
        reasons.append(f"{label} runtime clock must start at time zero.")
    if getattr(runtime, "cond_num_step", None) != 0:
        reasons.append(f"{label} runtime cond_num_step must be zero before restore.")
    if getattr(runtime, "i_save", None) != 0:
        reasons.append(f"{label} runtime i_save must be zero before restore.")

    runtime_inputs = getattr(runtime, "inputs", {})
    runtime_th_method = (
        runtime_inputs.get("METHOD")
        if isinstance(runtime_inputs, Mapping)
        else None
    )
    runtime_electric_method = (
        runtime_inputs.get("ELECTRIC_METHOD")
        if isinstance(runtime_inputs, Mapping)
        else None
    )
    if runtime_th_method != saved.th_method:
        reasons.append(
            f"{label} TH method differs: checkpoint={saved.th_method!r}, "
            f"runtime={runtime_th_method!r}."
        )
    if runtime_electric_method != saved.electric_method:
        reasons.append(
            f"{label} electric method differs: "
            f"checkpoint={saved.electric_method!r}, "
            f"runtime={runtime_electric_method!r}."
        )

    try:
        runtime_inventory = _component_inventory(runtime, label)
    except CheckpointValidationError as exc:
        reasons.append(str(exc))
        runtime_inventory = []
    runtime_components = {
        component.identifier: (component, kind)
        for component, kind in runtime_inventory
    }
    saved_ids = set(saved.components)
    runtime_ids = set(runtime_components)
    reasons.extend(
        f"{label} is missing component {identifier!r}."
        for identifier in sorted(saved_ids - runtime_ids)
    )
    reasons.extend(
        f"{label} contains unexpected component {identifier!r}."
        for identifier in sorted(runtime_ids - saved_ids)
    )
    for identifier in sorted(saved_ids & runtime_ids):
        saved_kind = saved.components[identifier]["kind"]
        runtime_component, runtime_kind = runtime_components[identifier]
        if saved_kind != runtime_kind:
            reasons.append(
                f"{label} component {identifier!r} kind differs: "
                f"checkpoint={saved_kind!r}, runtime={runtime_kind!r}."
            )
        elif saved_kind == "solid":
            runtime_node_state = getattr(runtime_component, "dict_node_pt", {})
            for key, saved_value in saved.components[identifier][
                "energy_history"
            ].items():
                if key not in runtime_node_state:
                    reasons.append(
                        f"{label} solid component {identifier!r} is missing "
                        f"runtime energy history {key!r}."
                    )
                else:
                    _validate_matching_shape(
                        saved_value,
                        runtime_node_state[key],
                        f"{label} solid component {identifier!r} {key}",
                        reasons,
                    )

    _validate_matching_mapping_shapes(
        saved.th_history,
        getattr(runtime, "dict_Step", None),
        f"{label} TH history",
        reasons,
    )
    _validate_matching_shape(
        saved.clock.get("EQTEIG"),
        getattr(runtime, "EQTEIG", None),
        f"{label} EQTEIG",
        reasons,
    )
    _validate_event_timeline(saved, runtime, label, reasons)
    _validate_electric_target(saved, runtime, label, reasons)
    for attribute in _BALANCE_ATTRIBUTES:
        if not hasattr(runtime, attribute):
            reasons.append(f"{label} is missing balance {attribute!r}.")
    _validate_output_target(
        saved,
        runtime,
        runtime_components,
        label,
        reasons,
    )


def _validate_matching_mapping_shapes(saved, runtime, label, reasons):
    if not isinstance(runtime, Mapping):
        reasons.append(f"{label} is missing from the runtime.")
        return
    saved_keys = set(saved)
    runtime_keys = set(runtime)
    reasons.extend(
        f"{label} is missing key {key!r}."
        for key in sorted(saved_keys - runtime_keys)
    )
    reasons.extend(
        f"{label} contains unexpected key {key!r}."
        for key in sorted(runtime_keys - saved_keys)
    )
    for key in sorted(saved_keys & runtime_keys):
        _validate_matching_shape(
            saved[key], runtime[key], f"{label}[{key!r}]", reasons
        )


def _validate_matching_shape(saved, runtime, label, reasons):
    if runtime is None:
        reasons.append(f"{label} is missing from the runtime.")
        return
    saved_shape = np.asarray(saved).shape
    runtime_shape = np.asarray(runtime).shape
    if saved_shape != runtime_shape:
        reasons.append(
            f"{label} shape differs: checkpoint={saved_shape!r}, "
            f"runtime={runtime_shape!r}."
        )


def _validate_event_timeline(saved, runtime, label, reasons):
    saved_events = np.asarray(saved.clock.get("events_time"))
    runtime_events = getattr(runtime, "events_time", None)
    if runtime_events is None:
        reasons.append(f"{label} event timeline is missing from the runtime.")
        return
    runtime_events = np.asarray(runtime_events)
    if (
        saved_events.shape != runtime_events.shape
        or not np.allclose(
            saved_events,
            runtime_events,
            rtol=1.0e-12,
            atol=1.0e-12,
        )
    ):
        reasons.append(f"{label} event timeline differs from the checkpoint.")
    saved_i_event = saved.clock.get("i_event")
    if not isinstance(saved_i_event, (int, np.integer)) or not (
        0 <= int(saved_i_event) < saved_events.size
    ):
        reasons.append(f"{label} saved i_event is outside the event timeline.")


def _validate_electric_target(saved, runtime, label, reasons):
    runtime_inputs = getattr(runtime, "inputs", {})
    runtime_active = (
        isinstance(runtime_inputs, Mapping)
        and runtime_inputs.get("I0_OP_MODE") is not None
    )
    saved_active = bool(saved.electric)
    if saved_active != runtime_active:
        reasons.append(
            f"{label} electric activation differs: "
            f"checkpoint={saved_active!r}, runtime={runtime_active!r}."
        )
        return
    if not saved_active:
        return
    for key in ("electric_solution", "electric_solution_steady"):
        _validate_matching_shape(
            saved.electric.get(key),
            getattr(runtime, key, None),
            f"{label} {key}",
            reasons,
        )


def _validate_output_target(saved, runtime, runtime_components, label, reasons):
    saved_output = saved.output_state
    runtime_num_step_save = getattr(runtime, "num_step_save", None)
    _validate_matching_shape(
        saved_output.get("num_step_save"),
        runtime_num_step_save,
        f"{label} num_step_save",
        reasons,
    )

    saved_i_save = saved_output.get("i_save")
    runtime_space_save = getattr(runtime, "Space_save", None)
    if runtime_space_save is None:
        reasons.append(f"{label} Space_save is missing from the runtime.")
    elif not isinstance(saved_i_save, (int, np.integer)) or not (
        0 <= int(saved_i_save) < np.asarray(runtime_space_save).size
    ):
        reasons.append(f"{label} saved i_save is outside runtime Space_save.")

    saved_buffers = saved_output.get("buffers", {})
    saved_conductor_buffer = saved_buffers.get("conductor")
    if not isinstance(saved_conductor_buffer, Mapping):
        reasons.append(f"{label} conductor output buffer is missing.")
    else:
        _validate_output_owner(
            saved_conductor_buffer, runtime, f"{label} output owner conductor", reasons
        )

    saved_component_buffers = saved_buffers.get("components", {})
    if not isinstance(saved_component_buffers, Mapping):
        reasons.append(f"{label} component output buffers are missing.")
        return
    for identifier in sorted(set(saved_component_buffers) & set(runtime_components)):
        saved_owners = saved_component_buffers[identifier]
        component = runtime_components[identifier][0]
        runtime_owners = {"object": component}
        if hasattr(component, "coolant"):
            runtime_owners["coolant"] = component.coolant
        if hasattr(component, "channel"):
            runtime_owners["channel"] = component.channel

        if not isinstance(saved_owners, Mapping):
            reasons.append(
                f"{label} component {identifier!r} output owners are invalid."
            )
            continue
        saved_owner_names = set(saved_owners)
        runtime_owner_names = set(runtime_owners)
        reasons.extend(
            f"{label} component {identifier!r} is missing output owner {name!r}."
            for name in sorted(saved_owner_names - runtime_owner_names)
        )
        reasons.extend(
            f"{label} component {identifier!r} has unexpected output owner {name!r}."
            for name in sorted(runtime_owner_names - saved_owner_names)
        )
        for owner_name in sorted(saved_owner_names & runtime_owner_names):
            _validate_output_owner(
                saved_owners[owner_name],
                runtime_owners[owner_name],
                f"{label} component {identifier!r} output owner {owner_name!r}",
                reasons,
            )


def _validate_output_owner(saved, runtime, label, reasons):
    if not isinstance(saved, Mapping):
        reasons.append(f"{label} buffer is invalid.")
        return
    for attribute in sorted(saved):
        if attribute not in _OUTPUT_BUFFER_ATTRIBUTES:
            reasons.append(f"{label} has unsupported buffer {attribute!r}.")
        elif not hasattr(runtime, attribute):
            reasons.append(f"{label} is missing buffer {attribute!r}.")
        else:
            _validate_restore_container_structure(
                saved[attribute],
                getattr(runtime, attribute),
                f"{label} buffer {attribute!r}",
                reasons,
            )


def _validate_restore_container_structure(saved, runtime, label, reasons):
    """Validate restorable mappings without requiring eager runtime keys.

    Output dictionaries are populated lazily during the transient.  A fresh
    runtime may therefore lack keys that legitimately exist in a checkpoint;
    those keys will be created by ``_copy_for_runtime``.  Runtime-only keys
    remain blocking because they indicate that the current runtime expects
    state that the checkpoint cannot supply.
    """

    saved_is_mapping = isinstance(saved, Mapping)
    runtime_is_mapping = isinstance(runtime, Mapping)
    if saved_is_mapping != runtime_is_mapping:
        reasons.append(f"{label} container type differs from the runtime.")
        return
    if not saved_is_mapping:
        return

    saved_keys = set(saved)
    runtime_keys = set(runtime)
    reasons.extend(
        f"{label} has unexpected runtime key {key!r}."
        for key in sorted(runtime_keys - saved_keys)
    )
    for key in sorted(saved_keys & runtime_keys):
        _validate_restore_container_structure(
            saved[key], runtime[key], f"{label}[{key!r}]", reasons
        )


def _read_checkpoint_file(h5file, checkpoint_path):
    _require_members(
        h5file,
        ("metadata", "simulation", "conductors"),
        "checkpoint root",
    )

    metadata = h5file["metadata"]
    if not isinstance(metadata, h5py.Group):
        raise CheckpointReadError("Checkpoint entry 'metadata' must be a group.")
    _require_attributes(
        metadata,
        ("schema_version", "created_utc", "trigger", "git_commit", "complete"),
        "metadata",
    )

    schema_version = _text(metadata.attrs["schema_version"], "schema_version")
    if schema_version != SCHEMA_VERSION:
        raise CheckpointReadError(
            f"Unsupported checkpoint schema {schema_version!r}; "
            f"expected {SCHEMA_VERSION!r}."
        )
    if not _strict_bool(metadata.attrs["complete"], "metadata.complete"):
        raise CheckpointReadError("Checkpoint is incomplete.")

    trigger = _text(metadata.attrs["trigger"], "trigger")
    if trigger not in VALID_TRIGGERS:
        raise CheckpointReadError(
            f"Invalid checkpoint trigger {trigger!r}."
        )
    created_utc = _text(metadata.attrs["created_utc"], "created_utc")
    git_commit = _text(metadata.attrs["git_commit"], "git_commit")
    try:
        parsed_created_utc = datetime.fromisoformat(created_utc)
    except ValueError as exc:
        raise CheckpointReadError(
            "Checkpoint creation time is not valid ISO-8601 text."
        ) from exc
    if parsed_created_utc.tzinfo is None:
        raise CheckpointReadError(
            "Checkpoint creation time must include a UTC offset."
        )
    if not git_commit:
        raise CheckpointReadError("Checkpoint git commit metadata is empty.")
    manifest = _read_input_manifest(metadata)

    simulation = h5file["simulation"]
    if not isinstance(simulation, h5py.Group):
        raise CheckpointReadError("Checkpoint entry 'simulation' must be a group.")
    _require_members(simulation, ("simulation_time", "num_step"), "simulation")
    simulation_time = _read_finite_time_dataset(
        simulation["simulation_time"], "simulation/simulation_time"
    )
    num_step = _read_non_negative_integer(
        simulation["num_step"], "simulation/num_step"
    )

    conductors_group = h5file["conductors"]
    if not isinstance(conductors_group, h5py.Group):
        raise CheckpointReadError("Checkpoint entry 'conductors' must be a group.")
    if len(conductors_group) == 0:
        raise CheckpointReadError("Checkpoint contains no conductors.")

    conductors = {}
    for group_name in conductors_group:
        conductor = _read_conductor(conductors_group[group_name], group_name)
        if conductor.identifier in conductors:
            raise CheckpointReadError(
                f"Duplicate conductor identifier {conductor.identifier!r}."
            )
        conductors[conductor.identifier] = conductor

    return CheckpointData(
        path=checkpoint_path,
        schema_version=schema_version,
        created_utc=created_utc,
        trigger=trigger,
        git_commit=git_commit,
        input_manifest=manifest,
        simulation_time=simulation_time,
        num_step=num_step,
        conductors=conductors,
    )


def _read_input_manifest(metadata):
    _require_members(metadata, ("input_manifest",), "metadata")
    group = metadata["input_manifest"]
    if not isinstance(group, h5py.Group):
        raise CheckpointReadError("metadata/input_manifest must be a group.")
    _require_members(group, ("paths", "sha256", "sizes"), "input manifest")

    paths = _read_string_vector(group["paths"], "input manifest paths")
    hashes = _read_string_vector(group["sha256"], "input manifest hashes")
    sizes = np.asarray(group["sizes"][()])
    if sizes.ndim != 1 or not np.issubdtype(sizes.dtype, np.integer):
        raise CheckpointReadError("Input manifest sizes must be a 1-D integer array.")
    if not (len(paths) == len(hashes) == len(sizes)):
        raise CheckpointReadError("Input manifest columns have different lengths.")

    entries = []
    seen_paths = set()
    for path, sha256, size in zip(paths, hashes, sizes):
        candidate = Path(path)
        if (
            not path
            or candidate.is_absolute()
            or ".." in candidate.parts
            or candidate.as_posix() != path
        ):
            raise CheckpointReadError(
                f"Invalid relative path in input manifest: {path!r}."
            )
        if path in seen_paths:
            raise CheckpointReadError(
                f"Duplicate path in input manifest: {path!r}."
            )
        if re.fullmatch(r"[0-9a-f]{64}", sha256) is None:
            raise CheckpointReadError(
                f"Invalid SHA-256 in input manifest for {path!r}."
            )
        if int(size) < 0:
            raise CheckpointReadError(
                f"Negative input size in manifest for {path!r}."
            )
        seen_paths.add(path)
        entries.append(InputManifestEntry(path, sha256, int(size)))
    return tuple(entries)


def _read_conductor(group, group_name):
    if not isinstance(group, h5py.Group):
        raise CheckpointReadError(
            f"conductors/{group_name} must be a group."
        )
    _require_attributes(group, ("identifier",), f"conductors/{group_name}")
    identifier = _text(group.attrs["identifier"], "conductor identifier")
    if not identifier or _safe_name(identifier) != group_name:
        raise CheckpointReadError(
            f"Conductor group {group_name!r} does not match identifier "
            f"{identifier!r}."
        )

    required_groups = (
        "metadata",
        "clock",
        "th_history",
        "electric",
        "balances",
        "components",
        "output_state",
    )
    _require_members(group, required_groups, f"conductor {identifier!r}")
    for name in required_groups:
        if not isinstance(group[name], h5py.Group):
            raise CheckpointReadError(
                f"Conductor {identifier!r} entry {name!r} must be a group."
            )

    metadata = group["metadata"]
    _require_attributes(
        metadata,
        ("th_method", "electric_method"),
        f"conductor {identifier!r} metadata",
    )
    clock = _read_mapping(group["clock"])
    _validate_read_clock(clock, identifier)

    th_history = _read_mapping(group["th_history"])
    for key in ("SYSVAR", "SYSLOD"):
        if key not in th_history:
            raise CheckpointReadError(
                f"Conductor {identifier!r} TH history is missing {key!r}."
            )

    electric = _read_mapping(group["electric"])
    if electric:
        for key in ("electric_solution", "electric_solution_steady"):
            if key not in electric:
                raise CheckpointReadError(
                    f"Conductor {identifier!r} electric state is missing {key!r}."
                )
        for key in ("cond_el_num_step", "electric_time"):
            if key not in clock:
                raise CheckpointReadError(
                    f"Conductor {identifier!r} clock is missing {key!r} for "
                    "an active electric state."
                )
        electric_time = np.asarray(clock["electric_time"])
        if electric_time.ndim != 0 or not np.isfinite(electric_time):
            raise CheckpointReadError(
                f"Conductor {identifier!r}: electric_time must be a finite scalar."
            )
        if clock["cond_num_step"] > 0 and not np.isclose(
            float(electric_time),
            np.asarray(clock["cond_time"])[-1],
            rtol=1.0e-12,
            atol=1.0e-12,
        ):
            raise CheckpointReadError(
                f"Conductor {identifier!r}: electric_time must coincide with "
                "cond_time[-1]."
            )

    balances = _read_mapping(group["balances"])
    for key in _BALANCE_ATTRIBUTES:
        if key not in balances:
            raise CheckpointReadError(
                f"Conductor {identifier!r} balances are missing {key!r}."
            )

    components = _read_components(group["components"], identifier)

    output_state = _read_mapping(group["output_state"])
    for key in ("i_save", "num_step_save", "buffers"):
        if key not in output_state:
            raise CheckpointReadError(
                f"Conductor {identifier!r} output state is missing {key!r}."
            )
    buffers = output_state["buffers"]
    if not isinstance(buffers, Mapping):
        raise CheckpointReadError(
            f"Conductor {identifier!r} output buffers must be a mapping."
        )
    component_buffers = buffers.get("components")
    if not isinstance(component_buffers, Mapping):
        raise CheckpointReadError(
            f"Conductor {identifier!r} output buffers are missing components."
        )
    if set(component_buffers) != set(components):
        raise CheckpointReadError(
            f"Conductor {identifier!r} component inventory does not match "
            "its component output buffers."
        )

    return ConductorCheckpointData(
        identifier=identifier,
        th_method=_text(metadata.attrs["th_method"], "TH method"),
        electric_method=_text(
            metadata.attrs["electric_method"], "electric method"
        ),
        clock=clock,
        th_history=th_history,
        electric=electric,
        balances=balances,
        components=components,
        output_state=output_state,
    )


def _read_components(group, conductor_identifier):
    if len(group) == 0:
        raise CheckpointReadError(
            f"Conductor {conductor_identifier!r} contains no components."
        )

    components = {}
    for group_name in group:
        component_group = group[group_name]
        if not isinstance(component_group, h5py.Group):
            raise CheckpointReadError(
                f"Component entry {group_name!r} must be a group."
            )
        _require_attributes(
            component_group,
            ("identifier", "kind"),
            f"component {group_name!r}",
        )
        identifier = _text(
            component_group.attrs["identifier"], "component identifier"
        )
        if not identifier or _safe_name(identifier) != group_name:
            raise CheckpointReadError(
                f"Component group {group_name!r} does not match identifier "
                f"{identifier!r}."
            )
        if identifier in components:
            raise CheckpointReadError(
                f"Duplicate component identifier {identifier!r}."
            )

        kind = _text(component_group.attrs["kind"], "component kind")
        if kind not in ("fluid", "solid"):
            raise CheckpointReadError(
                f"Component {identifier!r} has invalid kind {kind!r}."
            )
        state = _read_mapping(component_group)
        energy_history = state.get("energy_history")
        if kind == "solid":
            if not isinstance(energy_history, Mapping):
                raise CheckpointReadError(
                    f"Solid component {identifier!r} is missing energy_history."
                )
            for key in ("EEXT", "EJHT"):
                if key not in energy_history:
                    raise CheckpointReadError(
                        f"Solid component {identifier!r} energy history is "
                        f"missing {key!r}."
                    )
        elif energy_history is not None:
            raise CheckpointReadError(
                f"Fluid component {identifier!r} cannot contain energy_history."
            )

        components[identifier] = {"kind": kind, **state}
    return components


def _validate_read_clock(clock, identifier):
    for key in ("cond_time", "cond_num_step", "time_step", "EQTEIG", "i_event", "events_time"):
        if key not in clock:
            raise CheckpointReadError(
                f"Conductor {identifier!r} clock is missing {key!r}."
            )
    try:
        cond_time = _finite_time_history(
            clock["cond_time"], f"conductor {identifier!r}.cond_time"
        )
    except CheckpointValidationError as exc:
        raise CheckpointReadError(str(exc)) from exc
    cond_num_step = _python_non_negative_integer(
        clock["cond_num_step"],
        f"conductor {identifier!r}.cond_num_step",
    )
    if cond_num_step != len(cond_time) - 1:
        raise CheckpointReadError(
            f"Conductor {identifier!r}: cond_num_step ({cond_num_step}) must equal "
            f"len(cond_time) - 1 ({len(cond_time) - 1})."
        )
    time_step = np.asarray(clock["time_step"])
    if time_step.ndim != 0 or not np.isfinite(time_step) or float(time_step) < 0.0:
        raise CheckpointReadError(
            f"Conductor {identifier!r}: time_step must be a finite "
            "non-negative scalar."
        )


def _read_mapping(group):
    result = {}
    for name in group:
        node = group[name]
        original_name = _text(
            node.attrs.get("original_name", name),
            f"original name for {node.name}",
        )
        if original_name in result:
            raise CheckpointReadError(
                f"Duplicate reconstructed key {original_name!r} in {group.name}."
            )
        result[original_name] = _read_node(node)
    return result


def _read_node(node):
    if isinstance(node, h5py.Dataset):
        value = node[()]
        if isinstance(value, bytes):
            return value.decode("utf-8")
        if isinstance(value, np.ndarray):
            if value.dtype.kind in ("S", "O"):
                decoded = np.empty(value.shape, dtype=object)
                for index in np.ndindex(value.shape):
                    item = value[index]
                    decoded[index] = (
                        item.decode("utf-8") if isinstance(item, bytes) else item
                    )
                return decoded
            return value.copy()
        return value.item() if isinstance(value, np.generic) else value

    if not isinstance(node, h5py.Group):
        raise CheckpointReadError(f"Unsupported HDF5 object at {node.name}.")
    python_type = _text(
        node.attrs.get("python_type", "mapping"),
        f"python_type for {node.name}",
    )
    if python_type == "none":
        if len(node):
            raise CheckpointReadError(f"None node {node.name} must be empty.")
        return None
    values = _read_mapping(node)
    if python_type == "mapping":
        return values
    if python_type in ("list", "tuple"):
        expected_keys = [str(index) for index in range(len(values))]
        if set(values) != set(expected_keys):
            raise CheckpointReadError(
                f"Sequence node {node.name} has invalid indices."
            )
        sequence = [values[key] for key in expected_keys]
        return tuple(sequence) if python_type == "tuple" else sequence
    raise CheckpointReadError(
        f"Unsupported python_type {python_type!r} at {node.name}."
    )


def _require_members(group, names, label):
    missing = [name for name in names if name not in group]
    if missing:
        raise CheckpointReadError(
            f"{label} is missing required entries: {', '.join(missing)}."
        )


def _require_attributes(item, names, label):
    missing = [name for name in names if name not in item.attrs]
    if missing:
        raise CheckpointReadError(
            f"{label} is missing required attributes: {', '.join(missing)}."
        )


def _text(value, label):
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    if not isinstance(value, str):
        raise CheckpointReadError(f"{label} must be text.")
    return value


def _strict_bool(value, label):
    if not isinstance(value, (bool, np.bool_)):
        raise CheckpointReadError(f"{label} must be boolean.")
    return bool(value)


def _read_string_vector(dataset, label):
    if not isinstance(dataset, h5py.Dataset):
        raise CheckpointReadError(f"{label} must be a dataset.")
    values = np.asarray(dataset[()])
    if values.ndim != 1:
        raise CheckpointReadError(f"{label} must be a 1-D array.")
    return tuple(_text(value, label) for value in values)


def _read_finite_time_dataset(dataset, label):
    if not isinstance(dataset, h5py.Dataset):
        raise CheckpointReadError(f"{label} must be a dataset.")
    try:
        return _finite_time_history(dataset[()], label).copy()
    except CheckpointValidationError as exc:
        raise CheckpointReadError(str(exc)) from exc


def _read_non_negative_integer(dataset, label):
    if not isinstance(dataset, h5py.Dataset):
        raise CheckpointReadError(f"{label} must be a dataset.")
    return _python_non_negative_integer(dataset[()], label)


def _python_non_negative_integer(value, label):
    array = np.asarray(value)
    if array.ndim != 0 or not np.issubdtype(array.dtype, np.integer):
        raise CheckpointReadError(f"{label} must be an integer scalar.")
    result = int(array)
    if result < 0:
        raise CheckpointReadError(f"{label} cannot be negative.")
    return result


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
    return [
        {
            "path": entry.path,
            "sha256": entry.sha256,
            "size": entry.size,
        }
        for entry in _build_input_manifest_entries(
            base_path, excluded_dir=excluded_dir
        )
    ]


def _build_input_manifest_entries(input_directory, excluded_dir=None):
    """Build a deterministic, detached manifest for one input directory."""

    base_path = Path(input_directory)
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
            InputManifestEntry(
                path=path.relative_to(base_path).as_posix(),
                sha256=_sha256(path),
                size=path.stat().st_size,
            )
        )
    return tuple(manifest)


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

    component_inventory = _component_inventory(conductor, label)
    for component, kind in component_inventory:
        if kind != "solid":
            continue
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
    for component, kind in _component_inventory(
        conductor, f"conductor {conductor.identifier!r}"
    ):
        identifier = getattr(component, "identifier", component.__class__.__name__)
        component_group = components.create_group(_safe_name(identifier))
        component_group.attrs["identifier"] = identifier
        component_group.attrs["kind"] = kind
        if kind == "solid":
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
        child = parent.create_group(safe_name, track_order=True)
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


def _component_inventory(conductor, label):
    fluid_components = _collection(conductor, "FluidComponent")
    solid_components = _collection(conductor, "SolidComponent")
    all_components = _collection(conductor, "all_component")
    if not all_components:
        raise CheckpointValidationError(
            f"{label}: all_component cannot be empty."
        )

    fluid_ids = [id(component) for component in fluid_components]
    solid_ids = [id(component) for component in solid_components]
    all_ids = [id(component) for component in all_components]
    if len(set(fluid_ids)) != len(fluid_ids):
        raise CheckpointValidationError(
            f"{label}: FluidComponent contains duplicate objects."
        )
    if len(set(solid_ids)) != len(solid_ids):
        raise CheckpointValidationError(
            f"{label}: SolidComponent contains duplicate objects."
        )
    if set(fluid_ids) & set(solid_ids):
        raise CheckpointValidationError(
            f"{label}: a component cannot be both fluid and solid."
        )
    if len(set(all_ids)) != len(all_ids):
        raise CheckpointValidationError(
            f"{label}: all_component contains duplicate objects."
        )
    if set(all_ids) != set(fluid_ids) | set(solid_ids):
        raise CheckpointValidationError(
            f"{label}: all_component must contain exactly the fluid and solid "
            "components."
        )

    identifiers = [getattr(component, "identifier", None) for component in all_components]
    if any(not identifier for identifier in identifiers):
        raise CheckpointValidationError(
            f"{label}: every component needs an identifier."
        )
    if len(set(identifiers)) != len(identifiers):
        raise CheckpointValidationError(
            f"{label}: component identifiers must be unique."
        )

    fluid_id_set = set(fluid_ids)
    return [
        (component, "fluid" if id(component) in fluid_id_set else "solid")
        for component in all_components
    ]


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
