# OPENSC2 XLSX input-template finalization audit

Baseline: `benchmark-4c-e2e-verified-20260919`
(`9cbb8158b17161c7da09880f298ad51a834d7a67`).

Working branch: `input_template_xlsx_finalization`.

This branch finalizes the eight XLSX input templates together with the code
changes required to make the exposed controls consistent with their runtime
behavior. The runnable defaults use the ITER TF-like Case Study 1 as their
physical basis, with deliberately shortened time and space discretizations so
that the template can also serve as a quick installation smoke test.

## Reintroduced and implemented controls

- `EPS_INTERPOLATION` is present in `STACK` and `STR_MIX`. External strain
  profiles use this key; `IOP_INTERPOLATION` remains exclusive to current
  profiles.
- `MAXIMUM_ITERATION_NUMBER` is present in `CONDUCTOR_operation`. It is a
  positive-integer safety limit for the number of electric substeps inside one
  thermal-hydraulic step.
- If `ELECTRIC_TIME_STEP` is blank, the existing default of ten electric
  substeps is retained. If it is provided, OPENSC2 treats it as the largest
  allowed electric step, computes the required integer number of substeps, and
  partitions the thermal-hydraulic interval exactly. A request exceeding
  `MAXIMUM_ITERATION_NUMBER` is rejected explicitly.
- `trans_transp_multiplier` is restored to the coupling workbook. For each
  open fluid-fluid interface, it multiplies the transverse transport
  coefficients `K'`, `K''`, and `K'''`. It does not multiply the axial
  Darcy/Fanning friction factor of an individual channel. Values must be
  finite and non-negative and are entered in the upper triangular matrix;
  `1` leaves the transverse model unscaled.

## Archived or hidden controls

- `external_free_convection_correlation` is removed from `CONDUCTOR_input`.
  The incomplete selector is no longer exposed to users; the previous template
  default (`vertical_plate_churchill_chu_accurate`) is retained internally so
  that the established behavior does not change.
- `TIMEREF` and `TAUREF` are removed from `TRANSIENT`. Their only consumer was
  the placeholder `IADAPTIME=-2` function, which returned the previous time
  step without implementing the documented policy. Mode `-2` is no longer an
  accepted input value.
- `BITR` and `BOTR` are removed from `STACK`, `STR_MIX`, `STR_STAB`, and
  `Z_JACKET`. The old `IBIFUN=1` expression was not time dependent and could
  evaluate an undefined `0/0` ratio. Selecting it now raises an explicit
  `NotImplementedError`; supported magnetic-field definitions remain
  `IBIFUN=0` (`BISS`/`BOSS`) and `IBIFUN=-1` (`EXTERNAL_BFIELD`).

## Earlier cleanup retained

- Correct STACK keys: `N_tape`, `Stack_width`, and
  `superconducting_material`.
- Conditional CryoSoft fields: `STACK.RRR_Ag` and `Z_JACKET.RRR`.
- Removed unused or incomplete inputs: `TEMINI`, `PREINI`, `fixAlphaBvalue`,
  `QJFRACT`, `ISJOINT`, `XJBEG`, `XJBEIN`, `XJBEOUT`, and `XJENOUT`.
- Corrected `STR_MIX` example values and material documentation.
- Standalone workbook headers without external-workbook references.
- Checkpoint controls and the `CHECKPOINTS` sheet remain available.
- Diagnostic workbooks contain only the required object column and no trailing
  unnamed columns.

## ITER TF-like template basis

The active example represents a two-channel ITER TF-like conductor with one
Nb3Sn mixed strand component and one jacket component. Geometry, material,
magnetic-field, strain, hydraulic, coupling, and localized-heating inputs are
derived from Case Study 1. The following values are intentionally reduced for
the runnable template:

- final time: `0.01 s`;
- fixed thermal-hydraulic time step: `0.01 s`;
- spatial elements: `20`;
- active-component real-time plots disabled;
- diagnostic times aligned with the shortened run;
- transport current set to zero for a fast, numerically benign smoke case.

These reductions do not change the workbook schema or remove the Case Study 1
component topology.

## Verification results

### Workbook checks

- All eight workbooks reopen successfully.
- No spreadsheet formula-error token was detected.
- All sheets were rendered and visually inspected after export.
- OPENSC2-facing cached headers are available with `data_only=True`.
- The diagnostic sheets contain no surplus unnamed columns.
- The installed workbook SHA-256 hashes match the reviewed candidate package.

### Automated tests

- Modified Python modules pass `compileall`.
- The focused template and continuation group passes: `39 passed`.
- The complete automated suite passes: `232 passed`.
- `git diff --check` reports no whitespace error.

### Headless smoke test

The isolated template run was executed on Windows with Python 3.10.10 in the
`opensc2_mkl` environment.

- Exit code: `0`.
- OPENSC2 completion marker found.
- Elapsed process time: `20.2962747 s`.
- Simulated interval: one thermal-hydraulic step to `0.01 s`.
- Electric integration: ten substeps of `0.001 s`.
- Generated files: `176` (`101` SVG, `67` TSV, and `8` XLSX).
- No traceback or unhandled exception.
- No invalid scalar division caused by zero mass flow.
- No outlet-pressure overwrite warning.
- The copied input workbooks remained byte-for-byte unchanged.
- The repository state remained unchanged.

The remaining messages are expected initialization notices, small hydraulic
balancing corrections for the two parallel channels, and pre-existing pandas
`FutureWarning`/`PerformanceWarning` messages. They do not invalidate the
template run.

## Input-schema compatibility

This branch intentionally introduces a new input schema. Existing case-study
input sets must be reviewed before they are used with this code version. In
particular:

- every conductor definition requires `MAXIMUM_ITERATION_NUMBER`;
- every coupling workbook requires `trans_transp_multiplier`;
- active `STACK` and `STR_MIX` strain profiles require
  `EPS_INTERPOLATION` when `IEPS=-1`;
- `IADAPTIME=-2` and `IBIFUN=1` are no longer accepted.

The successful template smoke test does not by itself validate older case-study
input sets. Their migration and regression execution remain a separate task.

## Repository and release status

This audit validates the candidate on the feature branch only. It does not
constitute a merge into `development` or `main`, and it is not a public
versioned release. No release decision should be taken before the existing
case-study inputs have been migrated and exercised.
