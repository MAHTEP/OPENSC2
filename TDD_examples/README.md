# Running the OPENSC2 test cases

This directory contains four test cases that exercise different OPENSC2
models and cable configurations. The cases provide ready-to-run input files
for functional checks, examples, and numerical regression testing.

Run all commands in this guide from the `source_code` directory.

## Prerequisites

Use Python 3.10.10 and install the pinned dependencies in a virtual
environment as described in the main [README](../README.md#install-requirements).

On Windows PowerShell, activate the recommended environment with:

```powershell
.\opensc2\Scripts\Activate.ps1
python -m pip check
```

`pip check` should report `No broken requirements found`.

Keep simulation results outside the tracked test-case directories. This avoids
mixing generated files with the input data committed to the repository.

## Available cases

| Case | Input directory | Purpose | End time | Time step |
| --- | --- | --- | ---: | ---: |
| 1 | `CASE_1_ITER_like_LTS` | Heat-slug propagation in an ITER TF-like low-temperature-superconductor CICC | 100 s | 0.1 s |
| 2 | `CASE_2_ENEA_HTS_CICC` | Heat-slug propagation in a stacked-HTS slotted-core CICC for fusion applications | 50 s | 0.05 s |
| 3 | `CASE_3_HTS_HVDC` | Operation of a double-cryostat HTS HVDC cable for power transmission | 20000 s | 10 s |
| 4a | `CASE_4_4C_ELECTRIC_BENCHMARK/stabilization` | Stabilization stage of the 4C electrical benchmark | 610 s | 1 s |
| 4b | `CASE_4_4C_ELECTRIC_BENCHMARK/continuation` | Short electrical transient using the supplied component-current histories | 614 s | 0.00125 s |

All supplied cases use the fixed-time-step policy (`IADAPTIME=0`).

## Run a case with the graphical interface

Start OPENSC2 from `source_code`:

```powershell
python opensc2.py
```

In the graphical interface:

1. Select the directory containing the input workbooks.
2. Select **Add solution path** and choose a separate output location.
3. Load the input data.
4. Start the simulation.

For example, the Case 1 input directory is:

```text
..\TDD_examples\CASE_1_ITER_like_LTS
```

The simulation identifier stored in `transitory_input.xlsx` is used to name
the simulation-specific output directory.

## Run a case without the graphical interface

Headless execution uses a YAML file containing the input directory and the
output root:

```yaml
input_dir: "../TDD_examples/CASE_1_ITER_like_LTS"
output_dir: "C:/path/to/OPENSC2_results"
```

Save the YAML file outside the tracked test-case directories and run:

```powershell
python opensc2.py `
    --no-head `
    --io-path "C:\path\to\io_path.yaml"
```

Change `input_dir` to run another case. Use a different output root or an empty
output directory for each execution.

A successful run prints the input directory, output directory, simulation
identifier, progress information, and final elapsed time.

## Case 4: stabilization and continuation

Case 4 is distributed as two related input sets.

### Stabilization

The `stabilization` input set runs to 610 s with a 1 s thermal-hydraulic time
step. The conductor and electrically active components use
`I0_OP_MODE=0` and `IOP_MODE=0`, so OPENSC2 initializes the component currents
through its current-divider calculation.

### Continuation

The `continuation` input set runs to 614 s with a 0.00125 s time step. Selected
component currents are prescribed through:

```text
conductor_current_dump.xlsx
```

The conductor uses `I0_OP_MODE=-1`; the stack and stabilizer components use
`IOP_MODE=-1`. The supplied current workbook contains the time histories
required by the continuation case.

Because this workbook is already included, the distributed continuation case
can be run directly. Running the stabilization case first is not an operational
prerequisite for executing the supplied continuation input set.

`conductor_current_dump.xlsx` is an input current profile. It is not an OPENSC2
HDF5 checkpoint and it is not passed through the `--checkpoint` command-line
option.

Both Case 4 stages enable checkpoint creation in `transitory_input.xlsx`.
These recovery checkpoints are independent of the current-profile workbook.
Use separate output directories for stabilization and continuation.

## Checkpoint recovery and continuation

Checkpoint restoration is available only in headless mode. Recovery resumes
the original trajectory:

```powershell
python opensc2.py `
    --no-head `
    --io-path "C:\path\to\io_path.yaml" `
    --checkpoint "C:\path\to\checkpoint_step_000100.h5" `
    --restart-mode recovery
```

Continuation branches from a saved physical state while retaining permitted
time-policy and driver changes from the new input set:

```powershell
python opensc2.py `
    --no-head `
    --io-path "C:\path\to\io_path.yaml" `
    --checkpoint "C:\path\to\checkpoint_step_000100.h5" `
    --restart-mode continuation
```

The `--checkpoint` and `--restart-mode` options must be supplied together.
Refer to the main [README](../README.md#headless-simulations-and-checkpoint-restoration)
for the compatibility requirements.

## Outputs and validation

OPENSC2 writes the complete simulation history, spatial distributions, plots,
copied input files, and any enabled checkpoints below the selected output
root.

The tracked `Solution` directories included with Cases 1-3 contain only the
final spatial distributions of selected variables. They are useful for a
quick inspection, but they are not complete regression references because
they do not contain the full simulation history.

Complete reference outputs for Cases 1-4 are being prepared as a separate
versioned dataset. The dataset is intended to be published on Zenodo so that
a specific release can be cited through a DOI without adding large binary
archives to the Git repository. The DOI and dataset license will be added
here after publication.

When reference outputs are available, validate at least:

1. successful completion without an exception;
2. presence of the expected component-output files;
3. coverage of the expected final simulation time;
4. agreement of the numerical time histories and spatial distributions within
   the tolerances documented with the reference dataset.

Do not use byte-for-byte equality for generated figures or workbook
containers. Numerical regression checks should compare the underlying tabular
data with explicit absolute and relative tolerances.

## Automated tests

The automated test suite checks the core code paths, including headless
execution and checkpoint handling, but it does not replace the full reference
simulations:

```powershell
python -m pytest -q
```

Run the automated suite before integrating changes into a shared branch.
