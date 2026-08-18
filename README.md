# OPENSC2

## OPENSC2 in a nutshell

Object-oriented software for multiphysics simulations of Superconducting cables.

### Features

OPENSC2 is a software for the multi-physical analysis of thermal-hydraulic and electro-dynamic transients in Superconducting Cable-in-Conduit Conductors (CICC) for fusion magnets and power transmission.

Currently it is developed mainly in [Python](https://www.python.org/) but future versions will possibly take advantage of other programming languages such as [TypeScript](https://www.typescriptlang.org/) and [Rust](https://www.rust-lang.org/) as well as the [OpenModelica](https://www.openmodelica.org/) environment.

The software is built based on well-established numerical models and assumptions, re-arranged in an object-oriented framework to be user-friendly and easily manageable through a GUI. The input set is to be prescribed through self-explanatory excel files.
The developing team includes Prof. L. Savoldi[^1], Prof. F. Freschi, D. Placido[^2], S. Viarengo[^2] @ Dipartimento Energia “Galileo Ferraris” @ [Politecnico di Torino](https://www.polito.it/). Please, contact us at:

* laura.savoldi@polito.it
* fabio.freschi@polito.it
* daniele.placido@polito.it
* sofia.viarengo@polito.it.

[^1]: Head of the [**MAHTEP** research group](http://www.mahtep.polito.it/).
[^2]: PhD students @ the [**MAHTEP** research group](http://www.mahtep.polito.it/).

### Goals

The software is useful for steady state and transient analyses of CICC in operating conditions. It can deal with cables assemled with Low Temperatures (LTS) strands (both Nb3Sn and NbTi), and High Temperature Superconductors (HTS) tapes of different materials. Different coolants can be selected, together with very different cooling configurations. The software is useful to study the steady state operating conditions under environmental parasitic load, as well as transient operation such as: current variation in time, coolant flow variation in time, AC losses, quench, fast discharges, fault currents. The software is useful to assist the research for optimal configurations, subject to a set of constraints, and allows evaluating the temperature margin to current sharing along cables in any pre-defined operating scenarios.

A detailed description of the physics and of the first tests carried out for the initial phase of verification and validation of the software is available [here](https://doi.org/10.1016/j.cryogenics.2022.103457).

## Get started

Users can benefit from several test cases to check the software functionalities:

1. Heat slug propagation in an ITER TF-like CICC
2. Heat slug propagation in a stacked-HTS slotted-core CICC for fusion applications
3. Steady state operation for a double-cryostat HVDC cable for power transmission

To run a simulation with one of the above test cases, download the repository and install the requirements (more information is provided in [Install requirements](#install-requirements)). Then enter the `source_code` directory and start OPENSC2 with `python opensc2.py`. From the GUI, navigate to the `TDD_examples` directory and select one of the folders containing pre-compiled input files. Select **Add solution path** to choose where the results will be saved. By default, they are collected in the `Simulation_results` directory, which is created automatically when needed. Both tabular output and figures are written below the selected results directory.

### Headless simulations and checkpoint restoration

OPENSC2 can also run without opening the GUI. Headless execution reads the input and output directories from a YAML file such as:

```yaml
input_dir: 'C:\path\to\input_files'
output_dir: 'C:\path\to\simulation_results'
```

Run the following commands from the `source_code` directory. A normal headless simulation starts from the initial state:

```powershell
python opensc2.py `
    --no-head `
    --io-path "C:\path\to\io_path.yaml"
```

Checkpoint restoration is currently available only through this headless command-line interface. Strict recovery resumes the original trajectory and therefore requires an input directory whose files match the checkpoint manifest exactly. A new output directory may be selected in `io_path.yaml`:

```powershell
python opensc2.py `
    --no-head `
    --io-path "C:\path\to\io_path.yaml" `
    --checkpoint "C:\path\to\checkpoint_step_000100.h5" `
    --restart-mode recovery
```

Continuation branches from the saved physical state while retaining permitted definitions from a freshly initialized input configuration. The input directory should normally point to a copy of the original input set in which only the intended time-policy or driver values have been changed. Supported examples include the fixed/adaptive time-step policy and bounds, final time, checkpoint controls, conductor electric time step, current, magnetic field, axial magnetic-field gradient, and external heating. Immutable model changes are rejected before the checkpoint state is applied:

```powershell
python opensc2.py `
    --no-head `
    --io-path "C:\path\to\io_path.yaml" `
    --checkpoint "C:\path\to\checkpoint_step_000100.h5" `
    --restart-mode continuation
```

The `--checkpoint` and `--restart-mode` options must be provided together, and their use requires `--no-head`. In continuation mode, the new final time must be later than the checkpoint time. A separate output directory is recommended for every recovery or continuation branch.

The headless normal, recovery, and continuation paths are covered by automated unit and integration tests and by real end-to-end smoke tests. Developers can run the complete automated suite from `source_code` after installing `pytest`:

```powershell
python -m pytest -q
```

### Install requirements

The selected Python version is [3.10.10](https://www.python.org/downloads/release/python-31010/). From the `source_code` directory, create and activate a virtual environment before installing the dependencies. The suggested environment name is `opensc2_mkl`:

```powershell
python -m venv opensc2_mkl
.\opensc2_mkl\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m pip check
```

The electrical model solves sparse linear systems with [PyPardiso](https://github.com/haasad/PyPardiso), a Python interface to the PARDISO solver provided by the [Intel oneAPI Math Kernel Library (oneMKL)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html). The pinned PyPardiso and oneMKL runtime packages are included in `requirements.txt`; on supported Windows and Linux systems, the command above normally installs the required runtime libraries without a separate Intel oneAPI Toolkit installation.

After installation, the OPENSC2 PyPardiso/oneMKL path can be checked with a small functional solve:

```powershell
python -c "import numpy as np; from scipy.sparse import eye; from utility_functions.electric_auxiliary_functions import pardiso_spsolve; x = pardiso_spsolve(eye(2, format='csr'), np.array([1.0, 2.0])); np.testing.assert_allclose(x, [1.0, 2.0]); print('OPENSC2 PyPardiso/oneMKL: OK')"
```

If PyPardiso cannot locate or load the oneMKL runtime on the target platform, consult the [PyPardiso installation documentation](https://github.com/haasad/PyPardiso#installation) and the official [Intel oneMKL installation options](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html). Avoid copying individual DLL files into the OPENSC2 source tree; use the supported package or runtime installation procedures instead.

Among the dependencies there is [CoolProp](http://www.coolprop.org/), which may require additional packages depending on the operating system. Please refer to the [CoolProp Python documentation](http://www.coolprop.org/coolprop/wrappers/Python/index.html) and its [common prerequisites](http://www.coolprop.org/coolprop/wrappers/index.html#wrapper-common-prereqs).

## Help

Software documentation is under development, being the project at its initial stages. Detailed documentation will be provided as soon as an established version of the software is available.
For the time being feel free to send an e-mail to daniele.placido@polito.it if you need any help with your simulations.
Being currently an embryonic software, some of the possibilities provided in the input files may not yet be fully implemented or tested and you may get incorrect results or unexpected errors. A (not exhaustive) list of known issues is available in the [Issue](https://github.com/MAHTEP/OPENSC2/issues) section. To open a new issue, please [follow the procedure](https://github.com/MAHTEP/OPENSC2/blob/main/CONTRIBUTION.md).  
The development team apologizes for the inconvenience and is committed to fixing them as soon as possible.

## Contribution

The developing team wish to receive help form the users in the definition and test of new test cases, in the benchmark against other established software, in the inclusion of other functionalities.
To contribute please refer to [contribution](CONTRIBUTION.md).

## Code of Conduct

The developing team agreed to embrace the [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md) **Code of Conduct**.

## License

OPENSC2 is licensed under [![AGPL](https://www.gnu.org/graphics/agplv3-with-text-100x42.png)](LICENSE) or any other version of it.
