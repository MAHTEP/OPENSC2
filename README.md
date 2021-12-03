# OPENSC2 in a nutshell

Object-oriented software for multiphysics simulations of Superconducting cables. 

## Features

OPENSC2 is a software for the multi-physical analysis of thermal-hydraulic and electro-dynamic transients in Superconducting Cable-in-Conduit Conductors (CICC) for fusion magnets and power transmission.

Currently it is developed mainly in [Python](https://www.python.org/) but future versions will possibly take advantage of other programming languages such as [TypeScript](https://www.typescriptlang.org/) and [Rust](https://www.rust-lang.org/) as well as the [OpenModelica](https://www.openmodelica.org/) environment. 

The software is built based on well-established numerical models and assumptions, re-arranged in an object-oriented framework to be user-friendly and easily manageable through a GUI. The input set is to be prescribed through self-explanatory excel files.
The developing team includes Prof. L. Savoldi[^1], Prof. F. Freschi, D. Placido[^2], S. Viarengo[^2]  @ Dipartimento Energia “Galileo Ferraris” @ [Politecnico di Torino](https://www.polito.it/). Please, contact us at:
* laura.savoldi@polito.it
* fabio.freschi@polito.it 
* daniele.placido@polito.it 
* sofia.viarengo@polito.it.

[^1]: Head of the [**MAHTEP** research group](http://www.mahtep.polito.it/).
[^2]: PhD students @ the [**MAHTEP** research group](http://www.mahtep.polito.it/).
## Goals

The software is useful for steady state and transient analyses of CICC in operating conditions. It can deal with cables assemled with Low Temperatures (LTS) strands (both Nb3Sn and NbTi), and High Temperature Superconductors (HTS) tapes of different materials. Different coolants can be selected, together with very different cooling configurations. The software is useful to study the steady state operating conditions under environmental parasitic load, as well as transient operation such as: current variation in time, coolant flow variation in time, AC losses, quench, fast discharges, fault currents. The software is useful to assist the research for optimal configurations, subject to a set of constraints, and allows evaluating the temperature margin to current sharing along cables in any pre-defined operating scenarios.

# Get started

Users can benefit from several test cases to check the software functionalities:
1. Heat slug propagation in an ITER TF-like CICC
2. Heat slug propagation in a stacked-HTS slotted-core CICC for fusion applications
3. Steady state operation for a double-cryostat HVDC cable for power transmission

To run a simulation with one of the above test cases, download the repository and install the requirements (more informations in section [Install requirements](user-content-intall-requirements)). After that, you can run the software (the starting file is _simulation_starter.py_) and from the GUI you can navigate through the folder three until you enter directory _TDD_examples_ and then select one of the three folders contained with pre-compiled inpuput files. In the GUI window select **Add solution path** to select where to save the results (by default they are all collected in the directory _Simulation_results_, that is automatically created if does not already exist). User can create a new folder in this directory or open an existing one: the output (both .tsv files and .eps figures) will be saved in this folder.

## Install requirements

The selected Python version is [3.8.10](https://www.python.org/downloads/release/python-3810/). To install the requirements, create a virtual environment (suggested name _opensc2_) and activate it. In your terminal run the following command:

    python -m pip install --upgrade pip \\ to update pip to the last version  
    python -m pip install -r requirements.txt \\ to install the requirements  

Among the dependences there is [CoolProp](http://www.coolprop.org/) that, according to the operative system you use, may require some other dependences and/or packages. To deal with this, please follow the [documentation](http://www.coolprop.org/coolprop/wrappers/Python/index.html) and [prerequisites](http://www.coolprop.org/coolprop/wrappers/index.html#wrapper-common-prereqs).

# Help

Software documentation is under development, being the project at its initial stages. Detailed documentation will be provided as soon as we get an established version of the software is available.
For the time being feel free to send an e-mail to daniele.placido@polito.it if you need any help with your simulations.
Being currently an embryonic software, some of the possibilities provided in the input files may not yet be fully implemented or tested and you may get incorrect results or unexpected errors. 
The development team apologizes for the inconvenience and is committed to fixing them as soon as possible.

# Contribution

The developing team wish to receive help form the users in the definition and test of new test cases, in the benchmark against other established software, in the inclusion of other functionalities.
To contribute please refer to [contribution](CONTRIBUTION.md).

# Code of Conduct

The developing team agreed to embrace the [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md) **Code of Conduct**.
 
 # License
 OPENSC2 is licensed under [![AGPL](https://www.gnu.org/graphics/agplv3-with-text-100x42.png)](LICENSE) or any other version of it.
