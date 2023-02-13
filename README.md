# parametric_cgx_model

**Feb 2023 Update – Read the updated tutorials for the dapta app here: https://daptadocs.com**
This includes new files and step by step walk-throughs for [chaining components](https://daptadocs.com/Tutorials/Chaining%20component%20analyses.html), automating the parametric [design studies with Calculix](https://daptadocs.com/Tutorials/Automating%20parametric%20studies.html) and design [optimisation using OpenMDAO](https://daptadocs.com/Tutorials/Automating%20design%20optimisations.html). 

## Description

The aim is to demonstrate various engineering simulation processes using python3 and open source software. The simulation processes are developed around the example of a simple parametric composite wing model that is fixed at the root.

Core engineering simulation processes explored include:

* Geometry creation and meshing in CalculiX GraphiX
* FEA analysis (using Nastran or CalculiX CrunchiX) and analysis output recovery and processing (implemented for CalculiX CrunchiX only).
* Design optimisation using OpenMDAO
* Coupling of disciplinary analyses (here: structures and aero) for static analysis using preCICE  
* Calculate the aeroelastic flutter and divergence velocities from FEA exported stiffness and mass matrices and a simple unsteady aero model

The script assumes that CalculiX GraphiX / CrunchiX, Nastran (optional) and preCICE (optional) are installed and working.
Adjust the executable calls as needed, for example in [parametric_box.py](./parametric_box.py) global variable `LOCAL_EXECUTES`.

The examples are written for python3 on Windows and assuming access to Windows Systems for Linux 2 (`wsl`). Running on linux should be possible, but you will have to update any solver or subprocess calls that include the `wsl` command to local linux calls.

**Background and explainer videos for this script are available on our Dapta blog:**

* 16/12/2021: [Parametric FEM model creation with Python and CalculiX GraphiX (cgx)](https://www.dapta.com/parametric-fem-model-creation-with-python-and-calculix-graphix-cgx/);
Corresponding code release: [v0.0.1](https://github.com/daptablade/parametric_cgx_model/releases/tag/v0.0.1)
* 21/01/2022: [Automated FEM analysis and output processing with Python and CalculiX CrunchiX (ccx)](https://www.dapta.com/automated-fem-analysis-and-output-processing-with-python-and-calculix-crunchix-ccx/) Corresponding code release: [v0.0.2](https://github.com/daptablade/parametric_cgx_model/releases/tag/v0.0.2)
* 15/02/2022: [Design optimisation in practice with Python, OpenMDAO and Scipy](https://www.dapta.com/design-optimisation-in-practice-with-python-openmdao-and-scipy/)
Corresponding code release: [v0.0.3](https://github.com/daptablade/parametric_cgx_model/releases/tag/v0.0.3)
* 24/03/2022: [Multidisciplinary analysis with Python, PreCICE and Calculix Crunchix](https://www.dapta.com/multidisciplinary-analysis-with-python-precice-and-calculix-crunchix/) 
Corresponding code release: [v0.0.4](https://github.com/daptablade/parametric_cgx_model/releases/tag/v0.0.4)
* 30/06/2022: [Dynamic aeroelastic flutter and divergence analysis with Python and Calculix Crunchix](https://www.dapta.com/dynamic-aeroelastic-flutter-and-divergence-analysis-with-python-and-calculix-crunchix)
Corresponding code release: [v0.0.5](https://github.com/daptablade/parametric_cgx_model/releases/tag/v0.0.5.fix2)

## Quick Start

1. Make sure you have CalculiX GraphiX (CGX) installed. Nastran or CalculiX CrunchiX are only required if you want to run FEM analyses (see References section below).

3. Setup a virtual python environment with Venv (optional, but very much recommended - see basics below).

4. Execute the python script [parametric_box.py](./parametric_box.py) and inspect outputs: choose between `main(INPUT[0])` for a metallic Nastran model and `main(INPUT[1])` for a composite Calculix model.

5. Execute the python script [parametric_studies.py](./parametric_studies.py) and inspect outputs: this script executes multiple iterations of the main() function from `parametric_box.py` with incremental input changes.
Choose between `parametric_study_material(INPUTS[1])` for a composite material design study and `parametric_study_mesh(INPUTS[1])` for a mesh refinement study.

6. Execute the python script [design_optimisation.py](./design_optimisation.py) and inspect outputs:
this script executes a design optimisation using the [OpenMDAO](https://openmdao.org/) library. Warning: execution takes ~5-10min ...

7. Execute a coupled aero/structures analysis using preCICE, CalculiX CrunchiX and a simplified static strip theory aero model. First install preCICE - see separate instructions here: [preCICE README](./precice_modules/README.md). Then execute [aeroelastic_analysis.py](./aeroelastic_analysis.py). Warning: Again, this takes a bit of time as both solvers will execute several times until the loads are converged.

Note: FEM analysis output processing is only implemented for the Calculix analyses using INPUT[1] or INPUT[3].

## Venv virtual environment

We use a Venv to manage the imported python libraries - [Find out more about Venv here.](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/#creating-a-virtual-environment).

To create and configure the virtual environment in Windows execute the following in the command line:

```
> python3 -m venv venv
> venv\Scripts\activate
> python -m pip install -r requirements.txt
```

Where the requirements.txt lists the libraries this project is dependent on.

For preCICE, it is recommended you create a separate wsl virtual environment under the `precice_modules` folder. See separate preCICE installation instructions here: [preCICE README](./precice_modules/README.md).

## References

* CalculiX GraphiX and Calculix CrunchiX (V2.15 or later): [http://www.dhondt.de/](http://www.dhondt.de/)

* MSC Nastran (v2020 or later): [https://www.mscsoftware.com/product/msc-nastran](https://www.mscsoftware.com/product/msc-nastran)

* OpenMDAO: [https://openmdao.org/](https://openmdao.org/)

* preCICE: [https://precice.org/index.html](https://precice.org/index.html)

## Author and License

Dr Olivia Stodieck, olivia.stodieck@dapta.com

Dapta LTD, UK

Visit us at [www.dapta.com](https://www.dapta.com/)

This code is released under an MIT license - see LICENSE file.
