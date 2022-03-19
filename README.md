# parametric_cgx_model

## Description

Parametric shell half-wing model geometry creation and meshing in CalculiX GraphiX, followed by FEA analysis (using Nastran or CalculiX CrunchiX) and analysis output recovery and processing (implemented for CalculiX CrunchiX only).

The script assumes that CalculiX and Nastran are installed and working. Adjust the executable calls as needed in `LOCAL_EXECUTES`.

**Background and explainer video for this script are available on our Dapta blog:**

* 16/12/2021: [Parametric FEM model creation with Python and CalculiX GraphiX (cgx)](https://www.dapta.com/parametric-fem-model-creation-with-python-and-calculix-graphix-cgx/);
Corresponding code release: [v0.0.1](https://github.com/daptablade/parametric_cgx_model/releases/tag/v0.0.1)
* 21/01/2022: [Automated FEM analysis and output processing with Python and CalculiX CrunchiX (ccx)](https://www.dapta.com/automated-fem-analysis-and-output-processing-with-python-and-calculix-crunchix-ccx/) Corresponding code release: [v0.0.2](https://github.com/daptablade/parametric_cgx_model/releases/tag/v0.0.2)
* 15/02/2022: [Design optimisation in practice with Python, OpenMDAO and Scipy](https://www.dapta.com/design-optimisation-in-practice-with-python-openmdao-and-scipy/)
Corresponding code release: [v0.0.3](https://github.com/daptablade/parametric_cgx_model/releases/tag/v0.0.3)

## Quick Start

1. Make sure you have CalculiX GraphiX (CGX) installed.

2. Nastran or CalculiX CrunchiX are only required if you want to run FEM analyses (see References section below).

3. Setup a virtual python environment with Venv (optional, but very much recommended - see basics below).

4. Execute the python script 'parametric_box.py' and inspect outputs: choose between `main(INPUT[0])` for a metallic Nastran model and `main(INPUT[1])` for a composite Calculix model.

5. Execute the python script 'parametric_studies.py' and inspect outputs: this script executes multiple iterations of the main() function from parametric_box.py with incremental input changes.
Choose between `parametric_study_material(INPUTS[1])` for a composite material design study and `parametric_study_mesh(INPUTS[1])` for a mesh refinement study.

6. Execute the python script 'design_optimisation.py' and inspect outputs:
this script executes a design optimisation using the [OpenMDAO](https://openmdao.org/) library. See our Feb 2022 blog post on the topic: [Design optimisation in practice with Python, OpenMDAO and Scipy](https://www.dapta.com/design-optimisation-in-practice-with-python-openmdao-and-scipy/)

Note: analysis output processing is only implemented for the Calculix analysis using INPUT[1].

## Venv virtual environment

We use a Venv to manage the imported python libraries - [Find out more about Venv here.](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/#creating-a-virtual-environment).

To create and configure the virtual environment in Windows for example you can execute the following:

```
python3 -m venv venv
venv\Scripts\activate
python -m pip install -r requirements.txt
```

Where the requirements.txt lists the libraries this project is dependent on.

## References

* CalculiX GraphiX and Calculix CrunchiX (V2.15 or later): [http://www.dhondt.de/](http://www.dhondt.de/)

* MSC Nastran (v2020 or later): [https://www.mscsoftware.com/product/msc-nastran](https://www.mscsoftware.com/product/msc-nastran)

* OpenMDAO: [https://openmdao.org/](https://openmdao.org/)

## Author and License

Dr Olivia Stodieck, olivia.stodieck@dapta.com

Dapta LTD, UK

Visit us at [www.dapta.com](www.dapta.com)

This code is released under an MIT license - see LICENSE file.
