# parametric_cgx_model

## Description
Parametric shell half-wing model geometry creation and meshing in CalculiX GraphiX (cgx), followed by MSC Nastran normal modes (SOL103) analysis. 
The script assumes that cgx and Nastran are installed and working. Adjust the executable calls in the python script as needed.  

## Quick Start
1. Make sure you have CGX installed. 
2. Nastran is only required if you want to run the FEM analysis at the end. 
3. Execute the python script 'parametric_box.py' and inspect outputs. 

## References
CalculiX GraphiX (V2.15 or later): http://www.dhondt.de/
MSC Nastran (v2020 or later): https://www.mscsoftware.com/product/msc-nastran

## Author and License
Dr Olivia Stodieck, olivia.stodieck@daptablade.com
Daptablade LTD, UK
This code is released under an MIT license - see LICENSE file. 
