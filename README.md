# Alquimia

This code presents the implementation of the Alquimia project which aims to integrate metabolic models within a 
superstructure optimization framework. The code is written in Python and uses the [Pyomo](http://www.pyomo.org/)
to generate the superstructure and [GAMs](https://www.gams.com/) to solve the optimization problem. 

The publication of this work can be found in the following link: [Advancing Biorefinery Design through the Integration 
of Metabolic Models in Superstructure Optimization](https://www.sciencedirect.com/science/article/pii/S1369703X21002554)

The code is structured in the following way:
- file: *'case study propionate'* contains the code to generate the superstructure for the propionate case studies
  - 'case_study_v1.py' contains the code to generate and solve the superstructure of case study 1 
  - 'case_study_v2.gms' contains the code to generate and solve the superstructure of case study 2
- file: *'json models'* contains the surrogate models for the case studies
- file: *'SBML models'* contains the SBML models (GEMs) for the case studies
- file: *'SBML screening'* contains scripts for the curration of the sbml models
- file: *'excel files'* contains the Excel files with the data for the case studies
- scripts to build the superstructures and surrogate models
  - f_make_super_structure.py: contains the classes and functions to generate the superstructure
  - f_make_surrogate_models.py: contains the classes and functions to generate the surrogate models
  - f_screen_SBML.py: contains the classes and functions to screen the sbml models



# Software and Python Version

- Python Version: 3.10.1
- Packages see python_packages.md
- Gams Version: 37.1.0

To install the packages (OS: windos) run the script 'install_packages.py'.
Make sure you have pip installed and accessible in your system's PATH. If you are using a virtual environment, activate 
it before running the script to ensure that the packages are installed within the virtual environment.

## Acknolwedgements

This work was supported by project ALQUIMIA(PID2019-110993RJ-I00)
funded by the Agencia Estatal de Investigación Alquimia: Proyecto de I- D-i Programa Retos de la sociedad modalidad
Jovenes investigadores convocatoria.

A. Regueira would like to acknowledge the support of the Xunta de Galicia through a postdoctoral fellowship
(ED481B-2021-012). The authors belong to the Galician Competitive Research Group ED431C-2021/37, cofounded by ERRF (EU).