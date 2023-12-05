# Alquimia: Biorefinery Design and Optimisation

ALQUIMIA is a proof-of-concept software developed for the early-stage design of biorefineries. It leverages superstructure optimization to optimize the flow sheet of biorefinery processes. This approach is crucial in identifying efficient and cost-effective pathways in complex systems where multiple process routes are possible. A unique feature of ALQUIMIA is the integration of metabolic models within this optimization framework as surrogate models. These models represent the biochemical processes in bioreactors, a key component of biorefineries, allowing for a more holistic and accurate optimization of the entire process.

The code is written in Python, utilizing [Pyomo](http://www.pyomo.org/) to generate the superstructure and [GAMS](https://www.gams.com/) for solving the optimization problem.

## Content
The code is organized as follows:
- **File: 'Case Study Propionate'**: Contains the code to generate the superstructure for the propionate case studies.
  - `case_study_v1.py`: Code to generate and solve the superstructure of case study 1.
  - `case_study_v2.gms`: Code to generate and solve the superstructure of case study 2.
- **File: 'JSON Models'**: Contains the surrogate models for the case studies.
- **File: 'SBML Models'**: Includes the SBML models (GEMs) for the case studies.
- **File: 'SBML Screening'**: Scripts for the curation of the SBML models.
- **File: 'Excel Files'**: Excel files with data for the case studies.
- **Scripts for Building Superstructures and Surrogate Models**:
  - `f_make_super_structure.py`: Classes and functions to generate the superstructure.
  - `f_make_surrogate_models.py`: Classes and functions to generate the surrogate models.
  - `f_screen_SBML.py`: Classes and functions to screen the SBML models.

## Software and Python Version
- Python Version: 3.10.1
- For package details, see `python_packages.md`.
- Gams Version: 37.1.0

## Publication
DOI available upon publication.

## Acknowledgements
This work was supported by the ALQUIMIA project (PID2019-110993RJ-I00), funded by the Agencia Estatal de Investigación under the Programa Retos de la sociedad, modalidad Jovenes investigadores, convocatoria 2019.

A. Regueira acknowledges the support of the Xunta de Galicia through a postdoctoral fellowship (ED481B-2021-012). The authors are part of the Galician Competitive Research Group ED431C-2021/37, co-funded by the ERRF (EU).
