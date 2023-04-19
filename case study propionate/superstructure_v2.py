"""
@author: Lucas Van der Hauwaert
lucas.vanderhauwaert@usc.es
"""

from f_make_super_structure import make_super_structure, solve_model

excelFile = 'propionate_case_study_v2.xlsx'
superstructure = make_super_structure(excelFile= excelFile, printPyomoEq= False)

switchSolver = True
operatingDays= 1/24
if switchSolver:
    results = solve_model(superstructure,  operatingDays = operatingDays)
    #results = solve_model(superstructure, operatingDays = operatingDays ,solverType='BARON')
