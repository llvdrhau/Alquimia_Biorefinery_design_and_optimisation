"""
@author: Lucas Van der Hauwaert
lucas.vanderhauwaert@usc.es
"""

from f_make_super_structure import make_super_structure, solve_model

excelFile = 'propionate_case_study_v1.xlsx'
superstructure = make_super_structure(excelFile= excelFile, printPyomoEq= False)

switchSolver = True
operatingDays= 1 # days of operation
saveName = 'results_V1.xlsx'  #or something like 'results.xlsx'
if switchSolver:
    #results = solve_model(superstructure,  operatingDays = operatingDays, saveName = saveName)
    results = solve_model(superstructure, operatingDays = operatingDays ,solverType='BARON', saveName=saveName)
