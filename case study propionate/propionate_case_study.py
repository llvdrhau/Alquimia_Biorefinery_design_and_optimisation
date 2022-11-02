# -*- coding: utf-8 -*-
"""
@author: Lucas Van der Hauwaert lucas.vanderhauwaert@usc.es
"""
from Model_Constructor import makeModelWithCompositions
from f_makeIntervalObjects import check_excel_file
import pyomo.environ as pe
import pyomo.opt as po

# todo 1) make correct surrogate for the model of alberte in grams of C
excelFile = r'\data_propionibacteria.xlsx'
check_excel_file(excelFile)
model = makeModelWithCompositions(r'\data_propionibacteria.xlsx', "library_propionate")
model.pprint()

switchSolver = False
if switchSolver:

    solvername = 'gams'
    opt = po.SolverFactory(solvername)
    # could also introduce extra variable in opt.solve to specify solver eg: solver='cplex'

    # =============================================================================
    #     Possible solver are: 'BARON', 'ANTIGONE', 'CPLEX', 'DICOPT'
    # =============================================================================
    # results = opt.solve(model, solver='BARON', keepfiles=True, tee=True)
    results = opt.solve(model, keepfiles=True, tee=True)

    model.pprint()

    for v in model.component_objects(ctype=pe.Var):
        for index in v:
            print('{0} = {1}'.format(v[index], pe.value(v[index])))