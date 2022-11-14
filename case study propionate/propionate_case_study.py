# -*- coding: utf-8 -*-
"""
@author: Lucas Van der Hauwaert lucas.vanderhauwaert@usc.es
"""
#from Model_Constructor import makeModelWithCompositions
from f_make_interval_objects import check_excel_file
from f_make_super_structure import make_super_structure
import pyomo.environ as pe
import pyomo.opt as po

# todo 1) make correct surrogate for the model of alberte in grams of C
#excelFile = r'\data_propionibacteria.xlsx'
#excelFile = r'\play.xlsx' #just for testing
excelFile = r'\propionate_production.xlsx'
superstructure = make_super_structure(excelFile= excelFile)

#model = makeModelWithCompositions(r'\data_propionibacteria.xlsx', "library_propionate")
#model.pprint()

switchSolver = False
if switchSolver:
    solvername = 'gams'
    opt = po.SolverFactory(solvername)
    # could also introduce extra variable in opt.solve to specify solver eg: solver='cplex'

    # =============================================================================
    #     Possible solver are: 'BARON', 'ANTIGONE', 'CPLEX', 'DICOPT'
    # =============================================================================
    # results = opt.solve(model, solver='BARON', keepfiles=True, tee=True)
    results = opt.solve(superstructure, keepfiles=True, tee=True)

    superstructure.pprint()

    for v in superstructure.component_objects(ctype=pe.Var):
        for index in v:
            print('{0} = {1}'.format(v[index], pe.value(v[index])))