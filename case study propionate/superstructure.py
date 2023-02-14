# -*- coding: utf-8 -*-
"""
@author: Lucas Van der Hauwaert lucas.vanderhauwaert@usc.es
"""

from functions import make_super_structure
import time
import pyomo.environ as pe
import pyomo.opt as po

excelFile = r'\propionate_case_study.xlsx'

start_time = time.time()
superstructure = make_super_structure(excelFile= excelFile, printPyomoEq= False)
end_time = time.time()

run_time =  end_time - start_time

print('')
print('The run time is: {} seconds'.format(run_time))
print('')


switchSolver = True
if switchSolver:
    solvername = 'gams'
    opt = po.SolverFactory(solvername)
    # could also introduce extra variable in opt.solve to specify solver eg: solver='cplex'
    # =============================================================================
    #     Possible solver are: 'BARON', 'ANTIGONE', 'CPLEX', 'DICOPT'
    # =============================================================================
    #results = opt.solve(superstructure, solver='BARON', keepfiles=True, tee=True)
    results = opt.solve(superstructure, keepfiles=True, tee=True)
    for v in superstructure.component_objects(ctype=pe.Var):
        for index in v:
            if pe.value(v[index]) >= 0.1:
                print('{0} = {1}'.format(v[index], pe.value(v[index])))

    print('')
    for v2 in superstructure.component_objects(ctype=pe.Objective):
        for index2 in v2:
            a = pe.value(v2[index2])
            print('The objective value is:')
            print('{0} = {1}'.format(v2[index2], pe.value(v2[index2])))
