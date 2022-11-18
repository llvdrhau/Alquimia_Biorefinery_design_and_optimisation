# -*- coding: utf-8 -*-
"""
@author: Lucas Van der Hauwaert lucas.vanderhauwaert@usc.es
"""

from f_make_super_structure import make_super_structure
import pyomo.environ as pe
import pyomo.opt as po

# todo 1) make correct surrogate for the model of alberte in grams of C

excelFile = r'\propionate_production_pure_cultures.xlsx'
superstructure = make_super_structure(excelFile= excelFile)

switchSolver = True
if switchSolver:
    solvername = 'gams'
    opt = po.SolverFactory(solvername)
    # could also introduce extra variable in opt.solve to specify solver eg: solver='cplex'
    # =============================================================================
    #     Possible solver are: 'BARON', 'ANTIGONE', 'CPLEX', 'DICOPT'
    # =============================================================================
    results = opt.solve(superstructure, solver='BARON', keepfiles=True, tee=True)
    #results = opt.solve(superstructure, keepfiles=True, tee=True)
    for v in superstructure.component_objects(ctype=pe.Var):
        for index in v:
            if pe.value(v[index]) > 0:
                print('{0} = {1}'.format(v[index], pe.value(v[index])))