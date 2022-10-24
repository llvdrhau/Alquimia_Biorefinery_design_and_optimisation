# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import pyomo.environ as pyo
import pyomo.opt as po

sc = 3

model = pyo.ConcreteModel()
# input
model.glu = pyo.Var(bounds=(0, 15))
# output
model.prop = pyo.Var()
model.prop2 = pyo.Var()
# boolean vars
model.y1 = pyo.Var(domain=pyo.Boolean)
model.y2 = pyo.Var(domain=pyo.Boolean)

model.ConstraintsR = pyo.ConstraintList()
model.ConstraintsR.add(expr=model.prop == 0.5 * model.glu * model.y1)
model.ConstraintsR.add(expr=model.prop2 == 0.25 * model.glu * model.y2)
model.ConstraintsR.add(expr=model.y1 + model.y2 == 1)

objective = model.prop + model.prop2  # * 2 + model.p2 * 3
model.obj = pyo.Objective(sense=pyo.maximize, expr=objective)

solvername = 'gams'
opt = po.SolverFactory(solvername)
# could also introduce extra variable in opt.solve to specify solver eg: solver='cplex'

# =============================================================================
#     Possible solver are: 'BARON', 'ANTIGONE', 'CPLEX', 'DICOPT'
# =============================================================================
results = opt.solve(model, solver='BARON', keepfiles=True, tee=True)
# results = opt.solve(model, keepfiles=True, tee=True)

model.pprint()

for v in model.component_objects(ctype=pyo.Var):
    for index in v:
        print('{0} = {1}'.format(v[index], pyo.value(v[index])))

