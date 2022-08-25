import pyomo.opt as po
import pyomo.environ as pe
import keras
import numpy as np

from pyomo.core.base import constraint
from pyomo.core.base.constraint import Constraint
from pyomo.core.expr.numvalue import value

NN_gelatine = keras.models.load_model('NN_gelatine')
model = pe.ConcreteModel()

# create variables
# binary
model.y1 = pe.Var(domain=pe.Binary)  # path NN
model.y2 = pe.Var(domain=pe.Binary)  # path equation
# mass flow variables
model.gleatine = pe.Var(bounds=(0, 10))
model.out1 = pe.Var(domain=pe.PositiveReals)
model.out2 = pe.Var(domain=pe.PositiveReals)

# define constraints
def pathCon(model): # path constraint
    return model.y1 + model.y2 == 1

def eqWay(model): # path constraint
    return model.out1 == model.gleatine*0.01 + 0.05

def NNway(model): # path constraint
    gl = model.gleatine
    inputs = np.array([gl, 4, 6])
    inputs = inputs.reshape([1, -1])
    outputs = NN_gelatine.predict(inputs)
    a = outputs[0]
    return model.out2 ==  a

model.eqConstraint = pe.Constraint(rule=eqWay)
model.NNConstraint = pe.Constraint(rule=NNway)
model.pathConstraint = pe.Constraint(rule=pathCon)

# create objective
def objRule(model):
    return model.out1 * model.y1 + model.out2 * model.y2


model.obj = pe.Objective(sense=pe.maximize, rule=objRule)

model.pprint()

# run solver
solvername = 'gams'
opt = po.SolverFactory(solvername)
# could also introduce extra variable in opt.solve to specify solver eg: solver='cplex'
results = opt.solve(model, keepfiles=True, tee=True)

model.pprint()
# After solving, variable values may be accessed either by pe.value(model.myvar) or model.myvar.value.
# print resutls

print('glucose = ', pe.value(model.glu1))
print(pe.value(model.glu2))
print('propionate', pe.value(model.prop))
print('acetate ', pe.value(model.ace))
print('PHB', pe.value(model.phb))
print('path variables')
print(pe.value(model.y1))
print(pe.value(model.y2))
print('profit', pe.value(model.obj))
