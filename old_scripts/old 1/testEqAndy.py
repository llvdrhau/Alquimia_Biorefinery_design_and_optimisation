# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 10:16:00 2022

@author: 
"""


import pyomo.environ as pe
import pyomo.opt as po



model = pe.ConcreteModel()

# create bool vars 
model.y0 = pe.Var(domain=pe.Binary, initialize = 0) #
model.y1 = pe.Var(domain=pe.Binary, initialize = 0) # # 
model.y2 = pe.Var(domain=pe.Binary, initialize = 0) # # 
model.y3 = pe.Var(domain=pe.Binary, initialize = 0) #) #

# inputs 
model.in1 = pe.Var(domain= pe.PositiveReals, bounds=(0,101), initialize = 100)

# streams 
model.strm1 = pe.Var(domain= pe.PositiveReals)#, bounds=(0,300))
model.strm4 = pe.Var(domain= pe.PositiveReals) #, bounds=(0,300))    
        
# products 
model.px = pe.Var()#(bounds=(0,296))
model.py = pe.Var()#(bounds=(0,21))
model.p5 = pe.Var()#(domain = pe.PositiveReals)

#################################################################################################################  EQUATIONS 
#create constraints aka model library equartions

model.c1 = pe.Constraint(expr = model.strm1 == model.in1*0.8 )
model.c2 = pe.Constraint(expr = model.strm4 == model.strm1*0.9*model.y0 )
model.c3 = pe.Constraint(expr = model.p5 == 0.8*model.strm1*model.y1 )
model.c4 = pe.Constraint(expr = model.px == 0.1*model.strm4*model.y2 )
model.c5 = pe.Constraint(expr = model.py == 0.6*model.strm4*model.y3 )

# =============================================================================
# model.c1 = pe.Constraint(expr = model.strm1 == model.in1*0.8 )
# model.c2 = pe.Constraint(expr = model.strm4*model.y0 == model.strm1*0.9 )
# model.c3 = pe.Constraint(expr = model.p5*model.y1 == 0.8*model.strm1 )
# model.c4 = pe.Constraint(expr = model.px*model.y2 == 0.1*model.strm4 )
# model.c5 = pe.Constraint(expr = model.py*model.y3 == 0.6*model.strm4 )
# =============================================================================

# =============================================================================
# model.c1 = pe.Constraint(expr = model.strm1 == model.in1*0.8 )
# model.c2 = pe.Constraint(expr = model.strm4 == model.strm1*0.1 )
# model.c3 = pe.Constraint(expr = model.p5 == 0.99*model.strm1 )
# model.c4 = pe.Constraint(expr = model.px == 0.8*model.strm4 )
# model.c5 = pe.Constraint(expr = model.py == 0.2*model.strm4 )
# =============================================================================

################################################################################# boundry EQUATIONS 

# boundry constraints 
# =============================================================================

# 
model.b1a = pe.Constraint(expr = model.y2*0 <=  model.px )  
model.b1b = pe.Constraint(expr = model.px <= model.y2*100 )
#  
model.b2a = pe.Constraint(expr = model.y3*0 <=  model.py)
model.b2b = pe.Constraint(expr =  model.py<= model.y3*100)

model.b3a = pe.Constraint(expr = model.y1*0 <=  model.p5 )
model.b3b = pe.Constraint(expr =  model.p5 <= model.y1*100)

model.bstrm4 = pe.Constraint(expr = model.y0*0 <=  model.strm4 )
model.bstrm4b = pe.Constraint(expr = model.strm4 <= model.y0*100)
# 
# =============================================================================


#path constraints 
model.path1 = pe.Constraint(expr = model.y1 + model.y0 == 1)
model.path2 = pe.Constraint(expr = model.y3 + model.y2 == 1)

# create objective 
def objRule(model): 
    a =  4.0*model.p5*model.y1 + 2.8*model.px*model.y2 + 2.6*model.py*model.y3 - 2.0*model.in1
    b = 2.70*model.p5 + 2.8*model.px + 2.6*model.py - 0.05*model.in1
    return  b

model.obj = pe.Objective(sense=pe.maximize, rule=objRule)


# run solver 
solvername = 'gams'
opt = po.SolverFactory(solvername)
# could also introduce extra variable in opt.solve to specify solver eg: solver='cplex'
#results = opt.solve(model, solver = 'ANTIGONE',keepfiles=True, tee=True)
results = opt.solve(model,keepfiles=True, tee=True)


model.pprint()
#After solving, variable values may be accessed either by pe.value(model.myvar) or model.myvar.value.
# print resutls 

print('Y0 = ', pe.value(model.y0))    
print('Y1 = ', pe.value(model.y1)) 
print('Y2 = ', pe.value(model.y2)) 
print('Y3 = ', pe.value(model.y3)) 
print('IN = ', pe.value(model.in1)) 
print('STRM1 = ', pe.value(model.strm1)) 
print('STRM4 = ', pe.value(model.strm4)) 
print('px = ', pe.value(model.px)) 
print('py = ', pe.value(model.py)) 
print('p5 = ', pe.value(model.p5)) 


