# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 10:16:49 2022
this code contains the function to create the pyomo model specifide 
by an excel file 
@author: lucas van der hauwaert
"""

def makeModel(excelFile):
   # read excel file 
   import pandas as pd
   #import numpy as np
   import pyomo.environ as pe

   # create model instance
   model = pe.ConcreteModel()
   
   # retrive data 
   inputs = pd.read_excel(excelFile, sheet_name = 'inputs')  
   outputs = pd.read_excel(excelFile, sheet_name = 'outputs') 
   unitProcesses = pd.read_excel(excelFile, sheet_name = 'reactors') 
   
   # manipulate input 
   nIn = inputs.inputs.tolist()
   inputCost = inputs.price.to_numpy()
   #keyCostInput = list(range(1,len(inputCost)+1)) 
   #create dictionary so set can be made
   inputDic = {nIn[i]: inputCost[i] for i in range(len(inputCost))} 
   
   # manipulate output 
   nOut = outputs.outputs.tolist()
   outputCost = outputs.price.to_numpy()
   #keyCostOutput = list(range(1,len(outputCost)+1))
   #create dictionary so set can be made
   outputDic = {nOut[i]: outputCost[i] for i in range(len(outputCost))}
   
   # define sets (parameters that have fixed value eg: the cost)
   model.inSet = pe.Set(initialize = nIn)
   model.outSet = pe.Set(initialize = nOut)
   
   #define parameters 
   model.buyCte = pe.Param(model.inSet, initialize=inputDic)
   model.sellCte = pe.Param(model.outSet, initialize=outputDic)
   
   # define variables from input and output streams
   model.input = pe.Var(nIn, domain= pe.PositiveReals)
   model.output = pe.Var(nOut, domain= pe.PositiveReals)
   
   # objective function relating to profit and raw input costs  
   obj_expr_buy = sum(model.buyCte[i] * model.input[i] for i in model.inSet)
   obj_expr_sell = sum(model.sellCte[i] * model.output[i] for i in model.outSet)
   model.obj = pe.Objective(sense=pe.maximize, expr=obj_expr_sell- obj_expr_buy)






   return model 

# test to see if the varable/parameters 
if __name__ == '__main__' :
    model = makeModel('testModel.xlsx')
    model.pprint()

    
    







