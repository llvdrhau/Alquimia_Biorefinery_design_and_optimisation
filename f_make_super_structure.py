# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 10:43:34 2022

Builds further upone the script makeModelWithClassObjects
The goal is to make a superstructure where the streams have different
compositions plus intergrating the generic process interval
@author: Lucas Van der Hauwaert
"""
# imports
import pandas as pd
import numpy as np
import pyomo.environ as pe
import pyomo.opt as po
import importlib
import os
from f_makeIntervalObjects import make_reactor_intervals, make_input_intervals, update_reactor_intervals

def make_super_structure(excelFile):
    superStructure = pe.ConcreteModel()

    objectsInputDict = make_input_intervals(excelFile)
    objectsReactorDict = make_reactor_intervals(excelFile)

    allObjects = objectsInputDict | objectsReactorDict
    update_reactor_intervals(allObjects, excelFile)

    print(5)
    """
    Declare all interval variables (capital letters) and component variables (small letters)
    loop over all objects 
    declare all variables
    delare equations of each object 
    make the objective 
    RUN THE SOLVER 
    """
    # for i in objectsInputDict:
    #     try:
    #         interval_to_call = allObjects[i]
    #     except:  # create a warning if names in Excel don't match names in reator lib.
    #         raise Exception('nope, reactor interval object {} cannot be found'.format(i))
    #
    #     # head varible constrained by bounds
    #     intervalName = interval_to_call.inputName
    #     superStructure.intervalVar = pe.Var(intervalName, bounds=interval_to_call.boundryInput)
    #
    #     # free float variables
    #     compositionVarNames = interval_to_call.variables
    #     superStructure.variables = pe.Var(compositionVarNames, domain=pe.PositiveReals)
    #
    #     # introduce equations
    #     superStructure.InputConstraints = pe.ConstraintList()
    #     componentEquations = interval_to_call.componentEquations
    #     for eq in componentEquations:
    #         eq = eq.replace(intervalName, "superStructure.intervalVar['{}']".format(intervalName))
    #         for j in compositionVarNames:
    #             eq = eq.replace(j,"superStructure.variables['{}']".format(j))
    #         # print(eq) # for debugging
    #         constraintExpresion = eval(eq)
    #         superStructure.InputConstraints.add = pe.Constraint(expr=constraintExpresion)

    # Declare input, reactor and output blocks to the model structure
    ####################################################################################################################
    #nameInputs = list(objectsInputDict.keys())
    #superStructure.IntervalBlockInput = pe.Block(nameInputs, rule=input_block)
    # superStructure.IntervalBlockReactors = pe.Block(reactorNames, rule=IntervalBlockReactor)
    # superStructure.IntervalBlockOutput = pe.Block(nameOutputs, rule=IntervalBlockOutput)
    # superStructure.objective = po.Objective(sense=po.maximize, expr=obj_expresion)
    ####################################################################################################################
    # model.pprint() # debug check

    return superStructure