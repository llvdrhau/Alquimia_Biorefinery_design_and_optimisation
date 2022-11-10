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
from f_makeIntervalObjects import make_reactor_intervals, make_input_intervals, make_output_intervals, update_intervals, check_excel_file, get_vars_eqs_bounds, make_pyomo_equations

def make_super_structure(excelFile):
    model = pe.ConcreteModel()
    check_excel_file(excelName= excelFile)

    objectsInputDict = make_input_intervals(excelFile)
    objectsReactorDict = make_reactor_intervals(excelFile)
    objectsOutputDict  = make_output_intervals(excelFile)

    allObjects = objectsInputDict | objectsReactorDict | objectsOutputDict
    update_intervals(allObjects, excelFile)
    variables, equations, bounds = get_vars_eqs_bounds(allObjects)


    """
    Declare all interval variables (capital letters) and component variables (small letters)
    loop over all objects 
    declare all variables
    delare equations of each object 
    make the objective 
    RUN THE SOLVER 
    """
    def boundsRule(model,i):
        boudVar = bounds[i]
        lowerBound = 0
        upperBound = None
        if isinstance(boudVar,list):
            lowerBound = boudVar[0]
            upperBound = boudVar[1]
        elif isinstance(boudVar,str): # and 'bool' in boudVar or 'positiveReal'
            lowerBound = 0
            upperBound = None
        return (lowerBound,upperBound)

    model.var = pe.Var(variables['continuous'], domain=pe.PositiveReals, bounds=boundsRule)
    model.boolVar = pe.Var(variables['boolean'], domain=pe.Boolean)
    #model.fractionVar = pe.Var(variables['fractions'], domain=pe.PercentFraction)


    # introduce the equations to pyomo
    model.constraints = pe.ConstraintList()
    for eq in equations:
        print(eq)
        expresion = eval(eq)
        model.constraints.add(expresion)

    model.pprint()
    print(5)

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

    return model