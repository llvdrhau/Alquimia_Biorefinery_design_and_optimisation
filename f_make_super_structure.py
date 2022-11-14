# -*- coding: utf-8 -*-
"""
Created on Mon 29 nov  10:43:34 2022

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
        lowerBound = 0 # default
        upperBound = None # default
        if isinstance(boudVar,list):
            lowerBound = boudVar[0]
            upperBound = boudVar[1]
        elif isinstance(boudVar,str): #  'bool' or 'positiveReal' in boudVar
            lowerBound = 0
            upperBound = None
        return (lowerBound,upperBound)

    model.var = pe.Var(variables['continuous'], domain=pe.PositiveReals, bounds=boundsRule)
    if variables['boolean']:
        model.boolVar = pe.Var(variables['boolean'], domain=pe.Boolean)
    if variables['fraction']:
        model.fractionVar = pe.Var(variables['fraction'], domain=pe.PercentFraction)

    # introduce the equations to pyomo
    model.constraints = pe.ConstraintList()
    for eq in equations:
        print(eq)
        expresion = eval(eq)
        model.constraints.add(expresion)

    model.pprint() # debug check

    return model