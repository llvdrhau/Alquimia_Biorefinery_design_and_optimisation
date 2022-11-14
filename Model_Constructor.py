# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 10:43:34 2022

Builds further upone the script makeModelWithClassObjects
The goal is to make a superstructure where the streams have different
compositions plus intergrating the generic process interval
@author: Lucas Van der hauwaert
"""
# imports 
import pandas as pd
import numpy as np
import pyomo.environ as pe
import pyomo.opt as po
#import warnings
#import importlib
import os
from f_make_interval_objects import make_input_intervals, make_reactor_intervals

def removeSpacesInSeries(seriesOfStr):
    newList = []
    for i in seriesOfStr:
        i = i.replace(' ', '')
        newList.append(i)
    newSeries = pd.Series(newList)
    return newSeries
def makeModelWithCompositions(excelFile, reactorLib):
   # get the input and reactor objects
   objectsInputDict = make_input_intervals(excelFile)
   objectsReactorDict = make_reactor_intervals(excelFile)
   allObjects = objectsInputDict | objectsReactorDict

   # create model instance
   model = pe.ConcreteModel()
   # retrieve data
   #loc = r'C:\Users\Reboots\OneDrive - Universidade de Santiago de Compostela\Alquimia\PYOMO\pyomo_Alquimia\excel files'
   loc = os.getcwd()
   posAlquimia = loc.find('Alquimia')
   loc = loc[0:posAlquimia + 8]
   loc = loc + r'\excel files' + excelFile

   # read excel file
   components = pd.read_excel(loc, sheet_name = 'components')

   # inputs
   inPrice = components.input_price.to_numpy()
   posInputs = inPrice != 0 
   nameInputs = components.process_intervals[posInputs] # find input variable
   nameInputs = nameInputs.tolist()
   for i in range(len(nameInputs)):
       nameInputs[i] = nameInputs[i].replace(' ','') # remove space in excel
   
   inBoundsLow = components.lower_bound[posInputs].to_numpy() 
   inBoundsUpper = components.upper_bound[posInputs].to_numpy()
   inPrice= inPrice[posInputs]  
   
   # define fixed parameters cost raw material
   inputDic = {nameInputs[i]: inPrice[i] for i in range(len(inPrice))} # make dictionary
   inBoundDict = {nameInputs[i]: [inBoundsLow[i],inBoundsUpper[i]] for i in range(len(inPrice))} # make dictionary
   model.inSet = pe.Set(initialize = nameInputs) # define set inputs
   model.Input_Price  = pe.Param(model.inSet, initialize=inputDic) #define fixed params using dictoinary: price raw material 
   #model.Input_Price = pe.Param(initialize=inputDic)
   
   # outputs 
   outPrice = components.output_price.to_numpy()
   posOutputs = outPrice != 0
   nameOutputs = components.process_intervals[posOutputs].tolist() # find input variable
   for i in range(len(nameOutputs)): 
       nameOutputs[i] = nameOutputs[i].replace(' ','') # remove space in excel  
   
   outBoundLower = components.lower_bound[posOutputs].to_numpy()
   outBoundUpper = components.upper_bound[posOutputs].to_numpy()
   outPrice= outPrice[posOutputs]  # deleet nan in product  price
   #define fixed parameters markt price products
   outputDic = {nameOutputs[i]: outPrice[i] for i in range(len(outPrice))} # make dictionary 
   outBoundDict = {nameOutputs[i]: [outBoundLower[i], outBoundUpper[i]] for i in range(len(outPrice))} # make dictionary  
   
   model.outSet = pe.Set(initialize = nameOutputs)# define set outputs 
   model.Output_Price = pe.Param(model.outSet, initialize=outputDic) #define fixed params using dictoinary: price raw material 
   #model.sellCte = pe.Param(model.outSet, initialize=outPrice)
   
   # interval variables
   reactorBoundDict = {}
   processIntervalNameList = []
   for i in range(len(components.input_price.to_numpy())):
       if not components.input_price.to_numpy()[i] and not components.output_price.to_numpy()[i]:
           processIntervalName = components.process_intervals[i]  
           processIntervalName = processIntervalName.replace(' ','') # remove spaces
           processIntervalNameList.append(processIntervalName)
           lower_input_reactor = components.lower_bound[i]
           upper_input_reactor = components.upper_bound[i]
           reactorBoundDict.update({processIntervalName: [lower_input_reactor,upper_input_reactor]})
   
   def inBound(model, i) :
       bounds = inBoundDict[i]
       lb = bounds[0]
       ub = bounds[1]
       return (lb,ub)
   def outBound(model, i) :
       bounds = outBoundDict[i]
       lb = bounds[0]
       ub = bounds[1]
       return (lb,ub)
   def intervalBound(model, i) : # that is the max mass that can flow through the reactor 
       bounds = reactorBoundDict[i]
       lb = bounds[0]
       ub = bounds[1]
       return (lb,ub)

   ### define variables from input reactor and output streams
   model.input = pe.Var(nameInputs, domain= pe.PositiveReals, bounds = inBound)
   model.intervals = pe.Var(processIntervalNameList, domain= pe.PositiveReals, bounds = intervalBound)
   model.output = pe.Var(nameOutputs, domain= pe.PositiveReals, bounds = outBound)
   
# =============================================================================
#    ### objective function relating to profit and raw input costs
# =============================================================================
   obj_expr_buy = sum(model.Input_Price[i] * model.input[i] for i in model.inSet)
   obj_expr_sell = sum(model.Output_Price[i] * model.output[i] for i in model.outSet)
   obj_expresion = obj_expr_sell - obj_expr_buy
# =============================================================================
#    DECLARE all the component variables of the system
#    rename all components to original names and save to dict {old: new versions}
#     which are stored in the interval objects
# =============================================================================
   # find reactor names
   reactorData = pd.read_excel(loc, sheet_name = 'reactors')
   reactorNames = reactorData.reactor_name
   posReactors = reactorNames != 'none'
   reactorNames = reactorNames[posReactors]
   reactorNames = list(reactorNames)

   # preallocation
   NewComponentNames = []
   processIntervalNamesWithSpaces = components.process_intervals
   processIntervalNames = []

   for i,name in enumerate(processIntervalNamesWithSpaces): # get rid of spaces
       processIntervalNames.append(name.replace(' ', ''))  # remove space in excel

   for i in processIntervalNames :
       if i in nameInputs:  # if the interval is an input and not an output nor a reactor
           try:
               #interval_to_call = getattr(importlib.import_module(reactorLib), i)
               interval_to_call = allObjects[i]
           except:  # create a warning if names in excel don't match names in reator lib.
               errorStatement = 'ERROR make sure the reactor name or interval name %s in EXCEL is the same as the library.py file: %s ' % (i, reactorLib)
               raise NameError(errorStatement)

           ComponentNames = interval_to_call.compositionNames
           updateNamesInputObject = []

           for j in ComponentNames:
               toAdd = j + '_' + i
               updateNamesInputObject.append(toAdd)
               NewComponentNames.append(toAdd)

           interval_to_call.compositionNames = updateNamesInputObject
           interval_to_call.UpdateDict(updateNamesInputObject)
           interval_to_call.makeOldNewDict(ComponentNames, updateNamesInputObject)

       elif i in reactorNames:   # if the interval has a reactor
           try: # consider alternative: change the if statement it sucks -> look if the process interval has a reactor as the if stament and not it's names, horrible design
                interval_to_call = objectsReactorDict[i]
                #interval_to_call = getattr(importlib.import_module(reactorLib), i)
           except : # create a warning if names in excel don't match names in reator lib.
               raise Exception('nope, reactor interval object {} cannot be found'.format(i))

           ComponentNames = interval_to_call.outputs 
           
           if len(ComponentNames) == 1 and ComponentNames[0] in nameOutputs:
               # product names don't change, but need to right less code if 
               # dictionary which is a copy of each other (final products don't need to be renamed)
               interval_to_call.makeOldNewDictOutputs(ComponentNames,ComponentNames) 

           else: 
               updateNamesObjectOutputs = [] 
               
               for j in ComponentNames :
                   toAdd = j + '_' + i
                   updateNamesObjectOutputs.append(toAdd)
                   NewComponentNames.append(toAdd) 
                   
               interval_to_call.outputs = updateNamesObjectOutputs  # update the object OUTPUTS not imputs VERY IMPORTANT      
               interval_to_call.makeOldNewDictOutputs(ComponentNames,updateNamesObjectOutputs)

   model.allCompositionVariables =  pe.Var(NewComponentNames, domain= pe.PositiveReals)
   #model.pprint()
   # ==================================================================================================================
   #      # define equations of the input interval (mixing equations)
   # =================================================================================================================
   def IntervalBlockInput(b,i): # where i is the list of input interval
       model 
       try:
           interval_to_call = allObjects[i]
            #interval_to_call = getattr(importlib.import_module(reactorLib), i)
       except : # create a warning if names in excel don't match names in reator lib.
           raise Exception('nope, reactor interval object {} cannot be found'.format(i))
       
       inputName = interval_to_call.inputName
       nameComponents = interval_to_call.compositionNames
       ComponentCompositions  = interval_to_call.compositionDict
       
       b.setComponents = pe.Set(initialize = nameComponents)
       #b.componentVars = pe.Var(nameComponents, domain= pe.PositiveReals)
       b.componentFractions = pe.Param(b.setComponents, initialize = ComponentCompositions, within= pe.PercentFraction)
       b.massbalanceConstraints = pe.ConstraintList()
       for j in nameComponents: # TODO make a warning if input names are wrong
           b.massbalanceConstraints.add(expr =  model.allCompositionVariables[j] == model.input[inputName] * b.componentFractions[j])

       # handel boolean inputs
       connectionMatrix = pd.read_excel(loc, sheet_name='connectionMatrix')
       reactorRow = connectionMatrix.loc[connectionMatrix['process_intervals'] == i]
       nameBoolVariables = []
       reactorsWithBool = {}
       for j in reactorRow:
           a = str(reactorRow[j]) # uncomment for debugging
           if 'bool' in str(reactorRow[j]):  # TODO make sure a minimum of 2 bools are counted, display error if not
               boolVariableName = 'y_' + j
               nameBoolVariables.append(boolVariableName)
               reactorsWithBool.update({j: boolVariableName})

       # creat the bool variables and the sum constraint
       if nameBoolVariables:
           # initialise boolean values so sum is one
           initialBools = np.zeros((len(nameBoolVariables),), dtype=int)
           initialBools[0] = 1
           initialiseBoolVars = {nameBoolVariables[i]: initialBools[i] for i in range(len(nameBoolVariables))}

           # make set variables and constraints of the bools
           b.BooleanVariables = pe.Var(nameBoolVariables, domain=pe.Boolean, initialize=initialiseBoolVars)
           exprBoolean = sum(b.BooleanVariables[k] for k in b.BooleanVariables_index) == 1
           b.BooleanConstraints = pe.Constraint(expr=exprBoolean)

       # modify the expresions of the intervall that is affected by the booliean var
       for j in reactorsWithBool:
           try:
               interval_to_modify = allObjects[j]
               #interval_to_modify = getattr(importlib.import_module(reactorLib), j)
               originBool = i  # from which reactor it comes from
               interval_to_modify.isBool = (originBool, reactorsWithBool[j])
               equationsToModify = interval_to_modify.eq
               newEquations = []
               for k in equationsToModify:
                   eqCurrent = k
                   eqCurrent = eqCurrent.replace('==', '==(')
                   eqCurrent += ") * model.IntervalBlockReactors['%s'].BooleanVariables['%s']" % (
                   i, reactorsWithBool[j])
                   newEquations.append(eqCurrent)
                   # print(eqCurrent)
               interval_to_modify.eq = newEquations

           except:  # create a warning if names in excel don't match names in reator lib.
               print(
                   'ERROR make sure the reactor name or interval name %s in EXCEL is the same as the library.py file: %s ' % (
                   j, reactorLib))

   # ===================================================================================================================
   #      # define equations of the process (Reactor) intervals (mixing equations)
   # ===================================================================================================================
   # TODO as a fail safe make the list i (reactor intervals) in the right order to make sure every thing follows each other up sequentialy    
   def IntervalBlockReactor(b,i): # where I is the list intervals with reactors WHICH MUST BE IN THE CORRECT ORDER!!! inputs to outputs so follwing the sequence as described by the layer number
       model #  to carry model object over in function
       try:
           interval_to_call = allObjects[i]
            #interval_to_call = getattr(importlib.import_module(reactorLib), i)
       except : # create a warning if names in excel don't match names in reator lib.             
           print('ERROR make sure the reactor name or interval name %s in EXCEL is the same as the library.py file: %s ' % (i,reactorLib))
           
       # call equations, inputs and outputs of the reactor
       strExpr =   interval_to_call.eq
       inputList = interval_to_call.inputs # going to be original names expect in the case of mixing 
       outputList = interval_to_call.outputs # should be a dictionary with {origial names: new names}
       connectionMatrix = pd.read_excel(loc, sheet_name='connectionMatrix')
       # =============================================================================
       #        # Total of 4 options for a process interval
       #        # 1: there is mixing
       #        # 2: there is serparation
       #        # 3: there is a boolean stream
       #        # 4: there is a split stream
       #        # 5: the reactor
       # =============================================================================

    # =============================================================================
    #             # 1: if there is mixing
    # =============================================================================           
       mix = interval_to_call.mix
       if mix: # if is not empty
           AfterMixComponentsNames = []
           intervalsToMix = mix.keys()
           for j in inputList: 
               newName = j + '_mixed'
               AfterMixComponentsNames.append(newName)
               
           # TODO find out if it is necesary to update the input objects with the mixed input names      
           interval_to_call.makeOldNewDictInputs(inputList, AfterMixComponentsNames)  
           interval_to_call.inputMix = interval_to_call.oldNewDictInputs
           # add variables to block structure 
           b.afterMixVariables = pe.Var(AfterMixComponentsNames, domain= pe.PositiveReals)
           b.mixingConstraints = pe.ConstraintList() 

           for j in inputList: 
               afterMixName = interval_to_call.inputMix[j]
               expr = 0 #preallocate expresion 
               for k in intervalsToMix:
                   try:
                       PreviousProcessIntervalObject = allObjects[k]
                       #PreviousProcessIntervalObject = getattr(importlib.import_module(reactorLib), k)
                   except :
                       raise Exception('error in reading reactor object'.format(k))

                   connectInfo = mix[k]
                   # call output dictionary 
                   outputDictPreviousInterval = PreviousProcessIntervalObject.oldNewDictOutputs

                   if hasattr(PreviousProcessIntervalObject, 'outputSplit') and j in outputDictPreviousInterval:  
                      # compleet when you start putting splits in
                      pass
                           
                   elif PreviousProcessIntervalObject.separation and j in outputDictPreviousInterval:
                       # compleet when you start putting splits in
                       pass

                   elif j in outputDictPreviousInterval: 
                       nameToAdd = outputDictPreviousInterval[j]
                       expr += model.allCompositionVariables[nameToAdd]

               b.mixingConstraints.add(b.afterMixVariables[afterMixName] == expr)  
           

# =============================================================================
#            # 2: SEPERATION
# =============================================================================
       # if a dict is given for seoeration there is obviouly a seperation element
       streamToSeperate = interval_to_call.separation

       if streamToSeperate:
           # rename outputs of the process interval 
           seperationNames = streamToSeperate.keys()
           b.seperationVariables = pe.Var(seperationNames, domain= pe.PositiveReals)

           # define constraint lists 
           b.seperationConstraints = pe.ConstraintList()
           for j in seperationNames:
               seperationDictInterval = streamToSeperate[j]
               expr = 0
               for k in seperationDictInterval:
                   percentExtraction = seperationDictInterval[k]
                   compositionVar = interval_to_call.oldNewDictOutputs[k]
                   expr += model.allCompositionVariables[compositionVar] * percentExtraction
               expresion = b.seperationVariables[j] == expr
               b.seperationConstraints.add(expresion)

    # =============================================================================
    #                    #3 BOOLEAN VARIABLES/CONSTRAINTS
    #   Name and define constraints in the i'th reactor who's output is split by 
    #   boolean variables the boolean variables   
    # =============================================================================

       # look at the row of the reactor in the connenction matrix and make new bool variables and constraints
       reactorRow = connectionMatrix.loc[connectionMatrix['process_intervals'] == i]
       nameBoolVariables = []
       reactorsWithBool = {}
       for j in reactorRow:
           # print(str(reactorRow[j])) # uncomment for debugging
           if 'bool' in str(reactorRow[j]):  # TODO make sure a minimum of 2 bools are counted, display error if not
               boolVariableName = 'y_' + j
               nameBoolVariables.append(boolVariableName)
               reactorsWithBool.update({j:boolVariableName})

       # creat the bool variables and the sum constraint
       if nameBoolVariables:
           # initialise boolean values so sum is one
           initialBools = np.zeros((len(nameBoolVariables),), dtype=int)
           initialBools[0] = 1
           initialiseBoolVars = {nameBoolVariables[i]:initialBools[i] for i in range(len(nameBoolVariables))}

           # make set variables and constraints of the bools
           b.BooleanVariables = pe.Var(nameBoolVariables, domain = pe.Boolean, initialize = initialiseBoolVars)
           exprBoolean = sum(b.BooleanVariables[k] for k in b.BooleanVariables_index) == 1
           b.BooleanConstraints = pe.Constraint(expr = exprBoolean )

       # modify the expresions of the intervall that is affected by the booliean var
       for j in reactorsWithBool:
           try:
                interval_to_modify = allObjects[j]
                #interval_to_modify = getattr(importlib.import_module(reactorLib), j)
                originBool = i # from which reactor it comes from
                interval_to_modify.isBool = (originBool,reactorsWithBool[j])
                equationsToModify = interval_to_modify.eq
                newEquations = []
                for k in equationsToModify:
                    eqCurrent = k
                    eqCurrent = eqCurrent.replace('==', '==(')
                    eqCurrent += ") * model.IntervalBlockReactors['%s'].BooleanVariables['%s']" % (i,reactorsWithBool[j])
                    newEquations.append(eqCurrent)
                    #print(eqCurrent)
                interval_to_modify.eq = newEquations

           except : # create a warning if names in excel don't match names in reator lib.
               print('ERROR make sure the reactor name or interval name %s in EXCEL is the same as the library.py file: %s ' % (j,reactorLib))
           
       # if the current reactor is dependant on a bool, add constraints to the outputs of the reactor!!
       # I think this part is not necesary, what I was trying to do is force an equation to zero if the bool var is not chosen
       # but not necsary it is already in the equation
       ##uncomment for here should I be wrong

       # isDependentOnBool = interval_to_call.isBool
       # if isDependentOnBool:
       #     b.BooleanActivationConstraints = pe.ConstraintList()
       #     for j in interval_to_call.outputs:
       #         if j in nameOutputs:
       #             body = model.output[j] - reactorBoundDict[i][1] * model.IntervalBlockReactors[isDependentOnBool[0]].BooleanVariables[isDependentOnBool[1]]
       #         else:
       #             body = model.allCompositionVariables[j] -reactorBoundDict[i][1] * model.IntervalBlockReactors[isDependentOnBool[0]].BooleanVariables[isDependentOnBool[1]]
       #
       #         lowerBound = 0
       #         upperBound = 0 # reactorBoundDict[i][1] * model.IntervalBlockReactors[isDependentOnBool[0]].BooleanVariables[isDependentOnBool[1]] #
       #         exprBoolActivation = pe.inequality(lowerBound, body ,upperBound)   #TODO compleet the expresions by adding the bool variable
       #         b.BooleanActivationConstraints.add(exprBoolActivation)
           
    # =============================================================================
    #                      #4 SPLITTER VARIABLES/CONSTRAINTS
    #   Name and define constraints in the i'th reactor who's output is split by a splitter  
    # =============================================================================          
    # find the row of the reactor and make new split variables and constraints 
       
       # split2Reactor = {}
       # for j in reactorRow :
       #      #print(j)
       #      if 'split' in str(reactorRow[j]) and 'permeate' in str(reactorRow[j]):
       #          #splitVariableName = 'split_permeate_' + j
       #          split2Reactor.update({j:'permeate'})
       #
       #      elif 'split' in str(reactorRow[j]) and 'reject' in str(reactorRow[j]):
       #          # splitVariableName = 'split_reject_' + j
       #          split2Reactor.update({j:'reject'})
       #
       #
       #      elif 'split' in str(reactorRow[j]):
       #          # splitVariableName = 'split_fraction_' + j
       #          # nameSplitVariables.append(splitVariableName)
       #          split2Reactor.update({j:'single'})
       #
       #
       # OutSplitDict = {}
       # if  split2Reactor and len(split2Reactor) == 2:  #TODO if the reject AND perm are split then find a way to get the correct constraints
       #     b.SplitConstraints = pe.ConstraintList()
       #     b.SplitFraction = pe.Var(domain = pe.PercentFraction)
       #     counter = 0
       #     for j in split2Reactor:
       #          originalVarNames = list(interval_to_call.oldNewDictOutputs.keys())
       #          originalNameAndSplitNameDict = {}
       #          nameSplitVariables = []
       #          for k in originalVarNames:
       #              newOutputSplitName = k + '_split_' + j
       #              nameSplitVariables.append(newOutputSplitName)
       #              originalNameAndSplitNameDict.update({k:newOutputSplitName})
       #
       #
       #          # make variable names
       #          if counter == 0:
       #              b.SplitVariablesUp = pe.Var(nameSplitVariables, domain = pe.PositiveReals)
       #          else:
       #              b.SplitVariablesDown = pe.Var(nameSplitVariables, domain = pe.PositiveReals)
       #
       #          for k in originalVarNames:
       #              splitVar = originalNameAndSplitNameDict[k]
       #              if split2Reactor[j] == 'permeate': # interval_to_call.seperation:
       #                  var2split = interval_to_call.outputSeperation['permeate'][k]
       #                  if counter == 0:
       #                      exprSplit = b.SplitVariablesUp[splitVar] == b.SplitFraction * b.permeateVariables[var2split]
       #                  else:
       #                      exprSplit = b.SplitVariablesDown[splitVar] == (1 - b.SplitFraction) * b.permeateVariables[var2split]
       #
       #              elif split2Reactor[j] == 'reject':
       #                  var2split =   interval_to_call.outputSeperation['reject'][k]
       #                  if counter == 0:
       #                      exprSplit = b.SplitVariablesUp[splitVar] == b.SplitFraction * b.rejectVariables[var2split]
       #                  else:
       #                      exprSplit = b.SplitVariablesDown[splitVar]  ==(1 - b.SplitFraction) * b.rejectVariables[var2split]
       #
       #                  #"model.IntervalBlockReactors['%s'].rejectVariables['%s']"
       #              else:
       #                  var2split =  interval_to_call.oldNewDictOutputs[k]
       #                  if counter == 0:
       #                      exprSplit = b.SplitVariablesUp[splitVar] == b.SplitFraction * model.allCompositionVariables[var2split]
       #                  else:
       #                      exprSplit = b.SplitVariablesDown[splitVar] == (1 - b.SplitFraction) * model.allCompositionVariables[var2split]
       #
       #              b.SplitConstraints.add(exprSplit)
       #          counter += 1
       #          # update the output dict of the reactor
       #          OutSplitDict.update({j:originalNameAndSplitNameDict})
       #
       # # update currunt and next intervalls
       # if  split2Reactor and len(split2Reactor) == 2:
       #     # update the output of the reactor
       #     interval_to_call.outputSplit = OutSplitDict
       #     counter = 0
       #     for j in split2Reactor :
       #         nextReactor = allObjects[j]
       #         #nextReactor = getattr(importlib.import_module(reactorLib), j)
       #         if counter == 0 :
       #             nextReactor.split  ='up'
       #         else:
       #             nextReactor.split  ='down'
       #         counter += 1

       model.pprint() #todo get rid is this when done debugging
    # =============================================================================
    #            # Reactor phase   
    # =============================================================================
          
       # same for if their is mixing/seperation because you
       # find new input names by looking at it's connections i.e. the previous interval
       # only interested in the previous interval if it ends in a split or seperation (anything without a 0 or 1 in other words)
       pos = connectionMatrix[i] != 0
       nameConnectedInterval = connectionMatrix.process_intervals[pos].tolist()

       if len(nameConnectedInterval) > 1 :
           pos2 = nameConnectedInterval != 1
           nameConnectedInterval = nameConnectedInterval[pos2] # should now only have one element now
           nameConnectedInterval.replace(' ','') # remove spaces
           nameConnectedInterval = nameConnectedInterval
       else:
           nameConnectedInterval = nameConnectedInterval[0]

       pos = connectionMatrix.process_intervals == nameConnectedInterval
       connectInfo = str(connectionMatrix[i][pos]) #.tolist()
       try:
           #PreviousProcessIntervalObject = getattr(importlib.import_module(reactorLib), nameConnectedInterval)
           PreviousProcessIntervalObject = allObjects[nameConnectedInterval]
       except :
           raise Exception('mmmmmm')
       #find the inputs and outputs of the reactor
       # todo get rid of this if statement make sure the reator object oldNew dictionary is always reffering to what goes in or out
       if interval_to_call.mix :
           inputsReactor = interval_to_call.inputMix
           outputsReactor = interval_to_call.oldNewDictOutputs

       elif nameConnectedInterval not in nameInputs and 'split' in connectInfo:
           pass # get correct when spliting is compleet

       elif nameConnectedInterval not in nameInputs and 'sep' in connectInfo:
           pass # get correct when seperation is compleet

       else: # just simply connected without split/mix or seperation before the interval
           inputsReactor = PreviousProcessIntervalObject.oldNewDictOutputs
           outputsReactor = interval_to_call.oldNewDictOutputs

       # preallocate the constraints and massbalance equation
       b.reactorConstraints = pe.ConstraintList()
       massBalanceString = "model.intervals['%s'] == " % i
       # ajust the reactor expressions
       intervalsToMix = interval_to_call.mix
       for j in strExpr:
           strExprCurrent = j
           for k in inputsReactor:
                 if intervalsToMix:
                     replaceStr = "model.IntervalBlockReactors['%s'].afterMixVariables['%s']" % (i,inputsReactor[k]) # from the current reactor interval i
                     if replaceStr not in massBalanceString: # so you don't get duplicates in mass balance equation
                         massBalanceString += ' + ' + replaceStr

                 elif 'split' in  connectInfo:
                     pass

                 elif 'sep' in connectInfo:
                     pass
                     # replaceStr = "model.IntervalBlockReactors['%s'].permeateVariables['%s']" % (nameConnectedInterval,inputsReactor[k])
                     # if replaceStr not in massBalanceString: # avoid duplicates in expresion
                     #     massBalanceString += ' + ' + replaceStr
                 else:
                     replaceStr = "model.allCompositionVariables['%s']" % inputsReactor[k]
                     if replaceStr not in massBalanceString:
                         massBalanceString += ' + ' + replaceStr

                 strExprCurrent = strExprCurrent.replace(k,replaceStr)


           for k in outputsReactor:
                if k in nameOutputs:
                    replaceStr = "model.output['%s']" % outputsReactor[k]
                else:
                    replaceStr = "model.allCompositionVariables['%s']" %outputsReactor[k]
                strExprCurrent = strExprCurrent.replace(k,replaceStr)
           #model.pprint()
           print(strExprCurrent)
           expresion = eval(strExprCurrent)
           b.reactorConstraints.add(expresion)

       #print(massBalanceString)
       massBalanceExpresion = eval(massBalanceString)
       b.MassBalanceReactorConstraint = pe.Constraint(expr= massBalanceExpresion)

   # ==================================================================================================================
   #      # define equations of the output Intervals (mixing equations)
   #        i.e., the sum of all the intervalls that produce a certain end product 
   # =================================================================================================================
   # TODO check if the output block works
   def IntervalBlockOutput(b,i):
       model  # to carry model object over in function
       b.outputConstraints = pe.ConstraintList()
       connectionMatrix = pd.read_excel(loc, sheet_name = 'connectionMatrix')
       DFintervals = pd.read_excel(loc, sheet_name = 'components')
       # find connected intervals
       positionOfconectedIntervals = connectionMatrix[i].astype(dtype=bool)
       allIntervals = connectionMatrix['process_intervals']
       allIntervals = removeSpacesInSeries(allIntervals)
       connectedIntervalNames = allIntervals[positionOfconectedIntervals]
       equation = "model.output['{}']==".format(i)
       for j in connectedIntervalNames:
           try:
               interval_to_call = allObjects[j]
               #interval_to_call = getattr(importlib.import_module(reactorLib), j)
           except:  # create a warning if names in excel don't match names in reator lib.
               print('ERROR make sure the reactor name or interval name %s in EXCEL is the same as the library.py file: %s ' % (i, reactorLib))
           outputsInterval = interval_to_call.outputs
           # find the short name of the output in components sheet of the excel file
           posOutputInterval = allIntervals == i
           shortNameoutput = list(DFintervals.components[posOutputInterval])
           shortNameoutput = shortNameoutput[0] # there should only be ne element in the final output interval
           for k in outputsInterval :
               if shortNameoutput in k:
                   equation += "+ model.allCompositionVariables['{}']".format(k) # TODO check: carefull!!! check if all variables go to allCompositionVariables?
       outputMassBalanceExpresion = eval(equation)
       b.outputConstraints.add(expr= outputMassBalanceExpresion) #= pe.Constraint(expr= outputMassBalanceExpresion)
            # interesting warning when debuggung! block.del_component() and block.add_component().
       
   # Decalre input, reactor and output blocks to the model structure
   ####################################################################################################################
   model.IntervalBlockInput = pe.Block(nameInputs, rule=IntervalBlockInput)
   model.IntervalBlockReactors = pe.Block(reactorNames, rule=IntervalBlockReactor)
   model.IntervalBlockOutput = pe.Block(nameOutputs, rule=IntervalBlockOutput)
   model.objective = po.Objective(sense=po.maximize, expr=obj_expresion)
   ####################################################################################################################
   # model.pprint() # debug check
   return model 

# test to see if the varable/parameters 
if __name__ == '__main__' :
    
    model = makeModelWithCompositions(r'\testModelCompositions.xlsx', "testLibraryCompositions")
    #model = makeModelWithCompositions('testModelCompositionsOnlySeparation.xlsx',"testLibraryCompositions")
    model.pprint()
    
    switchSolver = False 
    if switchSolver: 
    
        solvername = 'gams'
        opt = po.SolverFactory(solvername)
        # could also introduce extra variable in opt.solve to specify solver eg: solver='cplex'
        
    # =============================================================================
    #     Possible solver are: 'BARON', 'ANTIGONE', 'CPLEX', 'DICOPT'
    # =============================================================================
        #results = opt.solve(model, solver='BARON', keepfiles=True, tee=True)
        results = opt.solve(model, keepfiles=True, tee=True)
        
        model.pprint()
        
        for v in model.component_objects(ctype=pe.Var):
            for index in v:
                print('{0} = {1}'.format(v[index], pe.value(v[index])))
    

