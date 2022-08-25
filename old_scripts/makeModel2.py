# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 15:42:07 2022

@author: lucas van der hauwaert
"""

# imports 
import pandas as pd
import numpy as np
import pyomo.environ as pe
import pyomo.opt as po
import warnings


def makeModel(excelFile, reactorLib):
   # create model instance
   model = pe.ConcreteModel()
   
   # retrive data 
   components = pd.read_excel(excelFile, sheet_name = 'components')  
   
   # inputs 
   inPrice = components.input_price.to_numpy()
   allInPrice = inPrice
   nIn = components.components[~np.isnan(inPrice)]# find input variable
  
   nIn = nIn.tolist()
   for i in range(len(nIn)): 
       nIn[i] = nIn[i].replace(' ','') # remove space in excel 
   
   inBoundsLow = components.lower_bound[~np.isnan(inPrice)].to_numpy() 
   inBoundsUpper = components.upper_bound[~np.isnan(inPrice)].to_numpy()
   inPrice= inPrice[~np.isnan(inPrice)]  # deleet nan in raw material price
   
   # define fixed parameters cost raw material
   inputDic = {nIn[i]: inPrice[i] for i in range(len(inPrice))} # make dictionary 
   inBoundDict = {nIn[i]: [inBoundsLow[i],inBoundsUpper[i]] for i in range(len(inPrice))} # make dictionary 
   model.inSet = pe.Set(initialize = nIn) # define set inputs 
   model.buyCte = pe.Param(model.inSet, initialize=inputDic) #define fixed params using dictoinary: price raw material 
   
   # outputs 
   outPrice = components.output_price.to_numpy()
   allOutPrice = outPrice
   nOut = components.components[~np.isnan(outPrice)].tolist() # find input variable
   for i in range(len(nOut)): 
       nOut[i] = nOut[i].replace(' ','') # remove space in excel  
   
   outBoundLower = components.lower_bound[~np.isnan(outPrice)].to_numpy()
   outBoundUpper = components.upper_bound[~np.isnan(outPrice)].to_numpy()
   outPrice= outPrice[~np.isnan(outPrice)]  # deleet nan in product  price
   #define fixed parameters markt price products
   outputDic = {nOut[i]: outPrice[i] for i in range(len(outPrice))} # make dictionary 
   outBoundDict = {nOut[i]: [outBoundLower[i], outBoundUpper[i]] for i in range(len(outPrice))} # make dictionary  
   
   model.outSet = pe.Set(initialize = nOut)# define set outputs 
   model.sellCte = pe.Param(model.outSet, initialize=outputDic) #define fixed params using dictoinary: price raw material 
   #model.sellCte = pe.Param(model.outSet, initialize=outPrice)
   
   # Stream variables
   streamDict = {}

   for i in range(len(allOutPrice)):
       if np.isnan(allOutPrice[i]) and np.isnan(allInPrice[i]):
           streamVariable = components.components[i]  
           streamVariable = streamVariable.replace(' ','') # remove spaces
           lower_stream = components.lower_bound[i]
           upper_stream = components.upper_bound[i]
           streamDict.update({streamVariable: [lower_stream,upper_stream]})
           
   
   def streamBound(model, i):
       bounds = streamDict[i]
       lb = bounds[0]
       ub = bounds[1]
       return (lb, ub)
   
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
   
   
     
   ### define variables from input and output streams
   model.streamVariables = pe.Var(streamDict, domain= pe.PositiveReals, bounds = streamBound)
   model.input = pe.Var(nIn, domain= pe.PositiveReals, bounds = inBound)
   model.output = pe.Var(nOut, domain= pe.PositiveReals, bounds = outBound)
   
# =============================================================================
#    ### objective function relating to profit and raw input costs  
# =============================================================================

# =============================================================================
# # TODO if primary inputs have boolean variables you have to define extra constratains 
# # ex m.IN == m.in_a * Y1 + m.in_b* Y2
# # with a for loop check if input var is boolean and add corret expersion eg
# expr = 0 
# for i in nIn :
#     if nIn in dictionaryBool :
#         expr += model.buyCte[i]* model.inpi[i]*model.bool[i]
#     else: 
#         expr += model.buyCte[i]* model.inpi[i]
# =============================================================================

   obj_expr_buy = sum(model.buyCte[i] * model.input[i] for i in model.inSet)
   obj_expr_sell = sum(model.sellCte[i] * model.output[i] for i in model.outSet)

   
# =============================================================================
#    ### boolean block function connenction/constraints == 2
# =============================================================================

   ### connections Bolean variables
   connectionMatrix = pd.read_excel(excelFile, sheet_name = 'connectionMatrix',index_col= 'Reactor') 
   # find streams with boolean variables 
   boolDic = {}
   boolCloneDict = {}
   counter = 0      
   BoleanVarAll = []
   cloneName = []
   
   for i in connectionMatrix: 
       a = connectionMatrix[i] 
       locateBolean =  a == 2
       nBolean = sum(locateBolean)  

       if nBolean > 1:
           
           BoleanVar = []
           cloneName = []

           startBolean = counter 
           counter += nBolean
           stopBolean = counter
           nBo = np.arange(startBolean,stopBolean) 

           for y in range(nBolean):
               cloneName.append('%s_%s' %(i,y))
           new_key_values_dict = {i: cloneName}
           boolCloneDict.update(new_key_values_dict)
           
           for x in nBo:    
               BoleanVar.append('Y%s' %x) 
               BoleanVarAll.append('Y%s' %x) 
           new_key_values_dict = {i: BoleanVar}
           boolDic.update(new_key_values_dict)
   
    
    
   def booleanBlock(b,itera):   # I think their is a way not to put this in a block but rather put everything in the main model     
       
   # make a mini dict to initilise the boolean variables with ones and zeros
       BoleanVar = boolDic.get(itera)    
       initialBools = np.zeros((len(BoleanVar),), dtype=int)
       initialBools[0] = 1
       miniDict = {BoleanVar[i]:initialBools[i] for i in range(len(BoleanVar))}

       
       b.boleanSet = pe.Set(initialize = BoleanVar) # define set inputs
       b.BoleanVariables = pe.Var(BoleanVar, domain = pe.Boolean, initialize = miniDict) 
       #b.BoleanVariables = pe.Var(BoleanVar,domain_type=pe.IntegerSet, lb=0, ub=1)
       boleanConstraint = sum(b.BoleanVariables[i] for i in b.boleanSet) == 1
       b.boleanConstraints = pe.Constraint(rule= boleanConstraint)
       
       # add second constraint?? 
# =============================================================================
#        boleanConstraint = np.prod(b.BoleanVariables[i] for i in b.boleanSet) == 1
#        b.boolList = pe.ConstraintList()
#        b.boolList.add(rule =boleanConstraint )
# =============================================================================
       
   def booleanConstrains(b,itera):   # I think their is a way not to put this in a block but rather put everything in the main model     
       BoleanVar = boolDic.get(itera)    
       b.boleanSet = pe.Set(initialize = BoleanVar) # define set inputs
       b.BoleanVariables = pe.Var(BoleanVar, domain = pe.Binary) 
       boleanConstraint = sum(b.BoleanVariables[i] for i in b.boleanSet) == 1
       
       return boleanConstraint
            
       
   # add the block to the main model  
   
   model.BoolBlock = pe.Block(boolDic,rule=booleanBlock)
   #model.Boolconstraint = pe.Constraint(boolDic,rule=booleanConstrains) 

# =============================================================================
#  block for massbalance connenction/constraints  == 3
# =============================================================================
   def massbalanceBlock(b,itera):        
       mbVar = mbDic.get(itera)    
       b.mbSet = pe.Set(initialize = mbVar) # define set inputs

       initialMb = np.ones((len(mbVar),), dtype=int)/len(mbVar)
       miniDict = {mbVar[i]:initialMb[i] for i in range(len(mbVar))}

       b.mbVariables = pe.Var(mbVar, domain = pe.PercentFraction, initialize = miniDict) 
       massaBalanceConstraint = sum(b.mbVariables[i] for i in b.mbSet) == 1
       b.massaBalanceConstraint = pe.Constraint(rule= massaBalanceConstraint)

   # find streams with massabalance variables 
   mbDic = {}
   mbCloneDict = {}
   counter = 0 
   mbVarAll = []
   
   for i in connectionMatrix: 
       a = connectionMatrix[i] 
       locateMb =  a == 3
       nMb = sum(locateMb)  
       if nMb > 1:
           mbVar = []
           cloneName = []
           startmb = counter 
           counter += nMb
           stopmb = counter
           nBo = np.arange(startmb,stopmb) 
           for y in range(nMb): 
               cloneName.append('%s_%s' %(i,y))
           new_key_values_dict = {i: cloneName}
           mbCloneDict.update(new_key_values_dict)
           
           for x in nBo:    
               mbVar.append('mb%s' %x) 
               mbVarAll.append('mb%s' %x)
           new_key_values_dict = {i: mbVar}
           mbDic.update(new_key_values_dict)
    
     # add the block to the main model   
   model.mbBlock = pe.Block(mbDic,rule=massbalanceBlock)
   
# =============================================================================
#    reactor dictionary 
# =============================================================================

########### needs to accept multiple expresions also carfull with input output which are both boolean!!!! 
   reactorDict = {}
   # reator data 
   reactorData = pd.read_excel(excelFile, sheet_name = 'reactors')   
   reactorNames = reactorData.reactors
   
   counterDict = {}
   namesCount = list(boolDic.keys()) + list(mbDic.keys())
   
   for i in namesCount:
        new_key = {i:0}
        counterDict.update(new_key)
           
   for i in reactorNames:  # loop over the reactor functions and compleet the dict
      import testLibrary # TODO change to official Library!!!
      try:
          reactor_to_call = getattr(testLibrary, i)
          reactorInfo = reactor_to_call() # first call is in_outdict second is reactor expresion 
      except NameError: # create a warning if names in excel don't match names in reator lib.
          print('ERROR make sure the reactor names in EXCEL are the same as the library')
      
      inOutDict = reactorInfo[0]
      strExpr =   reactorInfo[1]
      
      if isinstance(strExpr,str): # fail safe if you forget to code the string reactor expression as a list
          strExpr_help = [strExpr]
          strExpr = strExpr_help
          
      InputList = inOutDict['inputs']   #○ TODO make a warning to check if library names and excel names match
      
      CloneListBool = list(boolCloneDict.keys()) 
      commonSetBool = list(set(CloneListBool)&set(InputList))
      
      CloneListMb = list(mbCloneDict.keys())
      commonSetMb = list(set(CloneListMb)&set(InputList))
      
      ReactorExpression = []
      variablesList = inOutDict['inputs'] + inOutDict['outputs'] 

      for x in strExpr:
          strExpr2add = x
          if commonSetBool:
             for j in commonSetBool:
                 countKey = counterDict[j]
                 y = boolDic[j][countKey]
                 #counterDict[j]  +=  1 
                 strExpr2add = strExpr2add.replace('==', '==(')
                 strExpr2add += ')*%s' % y
                 
                 if y not in variablesList:
                     variablesList.append(y)
                 
                 #ReactorExpression.append(strExpr2add)
                 #reactorDict.update({i:[variablesList, ReactorExpression]})
                 
          if commonSetMb:
              for j in commonSetMb:
                  countKey = counterDict[j]
                  mb = mbDic[j][countKey]
                  #counterDict[j]  +=  1 
                  strExpr2add = strExpr2add.replace(j, '%s*%s' % (mb,j))
                  
                  if mb not in variablesList: # don't add it if thats aldeary been done by a previous equation
                      variablesList.append(mb) 
                      
          if commonSetMb or commonSetBool: 
              ReactorExpression.append(strExpr2add)
              reactorDict.update({i:[variablesList, ReactorExpression]})                                    
# =============================================================================
#                   try: 
#                       mb = mbDic[j][countKey]
#                       counterDict[j]  +=  1 
#                       strExpr2add = x.replace(j, '%s*%s' % (mb,j))
#                   except: 
#                       print('Error: check the inputs for the reactors, it seems like %s is connected to more reactors than defined in the excel sheet Connection data' % j )
#                   
# =============================================================================
        
          if len(commonSetBool)== 0 and len(commonSetMb) == 0:
              ReactorExpression.append(x)
              reactorDict.update({i:[variablesList, ReactorExpression]})
                 
      for u in commonSetBool: 
          counterDict[u]  +=  1 
      
      for u in commonSetMb: 
          counterDict[u]  +=  1               
                    
   

   # elements to predifine before the function reactorConstraints is used 
   #costUtility = 0 # predifine the utility cost 


   def reactorConstraints(b,i): # migth have to use indexed dictionary which is nested reactor name+ input output of said reactor 
       # turn to var  
       model
       components = reactorDict[i][0]
       strExpr = reactorDict[i][1]
       b.reactorConstrainsList = pe.ConstraintList()
       
       # call reator 
       #reactor_to_call = getattr(exec(reactorLib), i)  # TODO change name to offical library script 
       # getattr is great when using classes! might be super interesting if you change reactors into objects... 
       # link https://www.programiz.com/python-programming/methods/built-in/getattr
       reactor_to_call = getattr(testLibrary, i)       
       reactorInfo = reactor_to_call() # first call is in_outdict second is reactor expresion 

   
       
       # set an extra constraint on the reactor output  y*x_lo <= x <= y*x_up
       # coulb be done more elegantly by making a dict {output: accompaning Bool}
       if set(components)&set(BoleanVarAll): # TODO fix bug: if a reactor has multiple bools variables ur fucked 
           #outs = list(set(components)& set(nOut))
           outs = reactorInfo[0]['outputs']
           theBool = list(set(components)&set(BoleanVarAll))
           
           for u in boolDic:
               y = theBool[0]
               if y in boolDic[u]: # assuming one bool is posible for each output 
                   thestream = u    # TODO fix bug: if 2 reactor inputs are both boolen then their will be conflict!!! could be fixed if you put for o in outs in this loop!!! try example!!
                   
           for o in outs:  
               if o in nOut:
                   lower = outBoundDict[o][0]*model.BoolBlock[thestream].BoleanVariables[theBool[0]]
                   upper = outBoundDict[o][1]*model.BoolBlock[thestream].BoleanVariables[theBool[0]]
                   
                   b.reactorConstrainsList.add(expr = lower <= model.output[o])
                   b.reactorConstrainsList.add(expr =  model.output[o] <= upper ) 
               
               elif o in streamVariable: 
                   lower = streamDict[o][0]*model.BoolBlock[thestream].BoleanVariables[theBool[0]]
                   upper = streamDict[o][1]*model.BoolBlock[thestream].BoleanVariables[theBool[0]]
                   
                   b.reactorConstrainsList.add(expr = lower <= model.streamVariables[o])
                   b.reactorConstrainsList.add(expr =  model.streamVariables[o] <= upper ) 
                   
               
       ################################### mixing unit 
# =============================================================================
#        inputs = reactorInfo[0]['inputs']
#        totalMassIn = 0 
#        for z in inputs:
#            if z in nIn :
#                totalMassIn += model.input[z]
#            elif z in streamVariable: 
#                totalMassIn += model.streamVariables[z]
#                
#        #totaalMassIn =  sum(model.input[z] for z in inputs)   # in assumption that units are kg/h or smth like that
#        
#        #################################### utility unit 
#        costUtility =0 # TODO need to find smth that can accumulate the costs of the utility 
#        if reactorInfo[3]:  # activate if utility is used 
#            costUtility += (totalMassIn) * reactorInfo[3]
#          
# =============================================================================
       
       
       #################################### reactor units
       for k in strExpr:
           strExprCurrent = k 
           
           for j in components:
               
               if j in nIn: # 
                  replaceStr = "model.input['%s']" % j
                  strExprCurrent = strExprCurrent.replace(j,replaceStr)
                  
               elif j in nOut:#sum(nOut == j):
                  replaceStr = "model.output['%s']" % j
                  strExprCurrent = strExprCurrent.replace(j,replaceStr)
                  
                  
               elif j in list(streamDict.keys()):# may by use set() & set(stremDict:keys()) notation 
                  replaceStr = "model.streamVariables['%s']" % j
                  strExprCurrent = strExprCurrent.replace(j,replaceStr)
                  
               elif set([j])&set(BoleanVarAll):
                  for u in boolDic:
                      if j in boolDic[u]:
                          thestream = u 
                       
                  replaceStr = "model.BoolBlock['%s'].BoleanVariables['%s']" % (thestream,j)
                  # in this form: BoolBlock[str1].BoleanVariables[Y0]
                  strExprCurrent = strExprCurrent.replace(j,replaceStr)
                  
                  
               elif set([j]) & set(mbVarAll):
                   for u in mbDic: 
                       if j in mbDic[u]: 
                           thestream = u 
                           
                       replaceStr = "model.mbBlock['%s'].mbVariables['%s']" % (thestream,j)                            
                       strExprCurrent = strExprCurrent.replace(j,replaceStr)
                            
               else: 
                  warnings.warn("WARNING THE VARIABLE USED IN THE LIBRARY COULD NOT BE FOUND; CHECK IF THE NAMES IN EXCEL AND THE LIBRARY COINCIDE check reactor %s" % i)
       
       
           try: 
               # eg 
               #a = model.output['p1'] == model.streamVariables['str2'] #eval(strExprCurrent)
               a = eval(strExprCurrent)
               b.reactorConstrainsList.add(expr = a) 
      
           except : # create a warning if names in excel don't match names in reator lib.
               print('ERROR make sure the reactor names of the inputs outputs and other components in EXCEL are the same as the library the error is in the following eqation: %s' % strExprCurrent)
             
       #TODO 
        #################################### Waste unit 
   
   
   model.reactorBlock = pe.Block(reactorDict,rule=reactorConstraints)     

   # recall all objective elements and formulate the objective equation 
   model.obj = pe.Objective(sense=pe.maximize, expr=obj_expr_sell - obj_expr_buy)#  - costUtility  )




   return model 

# test to see if the varable/parameters 
if __name__ == '__main__' :
    
    #model = makeModel('testModel2.xlsx','testLibrary')
    model = makeModel('testModel3.xlsx','testLibrary.py')
    #model.pprint()
    
    
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
    

