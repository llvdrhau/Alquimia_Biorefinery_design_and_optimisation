# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 11:56:05 2022

@author: Lucas Van der hauwaert
"""

import os
import pandas as pd

class makeReactor:
    def __init__(self, inputs, outputs, eq, mix = [], utilities = [], isBool = [], split = [], separation = []):
        self.inputs  = inputs 
        self.outputs = outputs    # updated according to the modula (is there split/seperation/bool)
        self.mix = mix            # found by the excel file (don't need to specify in this script)
        self.utilities = utilities
        self.isBool = isBool       # is a tuple, 1) where the bool comes from and 2) the name of the bool affecting the outputs 
        self.split = split         # true or false # in that case of false just empty []
        self.separation = separation # dictionary defining seperation fractions of each component 
        
# =============================================================================
#           # introduced in the script in aplicable 
#         self.separationOuput
#         self.splitOuput 
#         self.inputMix 
# =============================================================================
        
        if isinstance(eq, str): # fail safe if you forget to code the string reactor expression as a list 
            self.eq = [eq]   # make it a list
        else: 
            self.eq = eq 
        
    def makeDict(self):
        in_outDict = { 
                'inputs': self.inputs,
                'outputs':self.outputs 
                }
        return in_outDict
    
    def makeOldNewDictOutputs(self,oldNames,newNames):
        oldNewDict = {oldNames[i]: newNames[i] for i in range(len(oldNames))} # make dictionary 
        self.oldNewDictOutputs = oldNewDict
        
    def makeOldNewDictInputs(self,oldNames,newNames):
        oldNewDict = {oldNames[i]: newNames[i] for i in range(len(oldNames))} # make dictionary 
        self.oldNewDictInputs = oldNewDict
    
    def makeSeperationDict(self, permeate, reject):  # define dictionary to see wat percentage of each stream component goes to the permeate and reject streams 
        percentExtractionVector =list(self.separation.values())
        seperationDict = {self.outputs[i]: [permeate[i], reject[i], percentExtractionVector[i]] for i in range(len(self.outputs))}
        return seperationDict
    
    # def updateSeperationDict(self):
    #     new_key = self.outputs
    #     old_key = 1 
    #     #dictionary[new_key] = dictionary.pop(old_key)
        
    # alternative class method 
    @classmethod
    def fromNeuralNetwork(cls,neuralNetwork):
        # so here comes the alternative constructor for neural networks wich will make the string equations from a machine learning model 
        # https://www.youtube.com/watch?v=rq8cL2XMM5M&ab_channel=CoreySchafer 
        pass 
   
  # TODO could you define a stream dependant on a boolean var before in enters a unit reactor
# think of sizing example   

def equations_from_FBA(excelName,reactorName,substrate):
    loc = os.getcwd()
    loc = loc + r'\excel files' + excelName
    conversionDF = pd.read_excel(loc, sheet_name='massFraction')


class InputCharaterisation:
    def __init__(self, inputName,compositionDict, isBool = [], split = [], separation = []):
        self.inputName  = inputName 
        self.compositionDict = compositionDict
        self.compositionNames = list(compositionDict.keys())
        self.isBool = isBool 
        self.split = split 
        self.separation = separation
        
    def makeOldNewDict(self, oldNames, newNames):
        oldNewDict = {oldNames[i]: newNames[i] for i in range(len(oldNames))} # make dictionary 
        self.oldNewDictOutputs = oldNewDict
        
        
    def UpdateDict(self,newKeys):
        oldDict = self.compositionDict
        oldKeys = list(oldDict.keys())
        rangeKeys = range(len(oldKeys))
        for i in rangeKeys: 
            new_key = newKeys[i]
            old_key = oldKeys[i]
            self.compositionDict[new_key] = self.compositionDict.pop(old_key) 
            
            
# =============================================================================
# # be very carefull with naming the compositions in the library, doen't use single letters!!!
# eg aa instead of a 
# =============================================================================



############################################################ reactor1

strExpr = ['xx == aa * 0.8 ',
           'yy == bb * 0.7' ,
           'zz == aa * 0.1 + bb * 0.3' ]   
          
ins = ['aa','bb','cc']
outs = ['xx','yy','zz']

reactor1 = makeReactor(ins, outs, strExpr)  


############################################################ reactor2

#mix = ['input2', 'reactor1']  # also defined by excel file now

ins = ['xx','yy','zz']
outs = ['uu','vv','ww']

strExpr = ['uu == xx *0.8 + zz*0.1 ',
           'vv == yy*0.3 + xx*0.1'  ,
           'ww == 0.2*zz+0.05*yy+0.09*xx'] 


separationPermeate = {'uu': 0.95,
                      'vv': 0.8,
                      'ww' : 0 }


is_bool = [] # 4 options or the bool is reject stream, permeate stream, reactor stream  or [] (no boolean variables)
is_split = []

reactor2 = makeReactor(ins, outs, strExpr,mix,[],is_bool,is_split,separationPermeate)  


############################################################ reactorSeperation

ins = ['xx','yy','zz']
outs = ['uu','ww','vv']

strExpr = ['uu == xx *0.8 + zz*0.1 ',
           'vv == yy*0.3 + xx*0.1'  ,
           'ww == 0.2*zz+0.05*yy+0.09*xx'] 

separationPermeate = {'uu': 0.95,
                      'vv': 0.8,
                      'ww' : 0 }


is_bool = []  # determined by excel file 
is_split = []

reactorSeperation = makeReactor(ins, outs, strExpr,[],[],is_bool,is_split,separationPermeate) 

############################################################ reactor3
strExpr = ['p1 == 0.5*uu + 0.05*vv' ]
ins = ['uu','vv']
outs = ['p1']
utilities  = 0
reactor3 = makeReactor(ins, outs, strExpr,utilities)  

############################################################ reactor4 
strExpr = ['p2 == 0.4*uu' ]
ins = ['uu']
outs = ['p2']
utilities  = 0

reactor4 = makeReactor(ins, outs, strExpr)  

############################################################ reactor5
strExpr = ['p3 == 0.2*ww']
            
ins = ['ww']
outs = ['p3']

reactor5 = makeReactor(ins, outs, strExpr)  

############################################################ reactor6

strExpr = ['p4 == ww * 0.6 ']
            
ins = ['ww']
outs = ['p4']
reactor6 = makeReactor(ins, outs, strExpr)  

############################################################ rExtra1 

strExpr = ['ww == ethl * 0.6 ']
            
ins = ['ethl']
outs = ['ww']
rExtra1 = makeReactor(ins, outs, strExpr)  

############################################################ rExtra2

strExpr = ['pExtra == ww * 0.6 ']
            
ins = ['ww']
outs = ['pExtra']
rExtra2 = makeReactor(ins, outs, strExpr)  


#############################################################################
 # Input charaterisation 
#############################################################################  

############################################################ input 1 
inputName = 'input1'
composition = {'aa': 0.3, 'bb': 0.7}

input1 = InputCharaterisation(inputName, composition)

############################################################ input 2 

inputName = 'input2'
composition = {'xx': 0.1, 'yy': 0.9}

input2 = InputCharaterisation(inputName, composition)

############################################################ input 3 

inputName = 'input3'
composition = {'ww': 1}

input3 = InputCharaterisation(inputName, composition)

############################################################ inExtra 

inputName = 'inExtra'
composition = {'ethl': 1}

inExtra = InputCharaterisation(inputName, composition)

