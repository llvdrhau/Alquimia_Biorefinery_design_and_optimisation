'''
Created on Tue sep 20 2022
Contains the classes and functions to make the process interval objects

@author: Lucas Van der hauwaert
email: lucas.vanderhauwaert@usc.es
'''

import os
import pandas as pd

class reactorIntervalClass:
    def __init__(self, inputs, outputs, eq, name, mix=[], utilities=[], isBool=[], split=[], separation=[]):
        self.inputs = inputs
        self.outputs = outputs  # updated according to the modula (is there split/seperation/bool)
        self.name = name
        self.mix = mix  # found by the excel file (don't need to specify in this script)
        self.utilities = utilities # consists of a dictionary {nameUtilty: [bounds]}
        self.isBool = isBool  # is a tuple, 1) where the bool comes from and 2) the name of the bool affecting the outputs
        self.split = split  # true or false # in that case of false just empty []
        self.separation = separation  # dictionary defining seperation fractions of each component

        # =============================================================================
        #           # introduced in the script in aplicable
        #         self.separationOuput
        #         self.splitOuput
        #         self.inputMix
        # =============================================================================

        if isinstance(eq, str):  # fail safe if you forget to code the string reactor expression as a list
            self.eq = [eq]  # make it a list
        else:
            self.eq = eq

    def makeDict(self):
        in_outDict = {
            'inputs': self.inputs,
            'outputs': self.outputs
        }
        return in_outDict

    def makeOldNewDictOutputs(self, oldNames, newNames):
        oldNewDict = {oldNames[i]: newNames[i] for i in range(len(oldNames))}  # make dictionary
        self.oldNewDictOutputs = oldNewDict

    def makeOldNewDictInputs(self, oldNames, newNames):
        oldNewDict = {oldNames[i]: newNames[i] for i in range(len(oldNames))}  # make dictionary
        self.oldNewDictInputs = oldNewDict

    def makeSeperationDict(self, permeate,
                           reject):  # define dictionary to see wat percentage of each stream component goes to the permeate and reject streams
        percentExtractionVector = list(self.separation.values())
        seperationDict = {self.outputs[i]: [permeate[i], reject[i], percentExtractionVector[i]] for i in
                          range(len(self.outputs))}
        return seperationDict

    # def updateSeperationDict(self):
    #     new_key = self.outputs
    #     old_key = 1
    #     #dictionary[new_key] = dictionary.pop(old_key)

    # alternative class method
    @classmethod
    def fromNeuralNetwork(cls, neuralNetwork):
        # so here comes the alternative constructor for neural networks wich will make the string equations from a machine learning model
        # https://www.youtube.com/watch?v=rq8cL2XMM5M&ab_channel=CoreySchafer
        pass

        #  could you define a stream dependant on a boolean var before in enters a unit reactor

class inputIntervalClass:
    def __init__(self, inputName, compositionDict, isBool=[], split=[], separation=[]):
        self.inputName = inputName
        self.compositionDict = compositionDict
        self.compositionNames = list(compositionDict.keys())
        self.isBool = isBool
        self.split = split
        self.separation = separation

    def makeOldNewDict(self, oldNames, newNames):
        oldNewDict = {oldNames[i]: newNames[i] for i in range(len(oldNames))}  # make dictionary
        self.oldNewDictOutputs = oldNewDict

    def UpdateDict(self, newKeys):
        oldDict = self.compositionDict
        oldKeys = list(oldDict.keys())
        rangeKeys = range(len(oldKeys))
        for i in rangeKeys:
            new_key = newKeys[i]
            old_key = oldKeys[i]
            self.compositionDict[new_key] = self.compositionDict.pop(old_key)

def splitAndRemoveSpaces(expr2split,splitCharacter):
    expresions = expr2split.split(splitCharacter)
    exprList = []
    for i in expresions:
        i = i.replace(' ', '')
        exprList.append(i)
    return exprList

def stringBounds2tupleBounds(stringBound):
    stringBound = stringBound.replace('[', '')
    stringBound = stringBound.replace(']', '')
    bounds = stringBound.split(',')
    upperbound = float(bounds[1])
    lowerbound = float(bounds[0])
    boundsArray = [lowerbound, upperbound]
    return boundsArray

def loadObjectesFromDictionary(dict):
    for i in dict:
        locals()[i] = dict[i]

# ============================================================================================================
# TODO see if you can incorporate these functions (makeInputIntervals, makeReactorInterval) in the classes?
# ============================================================================================================

# read function to automate making the interval classes
def makeInputIntervals(excelName):
    loc = os.getcwd()
    posAlquimia = loc.find('Alquimia')
    loc = loc[0:posAlquimia+8]
    loc = loc + r'\excel files' + excelName

    DFIntervals = pd.read_excel(loc, sheet_name='components')
    # inputs
    inputPrices = DFIntervals.input_price.to_numpy()
    posInputs = inputPrices != 0    #find where the input interval are (they have an input price)

    intervalNames = DFIntervals.process_intervals[posInputs]  # find names of input interval variable
    componentsList =  DFIntervals.components[posInputs]
    compositionsList =  DFIntervals.composition[posInputs]

    objectDictionary = {}
    #loop over all the inputs and make a class of each one
    for i, intervalName in enumerate(intervalNames):
        componentsOfInterval = componentsList[i].split(",")
        compositionsofInterval = compositionsList[i] # string or 1 depending if there are different components
        compsitionDictionary = {} # preallocate dictionary
        if compositionsofInterval == 1:  #if it is one no need to loop over the dictionary, there is only one compound
            component = componentsOfInterval[0].replace(' ','')
            fraction = compositionsofInterval # should allways be one i there is one component in the stream
            compsitionDictionary.update({component: fraction})
        else:
            compositionsofInterval = compositionsList[i].split(",")
            for j,component in enumerate(componentsOfInterval):
                component = component.replace(' ','')  #get rid of spaces
                fraction = compositionsofInterval[j]
                fraction = fraction.replace(' ','')
                fraction = float(fraction)
                compsitionDictionary.update({component:fraction})

        objectInput = inputIntervalClass(intervalName,compsitionDictionary)
        objectDictionary.update({intervalName:objectInput})

        # toExecute = '{0} = inputCharaterisation(intervalName,compsitionDictionary)'.format(intervalName)
        # exec(toExecute)
        return objectDictionary

    #return objectDictionary

def makeReactorIntervals(excelName):
    loc = os.getcwd()
    posAlquimia = loc.find('Alquimia')
    loc = loc[0:posAlquimia + 8]
    loc = loc + r'\excel files' + excelName

    DFreactors = pd.read_excel(loc, sheet_name='reactors')
    reactorIntervals = DFreactors.reactor_name

    objectDictionary = {} # preallcoate a dictionary with the interval names and the interval objects
    for i, intervalName in enumerate(reactorIntervals):
        #inputs of reactor
        inputsReactor = DFreactors.inputs[i]
        inputsReactor = splitAndRemoveSpaces(inputsReactor,',')
        #outputs of the reactor
        outputsReactor = DFreactors.outputs[i]
        outputsReactor = splitAndRemoveSpaces(outputsReactor,',')
        #find the equation of the reactor
        equations = DFreactors.equations[i]
        equations = splitAndRemoveSpaces(equations,',')
        # make initial object (with minimum requirements i.e., inputs outputs and reactor equations)
        objectReactor = reactorIntervalClass(inputsReactor,outputsReactor,equations, intervalName)

        if DFreactors.has_utility[i] != 0 :
            utilityVariableNames = DFreactors.has_utility[i]
            utilityVariableNames = splitAndRemoveSpaces(utilityVariableNames,';')
            utilityBounds = DFreactors.utility_bounds[i]
            utilityBounds = splitAndRemoveSpaces(utilityBounds,';')
            utilityDict = {}
            for j, unitName in enumerate(utilityVariableNames):
                unitBound = utilityBounds[j]
                tupleUnitBounds = stringBounds2tupleBounds(unitBound)
                utilityDict.update({unitName:tupleUnitBounds})
            objectReactor.utilities = utilityDict
        objectDictionary.update({intervalName:objectReactor})
        # Todo check utilities are made correctly.

        if DFreactors.has_seperation[i] != 0 and DFreactors.has_seperation[i] < 2 :  #lets not look at double serperation for now
            seperationDict = {}
            coefStr = DFreactors.seperation_coef[i]
            coefTuple = stringBounds2tupleBounds(coefStr)
            outputsStr = DFreactors.outputs[i]
            outputs = splitAndRemoveSpaces(outputsStr, ',' )
            for j, outputName in enumerate(outputs):
                seperationDict.update({outputName:coefTuple[j]})
            objectReactor.separation = seperationDict

    return objectDictionary

#def makeObjectDictionaries():

if __name__ == '__main__':
    testObje = makeReactorIntervals(r'\data_propionibacteria.xlsx')
    location = os.getcwd()
    print(location)
