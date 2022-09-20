'''
Created on Tue sep 20 2022
Contains the classes and functions to make the process interval objects

@author: Lucas Van der hauwaert
email: lucas.vanderhauwaert@usc.es
'''

import os
import pandas as pd

class makeReactor:
    def __init__(self, inputs, outputs, eq, mix=[], utilities=[], isBool=[], split=[], separation=[]):
        self.inputs = inputs
        self.outputs = outputs  # updated according to the modula (is there split/seperation/bool)
        self.mix = mix  # found by the excel file (don't need to specify in this script)
        self.utilities = utilities
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





class inputCharaterisation:
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


# read function to automate making the interval classes
def makeInputIntervals(excelName):
    loc = os.getcwd()
    loc = loc + r'\excel files' + excelName
    DFIntervals = pd.read_excel(loc, sheet_name='componets')
    # inputs
    inputPrices = DFIntervals.input_price.to_numpy()
    posInputs = inputPrices != 0    #find where the input interval are (they have an input price)

    intervalNames = DFIntervals.process_intervals[posInputs]  # find names of input interval variable
    componentsList =  DFIntervals.components[posInputs]
    compositionsList =  DFIntervals.composition[posInputs]

    #loop over all the inputs and make a class of each one
    for i, intervalName in enumerate(intervalNames):
        componentsOfInterval = componentsList[i].split(",")
        compositionsofInterval = compositionsList[i].split(",")

        compsitionDictionary = {}
        for j,component in enumerate(componentsOfInterval):
            component = component.replace(' ','')  #get rid of spaces
            fraction = compositionsofInterval[j]
            fraction = fraction.replace(' ','')
            compsitionDictionary.update({component:fraction})

        toExecute = '{0} = inputCharaterisation(intervalName,compsitionDictionary)'.format(intervalName)
        exec(toExecute)



