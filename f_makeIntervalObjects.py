'''
Created on Tue sep 20 2022
Contains the classes and functions to make the process interval objects

@author: Lucas Van der hauwaert
email: lucas.vanderhauwaert@usc.es
'''

import os
import pandas as pd
import numpy as np

# ============================================================================================================
# Input, reactor and out put Class
# ============================================================================================================

class ReactorIntervalClass:
    def __init__(self, inputs, outputs, eq, name, mix=[], utilities=[], isBool=[], split=[], separation=[]):
        self.inputs = inputs
        self.outputs = outputs  # updated according to the modula (is there split/seperation/bool)
        self.name = name
        self.mix = mix  # found by the Excel file (don't need to specify in this script)
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

        if isinstance(eq, str):  # failsafe if you forget to code the string reactor expression as a list
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

class InputIntervalClass:
    def __init__(self, inputName, compositionDict, inputPrice, boundry, Bool=None, split=None, separation=None):
        if separation is None:
            separation = []
        if split is None:
            split = []
        if Bool is None:
            Bool = []

        # declare input name
        self.inputName = inputName.upper()
        addOn4Variables = inputName.lower() + '_'

        # change the composition names in the dictionary
        compositionDictNew = {}
        for i in compositionDict:
            compositionDictNew.update({addOn4Variables + i: compositionDict[i]})
        self.compositionDict = compositionDictNew

        # make the component equations as string equations
        eqList = []
        for component in compositionDictNew:
            eq = "{} == {} * {}".format(component, self.compositionDict[component], self.inputName)
            eqList.append(eq)
        self.componentEquations = eqList
        self.variables =  list(compositionDictNew.keys()) # [self.inputName] not gona add this atm
        self.inputPrice = inputPrice

        # make boundry for all variables
        self.boundryInput = boundry

        # do the bounds in the f_make super_structure
        # variableBoundries = boundry
        # lowerB = boundry[0]
        # upper = boundry[1]
        # for i in self.variables:
        #     pass

        self.Bool = Bool
        self.split = split
        self.separation = separation
        self.compositionNames = list(compositionDict.keys())

    def rename_components(self):
        pass

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

# ============================================================================================================
# Validate the Excel file sheets, checks and error messages
# ============================================================================================================
"""
list of error 
* make sure all separations are acounted for 
* make sure all names of the intervals are the same in each excel sheet 
write all error you encounter here 
"""
def validate_seperation_coef(coef, intervalName, amountOfSep, connectionMatrix):
    coefList  = split_remove_spaces(coef, ';')
    arraysOfCoef = []
    for i in coefList:
        coefTuple = stringbounds_2_tuplebounds(i)
        coefArray = np.array(coefTuple)
        arraysOfCoef.append(coefArray)
    sumOfArrays = 0
    # check if there are as many separation processes as separation coefficient bounds
    if amountOfSep != len(arraysOfCoef):
        raise Exception("make sure all bounds are written for interval ".format(intervalName))
    # check if the sum of the separation coefficiets are one
    for i in arraysOfCoef:
        sumOfArrays += i
    nCol = len(sumOfArrays)  #get the amount of columns
    if sum(sumOfArrays == np.ones(nCol)):
        pass
    else:
        raise Exception("the sum of the seperation coefficients for {} do not add up to 1, check the excel file".format(intervalName))

    # check if all the seperation processes are accounted for in the connenction matrix
    rowId = list(connectionMatrix['process_intervals']).index(intervalName)
    rowMatrix = connectionMatrix.iloc[rowId].drop(['process_intervals'])
    counter = 0
    for i in rowMatrix:
        if 'sep' in str(i):
            counter += 1
    if counter != amountOfSep:
        raise Exception("Check interval {} there is a separation process missing in the connection matrix".format(intervalName))
def check_excel_file(excelName):
    loc = os.getcwd()
    posAlquimia = loc.find('Alquimia')
    loc = loc[0:posAlquimia+8]
    loc = loc + r'\excel files' + excelName

    DFIntervals = pd.read_excel(loc, sheet_name='Intervals')
    DFReactors =  pd.read_excel(loc, sheet_name='ConnectionMatrix')
    DFConnectionMatrix = pd.read_excel(loc, sheet_name='ConnectionMatrix')

    # check interval names in the connection matrix and interval list
    intervalNamesItervals = remove_spaces(DFIntervals['process_intervals'].to_list())
    intervalNamesConnenctionMatrixRow = remove_spaces(list(DFConnectionMatrix.columns))
    intervalNamesConnenctionMatrixRow.remove('process_intervals')
    intervalNamesConnenctionMatrixCol = remove_spaces(DFConnectionMatrix['process_intervals'].to_list())

    # check length
    if len(intervalNamesConnenctionMatrixCol) == len(intervalNamesConnenctionMatrixRow) == len(intervalNamesItervals):
        pass
    else:
        raise Exception('Interval name is missing in the connection matrix sheet or the interval sheet')
    # check names
    if intervalNamesItervals == intervalNamesConnenctionMatrixRow == intervalNamesConnenctionMatrixCol:
        pass
    else:
        positonError = [errorList for i, errorList in enumerate(intervalNamesItervals) if not intervalNamesItervals[i]
                                    == intervalNamesConnenctionMatrixRow[i] == intervalNamesConnenctionMatrixCol[i]]
        print(positonError)
        raise Exception('The names in the connection matrix sheet or the interval sheet are not the same')

# ============================================================================================================
# Usefull functions
# ============================================================================================================
def split_remove_spaces(expr2split,splitCharacter):
    expresions = expr2split.split(splitCharacter)
    exprList = []
    for i in expresions:
        i = i.replace(' ', '')
        exprList.append(i)
    return exprList

def stringbounds_2_tuplebounds(stringBound):
    stringBound = stringBound.replace('[', '')
    stringBound = stringBound.replace(']', '')
    bounds = stringBound.split(',')
    boundsList = []
    for i in bounds:
        boundsList.append(float(i))
    return boundsList

def load_objectes_from_dictionary(dict):
    for i in dict:
        locals()[i] = dict[i]

def remove_spaces(listOfInterest):
    exprList = []
    for i in listOfInterest:
        i = i.replace(' ','')
        exprList.append(i)
    return exprList
# ============================================================================================================
# Functions to make the interval objects
# ============================================================================================================

# read function to automate making the interval classes
def make_input_intervals(excelName):
    loc = os.getcwd()
    posAlquimia = loc.find('Alquimia')
    loc = loc[0:posAlquimia+8]
    loc = loc + r'\excel files' + excelName

    DFIntervals = pd.read_excel(loc, sheet_name='Intervals')
    # inputs
    inputPrices = DFIntervals.input_price.to_numpy()
    posInputs = inputPrices != 0    #find where the input interval are (they have an input price)

    intervalNames = DFIntervals.process_intervals[posInputs]  # find names of input interval variable
    componentsList =  DFIntervals.components[posInputs]
    compositionsList =  DFIntervals.composition[posInputs]

    ####
    inBoundsLow = DFIntervals.lower_bound[posInputs].to_numpy()
    inBoundsUpper = DFIntervals.upper_bound[posInputs].to_numpy()
    inputPrices = inputPrices[posInputs]

    # define fixed parameters cost raw material
    inputPriceDict = {intervalNames[i]: inputPrices[i] for i in range(len(inputPrices))}  # make dictionary
    boundryDict = {intervalNames[i]: [inBoundsLow[i], inBoundsUpper[i]] for i in range(len(inputPrices))}  # make dictionary
    ####
    objectDictionary = {}
    #loop over all the inputs and make a class of each one
    for i, intervalName in enumerate(intervalNames):
        inputPrice = inputPriceDict[intervalName]
        boundryInput = boundryDict[intervalName]
        componentsOfInterval = componentsList[i].split(",")
        compositionsofInterval = compositionsList[i] # string or 1, depending if there are different components
        compsitionDictionary = {} # preallocate dictionary
        if compositionsofInterval == 1:  # if it is one no need to loop over the dictionary, there is only one compound
            component = componentsOfInterval[0].replace(' ','')
            fraction = compositionsofInterval # should allways be one i there is one component in the stream
            compsitionDictionary.update({component: fraction})
        else:
            compositionsofInterval = compositionsList[i].split(",")
            for j,component in enumerate(componentsOfInterval):
                component = component.replace(' ','')  # get rid of spaces
                fraction = compositionsofInterval[j]
                fraction = fraction.replace(' ','')
                fraction = float(fraction)
                compsitionDictionary.update({component:fraction})

        objectInput = InputIntervalClass(intervalName,compsitionDictionary,inputPrice,boundryInput)
        objectDictionary.update({intervalName:objectInput})

        # toExecute = '{0} = inputCharaterisation(intervalName,compsitionDictionary)'.format(intervalName)
        # exec(toExecute)
    return objectDictionary

def make_reactor_intervals(excelName):
    # read Excel file
    loc = os.getcwd()
    posAlquimia = loc.find('Alquimia')
    loc = loc[0:posAlquimia + 8]
    loc = loc + r'\excel files' + excelName

    #read excel info
    DFreactors = pd.read_excel(loc, sheet_name='Reactors')
    connectionMatrix = pd.read_excel(loc, sheet_name='ConnectionMatrix')
    reactorIntervals = DFreactors.reactor_name
    objectDictionary = {} # preallcoate a dictionary with the interval names and the interval objects
    for i, intervalName in enumerate(reactorIntervals):
        #inputs of reactor
        inputsReactor = DFreactors.inputs[i]
        inputsReactor = split_remove_spaces(inputsReactor,',')
        #outputs of the reactor
        outputsReactor = DFreactors.outputs[i]
        outputsReactor = split_remove_spaces(outputsReactor,',')
        #find the equation of the reactor
        equations = DFreactors.equations[i]
        equations = split_remove_spaces(equations,',')
        # make initial object (with minimum requirements i.e., inputs outputs and reactor equations)
        objectReactor = ReactorIntervalClass(inputsReactor,outputsReactor,equations, intervalName)

        if DFreactors.has_utility[i] != 0 :
            utilityVariableNames = DFreactors.has_utility[i]
            utilityVariableNames = split_remove_spaces(utilityVariableNames,';')
            utilityBounds = DFreactors.utility_bounds[i]
            utilityBounds = split_remove_spaces(utilityBounds,';')
            utilityDict = {}
            for j, unitName in enumerate(utilityVariableNames):
                unitBound = utilityBounds[j]
                tupleUnitBounds = stringbounds_2_tuplebounds(unitBound)
                utilityDict.update({unitName:tupleUnitBounds})
            objectReactor.utilities = utilityDict


        if DFreactors.has_seperation[i] >= 2: #and DFreactors.has_seperation[i] < 2 :
            nrSeperations = DFreactors.has_seperation[i]
            outputsStr = DFreactors.outputs[i]
            amountOfSeperations = DFreactors.has_seperation[i]
            coefStr = DFreactors.seperation_coef[i]
            validate_seperation_coef(coefStr,intervalName,nrSeperations, connectionMatrix)
            coefList = split_remove_spaces(coefStr, ';')
            seperationDict = {}
            for j in range(amountOfSeperations):
                seperationName = intervalName + '_sep{}'.format(j+1)
                coefTuple = stringbounds_2_tuplebounds(coefList[j])
                outputs = split_remove_spaces(outputsStr, ',' )
                specificSeperationDict = {}
                for k, outputName in enumerate(outputs):
                    specificSeperationDict.update({outputName:coefTuple[k]})
                seperationDict.update({seperationName: specificSeperationDict})
            objectReactor.separation = seperationDict

            # check if it is mixed with other reactors
            processInvervalNames = connectionMatrix['process_intervals']
            #reactorRow = connectionMatrix.loc[processInvervalNames == i]
            reactorCol = connectionMatrix[intervalName]
            pos = reactorCol != 0 # find where mixing takes place # mixed streams are in the same colunm


            if sum(pos) >= 2:
                mixDict = {}
                intervalsToMix = list(processInvervalNames[pos])
                specifications = list(reactorCol[pos])
                for k,specs in enumerate(specifications):
                    if 'mix' in specs:
                        mixDict.update( {intervalsToMix[k]: specs} )
                #mixDict = {intervalsToMix[j]: specifications[j] for j in range(0, len(intervalsToMix))}
                objectReactor.mix = mixDict

        # put the object in the dictionary
        objectDictionary.update({intervalName:objectReactor})
    return objectDictionary

#def makeObjectDictionaries():

if __name__ == '__main__':
    testObje = make_reactor_intervals(r'\data_propionibacteria.xlsx')
    location = os.getcwd()
    print(location)
