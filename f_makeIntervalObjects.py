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
# Input, reactor and output Classes
# ============================================================================================================


class InputIntervalClass:
    def __init__(self, inputName, compositionDict, inputPrice, boundry, boolDict = None, split=None, separation=None):
        if separation is None:
            separation = []
        if split is None:
            split = []
        if boolDict is None:
            boolDict = []

        # declare input name
        self.inputName = inputName.upper()
        self.label = 'input'
        addOn4Variables = inputName.lower() + '_'

        # change the composition names in the dictionary
        compositionDictNew = {}
        initialCompositionNames = []
        for i in compositionDict:
            compositionDictNew.update({addOn4Variables + i: compositionDict[i]})
            initialCompositionNames.append(i)
        self.initialCompositionNames = initialCompositionNames
        self.compositionDict = compositionDictNew

        # make the component equations as string equations
        eqList = []
        for component in compositionDictNew:
            eq = "{} == {} * {}".format(component, self.compositionDict[component], self.inputName)
            eqList.append(eq)
        self.componentEquations = eqList
        self.componentVariables =  list(compositionDictNew.keys()) # [self.inputName] not gona add this atm
        self.allVariables = list(compositionDictNew.keys()) + [self.inputName]
        self.inputPrice = inputPrice

        # make boundary for all variables
        boundaryDict = {self.inputName:boundry} # intervals have bounds
        for i in self.componentVariables:
            boundaryDict.update({i:(0, None)})
        self.boundaries = boundaryDict
        self.boolDict = boolDict
        self.split = split
        self.separation = separation
        #self.compositionNames = list(compositionDict.keys())
        self.leavingInterval = self.componentVariables

        eqSumOfBoolsHelp = '1 == '
        eqSumOfBools = []
        if self.boolDict:
            for interval in boolDict:
                eqSumOfBoolsHelp += ' + ' + boolDict[interval]
            eqSumOfBools = [eqSumOfBoolsHelp]
        self.eqSumOfBools = eqSumOfBools

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

class ReactorIntervalClass:
    def __init__(self, inputs, outputs, reactionEquations, nameDict, mix=None, utilities=None, boolActivation = None ,boolDict=None, split=None, separation=None):
        if utilities is None:
            utilities = {}
        if separation is None:
            separation = {}
        if split is None:
            split = []
        if boolDict is None:
            boolDict = []
        if mix is None:
            mix = []
        if boolActivation is None:
            boolActivation = []

        self.label = 'reactor'
        self.name = nameDict
        self.inputs = inputs
        self.outputs = outputs  # updated according to the modula (is there split/seperation/bool)
        self.initialCompositionNames = inputs + outputs
        # error if you choose a wrong name
        for i in self.initialCompositionNames: #maybe not here the warning
            if i in list(nameDict.keys())[0]:
                raise Exception("the component {} is in interval name {}, change the component name of the reactor")

        self.mix = mix  # found by the Excel file (don't need to specify in this script)
        self.utilities = utilities # consists of a dictionary {nameUtilty: [bounds]}
        self.boolDict = boolDict  # is a tuple, 1) where the bool comes from and 2) the name of the bool affecting the outputs
        self.split = split  # true or false # in that case of false just empty []
        self.separation = separation  # dictionary defining separation fractions of each component

        if isinstance(reactionEquations, str):  # failsafe if you forget to code the string reactor expression as a list
            ReactionEquations = [reactionEquations]  # make it a list
        else:
            pass
        ouputs2change  = self.outputs
        if self.utilities:
            ouputs2change += list(self.utilities.keys())
        allEquations = []
        reactionVariablesOutput = []
        for eq in reactionEquations:
            for out in ouputs2change:
                newOutputName = out + '_{}'.format(self.name)
                reactionVariablesOutput.append(newOutputName)
                eq = eq.replace(out,newOutputName)
            allEquations.append(eq)
        self.reactionEquations = allEquations
        self.reactionVariablesOutput = reactionVariablesOutput
        self.intervalVariable = list(self.name.keys())[0]

        # mass equations (of the outputs from the reaction equations )
        eqMassInterval  = self.intervalVariable + " == "
        for out in reactionVariablesOutput:
            eqMassInterval += " + " + out
        self.totalMassEquation =  [eqMassInterval]

        # bool activation constraints
        boolActivationEquations = []
        if boolActivation: # if there is an activation constraint
            bounds = nameDict[self.intervalVariable]
            lowerActivationEq = "{} * {} <= {}".format(boolActivation[0],bounds[0],self.intervalVariable)
            upperActivationEq = "{} <= {} * {}".format(self.intervalVariable,boolActivation[0], bounds[1])
            boolActivationEquations.append(lowerActivationEq)
            boolActivationEquations.append(upperActivationEq)
        self.boolActivationEquations =  boolActivationEquations
        self.activationVariable =boolActivation[0]


        # sum of bool equations
        eqSumOfBoolsHelp = '1 == '
        eqSumOfBools = []
        if self.boolDict:
            for interval in boolDict:
                eqSumOfBoolsHelp += ' + ' + boolDict[interval]
            eqSumOfBools = [eqSumOfBoolsHelp]
        self.eqSumOfBools = eqSumOfBools

        # separation equations

        # mixing equations

        # spliting equations

        # define wat is leaving the reactor


    def makeDict(self):
        in_outDict = {
            'inputs': self.inputs,
            'outputs': self.outputs
        }
        return in_outDict

    def make_replacement_dict_output(self, oldNames, newNames):
        oldNewDict = {oldNames[i]: newNames[i] for i in range(len(oldNames))}  # make dictionary
        self.replacementDictOutput= oldNewDict
        return oldNewDict

    def make_replacement_dict_input(self, oldNames, newNames):
        oldNewDict = {oldNames[i]: newNames[i] for i in range(len(oldNames))}  # make dictionary
        self.replacementDictInput = oldNewDict
        return oldNewDict

    def makeSeperationDict(self, permeate,
                           reject):  # define dictionary to see wat percentage of each stream component goes to the permeate and reject streams
        percentExtractionVector = list(self.separation.values())
        seperationDict = {self.outputs[i]: [permeate[i], reject[i], percentExtractionVector[i]] for i in
                          range(len(self.outputs))}
        return seperationDict
    @classmethod
    def fromNeuralNetwork(cls, neuralNetwork):
        # so here comes the alternative constructor for neural networks wich will make the string equations from a machine learning model
        # https://www.youtube.com/watch?v=rq8cL2XMM5M&ab_channel=CoreySchafer
        pass

        #  could you define a stream dependant on a boolean var before in enters a unit reactor

class OutputIntervalClass:
    def __init__(self, outputName, outputPrice, Bool=None, split=None, separation=None):
        if separation is None:
            separation = []
        if split is None:
            split = []
        if Bool is None:
            Bool = []
        self.outputName = outputName
        self.outputPrice = outputPrice
        self.label = 'output'
# ============================================================================================================
# Validate the Excel file sheets, checks and error messages
# ============================================================================================================
"""
List of possible errors 
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

def get_connected_intervals(intervalName,conectionMatrix):
    conectionCol = conectionMatrix[intervalName]
    posConnect = conectionCol != 0
    nameConnectedIntervals = list(conectionMatrix['process_intervals'][posConnect])
    connectionInfo = list(conectionMatrix[intervalName][posConnect])
    connectionDict = {nameConnectedIntervals[i]:connectionInfo[i] for i in range(len(connectionInfo))}
    return  connectionDict

def get_replacement_dict(initialVars, newVars):
    replacementDict = {}
    for i in initialVars:
        for j in newVars:
            if i in j:  # the initial is always in the new name, thats how you can find wat belongs to where
                replacementDict.update({i:j})
    return replacementDict
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
    DFconnectionMatrix = pd.read_excel(loc, sheet_name='ConnectionMatrix')
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

        # check if it is an input to other intervals as a bool
        boolDict = {}
        processInvervalNames = DFconnectionMatrix['process_intervals'].to_list()
        rowIndex = processInvervalNames.index(intervalName)
        intervalRow = DFconnectionMatrix.iloc[rowIndex].to_list() # looking at the row will show to which intervals the curretn section is connencted to
        # intervalRow = DFconnectionMatrix.loc[processInvervalNames == intervalName].to_dict()
        for j, info in enumerate(intervalRow):
            if isinstance(info,str) and 'bool' in info:
                attachInterval = processInvervalNames[j]
                boolVar = 'y_' + attachInterval + '_' + intervalName
                boolDict.update({attachInterval: boolVar})


        # create object
        objectInput = InputIntervalClass(inputName=intervalName,compositionDict=compsitionDictionary,
                                         inputPrice=inputPrice,boundry=boundryInput,boolDict= boolDict)
        objectDictionary.update({intervalName:objectInput})
    return objectDictionary

def make_reactor_intervals(excelName):
    # read Excel file
    loc = os.getcwd()
    posAlquimia = loc.find('Alquimia')
    loc = loc[0:posAlquimia + 8]
    loc = loc + r'\excel files' + excelName

    #read excel info
    DFIntervals = pd.read_excel(loc, sheet_name='Intervals')
    DFreactors = pd.read_excel(loc, sheet_name='Reactors')
    DFconnectionMatrix = pd.read_excel(loc, sheet_name='ConnectionMatrix')

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
        # find the bounds of the interval
        listIntervals = list(DFIntervals.process_intervals)
        index = listIntervals.index(intervalName)
        lowerBound = DFIntervals.lower_bound[index]
        upperBound = DFIntervals.upper_bound[index]
        nameDict = {intervalName:[lowerBound,upperBound]}



        utilityDict = {}
        if DFreactors.has_utility[i] != 0 :
            utilityVariableNames = DFreactors.has_utility[i]
            utilityVariableNames = split_remove_spaces(utilityVariableNames,';')
            utilityBounds = DFreactors.utility_bounds[i]
            utilityBounds = split_remove_spaces(utilityBounds,';')
            for j, unitName in enumerate(utilityVariableNames):
                unitBound = utilityBounds[j]
                tupleUnitBounds = stringbounds_2_tuplebounds(unitBound)
                utilityDict.update({unitName:tupleUnitBounds})


        seperationDict = {}
        if DFreactors.has_seperation[i] >= 2: #and DFreactors.has_seperation[i] < 2 :
            nrSeperations = DFreactors.has_seperation[i]
            outputsStr = DFreactors.outputs[i]
            amountOfSeperations = DFreactors.has_seperation[i]
            coefStr = DFreactors.seperation_coef[i]
            validate_seperation_coef(coefStr,intervalName,nrSeperations, DFconnectionMatrix)
            coefList = split_remove_spaces(coefStr, ';')
            for j in range(amountOfSeperations):
                seperationName = intervalName + '_sep{}'.format(j+1)
                coefTuple = stringbounds_2_tuplebounds(coefList[j])
                outputs = split_remove_spaces(outputsStr, ',' )
                specificSeperationDict = {}
                for k, outputName in enumerate(outputs):
                    specificSeperationDict.update({outputName:coefTuple[k]})
                seperationDict.update({seperationName: specificSeperationDict})
            #objectReactor.separation = seperationDict

        # check if it is mixed with other reactors
        processInvervalNames = DFconnectionMatrix.process_intervals
        reactorCol = DFconnectionMatrix[intervalName]
        pos = reactorCol != 0 # find where mixing takes place # mixed streams are in the same colunm
        mixDict = {} #preallcoation
        if sum(pos) >= 2:
            intervalsToMix = list(processInvervalNames[pos])
            specifications = list(reactorCol[pos])
            for k,specs in enumerate(specifications):
                if 'mix' in specs:
                    mixDict.update( {intervalsToMix[k]: specs} )
            #mixDict = {intervalsToMix[j]: specifications[j] for j in range(0, len(intervalsToMix))}
            #objectReactor.mix = mixDict

        # check if it is an input to other intervals as a bool
        # check if it is an input to other intervals as a bool
        boolDict = {}
        processInvervalNames = DFconnectionMatrix['process_intervals'].to_list()
        rowIndex = processInvervalNames.index(intervalName)
        intervalRow = DFconnectionMatrix.iloc[rowIndex].to_list()  # looking at the row will show to which intervals the curretn section is connencted to
        # intervalRow = DFconnectionMatrix.loc[processInvervalNames == intervalName].to_dict()
        for j, info in enumerate(intervalRow):
            if isinstance(info, str) and 'bool' in info:
                attachInterval = processInvervalNames[j]
                boolVar = 'y_' + attachInterval + '_' + intervalName
                boolDict.update({attachInterval: boolVar})

        # get the boolean variable which the reactor is dependent on
        col = reactorCol.to_list()
        boolVariable = []
        for index, infoCol in enumerate(col):
            if isinstance(infoCol,str) and 'bool' in infoCol:
                connectingInterval = processInvervalNames[index]
                boolVariable.append('y_' + connectingInterval + '_' + intervalName)
        if len(boolVariable) > 1:
            raise Exception("Currently the iterval bloks can only except a bool stream from one location, not multiple")

        # make initial interval object
        objectReactor = ReactorIntervalClass(inputs = inputsReactor, outputs = outputsReactor,  reactionEquations= equations, nameDict =nameDict,
                             mix= mixDict, utilities=utilityDict, separation=seperationDict, boolActivation= boolVariable, boolDict= boolDict)
        # put the object in the dictionary
        objectDictionary.update({intervalName:objectReactor})
    return objectDictionary

# ============================================================================================================
# Functions to update the interval objects
# ============================================================================================================
def update_reactor_equations(intervalObject,allIntervalObjectsDict,connectedIntervals):
    previousIntervalName = list(connectedIntervals.keys())[0]
    previousIntervalObject = allIntervalObjectsDict[previousIntervalName]
    newInputs4Interval = previousIntervalObject.leavingInterval
    initialInputs4Interval = intervalObject.inputs
    replacementDict = get_replacement_dict(initialInputs4Interval, newInputs4Interval)
    equationsInterval = intervalObject.equations
    allEquations = []
    for eq in equationsInterval:
        for var in replacementDict:
            newVarName = replacementDict[var]
            eq = eq.replace(var, newVarName)
        allEquations.append(eq)
    intervalObject.equations = allEquations

def update_reactor_interval(allIntervalObjectsDict,excelName):
    loc = os.getcwd()
    posAlquimia = loc.find('Alquimia')
    loc = loc[0:posAlquimia + 8]
    loc = loc + r'\excel files' + excelName
    # read excel info
    DFreactors = pd.read_excel(loc, sheet_name='Reactors')
    connectionMatrix = pd.read_excel(loc, sheet_name='ConnectionMatrix')
    reactorIntervals = DFreactors.reactor_name
    for intervalName in allIntervalObjectsDict:
        intervalObject = allIntervalObjectsDict[intervalName]
        label = intervalObject.label
        if label == 'reactor':
            connectedIntervals = get_connected_intervals(intervalName= intervalName,conectionMatrix=connectionMatrix)
            #update_reactor_equations(intervalObject,allIntervalObjectsDict,connectedIntervals)
            if len(connectedIntervals) == 1 and list(connectedIntervals.values())[0] == 1:
                update_reactor_equations(intervalObject,allIntervalObjectsDict,connectedIntervals)
            if len(connectedIntervals) == 1 and 'bool' in list(connectedIntervals.values())[0]:
                update_reactor_equations(intervalObject,allIntervalObjectsDict,connectedIntervals)



            # if no mixing seperation or split




#def makeObjectDictionaries():
