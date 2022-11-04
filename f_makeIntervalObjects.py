'''
Created on Tue sep 20 2022
Contains the classes and functions to make the process interval objects

@author: Lucas Van der Hauwaert
email: lucas.vanderhauwaert@usc.es
'''

import os
import pandas as pd
import numpy as np
from classes_intervals import InputIntervalClass, ReactorIntervalClass #, OutputIntervalClass

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
    nCol = np.shape(sumOfArrays)  #get the amount of columns
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

# def get_replacement_dict(initialVars, newVars):
#     replacementDict = {}
#     for i in initialVars:
#         for j in newVars:
#             if i in j:  # the initial variable (of excel) is always in the new name, that's how you can find wat belongs to where
#                 replacementDict.update({i:j})
#     return replacementDict
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

def update_reactor_intervals(allIntervalObjectsDict,excelName):
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

            #update_reactor_equations: current interval connected by 1 interval
            if len(connectedIntervals) == 1: # and list(connectedIntervals.values())[0] == 1: or in the case of a bool
                previousIntervalName = list(connectedIntervals.keys())[0]
                previousIntervalObject = allIntervalObjectsDict[previousIntervalName]
                newReactorInputs4Interval = previousIntervalObject.leavingInterval
                intervalObject.update_reactor_equations(newReactorInputs4Interval)

            # update_reactor_equations: current interval connected by multiple intervals by MIXING
            if len(connectedIntervals) > 1: # so here is mixing
                objectDict2mix = {nameObjConect:allIntervalObjectsDict[nameObjConect] for nameObjConect in connectedIntervals}
                intervalObject.make_mix_equations(objectDict2mix)
                newReactorInputs4Interval = intervalObject.mixingVariables
                intervalObject.update_reactor_equations(newReactorInputs4Interval)

            # wat other
            # in the case of separation !!!!!!
            # if len(connectedIntervals) == 1 and list(connectedIntervals.values())[0] == 1:
            # if len(connectedIntervals) == 1 and 'bool' in list(connectedIntervals.values())[0]:
            #     update_reactor_equations(intervalObject,allIntervalObjectsDict,connectedIntervals)



            # if no mixing seperation or split


