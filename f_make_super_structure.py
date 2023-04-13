import pandas as pd
import cobra
import cobra.io
import numpy as np
from collections import OrderedDict
import json
import pyomo.environ as pe
from f_usefull_functions import *
from f_screen_SBML import count_atom_in_formula

# from typing import List

"""Created on the 30.01.2023
@author: Lucas Van der Hauwaert
@author contact: lucas.vanderhauwaert@usc.es

All functions used to create superstructure models
"""


########################################################################################################################
########################################################################################################################
# ============================================================================================================
# formulating the superstructure
# ============================================================================================================
########################################################################################################################
########################################################################################################################

# ============================================================================================================
# Generating surrogate models
# Functions to make the interval reaction equations from SBML models
# ============================================================================================================

def get_conversion_sbml(modelLocations, substrate_exchange_rnx, product_exchange_rnx,
                        substrate2zero='Ex_S_cpd00027_ext',
                        newObjectiveReaction=None, pFBA=None, printEq=False):
    allYields_pFBA = []
    allYields_FBA = []
    objectiveBiomass = []
    allEquations = []
    modelNames = []
    for i in modelLocations:
        model = cobra.io.read_sbml_model(i)
        # make sure the right objective is set
        if newObjectiveReaction:
            model.objective = newObjectiveReaction
        # change the glucose reaction to zero
        exchange_rnx_2_zero = substrate2zero
        model.reactions.get_by_id(exchange_rnx_2_zero).bounds = 0, 0

        # change bound of new substrate to -10 mol/h/gDW
        model.reactions.get_by_id(substrate_exchange_rnx).bounds = -10, 0
        # get names of the models
        modelName = i.split("\\")[-1]
        modelName = modelName.replace(".xml", "")
        modelNames.append(modelName)
        # run pFBA
        if pFBA:  # Todo fix flux to grams c for pFBA
            pfba_solution = cobra.flux_analysis.pfba(model)
            substrate_flux = pfba_solution.fluxes[substrate_exchange_rnx]
            product_flux = pfba_solution.fluxes[product_exchange_rnx]
            yield_pFBA = -product_flux / substrate_flux  # fixed
            allYields_pFBA.append(yield_pFBA)

        else:  # else do regular FBA
            solution = model.optimize()
            FBA_substrate_flux = solution.fluxes[substrate_exchange_rnx]
            substrateMet = model.reactions.get_by_id(substrate_exchange_rnx).reactants[0]
            substrateName = substrateMet.name
            # substrateFormula = substrateMet.formula
            Csub = count_atom_in_formula(substrateMet, 'C')
            strEqlist = []

            for j in product_exchange_rnx:
                productMet = model.reactions.get_by_id(j).reactants[0]
                productName = productMet.name
                productFormula = productMet.formula
                Cprod = count_atom_in_formula(productMet, 'C')
                FBA_product_flux = solution.fluxes[j]
                FBA_yield = abs((FBA_product_flux / FBA_substrate_flux) * (Cprod * 12) / (
                        Csub * 12))  # in gramsC / grams C: 12 gCarbon/mol

                allYields_FBA.append(FBA_yield)
                strEq = '{} == {} * {}'.format(productName, FBA_yield, substrateName)
                allEquations.append(strEq)
                if printEq:
                    print(modelName)
                    print(strEq)

    return allEquations, allYields_FBA


def make_str_eq_smbl(modelName, substrate_exchange_rnx, product_exchange_rnx, equationInfo):
    # modelLocations = modelName
    loc = os.getcwd()
    posAlquimia = loc.find('Alquimia')
    loc = loc[0:posAlquimia + 8]
    # modelLocations = loc + r'\excel files\' + modelName
    modelLocations = [loc + r"\SBML models\{}".format(modelName)]

    #  make extended reaction equtions for sbml models or maybe even make them into JSON files? idk
    # would run quicker instead of having to read the xlm files each time
    equations, yields = get_conversion_sbml(modelLocations, substrate_exchange_rnx, product_exchange_rnx)

    # make abbreviations dictionary
    abbrDict = {}
    inputSbmlName = equationInfo.input_name
    inputAbrr = equationInfo.input_abrr
    abbrDict.update({inputSbmlName: inputAbrr})

    outputSbmlName = split_remove_spaces(equationInfo.output_name, ',')
    outputAbrr = split_remove_spaces(equationInfo.output_abrr, ',')
    for i, out in enumerate(outputSbmlName):
        abbrDict.update({out: outputAbrr[i]})

    allEquations = []
    for eq in equations:
        for name in abbrDict:
            if name in eq:
                eq = eq.replace(name, abbrDict[name])
        allEquations.append(eq)

    return allEquations


def make_str_eq_json(modelObject, equationInfo):
    """Reconstructs the surrogate model which is in a .json file into a equation that can be read by Pyomo
    inputs:
    model object (dict): is the unpacked .json file
    equationInfo (DF): a pandas data frame containing the information from the Excel file about the surrogate model
                        i.e., the abbreviations for the model parameters
    """
    # extract all the information
    outputs = modelObject['outputs']
    coef = modelObject['coef']
    intercept = modelObject['intercept']
    name = modelObject['name']
    lable = modelObject['lable']

    if 'waterEq' in modelObject:
        waterEq = modelObject['waterEq']
    else:
        waterEq = ''  # make an empty equation

    # if the equation is to determin the yield of the reaction then we need to know from what component to take the
    # yield of! specified in the Excel sheet 'models', column  'yield_of'
    yieldOf = equationInfo['yield_of']

    # make abbreviations dictionary
    abbrDict = {}

    # get the variable names of the inputs and outputs in the reaction and their respective abbreviations
    if not equationInfo.input_name: # in other words the input names and abbreviations stay as their originals.
        # used for cluster inputs
        key1 = list(coef.keys())[0] # just take the first key, the substrates are the same for all the products
        varNameInputs = list(coef[key1].keys())
        AbrrInputs = varNameInputs
    else:
        varNameInputs = split_remove_spaces(equationInfo.input_name, ',')
        AbrrInputs = split_remove_spaces(equationInfo.input_abrr, ',')

    varNameOutputs = split_remove_spaces(equationInfo.output_name, ',')
    varNames = varNameInputs + varNameOutputs

    AbrrOutputs = split_remove_spaces(equationInfo.output_abrr, ',')
    Abrr = AbrrInputs + AbrrOutputs

    # make sure there are as many abbreviations as well as variable inputs
    if len(varNameInputs) != len(AbrrInputs):
        raise Exception(
            "Missing an abbreviation for one of the inputs of model {}. Check out the Excel sheet 'models'".format(
                name))
    if len(varNameOutputs) != len(AbrrOutputs):
        raise Exception(
            "Missing an abbreviation for one of the outputs of model {}. Check out the Excel sheet 'models'".format(
                name))

    for i, varN in enumerate(varNames):
        abbrDict.update({varN: Abrr[i]})

    equationList = []
    print(name)
    for out in outputs:
        outAbrr = abbrDict[out]
        coefOfOutputs = coef[out]

        # update the water Eq
        waterEq = waterEq.replace(out, outAbrr)

        # preallocate the right side of the reaction equation
        yieldEq = ''
        for feature in coefOfOutputs:
            ### this for loop to replace all the full variable names with the abbreviuations
            featureAbbr = ''
            for v in varNames:
                if v in feature:
                    featureAbbr = feature.replace(v, abbrDict[v])
            if not featureAbbr:
                raise Exception('the feature {} has no abbreviation check the JSON file and the sheet models, '
                                'are names and abrr correct?'.format(feature))
            ###
            featureCoef = coefOfOutputs[feature]
            yieldEq += ' + {} * {} '.format(featureAbbr, featureCoef)

        # if the eqation is from a sbml/GEM model the yield is already calculated
        if lable == 'SBML':
            equation = '{} == {}'.format(outAbrr, yieldEq)

        # if the yield of the reaction needs to be caluculated, we need to know what input to take the yield of
        elif lable == 'yield_equation':  # the lable is then "yield_equation"
            yieldEq += ' + {}'.format(intercept[out])
            equation = '{} == ({}) * {}'.format(outAbrr, yieldEq, yieldOf)
        else:
            raise Exception("the lable given for the surrogate model '{}' must be either 'SBML' or "
                            "'yield_equation' ".format(name))
        # print(equation)
        # print('')
        equationList.append(equation)
    if waterEq:  # if the string is not empty add it to the list of equations
        equationList.append(waterEq)
    return equationList

def make_str_eq_distilation_json(modelObject, intervalName):
    """Reconstructs the surrogate model which is in a .json file into a equation that can be read by Pyomo
    inputs:
    model object (dict): is the unpacked .json file
    equationInfo (DF): a pandas data frame containing the information from the Excel file about the surrogate model
                        i.e., the abbreviations for the model parameters
    """
    # extract all the information
    inputs = modelObject['inputs']
    outputs = modelObject['outputs']
    coef = modelObject['coef']
    intercept = modelObject['intercept']
    name = modelObject['name']
    lable = modelObject['lable']

    # get the variable and the features of the regresion
    varNameInputs = inputs
    varNameOutputs = outputs
    keysCoef = list(coef.keys())  # just take the first key, the substrates are the same for all the products
    varNames = varNameInputs + varNameOutputs

    replacementDict = {}
    for var in varNames:
        newVar = var.replace(var, '{}_{}'.format(var, intervalName))
        replacementDict.update({var:newVar})

    out = varNameOutputs[0] # there should only be one output name
    outVar = replacementDict[out]
    outVar = 'energy_consumption_{}'.format(intervalName)
    coefOfOutputs = coef[out]
    # preallocate the right side of the reaction equation
    energyRequiermentEqLeft = "model.var['{}'] == ".format(outVar)
    energyRequiermentEqRight = ""
    for feature in coefOfOutputs:
        featureCoef = coefOfOutputs[feature]
        energyRequiermentEqRight += " + {} * {} ".format(feature, featureCoef)
    # add intercept
    energyRequiermentEqRight += " + {}".format(intercept[out])
    energyRequiermentEqRight = '(' + energyRequiermentEqRight + ')'
    for var in replacementDict:
        rplc = "model.var['{}']".format(replacementDict[var])
        energyRequiermentEqRight = energyRequiermentEqRight.replace(var,rplc)
    eq = energyRequiermentEqLeft + energyRequiermentEqRight

    variableList = list(replacementDict.values())
    return eq, variableList


# Created on Tue Oct 04 2022
# Contains the classes to make the process interval objects

# ============================================================================================================
# Input, reactor and output Classes
# ============================================================================================================

def make_eqation_bool_dependent(equation, booleanVariable):
    """ Adds the boolean variable to the equation
    Params:
        equation (str): string equation in pyomo format
        booleanVariable (str): boolean variable to add
    """
    equationModified = equation.replace('==', '== ( ')
    equationModified += " ) * model.boolVar['{}'] ".format(booleanVariable)

    return equationModified


def define_connect_info(connectInfo):
    """
    determines the connection keys from the connectio matrix

    parameters:
    connectInfo (int or str): a cell from the connection matrix

    returns:
    connectKey (logical): if True the interval is connected one on one
    sepKey (str): sep1 or sep2: defines which seperation stream is entering
    splitKey (str): split1 or spilt2: defines which separation stream is entering
    boolKey (str): redundant but indicates if the entering stream is dependent on a 'boolean flow'
    """
    sepKey = ''
    splitKey = ''
    boolKey = ''
    connectKey = False  # connent key referes to if there is a simple conection without seperation or splitting

    if isinstance(connectInfo, str):
        if 'sep' in connectInfo:
            sepIndex = connectInfo.find('sep')
            sepKey = connectInfo[sepIndex: sepIndex + 4]
        if 'split' in connectInfo:
            splitIndex = connectInfo.find('split')
            splitKey = connectInfo[splitIndex: splitIndex + 6]
        if 'bool' in connectInfo:
            bIndex = connectInfo.find('split')
            boolKey = connectInfo[bIndex: bIndex + 4]

    if not splitKey and not sepKey:
        connectKey = True  # thus a simple connection

    return connectKey, sepKey, splitKey, boolKey


class BooleanClass:
    def __init__(self, ExcelDict):
        """ makes the equations that regulate if a certain process is chosen or not nl: 1 == sum(boolean variables)
        the function works as followed: the DFconnectionMatrix excludes the inputs!! important!! the input boolean variables
        are regulated in the first input objected.

        returns:
            boolean variables (list): list of boolean variables
            boolean equations (lsit): list of boolean equations
            """

        # create the label for the object
        self.label = 'bool object'

        # extract the necesary dataframes
        DFIntervals = ExcelDict['input_output_DF']
        DFconnectionMatrix = ExcelDict['connection_DF']

        # find the input intervals
        posInputs = DFIntervals.input_price.to_numpy() != 0
        inputIntervalNames = DFIntervals.loc[posInputs, 'process_intervals'].to_list()

        # ------------------------------ intput boolean equations---------------------------------------------------
        # We need to find out which inputs are bound to boolean variables amongst the input variable themselves i.e.,
        # if only certain inputs can be chosen amongst multiple possible inputs in other words we need to find the
        # boolean variables from the connenction matrix (diagonals of the matrix) so the activation equation sum(y)
        # == 1 can be made.

        # there are 2 possible ways, the specific inputs can be specified in the Excel file or selected from a list in
        # the case of "input clusters" clusters are made apparent when the column components in the Excel sheet
        # input_output_intervals contains a string to a .json file

        inputBooleanVariables = []  # prealloccate
        ClusterDict = {} # prealloccate
        inputBoolEquation = "1 == " # bool amoug different input intervals
        inputBoolEquationCluster = "1 == " #bool equation for selection of substrtate in cluster of 1 interval
        equationCheck = inputBoolEquation
        clusterSwitch = False
        for i, intervalName in enumerate(inputIntervalNames):
            componentSpecification = DFIntervals.loc[posInputs, 'components'][i]
            inputClusterDict = {}  # preallocate
            if '.json' in componentSpecification: # remember it's a pandas series still
                clusterSwitch = True
            if clusterSwitch:
                jsonFile = componentSpecification
                jsonLoc = get_location(file=jsonFile)
                with open(jsonLoc) as file:
                    inputs_prices = json.load(file)
                inputNames = list(inputs_prices.keys())
                inputBooleanVariables = []  # prealloccate

                # only one input from the cluster can be chosen
                for iName in inputNames:
                    inputBoolVar = 'y_{}_{}'.format(iName,inputIntervalNames[0])
                    inputBooleanVariables.append(inputBoolVar)
                    inputBoolEquationCluster += " + " + "model.boolVar['{}']".format(inputBoolVar)
                    price = inputs_prices[iName]
                    inputClusterDict.update({iName: {'price': price, 'bool': inputBoolVar}})
                # make the over arcing dictionary
                ClusterDict.update({intervalName:inputClusterDict})

            # look if there is a boolean amoung the input variables
            inputBoolVar = DFconnectionMatrix[intervalName][intervalName]  # diagonal position of the connection matrix
            if isinstance(inputBoolVar, str):
                inputBooleanVariables.append(inputBoolVar)
                inputBoolEquation += " + " + "model.boolVar['{}']".format(inputBoolVar)

        inputBoolEquation = [inputBoolEquation]
        inputBoolEquationCluster = [inputBoolEquationCluster]
        # if where is only a single input just make a boolean equation y = 1, so it is always chosen
        # (best to acctually avoid) how could we do this?
        if equationCheck == inputBoolEquation[0]:
            inputBoolEquation = [] # just empty
        if equationCheck == inputBoolEquationCluster[0]:
            inputBoolEquationCluster = []  # just empty


        # # find out if your working with clusters
        # componentSpecification = '' # preallocate to avoid errors
        # clusterSwitch = False
        # if len(inputIntervalNames) == 1:
        #     componentSpecification = DFIntervals.loc[posInputs, 'components'][0]
        #     if '.json' in componentSpecification: # remember it's a pandas series still
        #         clusterSwitch = True
        #
        # # if there is a cluster of inputs to chose from, make the boolean equation based on the input.json file given in
        # # the components column
        # inputClusterDict = {} # preallocate
        # if clusterSwitch:
        #     jsonFile = componentSpecification
        #     jsonLoc = get_location(file=jsonFile)
        #     with open(jsonLoc) as file:
        #         inputs_prices = json.load(file)
        #     inputNames = list(inputs_prices.keys())
        #     inputBooleanVariables = []  # prealloccate
        #
        #     inputBoolEquation = "1 == "
        #     for iName in inputNames:
        #         inputBoolVar = 'y_{}_{}'.format(iName,inputIntervalNames[0])
        #         inputBooleanVariables.append(inputBoolVar)
        #         inputBoolEquation += " + " + "model.boolVar['{}']".format(inputBoolVar)
        #         price = inputs_prices[iName]
        #         inputClusterDict.update({iName: {'price': price, 'bool': inputBoolVar}})
        #
        # else:
        #     inputBooleanVariables = []  # prealloccate
        #     inputBoolEquation = "1 == "
        #     equationCheck = inputBoolEquation
        #     for i, intervalName in enumerate(inputIntervalNames):
        #         inputBoolVar = DFconnectionMatrix[intervalName][intervalName]  # diagonal position of the connention matrix
        #         if isinstance(inputBoolVar, str):
        #             inputBooleanVariables.append(inputBoolVar)
        #             inputBoolEquation += " + " + "model.boolVar['{}']".format(inputBoolVar)

            # # if where is only a single input just make a boolean equation y = 1, so it is always chosen
            # # (best to acctually avoid) how could we do this?
            # if equationCheck == inputBoolEquation:
            #     inputBoolVar = 'y'
            #     inputBooleanVariables.append(inputBoolVar)
            #     inputBoolEquation += " + " + "model.boolVar['{}']".format(inputBoolVar)

        # if there exist bool inputs, make a unique list
        if len(inputBooleanVariables) != len(set(inputBooleanVariables)):
            # if not a unique list raise exception
            raise Exception("The boolean variables for the inputs are not unique, check the diagonals of the "
                            "input intervals of the connection matrix ")


        # ------------------------------ other interval boolean equations-----------------------------------------------
        # main idea:
        # 1) loop over the rows of the DF.
        # 2) count the columns of this row which are not zero in  sequence! (save the interval names in a list)
        # 3) the longest list is the list of intervals dependant of the boolean variable
        # 4) the boolean variable can be found on the diagonal (i.e. with the same row and column index)
        # 5) make the boolean equation
        # 6) delete the columns that contain the sequence
        # 7) restart at step 1 till the DF is empty

        # make a list of all the intervals excluding the inputs
        DFconnectionMatrix = DFconnectionMatrix.drop(labels=inputIntervalNames, axis=1)
        processIntervalnames = list(DFconnectionMatrix.index)[len(inputIntervalNames):]

        # prun the Dataframe, anything that does not have a bool label on the diagonal can be dropped
        toDrop = []
        for i in processIntervalnames:
            if not isinstance(DFconnectionMatrix[i][i], str):
                toDrop.append(i)
        DFconnectionMatrix = DFconnectionMatrix.drop(labels=toDrop, axis=1)

        # start the while loop, aslong as the dataframe is not empty continue
        switch = True
        equationsSumOfBools = []
        booleanVariables = []
        iteration = 0
        while switch:
            iteration += 1
            saveDict = {}
            for index, row in DFconnectionMatrix.iterrows():
                intervalNames = []
                for indexCol, element in row.items():
                    if isinstance(element, str) or element != 0:  # so comes from a separation or just connected by '1'
                        intervalNames.append(indexCol)
                    else:
                        break
                saveDict.update({index: intervalNames})

            # find the key in the saveDict that has the longest list
            key_max_sequential = max(saveDict, key=lambda k: len(saveDict[k]))

            # create the equation
            eq = '1 == '
            for interval in saveDict[key_max_sequential]:
                boolVar = DFconnectionMatrix[interval][interval]
                booleanVariables.append(boolVar)
                eq += "+ model.boolVar['{}'] ".format(boolVar)
            equationsSumOfBools.append(eq)

            # now we need to drop the colums in the original dataframe that already form 1 set of equations
            # the colunms that need to be droped are in saveDict[key_max_sequential]1
            cols2drop = saveDict[key_max_sequential]
            DFconnectionMatrix = DFconnectionMatrix.drop(labels=cols2drop, axis=1)

            # once the DF is empty we can stop the while loop
            if DFconnectionMatrix.empty:
                switch = False

            # if iteration 50 is hit, then stop the program, somthing wrong is going on
            if iteration > 50:
                raise Exception(
                    "50 iteration have passed to try and make the boolean equations, check to see if the "
                    "connection matrix is correctly formulated or check the function 'make_boolean_equations'")

        # while loop has stopped
        # check that all the boolean variables have unique values
        uniqueSet = set(booleanVariables)
        if len(uniqueSet) != len(booleanVariables):
            raise Exception("the boolean variables are not all unique, check the diagonal of the connection matrix")

        # add the cluster dictionary to the object, empty if there is none
        self.clusterDict = ClusterDict

        # add the boolean variables to the allVariable dictionary
        allBoolVariables = booleanVariables + inputBooleanVariables
        self.allVariables = {'continuous': [],
                             'boolean': allBoolVariables,
                             'fraction': []}

        # define it's boundries
        booleanBounds = {}
        for i in allBoolVariables:
            booleanBounds.update({i: 'bool'})
        self.boundaries = booleanBounds

        # add the equations to the pyomoEquations object
        self.pyomoEquations = equationsSumOfBools + inputBoolEquation + inputBoolEquationCluster


class InputIntervalClass:
    def __init__(self, inputName, compositionDict, inputPrice, boundryInputVar,clusterDict,boolDict=None,
                 split=None, separationDict=None, booleanVariable=None):
        if separationDict is None:
            separationDict = {}
        if split is None:
            split = []
        if boolDict is None:
            boolDict = {}
        if booleanVariable is None:
            booleanVariable = ''

        # declare (preallocate) empty pyomo equations list
        pyomoEq = []

        # declare input interval name
        self.label = 'input'
        self.inputName = inputName.upper()  # put in capitals
        self.inputPrice = inputPrice
        addOn4Variables = '_' + inputName.lower()

        # initial Boundry dictionary of all the variables (add where aprropriate )
        self.allBoundries = {}

        # change the composition names in the dictionary
        compositionDictNew = {}
        initialCompositionNames = []
        if len(compositionDict) == 1 and '.json' in list(compositionDict.keys())[0]:
            compositionDictNew = compositionDict # keep the json name in the dictionary
        else:
            for i in compositionDict:
                compositionDictNew.update({i + addOn4Variables: compositionDict[i]})
                initialCompositionNames.append(i)
        # self.initialCompositionNames = initialCompositionNames
        self.compositionDict = compositionDictNew

        # error if you choose a wrong name
        for i in initialCompositionNames:  # maybe don't place the warning here
            if i in inputName:
                raise Exception("the component {} is in interval name {}, change the component name of the reactor to "
                                "avoid conflict with the equations".format(inputName, i))

        # make the component equations as string equations
        # --------------------------------------------------------------------------- cluster section
        try:
            intervalCluster = clusterDict[inputName]
        except:
            intervalCluster = {}

        self.clusterDict = intervalCluster

        if intervalCluster: # i.e., working with clusters, if not clusterDict is empty
            for key in intervalCluster:
                # activation equations
                inputVariable = key
                booleanVariable = intervalCluster[key]['bool']
                activationEqPyoUB = "model.var['{}'] <= {} * model.boolVar['{}'] ".format(inputVariable,
                                                                                          boundryInputVar[1],
                                                                                          booleanVariable)

                activationEqPyoLB = "{} * model.boolVar['{}']  <= model.var['{}'] ".format(boundryInputVar[0],
                                                                                           booleanVariable,
                                                                                           inputVariable)
                # add to the list of equations
                pyomoEq.append(activationEqPyoLB)
                pyomoEq.append(activationEqPyoUB)

            # declare component variables
            componentVariables = list(intervalCluster.keys())
            continuousVariables = componentVariables
            self.leavingInterval = componentVariables
            #if compositionDictNew

        else: # --------------------------------------------------------------------  non- cluster section
            for component in compositionDictNew:
                eqPy = "model.var['{}'] == {} * model.var['{}']".format(component, self.compositionDict[component],
                                                                        self.inputName)
                pyomoEq.append(eqPy)

            componentVariables = list(compositionDictNew.keys())
            continuousVariables = componentVariables + [self.inputName]

            self.boolDict = boolDict  # necesary?
            self.split = split  # necesary?
            self.separationDict = separationDict  # necesary?
            self.leavingInterval = componentVariables


            if booleanVariable:
                activationEqPyoUB = "model.var['{}'] <= {} * model.boolVar['{}'] ".format(self.inputName,
                                                                                          boundryInputVar[1],
                                                                                          booleanVariable)
                activationEqPyoLB = "{} * model.boolVar['{}']  <= model.var['{}'] ".format(boundryInputVar[0],
                                                                                           booleanVariable, self.inputName)
                pyomoEq.append(activationEqPyoLB)
                pyomoEq.append(activationEqPyoUB)



        # ----------------- check ------------------------------------------------------------- point charly
        #  put all VARIABLES that pyomo needs to declare here
        self.allVariables = {'continuous': continuousVariables,  # continous variables
                             'boolean': [],  # boolean variables
                             'fraction': []}  # fraction variables

        # make BOUNDARIES for all variables
        # interval variables are bounded by inequality equations bounds (see above)
        boundaryDict = {self.inputName: [0, None]}
        # all other variables are positiveReals
        for i in componentVariables:
            boundaryDict.update({i: 'positiveReals'})
        self.boundaries = boundaryDict

        # put all EQUATIONS that pyomo needs to declare here
        self.pyomoEquations = pyomoEq


class ProcessIntervalClass:
    def __init__(self, inputs, outputs, reactionEquations, boundryInputVar, nameDict, mix=None,
                 utilities=None, energyUtility=None, booleanVariable=None, splitList=None, separationDict=None,
                 operationalVariablesDict=None, energyConsumption=0):

        if utilities is None:
            utilities = {}
        if separationDict is None:
            separationDict = {}
        if mix is None:
            mix = []
        if booleanVariable is None:
            booleanVariable = ''
        if operationalVariablesDict is None:
            operationalVariablesDict = {}
        if energyUtility is None:
            energyUtility = {}

        self.label = 'process_interval'
        self.booleanVariable = booleanVariable
        self.operationalVariablesDict = operationalVariablesDict

        # further the lable of the interval i.e.:reactor or seperator
        if reactionEquations:
            self.intervalType = 'reactor'  # lable
        elif reactionEquations is None and separationDict:
            self.intervalType = 'separator'  # lable

        # attributes you want to be able to call from the object (some are probably redundant)
        self.nameDict = nameDict
        self.inputs = inputs  # original unmodified names
        self.outputs = outputs  # original unmodified names
        self.initialCompositionNames = inputs + outputs
        self.mix = mix  # found by the Excel file (don't need to specify in this script)
        self.utilities = utilities  # consists of a dictionary {nameUtilty: [bounds]}
        self.separationDict = separationDict  # dictionary defining separation fractions of each component

        # error if you choose a wrong name
        intervalName = list(nameDict.keys())[0]
        self.intervalName = intervalName
        for i in self.initialCompositionNames:  # maybe don't place the warning here
            if i in intervalName:
                raise Exception("the component {} is in interval name {}, change the component name of the reactor to "
                                "avoid conflict with the equations")

        # interval Variable name
        originalIntervalName = list(self.nameDict.keys())[0]
        intervalVariable = originalIntervalName.upper()
        addOn4Variables = '_' + originalIntervalName

        # declare (preallocate) empty pyomo equations list
        pyomoEq = []

        # reactor equations
        reactionVariablesOutput = []  # preallocate to avoid error
        if reactionEquations != None:
            reactorEquations, reactionVariablesOutput, helpingDict = \
                self.make_reaction_equations(reactionEquations=reactionEquations,
                                             booleanVariable=booleanVariable,
                                             intervalVariable=intervalVariable)
            # pyomoEq += reactorEquations # added during the update fuction
        else:
            # If there is a dictionary for the separation equations but no reaction, then there is no helpingDict.
            # In other words now the separation equations need to be updated according to the incomming flow...
            # this is indicated the in function update_interval_equations by the self.intervalType label!!
            helpingDict = {}

        # separation equations
        separationVariables = []  # preallocate to avoid error
        if separationDict:  # if there is separation make the separation equations
            separationEquationsPyomo, separationVariables = self.make_separation_equations(separationDict, helpingDict,
                                                                                           booleanVariable=booleanVariable)
            if helpingDict:  # if the helping dict does noit exist the separation equations are added during the update function
                pyomoEq += separationEquationsPyomo  # otherwise they can be added strait away

        # spliting equations
        splitComponentVariables = []  # preallocate to avoid error
        splitFractionVariables = []  # preallocate to avoid error
        if splitList:
            splittingEquations, splitComponentVariables, splitFractionVariables = self.make_split_equations(splitList,
                                                                                                            addOn4Variables,
                                                                                                            booleanVariable=booleanVariable)
            pyomoEq += splittingEquations

        # mixing equations
        # see def make_mixing_equations(), these equations are made in the update when it is known
        # to where each stream is going to and hence which streams are mixed

        # utility (chemical) equations
        # are made in the update_function where reactor equations are also completed

        # utility (energy) equations
        # are made in the update_function because we need to know what is entering the interval
        # give the varibles needed to the object
        self.utilityEnergy = energyUtility

        # define wat is leaving the reactor
        # can either be from the reactor, the separation process or the spliting
        # maybe look at empty lists of the variables instead? check that out
        if not splitList and not separationDict:
            self.leavingInterval = reactionVariablesOutput
        elif separationDict and not splitList:
            self.leavingInterval = separationVariables
        elif splitList and not separationDict:
            self.leavingInterval = splitComponentVariables
        elif splitList and separationDict:
            VariablesGoingOutSep = separationVariables.copy()
            toRemove = []
            for i in splitList:
                separationStream = ''
                if 'sep' in i:
                    indexSep = i.find('sep')
                    separationStream = i[
                                       indexSep:indexSep + 4]  # the separation stream that is going te get knocked out
                for j in VariablesGoingOutSep:
                    if separationStream in j:
                        toRemove.append(j)

            # remove unwanted variables
            for i in toRemove:
                VariablesGoingOutSep.remove(i)

            self.leavingInterval = VariablesGoingOutSep + splitComponentVariables

        # put all self.VARIABLES that pyomo needs to declare here
        self.reactionVariablesOutput = reactionVariablesOutput
        # reactionVariablesInputs can be found in class function: update_reactor_equations
        self.intervalVariable = intervalVariable

        # define the bounds of the variables
        boundaryDict = {}  # intervals have bounds
        boundaryDict.update({intervalVariable: 'positiveReals'})
        for i in reactionVariablesOutput:
            boundaryDict.update({i: 'positiveReals'})

        for i in separationVariables:
            boundaryDict.update({i: 'positiveReals'})

        for i in boundryInputVar:  # this is for when you want to add a specific bound to a reaction variable SEE EXCEL
            boundaryDict[i] = boundryInputVar[i]
        self.boundryInputVar = boundryInputVar

        for i in splitComponentVariables:
            boundaryDict.update({i: 'positiveReals'})

        for i in splitFractionVariables:
            boundaryDict.update({i: 'fraction'})

        self.boundaries = boundaryDict

        # make a list with all the variables
        # continuousVariables = [self.intervalVariable] + self.reactionVariablesOutput + separationVariables + splitComponentVariables # self.separationVariables
        continuousVariables = self.reactionVariablesOutput + separationVariables + splitComponentVariables  # self.separationVariables

        # booleanVariables = self.boolVariables
        # + self.activationVariable don't need to add this, already in the previous interval under boolVariables
        self.allVariables = {'continuous': continuousVariables,
                             'boolean': [],
                             'fraction': splitFractionVariables}
        # add fraction variables

        # make a list with all the equations
        # self.allEquations = self.separationEquations + self.eqSumOfBools + self.boolActivationEquations + self.totalMassEquation
        self.pyomoEquations = pyomoEq

    def make_reaction_equations(self, reactionEquations, intervalVariable, booleanVariable=None):
        """ function that creates the (preliminary) equations of the reactions that take place in an interval.
        these equations do not have the input variable filled in. In the function update_interval_equations the
        reaction equations are updated with the correct inputs variables

        Parameters:
                reactionEquations (str or int) reference to where the equations are stored
                booleanVariable (str) the string referening to the boolean variable
                intervalVariable (str) name of the interval variable

        Returns:
               allEquationsPyomo (list) list of reqction equations
               reactionVariablesOut (list) list of variables after the reaction
               helpingDict (dict) a dictionary helping to make the separation equations

            """

        originalIntervalName = list(self.nameDict.keys())[0]
        addOn4Variables = '_' + originalIntervalName

        if isinstance(reactionEquations, str):  # failsafe if you forget to code the string reactor expression as a list
            reactionEquations = [reactionEquations]  # make it a list

        ouputs2change = self.outputs
        ReactorEquationsPyomo = []
        reactionVariablesOutput = []
        helpingDict = {}
        for eq in reactionEquations:
            eqPyo = eq
            for out in ouputs2change:
                newOutputName = out + '{}'.format(addOn4Variables)
                eq = eq.replace(out, newOutputName)
                # pyomo version
                newOutputNamePyo = "model.var['{}']".format(newOutputName)
                eqPyo = eqPyo.replace(out, newOutputNamePyo)

                if newOutputName not in reactionVariablesOutput:
                    reactionVariablesOutput.append(newOutputName)
                    helpingDict.update({out: newOutputName})  # helpìng dictionary for the serparation equations

            if booleanVariable:
                eqPyo = make_eqation_bool_dependent(equation=eqPyo, booleanVariable=booleanVariable)

            ReactorEquationsPyomo.append(eqPyo)

        # place all the equations in the object
        self.reactionEquations = ReactorEquationsPyomo

        # mass equations (of the outputs from the reaction equations )
        eqMassInterval = intervalVariable + " == "
        eqPyoMass = "model.var['{}'] == ".format(intervalVariable)
        for out in reactionVariablesOutput:
            eqMassInterval += " + " + out
            eqPyoMass += " + " + "model.var['{}']".format(out)  # pyomo version

        # decalre all equations and pass them on
        self.totalMassEquation = [eqMassInterval]  # TODO mass eq. not added to pyomo for the moment, necesary?
        self.reactionEquationsPyomo = ReactorEquationsPyomo
        allEquationsPyomo = [eqMassInterval] + ReactorEquationsPyomo

        # old activation equations (keep in comments for the moment)
        # bool activation constraints
        # has been de-activated. bool variables are now in the reaction equations (solver is more efficient that way)
        # boolActivationEquations = []
        # if boolActivation: # if there is an activation constraint
        #     bounds = nameDict[originalIntervalName]
        #     lowerActivationEq = "{} * {} <= {}".format(boolActivation[0],bounds[0],intervalVariable)
        #     upperActivationEq = "{} <= {} * {}".format(intervalVariable,boolActivation[0], bounds[1])
        #     boolActivationEquations.append(lowerActivationEq)
        #     boolActivationEquations.append(upperActivationEq)
        #
        #     # pyomo version
        #     lowerActivationEqPyo = "model.boolVar['{}'] * {} <= model.var['{}']".format(boolActivation[0], bounds[0], intervalVariable)
        #     upperActivationEqPyo = "model.var['{}'] <= model.boolVar['{}'] * {}".format(intervalVariable, boolActivation[0], bounds[1])
        #     pyomoEq.append(lowerActivationEqPyo)
        #     pyomoEq.append(upperActivationEqPyo)
        # self.boolActivationEquations = boolActivationEquations

        return allEquationsPyomo, reactionVariablesOutput, helpingDict

    def make_separation_equations(self, separationDict, helpingDict, booleanVariable=None):
        """ function that creates the separation equations of the process interval

        Parameters:
                separationDict (dict): list of components to separate
                helpingDict (dict): dictionary to help create the variables

        Returns:
               saparationEquations (list): list of separartion equations
               seperationVariables (list): list of variables
            """
        separationEquationsPyomo = []
        separationVariables = []
        for sep in separationDict:  # if it is empty it should not loop nmrly
            for componentSep in separationDict[sep]:
                # create the separation variable that is leaving
                sepVar = componentSep + '_' + sep
                separationVariables.append(sepVar)
                # sepVarPyo = "model.var['{}']".format(sepVar)

                if helpingDict:
                    var = helpingDict[componentSep]  # helpìng dictionary to get the right variable
                    # pyomo equations
                    eqSepPyo = "model.var['{}'] == {} * model.var['{}']".format(sepVar,
                                                                                separationDict[sep][componentSep], var)

                else:  # else if no helping dict, it will get updated in the update function
                    var = componentSep
                    eqSepPyo = "model.var['{}'] == {} * {}".format(sepVar, separationDict[sep][componentSep],
                                                                   var)
                if booleanVariable:  # add boolean variable if there are any
                    eqSepPyo = eqSepPyo.replace('==', '== ( ')
                    eqSepPyo += " ) * model.boolVar['{}'] ".format(booleanVariable)

                separationEquationsPyomo.append(eqSepPyo)

        self.separationEquations = separationEquationsPyomo
        self.separationVariables = separationVariables

        return separationEquationsPyomo, separationVariables

    def make_split_equations(self, splitList, addOn4Variables, booleanVariable=None):
        """ function that creates the equations of the split streams

        Params:
                splitList (list): list of components to split
                addOn4Variables (str): name of the add-on for the variables

        Returns:
               splittingEquations (list): list of reaction equations
               splitComponentVariables (list): list of variables
            """

        # Error messaging if there are not 2 stream to split to
        if len(splitList) != 1 and len(splitList) != 2:
            raise Exception('look at reactor row {} in the connection matrix, it looks like you are missing a split '
                            'statement'.format(addOn4Variables))

        # preallocate variables and equation lists
        splitFractionVariables = []
        splitComponentVariables = []
        splittingEquations = []

        # make split equations
        for splitStream in splitList:
            splitFractionVar = 'split_fraction_{}'.format(splitStream)
            splitFractionVariables.append(splitFractionVar)
            for component in self.outputs:  # the self.outputs are the original names given to the outputs of the reactor
                # make and get variables
                split1 = component + '_' + splitStream + '_' + 'split1'
                split2 = component + '_' + splitStream + '_' + 'split2'
                splitComponentVariables.append(split1)
                splitComponentVariables.append(split2)
                component2split = component + '_' + splitStream

                # make equations
                eqSplit1Pyo = "model.var['{}'] == model.fractionVar['{}'] * model.var['{}']".format(split1,
                                                                                                    splitFractionVar,
                                                                                                    component2split)
                eqSplit2Pyo = "model.var['{}'] == (1 - model.fractionVar['{}']) * model.var['{}']".format(split2,
                                                                                                          splitFractionVar,
                                                                                                          component2split)

                # add boolean variable if there are any
                if booleanVariable:
                    eqSplit1Pyo = eqSplit1Pyo.replace('==', '== ( ')
                    eqSplit2Pyo = eqSplit2Pyo.replace('==', '== ( ')
                    eqSplit1Pyo += " ) * model.boolVar['{}'] ".format(booleanVariable)
                    eqSplit2Pyo += " ) * model.boolVar['{}'] ".format(booleanVariable)

                # add equations to the lsit
                splittingEquations.append(eqSplit1Pyo)
                splittingEquations.append(eqSplit2Pyo)

        self.splitEquations = splittingEquations
        return splittingEquations, splitComponentVariables, splitFractionVariables

    def make_mix_equations(self, objects2mix):
        """
        makes the mixing equations based on all the intervals that are seen to be mixed in the connection interval

        parameters
        objects2mix (dict): dictionary containing the intervall classes which are to be mixed with the current one

        returns
        the pyomo equations
        """

        # depenent on bool?
        booleanVariable = self.booleanVariable

        # mixEquations: list[str] = []
        mixEquations = []
        initialInputNames = self.inputs

        # intervalNames2Mix = objects2mix.keys()

        # find the leaving variables
        leavingVars = []
        for objName in objects2mix:
            obj = objects2mix[objName][0]  # first element in tuple is the object
            connectInfo = objects2mix[objName][1]  # second element in tuple is connection info
            reactorKey, sepKey, splitKey, boolKey = define_connect_info(connectInfo)

            # if there is a separation process
            if isinstance(connectInfo, str) and sepKey or splitKey:
                allSepVars = obj.leavingInterval
                for var in allSepVars:
                    if sepKey and splitKey:
                        if sepKey in var and splitKey in var:
                            leavingVars.append(var)

                    elif sepKey and not splitKey:
                        if sepKey in var:
                            leavingVars.append(var)

                    elif splitKey and not sepKey:
                        if splitKey in var:
                            leavingVars.append(var)
            # otherwise just add the leaving variables
            else:
                leavingVars += obj.leavingInterval

        # make the equations
        mixingVariables = []
        eqMixPyo2Add = []
        intervalName = list(self.nameDict.keys())[0]
        for i, ins in enumerate(initialInputNames):
            mixVar = "{}_{}_mix".format(ins, intervalName)
            eqMix = mixVar + " == "
            startMixEq = eqMix

            # pyomo equation
            mixVarPyo = "model.var['{}']".format(mixVar)
            eqMixPyo = mixVarPyo + " == "

            # startMixEqPyo = eqMixPyo
            for lvar in leavingVars:
                if ins in lvar:
                    eqMix += " + " + lvar
                    eqMixPyo += " + " + "model.var['{}']".format(lvar)

            # For example in the case of pH this does not come from the previous interval!!
            # so the variable can stay as it is and no extra equations needs to be added, hence if eqMix != startMixEq:
            if eqMix != startMixEq:
                if booleanVariable:  # check if there is a dependance on a boolean variable
                    eqMixPyo = make_eqation_bool_dependent(equation=eqMixPyo, booleanVariable=booleanVariable)
                mixEquations.append(eqMix)
                mixingVariables.append(mixVar)
                eqMixPyo2Add.append(eqMixPyo)

        # # total flow going into the interval after mixing
        # totalMixVarible = "{}_total_mix".format(intervalName)
        # totalMixEqationPyomo = "model.var['{}'] == ".format(totalMixVarible)
        # for mixVars in mixingVariables:
        #     totalMixEqationPyomo += " + model.var['{}']".format(mixVars)
        #
        # mixingVariables.append(totalMixVarible)  # add to the list of variables
        # eqMixPyo2Add.append(totalMixEqationPyomo) # add to the list of equations

        # pass the varibales to the object and update the boundry dictionary
        self.mixingVariables = mixingVariables
        for i in mixingVariables:
            self.boundaries.update({i: (0, None)})

        ''' 
        # now the reaction equations need to be updated because
        # they're the mix variables are now the inputs of the reaction equations!
        # see function update_interval_equations & class function update_reactor_equations
        '''

        # add the variables and equations to the allVar/equations object
        # self.allEquations += mixEquations
        self.allVariables['continuous'] += mixingVariables
        self.pyomoEquations += eqMixPyo2Add
        self.mixEquations = eqMixPyo2Add

    def make_incoming_massbalance_equation(self, enteringVariables):
        """
        make the mass balance of the in coming components of an interval before the reaction phase.
        sum of the mixed components or sum of in coming components

        parameters:
        incoming_components (list):

        return: updated interval object
        """

        intervalName = self.intervalName
        # total flow going into the reactor stage of the interval (after mixing if there is any)
        enteringMassVarible = "{}_mass_entering".format(intervalName)
        enteringMassEqationPyomo = "model.var['{}'] == ".format(enteringMassVarible)

        for enteringVars in enteringVariables:
            enteringMassEqationPyomo += " + model.var['{}']".format(enteringVars)

        if self.booleanVariable:
            enteringMassEqationPyomo = make_eqation_bool_dependent(equation=enteringMassEqationPyomo,
                                                                   booleanVariable=self.booleanVariable)

        # add to the list of equations + variables and update the boundry dictionary
        self.pyomoEquations += [enteringMassEqationPyomo]
        self.incomingFlowEquation = [enteringMassEqationPyomo]
        self.incomingFlowVariable = enteringMassVarible
        self.allVariables['continuous'] += [enteringMassVarible]
        self.boundaries.update({enteringMassVarible: (0, None)})
        self.enteringVariables = enteringVariables

    def update_interval_equations(self, newInputs4Interval):
        """
        fills in the uncompleted reaction equations or separation equations
        I.e., the inputs of the reaction or separation processes

        REMARK: # the operational variables are not dependent on the previous interval, so you could replace these
        variable in the make_reactor_equation function as well. Anyway we'll replace the operational variable in this
        function (see also make_replacement_dict function)

        Params:
            newInputs4Interval (list): is a list of the input variables entering the process interval
            i.e., after transformation from the previous interval

        Returns: Updated reaction or separation equations
        """

        intervalType = self.intervalType
        name = list(self.nameDict.keys())[0]  # for debugging porposes
        # print(name)

        replacementDict, boundsDict = self.get_replacement_dict(newInputs4Interval)

        # determine which equations need to be updated
        if intervalType == 'reactor':
            equationsInterval = self.reactionEquations  # the reactor equations

        else:  # so intervalType == 'separator':
            equationsInterval = self.separationEquations  # the separation equations

        allEquations = []
        for eq in equationsInterval:
            for var in replacementDict:
                newVarName = replacementDict[var]
                # the output of the equation is already in pyomo format, to avoid conflict split
                # the eqution left and right of the '==', The right side is the part that needs to be updated
                posEqualSign = eq.find('==')
                leftEquation = eq[0:posEqualSign]
                rightEquation = eq[posEqualSign:]  # slice till the end [position:]
                rightEquation = rightEquation.replace(var, "model.var['{}']".format(newVarName))
                eq = leftEquation + rightEquation
            allEquations.append(eq)

        # update the reactor/separation equations to the object
        if intervalType == 'reactor':
            self.reactionEquations = allEquations
        else:  # so intervalType == 'separator':
            self.separationEquations = allEquations

        # add the variables and equations to the allVariables/pyomoEquations object
        reactionVariablesInputs = list(replacementDict.values())
        self.reactionVariablesInputs = reactionVariablesInputs
        self.allVariables[
            'continuous'] += self.reactionVariablesInputs  # + [enteringMassVarible] # add to the list of variables
        self.pyomoEquations += allEquations

        # add variables to boundary dictionary
        self.boundaries.update(boundsDict)

        # replace mass variables which have specific bounds
        boundryInputVar = self.boundryInputVar
        for massVar in boundryInputVar:  # this is for when you want to add a specifice bound to a reaction variable SEE EXCEL
            self.boundaries[massVar] = boundryInputVar[massVar]

    def make_utility_equations(self):
        """
        makes the equations for the flow of the chemical utility and the cost of said utlity
        """
        booleanVariable = self.booleanVariable
        reactorName = list(self.nameDict.keys())[0]  # get reactor name
        utilities = self.utilities  # get the specified utilities
        utlityCostVariable = 'cost_utility_chemicals_{}'.format(
            reactorName)  # cost variable (can be defined before the loop)

        # add variable to list
        self.allVariables['continuous'] += [utlityCostVariable]  # add to the lsit of variables
        # specify the bounds of the variable
        self.boundaries.update({utlityCostVariable: (0, None)})

        # preallocate the cost equation
        utilityCostEqPyomo = "model.var['{}'] ==  ".format(utlityCostVariable)
        # first get incoming mass from the variable 'total_mix' variable (needed to calculate incoming utility stream)
        incomingFlowVariable = self.incomingFlowVariable

        # loop over al the utilities
        massEquations = []
        for ut in utilities:
            utilityName = ut
            utilityParameter = utilities[ut]['parameter']
            utilityCost = utilities[ut]['cost']
            allUtilityVariables = []  # preallocation
            # make the variable name
            utilityVariable = '{}_{}'.format(utilityName, reactorName)
            allUtilityVariables.append(utilityVariable)

            # decalre variable in list
            self.allVariables['continuous'] += allUtilityVariables
            self.boundaries.update({utilityVariable: (0, None)})

            # declare the utility equations i.e., mass utility == parameter_ut * incoming flow
            utilityMassEqPyomo = "model.var['{}'] == {} * model.var['{}']".format(utilityVariable, utilityParameter,
                                                                                  incomingFlowVariable)

            # should be zero if the interval is not chosen
            if booleanVariable:
                utilityMassEqPyomo = utilityMassEqPyomo.replace('==', '== ( ')
                utilityMassEqPyomo += " ) * model.boolVar['{}']".format(booleanVariable)

            # add the equation to the list
            massEquations.append(utilityMassEqPyomo)

            # cost equation
            utilityCostEqPyomo += "+ (model.var['{}'] * {})".format(utilityVariable, utilityCost)

        # add the cost equation to pyomo equation list
        self.pyomoEquations += massEquations + [utilityCostEqPyomo]
        self.utilityEquations = massEquations + [utilityCostEqPyomo]
        self.utilityCostVariable = utlityCostVariable
        # add the inflow equation as well as it relates to the utility equations

    def make_energy_utility_equations(self):
        """ Makes the equations for the energy consumption of an interval

        Either the energyConsumptionParameter is a parameter (float) in kJ/kg_feed
        or energyConsumptionParameter is a json file containing the equations for a specific separation process

        """
        energyConsumptionParameter = self.utilityEnergy['consumption_parameter']
        energyPrice = self.utilityEnergy['price']

        allEnergyEquation = []
        allEnergyVars = []
        energyVar = "energy_consumption_{}".format(self.intervalName)
        if isinstance(energyConsumptionParameter, str) and 'json' in energyConsumptionParameter:
            jsonFile = energyConsumptionParameter
            jsonLoc = get_location(file=jsonFile)
            try:
                with open(jsonLoc) as file:
                    energyConsumptionObject = json.load(file)
            except:
                raise Exception('check the spelling of the json file, of the energy model for the '
                                'interval {}'.format(self.intervalName))

            # if the separation interval is a distilation
            if energyConsumptionObject["lable"] == "short_cut":
                eqList = energyConsumptionObject['equations']
                varList = energyConsumptionObject['variables']
                name = energyConsumptionObject['name']

                # get the variables you need
                # light key var
                lightKey = energyConsumptionObject['light_key']
                enteringVars = self.enteringVariables
                lightKeyVar = [s for s in enteringVars if lightKey in s][0]

                # entering feed var
                Feed_Var = self.incomingFlowVariable  # in kg/h

                # get the seperation dict to find the desired composition in the distillate and bottom streams
                separationDict = self.separationDict
                fractionArry = []
                for key in separationDict:
                    streamComposition = separationDict[key]  # [s for s in enteringVars if lightKey in s][0]
                    fractionArry.append(streamComposition[lightKey])

                x_D_var = max(fractionArry)  # fraction of the light key in the distillate is by convention the maximum
                x_B_var = min(fractionArry)  # fraction of the light key in the bottom is by convention the minimum

                # build equations for the fraction of the entering light key in the feed
                x_F_var = 'x_F_{}'.format(self.intervalName)  # feed variable
                x_F_eq = " model.var['{}'] == (model.var['{}'] * 0 / model.var['{}']) ".format(x_F_var, lightKeyVar,
                                                                                               Feed_Var)  # in mass %

                # build the equations of the shortcut method
                equationsShortcut = []
                for eq in eqList:
                    eq = eq.replace('Feed', "model.var['{}']".format(Feed_Var))
                    eq = eq.replace('x_F', "model.var['{}']".format(x_F_var))
                    eq = eq.replace('x_D', str(x_D_var))
                    eq = eq.replace('x_B', str(x_B_var))
                    eq = eq.replace('Q_tot_{}'.format(name), "energy_consumption_{}".format(self.intervalName))
                    equationsShortcut.append(eq)

                # add all the equations and variable to the 'collecting list'
                allEnergyEquation += [x_F_eq] + equationsShortcut
                allEnergyVars += [x_F_var] + [energyVar] + varList
                for var in allEnergyVars:
                    self.boundaries.update({var: (1e-6, None)})  # make it so that these variables can not hit zero
                    # change the feed var so there is no possiblity to dived by zero in the disitillation equations
                self.boundaries.update({Feed_Var: (1e-6, None)})

            elif energyConsumptionObject["lable"] == "Distillation_Regresion":
                # infoDict = {'input_name': [], 'yield_of':[], }
                # DFinfo = pd.DataFrame(infoDict)
                eq, variables = make_str_eq_distilation_json(modelObject=energyConsumptionObject, intervalName=self.intervalName)
                eq += " * model.var['{}']".format(self.incomingFlowVariable)
                allEnergyEquation.append(eq)
                # make the equation that defines the feed composition of the light key
                lightKey = energyConsumptionObject['lightKey']
                lightKeyVar = ''
                for var in self.enteringVariables:
                    if lightKey in var:
                        lightKeyVar = var
                if not lightKeyVar:
                    raise Exception("The light key variable '{}' used in the surrogate Distillation model does not match"
                                    " any abrreviations\n of the inputs to the process interval {}. Plz check the "
                                    "distilation model to see if the name of the lightKey is correct "
                                    "".format(lightKey, self.intervalName))

                compositionEq = "model.var['{}'] == model.var['{}'] / (model.var['{}'] + 1e-12)" \
                                "".format(variables[0],lightKeyVar,self.incomingFlowVariable)
                allEnergyEquation.append(compositionEq)
                # add the variable to the list
                allEnergyVars += variables
                for var in variables:
                    self.boundaries.update({var: (0, None)})  # make it so that these variables can not hit zero
        else:
            energyConsumptionEq = "model.var['{}'] == {} * model.var['{}']".format(energyVar,
                                                                                   energyConsumptionParameter,
                                                                                   self.incomingFlowVariable)
            allEnergyEquation.append(energyConsumptionEq)
            allEnergyVars.append(energyVar)
            for var in allEnergyVars:
                self.boundaries.update({var: (0, None)})  # make it so that these variables can hit zero

        # -------- cost eqution for energy consumption
        costVariable = "cost_utility_energy_{}".format(self.intervalName)
        self.boundaries.update({costVariable: (0, None)})
        costEq = "model.var['{}'] == {} * model.var['{}'] ".format(costVariable, energyPrice, energyVar)

        # add all the varibles and equations to the pyomo list, and it's individual list
        self.pyomoEquations += allEnergyEquation + [costEq]
        self.allVariables['continuous'] += allEnergyVars + [costVariable]
        self.utilityEnergyEquations = allEnergyEquation + [costEq]
        self.utilityEnergyCostVariable = costVariable

    # helping functions not related to making equations
    def get_replacement_dict(self, newVars):
        """ This function makes the dictionary to facilitate switching out new variable in an interval with the new ones

        Params:
            newVars (list): list of the new incoming variables

        Returns:
            replacementDict (Dict): dictionary of the old names and their new variable names
            boundsDict (Dict): dictionary to update the boundries of the new variables!!
        """

        # get the name of the interval
        itervalName = self.intervalName

        # get the original input and operational variable names
        initialVars = self.inputs
        operationalVarsDict = self.operationalVariablesDict
        utilityVars = self.utilities

        boundsDict = {}  # preallocate
        replacementDict = {}  # preallocate
        for i in initialVars:
            for j in newVars:
                if i in j:  # the initial variable (of Excel) is always in the new name, that's how you can find wat belongs to where
                    replacementDict.update({i: j})
                    boundsDict.update({j: 'positiveReals'})

        # ideally you should call the new varible that is missing from the object self. This way you won't come into
        # conflict if you change the naming convention...
        for var in operationalVarsDict:
            newOperationalVar = '{}_{}'.format(var, itervalName)
            replacementDict.update({var: newOperationalVar})
            boundsDict.update({newOperationalVar: operationalVarsDict[var]})

        for var in utilityVars:
            newOperationalVar = '{}_{}'.format(var, itervalName)
            replacementDict.update({var: newOperationalVar})

        return replacementDict, boundsDict


class OutputIntervalClass:
    def __init__(self, outputName, outputBound, outputPrice, outputVariable, mixDict=None):
        if mixDict is None:
            mixDict = {}
        self.label = 'output'
        outputName = outputName.upper()
        self.outputName = outputName
        self.outputPrice = outputPrice
        if outputBound[1] == 'None' or outputBound[1] == 'none':
            self.outputBound = 'positiveReal'
        else:
            self.outputBound = outputBound
        self.outputVariable = outputVariable  # eg: acetate is output name but is refered to as ace in all the reactions
        self.allVariables = {'continuous': [outputName],
                             'boolean': [],  # there are no boolean variables for output intervals
                             'fraction': []}  # there are no fraction variables for output intervals
        self.boundaries = {outputName: self.outputBound}

    def make_output_equations(self, objects2mix):

        # find the leaving variables of the object(s) that go into the output
        leavingVars = []
        for objName in objects2mix:
            obj = objects2mix[objName][0]  # first element in tuple is the object
            connectInfo = objects2mix[objName][1]  # first element in tuple is connection info

            # if there is a separation process
            if isinstance(connectInfo, str) and 'sep' in connectInfo:
                allSepVars = obj.leavingInterval
                for var in allSepVars:
                    if connectInfo in var and self.outputVariable in var:
                        leavingVars.append(var)
            # otherwise just add the leaving variables
            else:
                leavingVars += obj.leavingInterval

        # make the equations
        endVar = self.outputName
        eqEnd = endVar + " == "
        pyomoEqEnd = "model.var['{}']".format(endVar) + " == "
        # startMixEq = eqEnd
        for i, lvar in enumerate(leavingVars):
            eqEnd += " + " + lvar
            pyomoEqEnd += " + " + "model.var['{}']".format(lvar)

        self.endEquations = eqEnd
        self.allEquations = [eqEnd]
        self.pyomoEquations = [pyomoEqEnd]
        # self.endVariables = endVar


class WastIntervalClass:
    def __init__(self, wasteVariables, wastePrice, mixDict=None):
        """ initiates the class object for the waste interval

        Parameters:
            wasteVariables (DF): variable that gives the components that go to the waste stream for each interval
            wastePrice (DF): data frame that gives the price of the waste stream per kg of waste generated
            mixDict (Dict):

        Retruns:
            initital waste oject (gets updated in 'update_intervals')

            """
        if mixDict is None:
            mixDict = {}
        self.label = 'waste'
        outputName = 'WASTE'
        self.outputName = outputName
        self.wastePrice = wastePrice
        self.wasteVariables = wasteVariables

        # iniciate equation list, varibale list and boundry dictionary
        self.allVariables = {'continuous': [], 'boolean': [], 'fraction': []}
        self.pyomoEquations = []
        self.boundaries = {}

    def make_waste_equations(self, objects2mix):
        """ makes the equations relating to the waste of all the intervals
         Parameters:
             objects2mix(dict)

         Return:
             pyomoEquations (list)
             """

        # find the leaving variables of the object(s) that go into the waste interval
        wasteVariableList = []
        equationList = []
        variableList = []

        for objName in objects2mix:
            wasteComponents = []
            obj = objects2mix[objName][0]  # first element in tuple is the object
            connectInfo = objects2mix[objName][1]  # first element in tuple is connection info

            # make waste variable
            wasteVariableMass = "waste_{}".format(objName)
            wasteVariableList.append(wasteVariableMass)

            # iniciate waste equation
            wasteEqPyomo = "model.var['{}'] == ".format(wasteVariableMass)

            # all the objects have to come from a seperation process otherwise it isn't waste: if so give error message
            if 'sep' not in connectInfo:
                raise Exception("The waste generated in interval {} does not come into the waste interval through a "
                                "seperation process: check the connenction matrix in the column 'waste'".format(
                    objName))

            # find the variables (of which seperated stream) that go into the waste equation
            allSepVars = obj.leavingInterval
            for var in allSepVars:
                if connectInfo in var:  # and outputVars in var:
                    wasteEqPyomo += "+ model.var['{}'] ".format(var)
                    # wasteComponents.append(var)
            equationList.append(wasteEqPyomo)

        # iniciate cost equation for waste
        wasteVariableCost = "cost_waste"  # make waste cost variable
        wastePrices = self.wastePrice  # get the prices of waste for each interval
        wasteCostEquation = "model.var['{}'] == ".format(wasteVariableCost)

        # change the index of the DF wastePrices to that of the waste variables
        # (should always be in the same order so should be ok)
        wastePrices = wastePrices.reset_index()
        wastePrices['waste_variables'] = wasteVariableList
        wastePrices = wastePrices.set_index('waste_variables')
        # wastePrices = wastePrices.reindex(wasteVariableList)

        # iterate over the waste variuables to calculate the cost of disposed waste
        for wasteVar, row in wastePrices.iterrows():
            wasteCostEquation += "+ model.var['{}'] * {}".format(wasteVar, row.waste_price)

        # make/ add a list of all equations and variable to pass on
        equationList.append(wasteCostEquation)
        variableList += wasteVariableList + [wasteVariableCost]

        # add eqautions and variable top the object, so they can be read later on
        self.allVariables['continuous'] += variableList
        self.costVariable = wasteVariableCost
        self.pyomoEquations += equationList
        for var in variableList:
            self.boundaries.update({var: (0, None)})  # positive reals


class MakeObjectiveClass:
    def __init__(self, inputObjects, processObjects, outputObjects, wasteObject):
        """ makes all  equations and variables related to economic parameters of each interval
        Params:
            inputObjects (Dict): a dictionary containing all the input interval objects
            outputObjects (Dict): a dictionary containing all the output interval objects
        """
        self.label = 'cost_model'
        self.allVariables = {'continuous': [],
                             'boolean': [],
                             'fraction': []}

        GREV_var, GREV_eq = self.make_GREV_equation(objectsOutputDict=outputObjects)
        OPEX_var, variables, equations = self.make_OPEX_equations(inputObjects=inputObjects,
                                                                  processObjects=processObjects,
                                                                  wasteObject=wasteObject)

        variableList = [GREV_var] + variables
        self.pyomoEquations = [GREV_eq] + equations
        self.allVariables['continuous'] += variableList

        # initiate the dictionary for boundries of the variables
        self.boundaries = {}

        # all the opex variable can be any real number
        for var in variables:
            self.boundaries.update({var: 'Reals'})  # all varible are Real numbers (so can be negative)

        # the GREV_var (Gross revenue) MUST be positive and bigger than 0!!!
        # : that way you ensure the production of products
        self.boundaries.update({GREV_var: [0.001, None]})

        # now make the objective that you want to maximise
        # The EBIT represents the objective function which is the (yearly? or hourly?) profit from
        # a given operation
        EBIT = "model.var['{}'] - model.var['{}']".format(GREV_var, OPEX_var)
        self.EBIT = EBIT

    def make_GREV_equation(self, objectsOutputDict):
        """ make the equation for the Gross revenue """
        GREV_var = "GREV"
        GREV_eq = "model.var['{}'] == ".format(GREV_var)
        for nameObj in objectsOutputDict:
            outObj = objectsOutputDict[nameObj]
            outputPrice = outObj.outputPrice
            outVar = outObj.outputName
            if not outputPrice:
                raise Exception('Hey you forgot to give a price for the output interval {}'.format(outVar))
            else:
                GREV_eq += "model.var['{}'] * {} + ".format(outVar, outputPrice)

        posPlus = GREV_eq.rfind('+')
        GREV_eq = GREV_eq[0:posPlus]
        # GREV_eq = '(' + GREV_eq + ')'

        return GREV_var, GREV_eq

    def make_OPEX_equations(self, inputObjects, processObjects, wasteObject):
        """ calculates the OPerating EXpenses of the superstructure: calulations for
        1) cost of raw material
        2) cost of utilities (chemicals)
        3) cost of utilities (energy)
        4) cost of waste
        """
        # preallocate the cost variables and equations that need to be added
        allCostVariables = []
        allCostEquations = []

        # -------- cost raw materials

        CostRawMaterialVar = "Raw_material_cost"
        CostRawMaterialEq = "model.var['{}'] == ".format(CostRawMaterialVar)
        for nameObj in inputObjects:
            # retrive the object from the dictionary
            inObj = inputObjects[nameObj]

            # retrive info from the object
            inputPrice = inObj.inputPrice
            inputVar = inObj.inputName
            inputCluster = inObj.clusterDict

            if inputCluster:
              for substrate in  inputCluster:
                  inputPrice = inputCluster[substrate]['price']
                  CostRawMaterialEq += "model.var['{}'] * {} + ".format(substrate, inputPrice)

            elif isinstance(inputPrice, float) or isinstance(inputPrice, int):
                CostRawMaterialEq += "model.var['{}'] * {} + ".format(inputVar, inputPrice)
            else:
                raise Exception('Hey you forgot to give a price or list of substrates (clusterDict) for the input '
                                'interval {}'.format(inputVar))

        posPlus = CostRawMaterialEq.rfind('+')  # finds the last '+' in the equation
        CostRawMaterialEq = CostRawMaterialEq[0:posPlus]  # delete the last plus

        # add to the lists of variables and equations
        allCostEquations.append(CostRawMaterialEq)
        allCostVariables.append(CostRawMaterialVar)
        self.CostRawMaterialEq = CostRawMaterialEq

        # -------- cost utilities (chemicals & energy)
        costUtilitiesVar = 'Utility_cost'
        costUtilitiesEqLeft = "model.var['{}'] == ".format(costUtilitiesVar)
        costUtilitiesEqRight = ''
        for nameObj, obj in processObjects.items():
            if hasattr(obj, 'utilityCostVariable'):
                utCostVar = obj.utilityCostVariable  # cost of chemicals
                costUtilitiesEqRight += " model.var['{}'] +".format(utCostVar)
            if hasattr(obj, 'utilityEnergyCostVariable'):
                utCostVar = obj.utilityEnergyCostVariable  # cost of energy
                costUtilitiesEqRight += " model.var['{}'] +".format(utCostVar)

        posPlus = costUtilitiesEqRight.rfind('+')  # finds the last '+' in the equation
        costUtilitiesEqRight = costUtilitiesEqRight[0:posPlus]  # delete the last plus

        # if there are cost associated to the use of utilities (that is the right hand side of the equation is not empty)
        # add the variable and equation to the list of variable and equations to declare
        if costUtilitiesEqRight:
            costUtilitiesEq = costUtilitiesEqLeft + costUtilitiesEqRight
            allCostVariables.append(costUtilitiesVar)
            allCostEquations.append(costUtilitiesEq)
            self.CostUtilitiesEq = costUtilitiesEq

        # -------- cost of waste
        wasteObj = wasteObject['waste']  # the waste object is in a dictionary with length 1 and key: 'waste'
        if hasattr(wasteObj, 'costVariable'):
            costWasteVar = wasteObj.costVariable
            allCostVariables.append(costWasteVar)
            self.costWasteEq = wasteObj.pyomoEquations[-1]  # the last element in the list is the cost equation

        # -------- make the OPEX equation
        OPEX_var = "OPEX"
        OPEX_eq = "model.var['{}'] == ".format(OPEX_var)
        for var in allCostVariables:
            OPEX_eq += "model.var['{}'] + ".format(var)

        posPlus = OPEX_eq.rfind('+')  # finds the last '+' in the equation
        OPEX_eq = OPEX_eq[0:posPlus]  # delete the last plus

        self.OPEX_eq = OPEX_eq

        # add to lists
        allCostVariables.append(OPEX_var)
        allCostEquations.append(OPEX_eq)

        return OPEX_var, allCostVariables, allCostEquations


# ============================================================================================================
# Validate the Excel file sheets, checks and error messages
# ============================================================================================================
# Write all error you encounter here, so you can write error messages later
# List of possible errors:
# 1) Make sure only split, sep bool ands mix are the only words in the connection matrix

def check_seperation_coef(coef, intervalName, amountOfSep, connectionMatrix):
    """
    makes sure the mass balances of the seperation processes are respected
    not sure if this really makes sens for distilation where you have different compositions in top and bottom however...
    revise is function

    Parameters:
        coef (str): string of all the seperation coef as read from the Excel file in the column speration_coef
        intervalName (str): name of the interval
        amountOfSep(int): the amount of separations by counting the number intervals in coef
        connectionMatrix (DF): a dataframe with the entire connecntion matrix

    Returns:
        erros made in the Excel file

    """
    # Set the index of the connenction matrix as the process inbnterval names
    # connectionMatrix = connectionMatrix.set_index('process_intervals') #already done when reading the Excel File

    # check if you have not forgotten to define the bounds if in the connection matrix, seperated streams have been defined
    if coef == 0:
        rowConnenctionMatrix = connectionMatrix.loc[intervalName]
        countSeperationInMatrix = 0
        for i in rowConnenctionMatrix:
            if isinstance(i, str) and 'sep' in i:
                raise Exception("No bounds where giving for the separation of the the interval '{}' but separated "
                                "streams are found in the connnection matrix. Plz "
                                "define seperation bounds in the Excel sheet 'process_intervals'".format(intervalName))

    coefList = split_remove_spaces(coef, ';')
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
    nCol = np.shape(sumOfArrays)  # get the amount of columns
    if sum(sumOfArrays != np.ones(nCol)):
        raise Exception("the sum of the seperation coefficients for {} do not add up to 1, check the excel file".format(
            intervalName))

    # check if all the seperation processes are accounted for in the connenction matrix
    try:
        rowMatrix = connectionMatrix.loc[intervalName]
    except:
        raise Exception("The interval name '{}' could not be found make sure the interval name does not contain spaces "
                        "and is written correctly".format(intervalName))

    counterSeparation = []
    counterSplit = []
    for connectCell in rowMatrix:
        keys = define_connect_info(connectInfo=connectCell)
        sepKey = keys[1]
        splitKey = keys[2]

        if 'sep' in sepKey:
            counterSeparation.append(sepKey)
        # if 'split' in splitKey:
        #     counterSplit.append(splitKey) # uncomment later to make checks on split streams

    # make a set so you can determin the unique seperation streams
    list_set = set(counterSeparation)
    nUniqueSeparation = len(list_set)

    #  nUniqueSeparation counted from the connection matrix
    # amountOfSep counted from the number of bounds found in the column separation_coef (sheet process_intervals)
    if nUniqueSeparation != amountOfSep:
        if nUniqueSeparation > amountOfSep:
            raise Exception(
                "There are more separation processes defined in the connection matrix then the amount of separation arrays defined in the column separation_coef"
                ""
                "Check interval {} there is a separation array missing in the column sepqrqtion_coef".format(
                    intervalName))
        else:
            raise Exception(
                "There are more separation processes defined in the amount of separation arrays (defined in the column separation_coef)"
                "then the connection matrix"
                "Check interval {} there is a separation process missing in the connection matrix".format(intervalName))


def check_excel_file(excelName):
    """ checks if the Excel file does not contain fatal errors for the generation of the super structure
    """

    loc = get_location(excelName)

    DFIntervals = pd.read_excel(loc, sheet_name='input_output_intervals')
    DFprocessIntervals = pd.read_excel(loc, sheet_name='process_intervals')
    DFeconomicParameters = pd.read_excel(loc, sheet_name='economic_parameters')
    DFConnectionMatrix = pd.read_excel(loc, sheet_name='connection_matrix')
    DFAbbr = pd.read_excel(loc, sheet_name='abbreviations')

    # check interval names in the connection matrix and interval list
    intervalNamesIn = remove_spaces(DFIntervals.process_intervals[DFIntervals.input_price != 0].to_list())
    intervalNamesReactors = remove_spaces(DFprocessIntervals['process_intervals'].to_list())
    intervalNamesOut = remove_spaces(DFIntervals.process_intervals[DFIntervals.output_price != 0].to_list())
    intervalNames = intervalNamesIn + intervalNamesReactors + intervalNamesOut

    intervalNamesConnenctionMatrixRow = remove_spaces(list(DFConnectionMatrix.columns))
    intervalNamesConnenctionMatrixRow.remove('process_intervals')
    intervalNamesConnenctionMatrixRow.remove('waste')
    intervalNamesConnenctionMatrixCol = remove_spaces(DFConnectionMatrix['process_intervals'].to_list())
    intervalNamesConnenctionMatrixCol.remove('waste')

    # check length
    if len(intervalNamesConnenctionMatrixCol) == len(intervalNamesConnenctionMatrixRow) == len(intervalNames):
        pass
    else:
        raise Exception('Interval name is missing in the connection matrix sheet or the interval sheet')
    # check names
    if intervalNames == intervalNamesConnenctionMatrixRow == intervalNamesConnenctionMatrixCol:
        pass
    else:
        positonError = [errorList for i, errorList in enumerate(intervalNames) if not intervalNames[i]
                                                                                      ==
                                                                                      intervalNamesConnenctionMatrixRow[
                                                                                          i] ==
                                                                                      intervalNamesConnenctionMatrixCol[
                                                                                          i]]
        print(positonError)
        raise Exception('The names in the connection matrix sheet or the interval sheet are not the same')

    # check if all abbreviations are defined
    abbreviations = split_remove_spaces(list(DFprocessIntervals.inputs), ',') \
                    + split_remove_spaces(list(DFprocessIntervals.outputs), ',') \
                    + split_remove_spaces(list(DFIntervals.components), ',')

    # uniqueListAbr = list(OrderedDict.fromkeys(abbreviations))
    abbrSet = set(abbreviations)
    abbrExcel = set(split_remove_spaces(list(DFAbbr.abbreviation), ','))

    missingAbbr = abbrSet - abbrExcel
    if missingAbbr:
        for abv in missingAbbr:
            if ".json" in abv:
                pass
            else:
                raise Exception('You are missing a definition for the following abbreviations: {}'.format(missingAbbr))
    else:
        pass

    # check if the process interval names are the same in the sheet economic_parameters and process_intervals
    if not DFeconomicParameters['process_intervals'].equals(DFprocessIntervals['process_intervals']):
        raise ValueError("Interval names in 'economic_parameters' and 'process_intervals' do not match")


# ============================================================================================================
# Functions to make the interval objects
# ============================================================================================================
def read_excel_sheets4_superstructure(excelName):
    '''
    Reads the Excel file containing the data for superstructure generation and returns each sheet as a dataframes
    The dataframes are stored in a dictionary

    parameters:
    excelName (str): name of the Excel file saved in the file directory 'excel files'

    returns:
    ExcelDict (Dict): a dictionary containing DF
    '''
    # read Excel file
    loc = get_location(file=excelName)

    # read excel info
    DFInOutIntervals = pd.read_excel(loc, sheet_name='input_output_intervals')
    DFconnectionMatrix = pd.read_excel(loc, sheet_name='connection_matrix', index_col='process_intervals')
    DFprocessIntervals = pd.read_excel(loc, sheet_name='process_intervals', index_col='process_intervals')
    DFeconomicParameters = pd.read_excel(loc, sheet_name='economic_parameters', index_col='process_intervals')
    DFmodels = pd.read_excel(loc, sheet_name='models', index_col='model_name')
    DFabrriviations = pd.read_excel(loc, sheet_name='abbreviations')

    ExcelDict = {
        'input_output_DF': DFInOutIntervals,
        'connection_DF': DFconnectionMatrix,
        'process_interval_DF': DFprocessIntervals,
        'economic_parameters_DF': DFeconomicParameters,
        'models_DF': DFmodels,
        'abbreviations_DF': DFabrriviations
    }

    return ExcelDict


def make_mix_dictionary(intervalName, DFconnectionMatrix):
    """ Returns the mixing dictionary based on the connection matrix

    Parameters:
        intervalName (str): name of the interval where the mixed streams arrive
        DFconnectionMatrix (DF): dataframe with the connection matrix

    Return:
        MixDict (dict): dictionary of interval objects that needs to be mixed
    """
    # check if it is mixed with other reactors
    # processInvervalNames = list(DFconnectionMatrix.index)

    reactorCol = DFconnectionMatrix[intervalName]
    specifications = reactorCol.where(
        reactorCol != 0).dropna()  # find where mixing takes place, mixed streams are in the same colunm
    intervals2Mix = list(specifications.index)  # the indexs are the names of the process interval to mix

    # get rid of the diagonal element refering to itself (if it is there due to defining the boolean variable)
    try:
        # drop the boolean variable which has as index the current interval name
        intervals2Mix.remove(intervalName)
        specifications = specifications.drop(intervalName)
    except:
        pass

    mixDict = {}  # preallcoation
    if len(specifications) >= 2:
        for k, specs in enumerate(specifications):
            mixDict.update({intervals2Mix[k]: specs})
    return mixDict


def make_boolean_equations(DFconnectionMatrix, processIntervalnames):
    """ makes the equations that regulate if a certain reactor is chosen or not nl: 1 == sum(boolean variables)
    the function works as followed: the DFconnectionMatrix excludes the inputs!! important!! the input boolean variables
    are regulated in the first input objected.
    main idea:
    1) loop over the rows of the DF.
    2) count the colums of this row which are not zero in  sequence! (save the interval names in a list)
    3) the longest list is the list of intervals dependant of the boolean variable
    4) the boolean variable can be found on the diagonal (i.e. with the same row and column index)
    5) make the boolean equation
    6) deleet the columns that contain the sequence
    7) restart at step 1 till the DF is empty

    returns:
        boolean variables (list): list of boolean variables
        boolean equations (lsit): list of boolean equations
    """

    # prun the Dataframe, anything that does not have a bool label on the diagonal can be droped
    toDrop = []
    for i in processIntervalnames:
        if not isinstance(DFconnectionMatrix[i][i], str):
            toDrop.append(i)

    DFconnectionMatrix = DFconnectionMatrix.drop(labels=toDrop, axis=1)

    switch = True
    equationsSumOfBools = []
    booleanVariables = []
    iteration = 0
    while switch:
        iteration += 1
        saveDict = {}
        for index, row in DFconnectionMatrix.iterrows():
            intervalNames = []
            for indexCol, element in row.items():
                # element = row[colName]
                if isinstance(element, str) or element != 0:  # so comes from a separation or just connected by '1'
                    intervalNames.append(indexCol)
                else:
                    break

            saveDict.update({index: intervalNames})

        # find the key in the saveDict that has the longest list
        key_max_sequential = max(saveDict, key=lambda k: len(saveDict[k]))

        # create the equation
        eq = '1 == '
        for interval in saveDict[key_max_sequential]:
            boolVar = DFconnectionMatrix[interval][interval]
            booleanVariables.append(boolVar)
            eq += "+ model.boolVar['{}'] ".format(boolVar)
        equationsSumOfBools.append(eq)

        # now we need to drop the colums in the original dataframe that already form 1 set of equations
        # the colunms that need to be droped are in saveDict[key_max_sequential]1
        cols2drop = saveDict[key_max_sequential]
        DFconnectionMatrix = DFconnectionMatrix.drop(labels=cols2drop, axis=1)

        if DFconnectionMatrix.empty:  # once the DF is empty we can stop the while loop
            switch = False

        if iteration > 50:
            raise Exception("50 iteration have passed to try and make the boolean equations, check to see if the "
                            "connection matrix is correctly formulated or check the function 'make_boolean_equations'")

    # check that all the boolean variables have unique values
    uniqueSet = set(booleanVariables)
    if len(uniqueSet) != len(booleanVariables):
        raise Exception("the boolean variables are not all unique, check the diagonal of the connection matrix")

    return booleanVariables, equationsSumOfBools


# functions to automate making the interval class objects
def make_input_intervals(ExcelDict, clusterDict):
    """ Makes the process intervals of inputs.

    parameters:
    ExcelDict (Dict): a dictionary with all the data to make the superstructure

    return:
    objectDict (Dict): a dictionary holding all the class objects for input intervals
    """

    DFIntervals = ExcelDict['input_output_DF']
    DFconnectionMatrix = ExcelDict['connection_DF']

    # inputs
    inputPrices = DFIntervals.input_price.to_numpy()
    posInputs = inputPrices != 0  # find where the input interval are (they have an input price)

    # find names of input interval variable
    InputIntervalNames = DFIntervals.process_intervals[posInputs]
    componentsList = DFIntervals.components[posInputs]
    compositionsList = DFIntervals.composition[posInputs]

    # find upper and lower bounds on the mass
    inBoundsLow = DFIntervals.lower_bound[posInputs].to_numpy()
    inBoundsUpper = DFIntervals.upper_bound[posInputs].to_numpy()

    ### the input price can be either specifed or can be read from the json file incase it's a cluster
    inputPrices = inputPrices[posInputs]

    # define fixed parameters cost raw material
    inputPriceDict = {InputIntervalNames[i]: inputPrices[i] for i in range(len(inputPrices))}  # make dictionary
    boundryDict = {InputIntervalNames[i]: [inBoundsLow[i], inBoundsUpper[i]] for i in
                   range(len(inputPrices))}  # make dictionary

    # loop over all the inputs and make a class of each one
    objectDictionary = {}
    for i, intervalName in enumerate(InputIntervalNames):
        # pass the bool variable responsible for activating an input if present
        diagonalInputInfo = DFconnectionMatrix[intervalName][intervalName]
        booleanVar = None
        if isinstance(diagonalInputInfo, str):
            booleanVar = diagonalInputInfo

        inputPrice = inputPriceDict[intervalName]
        boundryInput = boundryDict[intervalName]
        componentsOfInterval = componentsList[i].split(",")
        compositionsofInterval = compositionsList[i]  # string or 1, depends on if there are different components

        # preallocate lists and dictionaries
        compsitionDictionary = {}  # preallocate dictionary

        # make the composition dictionary
        # Todo for future me, may be in a cluster you would like to select mixtures... adapt the code for this necessity
        #  if required in the future
        # compositionsofInterval would be the variable specifying how many substrates you would want to select

        if not isinstance(compositionsofInterval, str) and compositionsofInterval == 1:  # if it is one, no need to loop over the dictionary, there is only one compound
            component = componentsOfInterval[0].replace(' ', '')
            fraction = compositionsofInterval  # should always be one component in the stream
            compsitionDictionary.update({component: fraction})
        else:
            compositionsofInterval = compositionsList[i].split(",")
            for j, component in enumerate(componentsOfInterval):
                component = component.replace(' ', '')  # get rid of spaces
                fraction = compositionsofInterval[j]
                fraction = fraction.replace(' ', '')
                fraction = float(fraction)
                compsitionDictionary.update({component: fraction})

        # create object
        objectInput = InputIntervalClass(inputName=intervalName, compositionDict=compsitionDictionary,
                                         inputPrice=inputPrice, boundryInputVar=boundryInput,
                                         booleanVariable=booleanVar, clusterDict=clusterDict )
        objectDictionary.update({intervalName: objectInput})

    return objectDictionary


def make_process_intervals(ExcelDict):
    """ Makes the process intervals of the process intervals (excluding inputs and outputs).

        parameters:
        ExcelDict (Dict): a dictionary with all the data to make the superstructure

        return:
        objectDict (Dict): a dictionary holding all the class objects for process intervals
        """

    # DFInOutIntervals = ExcelDict['input_output_DF']
    DFconnectionMatrix = ExcelDict['connection_DF']
    DFprocessIntervals = ExcelDict['process_interval_DF']
    DFeconomicParameters = ExcelDict['economic_parameters_DF']
    DFmodels = ExcelDict['models_DF']

    # loop over all the process intervals to make their respective objects
    listProcessIntervals = list(DFprocessIntervals.index)
    objectDictionary = {}  # preallcoate a dictionary with the interval names and the interval objects
    for intervalName in listProcessIntervals:
        intervalName = intervalName.replace(' ', '')  # remove annoying spaces
        # inputs of reactor
        inputsReactor = DFprocessIntervals.inputs[intervalName]
        if ".json" in inputsReactor:
            jsonLoc = get_location(file=inputsReactor)
            with open(jsonLoc) as file:
                inputsObject = json.load(file)
                inputsReactor = list(inputsObject.keys())
        else:
            inputsReactor = split_remove_spaces(inputsReactor, ',')
        # outputs of the reactor
        outputsReactor = DFprocessIntervals.outputs[intervalName]
        outputsReactor = split_remove_spaces(outputsReactor, ',')
        # find the bounds of the interval
        boundsReactor = eval(DFprocessIntervals.interval_bounds[intervalName])
        nameDict = {intervalName: boundsReactor}

        # find the correct equation of the reactor
        reactionIndicator = DFprocessIntervals.reaction_model[
            intervalName]  # indicates in wat format to read the equations

        if isinstance(reactionIndicator,
                      int) and reactionIndicator == 0:  # if it is zero there is no reaction, only seperation, e.g., for a distilation process
            equations = None

        elif 'xml' in reactionIndicator:
            modelName = DFprocessIntervals.reaction_model[intervalName]
            try:
                equationInfo = DFmodels.loc[modelName]  # check if said model is present in the Excel sheet models
            except:
                raise Exception(
                    'make sure the name {} is identical in the sheet sbml_models and process_intervals'.format(
                        modelName))

            inputID = equationInfo.SBML_input_ID
            outputID = equationInfo.SBML_output_ID
            outputID = equationInfo(outputID, ',')

            equations = make_str_eq_smbl(modelName=modelName, substrate_exchange_rnx=inputID,
                                         product_exchange_rnx=outputID, equationInfo=equationInfo)
            # print('s')
        elif 'json' in reactionIndicator:
            jsonFile = DFprocessIntervals.reaction_model[intervalName]
            try:
                equationInfo = DFmodels.loc[jsonFile]  # check if said model is present in the execel sheet models
            except:
                raise Exception(
                    "make sure the name for the json file {} is identical in the sheet 'models' and 'reactor_interval'"
                    "s__".format(jsonFile))
            # find the save location
            jsonLoc = get_location(file=jsonFile)
            with open(jsonLoc) as file:
                reactionObject = json.load(file)

            equations = make_str_eq_json(reactionObject, equationInfo)

        else:  # the equation is written as a string in the Excel file handy for testing porposes
            equations = DFprocessIntervals.reaction_model[intervalName]
            equations = split_remove_spaces(equations, ',')
            for eq in equations:
                if '==' not in eq:
                    raise Exception(
                        'Take a look at the reaction model mate, there is no json, xml or correct reaction given'
                        ' for reactor {}'.format(intervalName))

        # find if the interval is dependent on a boolean variable
        boolVar = None
        boolInformation = DFconnectionMatrix[intervalName][intervalName]
        if isinstance(boolInformation, str):
            boolVar = boolInformation

        # find special component bounds like that for pH
        boundsComponentStr = DFprocessIntervals.input_bounds[intervalName]
        if not isinstance(boundsComponentStr, str):
            boundsComponent = {}
        else:
            boundsComponent = str_2_dict(boundsComponentStr, intervalName)

        # find if there are any utilities that are used
        utilityDict = {}  # only one utilty per process interv al is permitted for the moment
        if DFprocessIntervals.ut_chemical[intervalName] != 0:
            # get utility name
            utilityVariableNames = DFprocessIntervals.ut_chemical[intervalName]
            utilityVariableNames = split_remove_spaces(utilityVariableNames, ',')
            # get utility parameter
            utilityParameter = [DFprocessIntervals.mu_ut[
                                    intervalName]]  # made as a list so in the future multiple utility componets can be added
            # get price utility
            utilityPrice = [DFeconomicParameters.ut_chem_price[
                                intervalName]]  # made as a list so in the future multiple utility componets can be added
            # utilityPrice = split_remove_spaces(utilityPrice,',')
            for j, utilityName in enumerate(utilityVariableNames):
                utilityDict.update({utilityName: {'cost': utilityPrice[j], 'parameter': utilityParameter[j]}})

        # find if the separated streams and where they go to/ the components to separate
        seperationDict = {}
        outputsStr = DFprocessIntervals.outputs[intervalName]
        coefStr = DFprocessIntervals.seperation_coef[intervalName]
        coefList = split_remove_spaces(coefStr, ';')
        amountOfSeperations = len(coefList)
        check_seperation_coef(coefStr, intervalName, amountOfSeperations, DFconnectionMatrix)

        if DFprocessIntervals.seperation_coef[intervalName] != 0:  # and DFprocessIntervals.has_seperation[i] < 2 :
            for j in range(amountOfSeperations):
                seperationName = intervalName + '_sep{}'.format(j + 1)
                coefTuple = stringbounds_2_tuplebounds(coefList[j])
                outputs = split_remove_spaces(outputsStr, ',')
                specificSeperationDict = {}
                for k, outputName in enumerate(outputs):
                    specificSeperationDict.update({outputName: coefTuple[k]})
                seperationDict.update({seperationName: specificSeperationDict})
            # objectReactor.separation = seperationDict

        # check if it is mixed with other reactors
        mixDict = make_mix_dictionary(intervalName=intervalName, DFconnectionMatrix=DFconnectionMatrix)

        # find to which interval the stream is split to (indicated in the connection matrix)
        intervalRow = DFconnectionMatrix.loc[intervalName]
        splitList = []  # find the reactor or separation stream to split
        for j, info in enumerate(intervalRow):
            if isinstance(info, str) and 'split' in info and 'sep' in info:
                indexSep = info.find('sep')
                separationStream = info[indexSep:indexSep + 4]
                splitList.append('{}_{}'.format(intervalName, separationStream))

            elif isinstance(info, str) and 'split' in info and not 'sep' in info:
                splitList.append(intervalName)

        # trick to get unique values
        setSplits = set(splitList)
        listSplits = list(setSplits)

        # pass on the connection information
        connectedIntervals = get_connected_intervals(intervalName=intervalName, conectionMatrix=DFconnectionMatrix)

        # pass on the operational variables if there are any
        operationalVars = DFprocessIntervals.operation_bounds.loc[intervalName]
        operationalVarDict = {}
        if isinstance(operationalVars, str):
            try:
                operationalVarDict = eval(operationalVars)
            except:  # raise an exception if the dictionary is not well writen
                raise Exception("The operational variables {} for interval {} are not in dictionary "
                                "format: \n {} the format is {'var1': [lb, ub], 'var2': [lb, ub] ... }")

        # pass on parameters for the energy consumption of the interval
        energyConsumptionParameter = DFprocessIntervals.energy_consumption.loc[intervalName]
        energyPrice = DFeconomicParameters.ut_energy_price.loc[intervalName]
        utilityEnergyDict = {'price': energyPrice, 'consumption_parameter': energyConsumptionParameter}

        # make initial interval object
        objectReactor = ProcessIntervalClass(inputs=inputsReactor, boundryInputVar=boundsComponent,
                                             outputs=outputsReactor, reactionEquations=equations, nameDict=nameDict,
                                             mix=mixDict, utilities=utilityDict, energyUtility=utilityEnergyDict,
                                             separationDict=seperationDict,
                                             splitList=listSplits, booleanVariable=boolVar,
                                             operationalVariablesDict=operationalVarDict)
        # put the object in the dictionary
        objectDictionary.update({intervalName: objectReactor})
    return objectDictionary


def make_output_intervals(ExcelDict):
    """ Makes the process intervals of outputs.

        parameters:
        ExcelDict (Dict): a dictionary with all the data to make the superstructure

        return:
        objectDict (Dict): a dictionary holding all the class objects for output intervals
        """

    # read excel info
    DFconnectionMatrix = ExcelDict['connection_DF']
    DFIntervals = ExcelDict['input_output_DF']

    # find the output interval names and information in one step
    output_data = DFIntervals[DFIntervals.output_price != 0].drop(['output_price'], axis=1)
    objectDictionary = {}  # preallcoate a dictionary with the interval names and the interval objects
    for i, row in output_data.iterrows():
        intervalName = row.process_intervals.replace(' ', '')

        # find the bounds of the interval
        lowerBound = row.lower_bound
        upperBound = row.upper_bound
        outputBound = [lowerBound, upperBound]

        outputPrice = DFIntervals.output_price[i]

        # find the variable name acossiated with the output
        outputVariable = row.components.replace(' ', '')

        # check if it is mixed with other reactors
        mixDict = make_mix_dictionary(intervalName=intervalName, DFconnectionMatrix=DFconnectionMatrix)

        # make initial interval object
        objectReactor = OutputIntervalClass(outputName=intervalName, outputBound=outputBound,
                                            outputPrice=outputPrice, outputVariable=outputVariable, mixDict=mixDict)
        # put the object in the dictionary
        objectDictionary[intervalName] = objectReactor

    return objectDictionary


def make_waste_interval(ExcelDict):
    """ Makes the process intervals of the waste interval.

            parameters:
            ExcelDict (Dict): a dictionary with all the data to make the superstructure

            return:
            objectDict (Dict): a dictionary holding the waste objects for the waste interval
            """

    # retrive data
    DFconnectionMatrix = ExcelDict['connection_DF']
    DFprocessIntervals = ExcelDict['process_interval_DF']
    DFeconomicParameters = ExcelDict['economic_parameters_DF']
    intervalName = 'waste'

    # find the prices of waste per interval
    priceWasteDF, wasteVariables = DFeconomicParameters.loc[DFconnectionMatrix.waste[
                                                                DFconnectionMatrix.waste != 0].dropna().index, "waste_price"], \
                                   DFprocessIntervals.loc[DFconnectionMatrix.waste[
                                                              DFconnectionMatrix.waste != 0].dropna().index, "outputs"]

    # check if it is mixed with other reactors
    mixDict = make_mix_dictionary(intervalName=intervalName, DFconnectionMatrix=DFconnectionMatrix)

    # make initial interval object
    objectWaste = WastIntervalClass(mixDict=mixDict, wastePrice=priceWasteDF, wasteVariables=wasteVariables)

    # put the object in the dictionary
    objectDictionary = {intervalName: objectWaste}

    return objectDictionary


# ============================================================================================================
# Functions to update the interval objects
# ============================================================================================================
def update_intervals(allIntervalObjectsDict, ExcelDict):
    """ Updated the equations of all interval objects. the mixing equations are added, the reaction equations are updated
    and the utlity equations are added.
    # check if this function is really necesary, can't I just get the equations right in one loop? i.e., in the function
    make_process_intervals

    Parameters:
        allIntervalObjectsDict (Dict): dictionary containing all objects and their interval names
        ExcelDict (Dict): Dictionary containing all info on the superstructure in the form of dataframes

        """

    connectionMatrix = ExcelDict['connection_DF']
    # DFprocessIntervals = ExcelDict['process_interval_DF']
    for intervalName in allIntervalObjectsDict:
        intervalObject = allIntervalObjectsDict[intervalName]
        label = intervalObject.label
        connectedIntervals = get_connected_intervals(intervalName=intervalName, conectionMatrix=connectionMatrix)

        if label == 'process_interval':
            # get the connection info if there is only one connecting interval (is there mixing or not)
            try:
                connectInfo = list(connectedIntervals.values())[0]
            except:
                raise Exception('The interval {} is not connected to any previous interval, '
                                'check the connection matrix'.format(intervalName))

            simpleConcention, sepKey, splitKey, boolKey = define_connect_info(connectInfo)
            # simpleConcention = True  # just connecting from one reactor to the next with the connection possibly being a bool

            # get the previous interval object (if there is mixing these variables are ignored)
            previousIntervalName = list(connectedIntervals.keys())[0]
            previousIntervalObject = allIntervalObjectsDict[previousIntervalName]
            enteringVariables = previousIntervalObject.leavingInterval

            if len(connectedIntervals) == 1:  # in other words no mixing
                # update_reactor_equations: current interval connected by 1 interval
                if simpleConcention:  # and connectInfo == 1 or if it is 'bool' (does not matter)
                    intervalObject.update_interval_equations(enteringVariables)
                    intervalObject.make_incoming_massbalance_equation(enteringVariables)

                # update_reactor_equations: current interval is connected by a separation stream and/or split stream
                elif sepKey and splitKey:
                    newReactorInputs4Interval = []
                    for var in enteringVariables:
                        if sepKey in var and splitKey in var:
                            newReactorInputs4Interval.append(var)
                    intervalObject.update_interval_equations(newReactorInputs4Interval)
                    intervalObject.make_incoming_massbalance_equation(newReactorInputs4Interval)

                elif sepKey and not splitKey:
                    newReactorInputs4Interval = []
                    for var in enteringVariables:
                        if sepKey in var:
                            newReactorInputs4Interval.append(var)
                    intervalObject.update_interval_equations(newReactorInputs4Interval)
                    intervalObject.make_incoming_massbalance_equation(newReactorInputs4Interval)

                elif splitKey and not sepKey:  # only splitting remains
                    newReactorInputs4Interval = []
                    for var in enteringVariables:
                        if splitKey in var:
                            newReactorInputs4Interval.append(var)
                    intervalObject.update_interval_equations(newReactorInputs4Interval)
                    intervalObject.make_incoming_massbalance_equation(newReactorInputs4Interval)

            # update_reactor_equations:
            # current interval is connected by multiple intervals by MIXING (including mixing separated streams)
            elif len(connectedIntervals) > 1:  # so here is mixing
                objectDict2mix = {
                    nameObjConect: (allIntervalObjectsDict[nameObjConect], connectedIntervals[nameObjConect]) for
                    nameObjConect in connectedIntervals}
                intervalObject.make_mix_equations(objectDict2mix)
                newReactorInputs4Interval = intervalObject.mixingVariables

                intervalObject.update_interval_equations(newReactorInputs4Interval)
                intervalObject.make_incoming_massbalance_equation(newReactorInputs4Interval)

            if intervalObject.utilities:
                # if the utilities dictionary is not empty, there is a utility to be added to the interval
                intervalObject.make_utility_equations()

            energyConsumptionParam = intervalObject.utilityEnergy['consumption_parameter']
            if isinstance(energyConsumptionParam, str) or energyConsumptionParam != 0:
                intervalObject.make_energy_utility_equations()

        elif label == 'output':
            objectDict2mix = {nameObjConect: (allIntervalObjectsDict[nameObjConect], connectedIntervals[nameObjConect])
                              for nameObjConect in connectedIntervals}
            intervalObject.make_output_equations(objectDict2mix)

        elif label == 'waste':
            objectDict2mix = {nameObjConect: (allIntervalObjectsDict[nameObjConect], connectedIntervals[nameObjConect])
                              for nameObjConect in connectedIntervals}
            intervalObject.make_waste_equations(objectDict2mix)


def get_vars_eqs_bounds(objectDict):
    """ Returns all the varible that pyomo needs to declare, the equations (in pyomo format) and a dictionary with the
    bounds of the variable (if there are any). Also prints the equations in an orderly fashion to help with debuging

    Params:
        objectDict (dict) : a dictionary holding all the objects of the superstructure

    Returns:
        variables (list) : all variable of the whole super structure
        equations (list) : all equaitions of the superstructure
        boundsContinousVars (dict) : all bounds of the variabels

    """
    variables = {}
    continuousVariables = []
    booleanVariables = []
    fractionVariables = []
    equations = []
    boundsContinousVars = {}
    for objName in objectDict:
        obj = objectDict[objName]
        equations += obj.pyomoEquations
        pluses = '++'
        # print(pluses*10 + objName + pluses*10)
        print('{}__{}__{}'.format(pluses * 10, objName, pluses * 10))
        if obj.label == 'process_interval':

            # print the equations per classification
            if hasattr(obj, 'mixEquations'):
                print('------ mixing equations ------')
                mixEq = obj.mixEquations
                for e in mixEq:
                    print(e)
                print('')

            if hasattr(obj, 'incomingFlowEquation'):
                print('------ Incoming mass equation ------')
                enteringEq = obj.incomingFlowEquation[0]
                print(enteringEq)
                print('')

            if hasattr(obj, 'utilityEquations'):
                print('------ chemical utility equations ------')
                utEq = obj.utilityEquations
                for e in utEq:
                    print(e)
                print('')

            if hasattr(obj, 'reactionEquations'):
                print('------ reaction equations ------')
                rxnEq = obj.reactionEquations
                for e in rxnEq:
                    print(e)
                print('')

            if hasattr(obj, 'separationEquations'):
                print('------ separation equations ------')
                sepEq = obj.separationEquations
                for e in sepEq:
                    print(e)
                print('')

            if hasattr(obj, 'utilityEnergyEquations'):
                print('------ Energy utility equations ------')
                utEnEq = obj.utilityEnergyEquations
                for e in utEnEq:
                    print(e)
                print('')

            if hasattr(obj, 'splitEquations'):
                print('------ separation equations ------')
                splitEq = obj.splitEquations
                for e in splitEq:
                    print(e)
                print('')

            print('')

        else:  # obj.label == 'output':
            for eq_interval in obj.pyomoEquations:
                print(eq_interval)
            print('')  # print a space for readability

        # collect all the variables
        continuousVariables += obj.allVariables['continuous']
        booleanVariables += obj.allVariables['boolean']
        fractionVariables += obj.allVariables['fraction']
        boundsContinousVars = boundsContinousVars | obj.boundaries

    # remove double variables in the list of continuous variables
    unique_list_var_continous = list(OrderedDict.fromkeys(
        continuousVariables))  # preserves order (easier to group the equations per interval this way)
    unique_list_var_bool = list(OrderedDict.fromkeys(booleanVariables))

    # dictionary to bundel  all the varibles
    variables = {'continuous': unique_list_var_continous,
                 'boolean': unique_list_var_bool,
                 'fraction': fractionVariables}

    return variables, equations, boundsContinousVars


# ============================================================================================================
# Master function: generates the superstructure
# ============================================================================================================

def make_super_structure(excelFile, printPyomoEq=False):
    """ Master function: calls all other functions to make the superstructure

    Declare all interval variables (capital letters) and component variables (small letters)
    loop over all objects
    declare all variables
    declare equations of each object
    make the objective

    params:
        excelFile (str): name of the Excel file saved in the file location 'excel files'

    returns:
        model (pyomo structure): the model of the super structure
    """

    model = pe.ConcreteModel()
    check_excel_file(excelName=excelFile)
    excelDict = read_excel_sheets4_superstructure(excelName=excelFile)

    boolObject = BooleanClass(ExcelDict=excelDict)
    clusterDict = boolObject.clusterDict # the cluster dictionary is already made in the boolean object

    boolObjectDict = {'boolean_object': boolObject}
    objectsInputDict = make_input_intervals(ExcelDict=excelDict,clusterDict=clusterDict)
    objectsProcessDict = make_process_intervals(ExcelDict=excelDict)
    objectsOutputDict = make_output_intervals(ExcelDict=excelDict)
    objectsWasteDict = make_waste_interval(ExcelDict=excelDict)

    allObjectsDict = objectsInputDict | objectsProcessDict | objectsOutputDict | objectsWasteDict
    # update the intervals, so they are conected to the right interval
    update_intervals(allObjectsDict, excelDict)

    # make an object with all the equations of the cost models
    CostModelObj = MakeObjectiveClass(inputObjects=objectsInputDict,
                                      outputObjects=objectsOutputDict,
                                      processObjects=objectsProcessDict,
                                      wasteObject=objectsWasteDict)

    CostModelDict = {'cost_model': CostModelObj}

    allObjects = boolObjectDict | allObjectsDict | CostModelDict  # add the logic model and the cost model to the object list
    variables, equations, bounds = get_vars_eqs_bounds(allObjects)

    def boundsRule(model, i):
        boudVar = bounds[i]
        if isinstance(boudVar, list) or isinstance(boudVar, tuple):
            lowerBound = boudVar[0]
            upperBound = boudVar[1]

        elif boudVar == 'positiveReals':  # including 0
            lowerBound = 0  # 0.0001
            upperBound = None

        elif boudVar == 'Reals':
            lowerBound = None
            upperBound = None

        else:  # elif isinstance(boudVar, str):  # 'bool' or 'positiveReal' in boudVar
            lowerBound = 0
            upperBound = None
        return (lowerBound, upperBound)

    model.var = pe.Var(variables['continuous'], domain=pe.Reals, bounds=boundsRule)
    if variables['boolean']:
        # noinspection PyUnresolvedReferences
        model.boolVar = pe.Var(variables['boolean'], domain=pe.Boolean)
    if variables['fraction']:
        # noinspection PyUnresolvedReferences
        model.fractionVar = pe.Var(variables['fraction'], domain=pe.PercentFraction)

    # introduce the equations to pyomo
    model.constraints = pe.ConstraintList()
    for eq in equations:
        # print(eq) # printing of equation is now in the function get_vars_eq_bounds
        # expresion = eval(eq)
        # model.constraints.add(expresion)
        try:
            expresion = eval(eq)
        except:
            raise Exception('The following equation can not be read by eval: {}'.format(eq))

        try:
            model.constraints.add(expresion)
        except:
            raise Exception('The following equation can not be read by pyomo: {}'.format(eq))

    # define the objective
    # EBIT
    objectiveExpr = CostModelObj.EBIT
    print("the objective of the model is to maximise the EBIT: ")
    print(objectiveExpr)
    model.objectiveValue = pe.Objective(expr=eval(objectiveExpr), sense=pe.maximize)  # sense=pe.maximize

    # model.pprint()
    if printPyomoEq:
        model.pprint()  # debug check

    return model
