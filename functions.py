import pandas as pd
import cobra
import os
import cobra.io
import re
import numpy as np
import warnings

"""Created on the 30.01.2023
@author: Lucas Van der Hauwaert
@author contact: lucas.vanderhauwaert@usc.es

All functions used to create superstructures 
All functions used to analise GEMs 
"""

""" TODO
add all the functions used into one file
separate the different function using comments
e.g.;
"""

########################################################################################################################
# ============================================================================================================
# Usefull functions
# ============================================================================================================
########################################################################################################################

def split_remove_spaces(expr2split,splitCharacter):
    """ Splits a given string according to a specified character e.g., ','
    also removes the spaces ' '
    """
    exprList = []
    if not isinstance(expr2split, str) and not isinstance(expr2split,list): # so a float or int just one number
        return [float(expr2split)]   # if mulptiple prices for utilty you're gona have to change str to floats

    elif isinstance(expr2split,list):
        for exp in expr2split:
            expresions = exp.split(splitCharacter)
            for i in expresions:
                i = i.replace(' ', '')  # remove annoying spaces
                exprList.append(i)
        return exprList

    elif isinstance(expr2split, str):
        expresions = expr2split.split(splitCharacter)
        for i in expresions:
            i = i.replace(' ', '') # remove annoying spaces
            exprList.append(i)
        return exprList

def stringbounds_2_tuplebounds(stringBound):
    """ transforms a bound that is writen as a string (e.g., '[20, 50]') to a tuple
    """
    stringBound = stringBound.replace('[', '')
    stringBound = stringBound.replace(']', '')
    bounds = stringBound.split(',')
    boundsList = []
    for i in bounds:
        boundsList.append(float(i))
    return boundsList

def load_objectes_from_dictionary(dict):
    """ deleet I think, not used"""
    for i in dict:
        locals()[i] = dict[i]

def remove_spaces(listOfInterest):
    """ removes spaces """
    exprList = []
    for i in listOfInterest:
        i = i.replace(' ','')
        exprList.append(i)
    return exprList

def get_connected_intervals(intervalName,conectionMatrix):
    """From the conection matrix (see Excel file) the connected intervals are found which go into the specified
    interval (intervalName) is found
    """
    conectionCol = conectionMatrix[intervalName]
    posConnect = conectionCol != 0
    nameConnectedIntervals = list(conectionMatrix['process_intervals'][posConnect])
    connectionInfo = list(conectionMatrix[intervalName][posConnect])
    connectionDict = {nameConnectedIntervals[i]:connectionInfo[i] for i in range(len(connectionInfo))}
    return  connectionDict

def str_2_dict(string,intervalname):
    D = eval(string)
    inputBoundsDict = {}
    for i in D:
        inputBoundsDict.update({i+'_'+intervalname : D[i]})
    return inputBoundsDict

def get_location(file):
    """ gets the file location from the Directory 'excel files'
    """
    loc = os.getcwd()
    posAlquimia = loc.find('Alquimia')
    loc = loc[0:posAlquimia + 8]

    if '/' in loc: # in the case of Mac OS
        file = r"/{}".format(file)
        if '.xslx' in file:
            loc = loc + r'/excel files' + file
        elif '.xml' in file:
            loc = loc + r'/SBML models' + file

    elif "\\" in loc : # in the case of Windows OS
        file = r"\{}".format(file)
        if '.xslx' in file:
            loc = loc + r'\excel files' + file
        elif '.xml' in file:
            loc = loc + r'\SBML models' + file
    return loc

########################################################################################################################
# ============================================================================================================
# SBML screening
# ============================================================================================================
########################################################################################################################

# originaly from the file f_find_carbons
def findCarbonsMissingMetabolite(model, metID):
    met = model.metabolites.get_by_id(metID)
    ProducingRct = getProducingReactions(model,metID)
    reaction = ProducingRct[0] # just need one producing reaction so you can stop at the first one

    stoiCoef_rct = abs(reaction.metabolites[met])
    reaction_reactants = reaction.reactants
    reaction_products = reaction.products
    reaction_products.remove(met)
    #remove the metabolite which we are looking that way you
    # can check with a reaction that does have the correct carbons if it works
    carbonInReactants = countCarbonInList(reaction,  reactionList= reaction_reactants)
    carbonInProducts = countCarbonInList(reaction, reactionList = reaction_products)

    CarbonsMissingMet = (carbonInReactants-carbonInProducts)/stoiCoef_rct
    return  CarbonsMissingMet

def count_carbon_in_formula(metFormula):
    if not metFormula or 'C' not in metFormula:  # if there is no carbon in the formula or the formula is of type NONE
        allCarbons = 0

    #nrOfCarbons = 0  # just in case something wierd happens
    else:
        splitFormula = re.split('(\d+)', metFormula)
        allCarbons = 0
        for j, element in enumerate(splitFormula):
            element = element.replace(' ','')
            if 'C' in element and len(element) == 1:
                allCarbons += int(splitFormula[j + 1])
            elif 'C' in element and len(element) > 1:
                posCarbon = element.index('C')
                if (posCarbon+1) == len(element): # e.g., ['RC' '11'] then there are 11 carbons to be counted
                    if splitFormula[j+1].isnumeric():
                        allCarbons += int(splitFormula[j + 1])
                    else:
                        allCarbons += 1 # the next element is the next atom so only add one carbon
                elif (posCarbon+1) != len(element):
                    if element[posCarbon + 1].isupper():  # for case like CH4 there is no 1 next to the C
                        allCarbons += 1  # only one carbon
                    elif element[posCarbon + 1].islower(): # for cases like Co (cobalt) or Cu
                        continue
            else:
                continue
    return allCarbons

def countCarbonInList(reaction, reactionList):
    carbonNrAll = []
    for met in reactionList:
        stoiCoef_i = abs(reaction.metabolites[met])
        metFormula = met.formula
        nrOfCarbons = count_carbon_in_formula(metFormula)
        carbonNrAll.append(nrOfCarbons*stoiCoef_i)
    totaalCarbonAtoms = sum(carbonNrAll)
    return  totaalCarbonAtoms

def getIDlistOfProducingMetabolties(model, metID):
    reactions = getProducingReactions(model , metID)
    ids = []
    coef = []
    r = reactions[0] # only interested in one reaction(doesn't matter which)
    metCoef = r.metabolites
    coefProduct = abs(metCoef[model.metabolites.get_by_id(metID)])
    reactants = r.reactants
    for met in reactants:
        ids.append(met.id)
        coef.append(abs(metCoef[met]))

    return ids, coef, coefProduct

def getProducingReactions(model, metID):
    met = model.metabolites.get_by_id(metID)
    frozenRct = met.reactions

    # get the reaction that produces the metabolite
    rcts = list(frozenRct)
    ProducingRct = []
    for rct in rcts:
        products = rct.products
        for prod in products:
            if prod == model.metabolites.get_by_id(metID):
                ProducingRct.append(rct)

    reaction = ProducingRct[0]
    return ProducingRct

def findCarbonsOfReaction(model,reactionID):
    reaction = model.reactions.get_by_id(reactionID)
    reactantList = reaction.reactants
    product =  reaction.products
    if len(product) > 1:
        raise ValueError("there should only be one product in the reaction, check")
    coefOfProduct = [abs(reaction.metabolites[i]) for i in product]
    coefOfReactants = [abs(reaction.metabolites[i]) for i in reactantList]
    carbonOfEachMolecule = []
    for met in reactantList:
        if met.formula:  # if it has formula count the carbons
            formula = met.formula
            nCarbon = count_carbon_in_formula(formula)
            carbonOfEachMolecule.append(nCarbon)
        else: #else go one reaction deeper to find the amount of carbons
            metID = met.id
            nCarbon = findCarbonsMissingMetabolite(model=model, metID=metID)
            if nCarbon > 0:
                carbonOfEachMolecule.append(nCarbon)
            else:
                switch = True
                # get missing metabolites and run
                name = model.metabolites.get_by_id(metID).name
                subMetabolites, coef, coefMet = getIDlistOfProducingMetabolties(model, metID)
                carbonSubReactions = []
                for subMetID in subMetabolites:
                    namesubMet = model.metabolites.get_by_id(subMetID).name
                    nSubCarbon = findCarbonsMissingMetabolite(model=model, metID=subMetID)
                    carbonSubReactions.append(nSubCarbon)

                control = [carbonSubReactions[i] >= 0 for i in range(len(carbonSubReactions))]
                if sum(control) == len(carbonSubReactions): #all reactions have carbon that makes them
                    coefficients = np.array(coef)
                    carbons = np.transpose(np.array(carbonSubReactions))
                    sumOfCarbons = np.matmul(coefficients, carbons)
                    nCarbon = sumOfCarbons/coefMet
                    carbonOfEachMolecule.append(nCarbon)
                else:
                    carbonOfEachMolecule.append(0)
                    # program to display warning a message
                    # displaying the warning message
                    warnings.warn('the carbons in metabolite {} could not be found, check manually. Metabolite ID is {}'.format(name, metID))


    coefficientsHeadReaction = np.array(coefOfReactants)
    carbonsHeadReaction = np.transpose(np.array(carbonOfEachMolecule))
    sumOfCarbons = np.matmul(coefficientsHeadReaction, carbonsHeadReaction)
    carbonProduct = sumOfCarbons / coefOfProduct
    return carbonProduct[0], carbonOfEachMolecule, coefOfReactants

def carbon_balance_in_out(modelLocation, metIDsMissingCarbon=None, tol = 0):
    if isinstance(modelLocation,str):
        model = cobra.io.read_sbml_model(modelLocation)
        modelName = modelLocation.split("\\")[-1]
        modelName = modelName.replace(".xml", "")
    else:
        model = modelLocation #then it will be the passed on model
        modelName = 'derp idk look if the model has a name in its struture'

    if metIDsMissingCarbon is None:
        metIDsMissingCarbon = []

    allRctIDMissingCarbon = []
    missingCarbonDict = {}
    if not isinstance(metIDsMissingCarbon, list):
        metIDsMissingCarbon = [metIDsMissingCarbon] #change ito a list if it is not
        if metIDsMissingCarbon: #if it is not empty
            for metID in metIDsMissingCarbon: #write a for loop to go over all the missing metabolites and find the producing reaction
                reactions = getProducingReactions(model= model, metID=metID)
                rctIDMissingCarbon = reactions[0] #only want the first reaction
                allRctIDMissingCarbon.append(rctIDMissingCarbon)
                missingCarbonDict.update({metID:rctIDMissingCarbon})

    df = model.summary()
    #exchangeRxn = model.exchanges
    uptake = df.uptake_flux
    secretion = df.secretion_flux

    dfUptake = carbonBalance(model= model , reactionDF= uptake, missingCarbonDict= missingCarbonDict, tol = tol)
    dfSecretion = carbonBalance(model= model , reactionDF= secretion, missingCarbonDict= missingCarbonDict, tol = tol)
    totalCarbonIn = sum(dfUptake['flux (gram-C/g-DW/h)'])
    CgramsOut = sum(dfSecretion['flux (gram-C/g-DW/h)'])
    CarbonBalance = abs(CgramsOut / totalCarbonIn)*100
    print('')
    print("{:2f} of the carbon is accounted for".format(CarbonBalance))
    fluxSecretion = dfSecretion['flux (gram-C/g-DW/h)']
    percentListSecretion = [round(abs(i/totalCarbonIn)*100,2) for i in fluxSecretion]
    dfSecretion['% carbon'] = percentListSecretion

    print('')
    print(modelName)
    print('')
    print(dfUptake)
    print('')
    print('')
    print(dfSecretion)
    return dfUptake, dfSecretion

def carbonBalance(model, reactionDF, missingCarbonDict, tol=0):
    metNamesAll = []
    carbonNrAll = []
    gramsCAll = []
    fluxAll = []
    for i, metID in enumerate(reactionDF.metabolite):
        # tttt = abs(reactionDF.flux[i])
        # t = tol
        if abs(reactionDF.flux[i]) > tol:
            fluxAll.append(reactionDF.flux[i])
            met = model.metabolites.get_by_id(metID)
            metName = met.name
            metNamesAll.append(metName)
            metFormula = met.formula
            if metID in missingCarbonDict.keys():
                rct = missingCarbonDict[metID]
                rctID = rct.id
                c = findCarbonsOfReaction(model=model, reactionID=rctID)
                nrOfCarbons = c[0]
            else:
                nrOfCarbons  = count_carbon_in_formula(metFormula)

            carbonNrAll.append(nrOfCarbons)
            gramsC = nrOfCarbons * 12 * reactionDF.flux[i]
            gramsCAll.append(gramsC)


    dictReactions = {'Metabolite': metNamesAll,
                     '# of C': carbonNrAll,
                     'flux (mmol/g-DW/h)':fluxAll ,
                     'flux (gram-C/g-DW/h)': gramsCAll}

    dataFrameReactions = pd.DataFrame(dictReactions)
    return dataFrameReactions
