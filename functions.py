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
def find_Carbons_Missing_Metabolite(model, metID):
    met = model.metabolites.get_by_id(metID)
    ProducingRct = getProducingReactions(model,metID)
    reaction = ProducingRct[0] # just need one producing reaction so you can stop at the first one

    stoiCoef_rct = abs(reaction.metabolites[met])
    reaction_reactants = reaction.reactants
    reaction_products = reaction.products
    reaction_products.remove(met)
    #remove the metabolite which we are looking that way you
    # can check with a reaction that does have the correct carbons if it works
    carbonInReactants = count_carbon_in_list(reaction,  reactionList= reaction_reactants)
    carbonInProducts = count_carbon_in_list(reaction, reactionList = reaction_products)

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

def count_carbon_in_list(reaction, reactionList):
    carbonNrAll = []
    for met in reactionList:
        stoiCoef_i = abs(reaction.metabolites[met])
        metFormula = met.formula
        nrOfCarbons = count_carbon_in_formula(metFormula)
        carbonNrAll.append(nrOfCarbons*stoiCoef_i)
    totaalCarbonAtoms = sum(carbonNrAll)
    return  totaalCarbonAtoms

def get_IDlist_Of_Producing_Metabolties(model, metID):
    reactions = get_Producing_Reactions(model , metID)
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

def get_Producing_Reactions(model, metID):
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

def find_carbons_of_reaction(model,reactionID):
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
            nCarbon = find_Carbons_Missing_Metabolite(model=model, metID=metID)
            if nCarbon > 0:
                carbonOfEachMolecule.append(nCarbon)
            else:
                switch = True
                # get missing metabolites and run
                name = model.metabolites.get_by_id(metID).name
                subMetabolites, coef, coefMet = get_IDlist_Of_Producing_Metabolties(model, metID)
                carbonSubReactions = []
                for subMetID in subMetabolites:
                    namesubMet = model.metabolites.get_by_id(subMetID).name
                    nSubCarbon = find_Carbons_Missing_Metabolite(model=model, metID=subMetID)
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

def carbon_balance(model, reactionDF, missingCarbonDict, tol=0):
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
                c = find_carbons_of_reaction(model=model, reactionID=rctID)
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
                reactions = get_Producing_Reactions(model= model, metID=metID)
                rctIDMissingCarbon = reactions[0] #only want the first reaction
                allRctIDMissingCarbon.append(rctIDMissingCarbon)
                missingCarbonDict.update({metID:rctIDMissingCarbon})

    df = model.summary()
    #exchangeRxn = model.exchanges
    uptake = df.uptake_flux
    secretion = df.secretion_flux

    dfUptake = carbon_balance(model= model , reactionDF= uptake, missingCarbonDict= missingCarbonDict, tol = tol)
    dfSecretion = carbon_balance(model= model , reactionDF= secretion, missingCarbonDict= missingCarbonDict, tol = tol)
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

# originaly from the file print_model_2_excel

def string_reactions(reaction, case='names'):
    """ returns the reaction as a string with the stoichiometry and prints it to the terminal
    """

    rxn = reaction
    reactants = rxn.reactants
    reactantStr = ''
    stoiFactor = rxn.metabolites

    for metReac in reactants:
        if case == 'names':
            reactName = metReac.name
        else:
            reactName = metReac.formula
            if reactName == '':
                reactName = "###"

        stoi = stoiFactor[metReac]
        reactantStr += '{} {} + '.format(stoi, reactName)

    productStr = ""
    products = rxn.products
    for metProd in products:
        if case == 'names':
            prodName = metProd.name
        else:
            prodName = metProd.formula
            if prodName == '':
                prodName = "###"  # so it is easier to see if it is missing

        stoi = stoiFactor[metProd]
        productStr += '{} {} + '.format(stoi, prodName)

    reactionStr = "{} = {}".format(reactantStr, productStr)
    return reactionStr

def print_SBML_info_2_excel(modelName, idMissingCarbon=None):
    """
    This function imports information of the model to an Excel file to check the missing carbon in the metabolic reactions
    the percent of missing carbon can also  be attributed to the reactions to check on their importants of the reactions to
    the outcome of the FBA

    Inputs:
    modelName: str name or location of the model, optional

    idMissingCarbon: if you know a metabolite does not have carbon in it's formula you can specify it as a string and
    the code will estimate it for you. default = None

    blanace elements: list of elements to balance. default = ['C', 'O', 'H']

    output:
    Excel files with info
    """
    if isinstance(modelName, str):
        modelLocation = get_location(modelName)
        model = cobra.io.read_sbml_model(modelLocation)
        saveName = modelName
        saveName = saveName.replace('.xml', '')
        saveName = '{}_analysis.xlsx'.format(saveName)
    else:
        model = modelName
        saveName = model.name
        saveName = saveName.replace('.xml', '')
        saveName = '{}_analysis.xlsx'.format(saveName)

    StoiMatrixDF = cobra.util.create_stoichiometric_matrix(model, array_type='DataFrame')  # [:,0:posExchangeRxn]
    # print(StoiMatrixDF)

    FBA = model.optimize()
    fluxArray = FBA.fluxes  # [0:posExchangeRxn] #.to_numpy()

    # fluxArray.drop(ExRxn, inplace=True)
    # rxnFluxes = np.dot(StoiMatrixDF, fluxArray)
    # metabolicFluxDF = pd.DataFrame(metabolicFlux, index=StoiMatrixDF.index)
    # print(metabolicFluxDF)

    metabolites = model.metabolites
    metID = []
    metName = []
    metNameWithSpaces = []
    strRxns = []
    fluxRxn = []
    stoiMetMissingFormula = []
    carbonCount = []
    for met in metabolites:
        formula = met.formula
        carbon = countCarbonInFormula(metFormula=formula)
        carbonCount.append(carbon)
        if not formula:
            metID.append(met.id)
            metName.append(met.name)
            metNameWithSpaces.append(met.name)
            frozenRct = met.reactions
            # get the reaction that produces the metabolite
            allRxn = list(frozenRct)

            extraNames = [''] * (len(allRxn) - 1)
            metNameWithSpaces = metNameWithSpaces + extraNames
            # metName = metName + extraNames

            # helpNames = metName + extraNames
            # metNameWithSpaces = metNameWithSpaces + helpNames
            # metName.append(extraNames)

            for rxn in allRxn:
                reactionStr = string_reactions(reaction=rxn)
                strRxns.append(reactionStr)

                fluxRxn.append(rxn.flux)
                stoiFactor = rxn.metabolites
                stoiMetMissingFormula.append(stoiFactor[met])

    DictCarbons = {'ID metabolite': list(StoiMatrixDF.index), '# Carbons': carbonCount}
    CarbonsDF = pd.DataFrame(DictCarbons)

    DictMetabolites = {'ID': metID, 'Name': metName}
    DFmetabolites = pd.DataFrame(DictMetabolites)
    # print(DFmetabolites)

    DictMetRnx = {'Name': metNameWithSpaces, 'Reactions': strRxns, 'Flux': fluxRxn,
                  'Stoichiometry': stoiMetMissingFormula}
    DFmetRnx = pd.DataFrame(DictMetRnx)
    # print(DFmetRnx)

    # find the ingoing and outgoing fluxes
    inputDF, outputDF = carbon_balance_in_out(modelLocation=model, metIDsMissingCarbon=idMissingCarbon, tol=0.0001)

    # calculate the percentage of carbon at goes missing in each reaction (exculed the exchage reactions?)
    # exclude the transfer (exchange reactions) reactions
    exchRxn = model.exchanges
    # transform id's to a list
    keysListExRxn = [rx.id for rx in exchRxn]

    # to drop or not?
    fluxArray.drop(keysListExRxn, inplace=True)
    fluxArrayNp = np.array(fluxArray).reshape(-1, 1)
    st = StoiMatrixDF.drop(keysListExRxn, axis='columns')  # ,inplace=False)
    # st = StoiMatrixDF

    carbonArray = np.array(carbonCount)
    carbonArray = carbonArray.reshape(-1, 1)
    StoiMatrixTranspose = np.transpose(st)

    carbonBalance = np.dot(StoiMatrixTranspose, carbonArray)  # in mols of Carbon
    # multiply by 12 gC/mol
    carbonMassBalance = carbonBalance * 12
    carbonMissing = carbonMassBalance * fluxArrayNp
    rxnIDreduced = list(fluxArray.index)

    DictMissingCarbonALL = {'Reaction Id': rxnIDreduced,
                            'carbon (gCarbon)': carbonMissing.reshape(len(carbonMissing), )}
    DFMissingCarbonALL = pd.DataFrame(DictMissingCarbonALL)

    unbalancedId = []
    unbalanced = []
    rxnUnblanaced = []
    for i, c in enumerate((carbonMissing)):
        if abs(c) > 0:
            unbalancedId.append(rxnIDreduced[i])
            unbalanced.append(c[0])
            reactionStr = string_reactions(model.reactions.get_by_id(rxnIDreduced[i]))
            rxnUnblanaced.append(reactionStr)
    DictMissingCarbon = {'Reaction Id': unbalancedId, 'carbon (gCarbon)': unbalanced, 'reaction': rxnUnblanaced}
    DFUnbalancedCarbonReactions = pd.DataFrame(DictMissingCarbon)

    # UnbalancedCarbonReactions
    # carbonStoiMatrix
    # carbonStoiMatrix = np.dot(np.array(carbonCount),np.transpose(StoiMatrixDF))
    # numpy.array(carbonCount)

    with pd.ExcelWriter(saveName) as writer:
        StoiMatrixDF.to_excel(writer, sheet_name='StoichiometricMatrix')
        fluxArray.to_excel(writer, sheet_name='ReactionFluxes')
        CarbonsDF.to_excel(writer, sheet_name='CarbonsPerMetbolite')
        DFmetRnx.to_excel(writer, sheet_name='MissingFormulaReactions')
        DFUnbalancedCarbonReactions.to_excel(writer, sheet_name='Carbon Balance')
        inputDF.to_excel(writer, sheet_name='Carbon Input')
        outputDF.to_excel(writer, sheet_name='Carbon output')


########################################################################################################################
# ============================================================================================================
# formulating the superstructure
# ============================================================================================================
########################################################################################################################