import re
import warnings

import cobra
import cobra.io
import numpy as np
import pandas as pd

from f_usefull_functions import *


########################################################################################################################
# ============================================================================================================
# SBML screening
# ============================================================================================================
########################################################################################################################

# originaly from the file f_find_carbons
def find_Carbons_Missing_Metabolite(model, metID):
    met = model.metabolites.get_by_id(metID)
    ProducingRct = get_Producing_Reactions(model, metID)
    reaction = ProducingRct[0]  # just need one producing reaction so you can stop at the first one

    stoiCoef_rct = abs(reaction.metabolites[met])
    reaction_reactants = reaction.reactants
    reaction_products = reaction.products
    reaction_products.remove(met)
    # remove the metabolite which we are looking that way you
    # can check with a reaction that does have the correct carbons if it works
    carbonInReactants = count_element_in_list(reaction, reactionList=reaction_reactants, element='C')
    carbonInProducts = count_element_in_list(reaction, reactionList=reaction_products, element='C')

    CarbonsMissingMet = (carbonInReactants - carbonInProducts) / stoiCoef_rct
    return CarbonsMissingMet


def count_carbon_in_formula(metFormula):
    if not metFormula or 'C' not in metFormula:  # if there is no carbon in the formula or the formula is of type NONE
        allCarbons = 0

    # nrOfCarbons = 0  # just in case something wierd happens
    else:
        splitFormula = re.split('(\d+)', metFormula)
        allCarbons = 0
        for j, element in enumerate(splitFormula):
            element = element.replace(' ', '')
            if 'C' in element and len(element) == 1:
                allCarbons += int(splitFormula[j + 1])
            elif 'C' in element and len(element) > 1:
                posCarbon = element.index('C')
                if (posCarbon + 1) == len(element):  # e.g., ['RC' '11'] then there are 11 carbons to be counted
                    if splitFormula[j + 1].isnumeric():
                        allCarbons += int(splitFormula[j + 1])
                    else:
                        allCarbons += 1  # the next element is the next atom so only add one carbon
                elif (posCarbon + 1) != len(element):
                    if element[posCarbon + 1].isupper():  # for case like CH4 there is no 1 next to the C
                        allCarbons += 1  # only one carbon
                    elif element[posCarbon + 1].islower():  # for cases like Co (cobalt) or Cu
                        continue
            else:
                continue
    return allCarbons


def count_atom_in_formula(metabolite, atom):
    if atom == 'e-':
        count = metabolite.charge
    else:
        try:
            count = metabolite.elements[atom]
        except:
            count = 0  # if the atom is not a key, it is not in the formula and therefore zero
    return count


def count_element_in_list(reaction, reactionList, element):
    elementNrAll = []
    for met in reactionList:
        stoiCoef_i = abs(reaction.metabolites[met])
        # metFormula = met.formula
        nrOfelements = count_atom_in_formula(metabolite=met, atom=element)
        elementNrAll.append(nrOfelements * stoiCoef_i)
    totaalElementAtoms = sum(elementNrAll)
    return totaalElementAtoms


def get_IDlist_Of_Producing_Metabolties(model, metID):
    reactions = get_Producing_Reactions(model, metID)
    ids = []
    coef = []
    r = reactions[0]  # only interested in one reaction(doesn't matter which)
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


def find_carbons_of_reaction(model, reactionID):
    reaction = model.reactions.get_by_id(reactionID)
    reactantList = reaction.reactants
    product = reaction.products
    if len(product) > 1:
        raise ValueError("there should only be one product in the reaction, check")
    coefOfProduct = [abs(reaction.metabolites[i]) for i in product]
    coefOfReactants = [abs(reaction.metabolites[i]) for i in reactantList]
    carbonOfEachMolecule = []
    for met in reactantList:
        if met.formula:  # if it has formula count the carbons
            formula = met.formula
            nCarbon = count_atom_in_formula(metabolite=met, atom='C')
            carbonOfEachMolecule.append(nCarbon)
        else:  # else go one reaction deeper to find the amount of carbons
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
                if sum(control) == len(carbonSubReactions):  # all reactions have carbon that makes them
                    coefficients = np.array(coef)
                    carbons = np.transpose(np.array(carbonSubReactions))
                    sumOfCarbons = np.matmul(coefficients, carbons)
                    nCarbon = sumOfCarbons / coefMet
                    carbonOfEachMolecule.append(nCarbon)
                else:
                    carbonOfEachMolecule.append(0)
                    # program to display warning a message
                    # displaying the warning message
                    warnings.warn(
                        'the carbons in metabolite {} could not be found, check manually. Metabolite ID is {}'.format(
                            name, metID))

    coefficientsHeadReaction = np.array(coefOfReactants)
    carbonsHeadReaction = np.transpose(np.array(carbonOfEachMolecule))
    sumOfCarbons = np.matmul(coefficientsHeadReaction, carbonsHeadReaction)
    carbonProduct = sumOfCarbons / coefOfProduct
    return carbonProduct[0], carbonOfEachMolecule, coefOfReactants


def carbon_balance(model, reactionDF, missingCarbonDict, tol=0.0001):
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
                nrOfCarbons = count_atom_in_formula(metabolite=met, atom='C')

            carbonNrAll.append(nrOfCarbons)
            gramsC = nrOfCarbons * 12 * reactionDF.flux[i]
            gramsCAll.append(gramsC)

    dictReactions = {'Metabolite': metNamesAll,
                     '# of C': carbonNrAll,
                     'flux (mmol/g-DW/h)': fluxAll,
                     'flux (gram-C/g-DW/h)': gramsCAll}

    dataFrameReactions = pd.DataFrame(dictReactions)
    return dataFrameReactions


def carbon_balance_in_out(modelLocation, metIDsMissingCarbon=None, tol=0.0001):
    if isinstance(modelLocation, str):
        model = cobra.io.read_sbml_model(modelLocation)
        modelName = modelLocation.split("\\")[-1]
        modelName = modelName.replace(".xml", "")
    else:
        model = modelLocation  # then it will be the passed on as a COBRA model
        modelName = 'derp idk look if the model has a name in its struture'

    if metIDsMissingCarbon is None:
        metIDsMissingCarbon = []

    allRctIDMissingCarbon = []
    missingCarbonDict = {}
    if not isinstance(metIDsMissingCarbon, list):
        metIDsMissingCarbon = [metIDsMissingCarbon]  # change ito a list if it is not
        if metIDsMissingCarbon:  # if it is not empty
            for metID in metIDsMissingCarbon:  # write a for loop to go over all the missing metabolites and find the producing reaction
                reactions = get_Producing_Reactions(model=model, metID=metID)
                rctIDMissingCarbon = reactions[0]  # only want the first reaction
                allRctIDMissingCarbon.append(rctIDMissingCarbon)
                missingCarbonDict.update({metID: rctIDMissingCarbon})

    df = model.summary()
    # exchangeRxn = model.exchanges
    uptake = df.uptake_flux
    secretion = df.secretion_flux

    dfUptake = carbon_balance(model=model, reactionDF=uptake, missingCarbonDict=missingCarbonDict, tol=tol)
    dfSecretion = carbon_balance(model=model, reactionDF=secretion, missingCarbonDict=missingCarbonDict, tol=tol)
    totalCarbonIn = sum(dfUptake['flux (gram-C/g-DW/h)'])
    CgramsOut = sum(dfSecretion['flux (gram-C/g-DW/h)'])
    CarbonBalance = abs(CgramsOut / totalCarbonIn) * 100
    print('')
    print("{:2f} of the carbon is accounted for".format(CarbonBalance))
    fluxSecretion = dfSecretion['flux (gram-C/g-DW/h)']
    percentListSecretion = [round(abs(i / totalCarbonIn) * 100, 2) for i in fluxSecretion]
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
def string_reactions(reaction, case='names', printFlux=False):
    """ returns the reaction as a string with the stoichiometry and prints it to the terminal
    """

    rxn = reaction
    if printFlux:
        try:
            flux = rxn.flux
        except:
            flux = 'Undetermined'
            print('use model.optimize() before calling this function to get the fluxes')
    else:
        flux = 'Undetermined'

    reactants = rxn.reactants
    reactantStr = ''
    stoiFactor = rxn.metabolites

    for metReac in reactants:
        if case == 'names':
            reactName = metReac.name
        elif case == 'formulas':
            reactName = metReac.formula
            if reactName == '':
                reactName = "###"
        elif case == 'id':
            reactName = metReac.id
        else:
            raise Exception("the option 'case' can only be 'names' or 'formulas' ")

        stoi = stoiFactor[metReac]
        reactantStr += '{} {} + '.format(stoi, reactName)

    productStr = ""
    products = rxn.products
    for metProd in products:
        if case == 'names':
            prodName = metProd.name

        elif case == 'formulas':
            prodName = metProd.formula
            if prodName == '':
                prodName = "###"  # so it is easier to see if it is missing

        elif case == 'id':
            prodName = metProd.id

        stoi = stoiFactor[metProd]
        productStr += '{} {} + '.format(stoi, prodName)

    reactionStr = "{} = {}".format(reactantStr, productStr)
    return reactionStr, flux


def print_all_rxn_of_metabolite(metabolite, case='names', printFlux=False):
    allRxn = list(metabolite.reactions)
    listRxn = []
    for rxn in allRxn:
        strRnx, flux = string_reactions(reaction=rxn, case=case, printFlux=printFlux)
        listRxn.append(strRnx)
        print(strRnx)
        if printFlux:
            print('the flux of reaction {} is : {} mmol/gDW/h \n '
                  'the bounds are {} mmol/gDW/h \n'
                  'the compartment(s) of the reactions are {} \n'.format(rxn.id, flux, rxn.bounds, rxn.compartments))
    return  listRxn

def find_maintanace_reaction(metabolite, case='names', printFlux=False):
    allRxn = list(metabolite.reactions)
    listRxn = []
    for rxn in allRxn:
        len_reaction = len(rxn.products) + len(rxn.reactants)
        if len_reaction == 4:
            strRnx, flux = string_reactions(reaction=rxn, case=case, printFlux=printFlux)
            listRxn.append(strRnx)
            print(strRnx)
            if printFlux:
                print('the flux of reaction {} is : {} mmol/gDW/h \n '
                      'the bounds are {} mmol/gDW/h \n'
                      'the compartment(s) of the reactions are {} \n'.format(rxn.id, flux, rxn.bounds, rxn.compartments))
    return  listRxn

def find_unbalanced_rxn_of_element(model, stoiMatrix, fluxArray, element, elementCount):
    """ finds the reactions where the given element is unbalanced
    Params:
        model: COBRA model
        stoiMatrix (DF): Pandas Dataframe of the stoichiometric matrix
        fluxArray (Series): Pandas series with the flux of all reactions'
        element (str): can either be 'C', 'O', 'H' or 'e-'
        elementCount (list): list of # carbons per metabolite of the model

    Returns:
          DFUnbalancedElementReactions (DF): A dataframe with all the unbalanced reactions
    """

    MM = {'C': 12, 'O': 15.99, 'H': 1, 'e-': 1}
    fluxArrayNp = np.array(fluxArray).reshape(-1, 1)  # reshape the flux array
    elementArray = np.array(elementCount)
    elementArray = elementArray.reshape(-1, 1)
    StoiMatrixTranspose = np.transpose(stoiMatrix)

    elementBalance = np.dot(StoiMatrixTranspose, elementArray)  # in mols of element
    # multiply by MM g[C,O,H]/mol
    elementMassBalance = elementBalance * MM[element]
    elementMissing = elementMassBalance * fluxArrayNp
    rxnIDreduced = list(fluxArray.index)

    unbalancedId = []
    unbalanced = []
    rxnUnblanaced = []
    for i, c in enumerate((elementMissing)):
        if abs(c) > 0:
            unbalancedId.append(rxnIDreduced[i])
            unbalanced.append(c[0])
            reactionStr = string_reactions(model.reactions.get_by_id(rxnIDreduced[i]))
            rxnUnblanaced.append(reactionStr)
    DictMissingElementInRxn = {'Reaction Id': unbalancedId, '{}(g{})'.format(element, element): unbalanced,
                               'reaction': rxnUnblanaced}
    DFUnbalancedElementReactions = pd.DataFrame(DictMissingElementInRxn)

    return DFUnbalancedElementReactions


def get_list_metabolite_ids_names(model):
    """returns a nx2 arrary whith the first column the IDs of the metaboliets and the second with the names of the
    metabolites

    Input = model (COBRA)
    """

    metabolites = model.metabolites
    names = []
    ids = []
    for met in metabolites:
        names.append(met.name)
        ids.append(met.id)

    DFidName = pd.DataFrame({'IDs': ids, 'Names': names})
    return DFidName


def print_SBML_info_2_excel(modelName, idMissingCarbon=None, saveName=None, tolerance=0.0001,
                            print2Excel=True):
    """
    This function imports information of the model to an Excel file to check the missing carbon in the metabolic reactions
    the percent of missing carbon can also  be attributed to the reactions to check on their importants of the reactions to
    the outcome of the FBA

    Inputs:
    modelName: str name or location of the model, optional

    idMissingCarbon: if you know a metabolite does not have carbon in its formula you can specify the ID as a string and
    the code will estimate it for you. default = None (a bit hacky however)

    blanace elements: list of elements to balance. default = ['C', 'O', 'H']

    output:
    Excel files with info
    """
    if isinstance(modelName, str):  # the model still needs to be retrived from the file
        modelLocation = get_location(modelName)
        model = cobra.io.read_sbml_model(modelLocation)
        name = modelName
    else:  # the COBRA model is given as the input
        model = modelName
        name = model.name

    if saveName is None:
        saveName = name
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
    metIdWithSpaces = []
    strRxns = []
    idRxns = []
    fluxRxn = []
    stoiMetMissingFormula = []
    carbonCount = []
    oxygenCount = []
    hydrogenCount = []
    chargeCount = []

    for met in metabolites:
        formula = met.formula
        carbon = count_atom_in_formula(metabolite=met, atom='C')
        oxygen = count_atom_in_formula(metabolite=met, atom='O')
        hydrogen = count_atom_in_formula(metabolite=met, atom='H')
        charge = count_atom_in_formula(metabolite= met, atom='e-')

        carbonCount.append(carbon)
        oxygenCount.append(oxygen)
        hydrogenCount.append(hydrogen)
        chargeCount.append(charge)

        if not formula:
            metID.append(met.id)
            metName.append(met.name)
            metNameWithSpaces.append(met.name)
            metIdWithSpaces.append(met.id)
            frozenRct = met.reactions
            # get the reaction that produces the metabolite
            allRxn = list(frozenRct)

            extraNames = [''] * (len(allRxn) - 1)
            metNameWithSpaces = metNameWithSpaces + extraNames
            metIdWithSpaces = metIdWithSpaces + extraNames
            # metName = metName + extraNames

            # helpNames = metName + extraNames
            # metNameWithSpaces = metNameWithSpaces + helpNames
            # metName.append(extraNames)

            for rxn in allRxn:
                reactionStr = string_reactions(reaction=rxn)
                strRxns.append(reactionStr)
                idRxns.append(rxn.id)
                fluxRxn.append(rxn.flux)
                stoiFactor = rxn.metabolites
                stoiMetMissingFormula.append(stoiFactor[met])

    DictCarbons = {'ID metabolite': list(StoiMatrixDF.index), '# Carbons': carbonCount}
    CarbonsDF = pd.DataFrame(DictCarbons)

    # DictMetabolites = {'ID': metID, 'Name': metName}
    # DFmetabolites = pd.DataFrame(DictMetabolites)
    # print(DFmetabolites)

    DictMetRnx = {'Name': metNameWithSpaces, 'ID': metIdWithSpaces, 'Reactions': strRxns, 'Rxn ID': idRxns,
                  'Flux': fluxRxn,
                  'Stoichiometry': stoiMetMissingFormula}
    DFmetRnx = pd.DataFrame(DictMetRnx)
    # print(DFmetRnx)

    # find the ingoing and outgoing fluxes
    inputDF, outputDF = carbon_balance_in_out(modelLocation=model, metIDsMissingCarbon=idMissingCarbon, tol=tolerance)

    # calculate the mass of carbon at goes missing in each reaction (excluded the exchange reactions?)
    # exclude the transfer (exchange reactions) reactions
    exchRxn = model.exchanges
    # transform id's to a list
    keysListExRxn = [rx.id for rx in exchRxn]

    #  drop the exchange reactions, they are never balanced so don't bother looking at them
    fluxArray.drop(keysListExRxn, inplace=True)
    st = StoiMatrixDF.drop(keysListExRxn, axis='columns')  # ,inplace=False)

    DFcarbon = find_unbalanced_rxn_of_element(model=model, stoiMatrix=st, fluxArray=fluxArray,
                                              element='C', elementCount=carbonCount)

    DFoxygen = find_unbalanced_rxn_of_element(model=model, stoiMatrix=st, fluxArray=fluxArray,
                                              element='O', elementCount=oxygenCount)

    DFhydrogen = find_unbalanced_rxn_of_element(model=model, stoiMatrix=st, fluxArray=fluxArray,
                                                element='H', elementCount=hydrogenCount)

    DFcharge =  find_unbalanced_rxn_of_element(model=model, stoiMatrix=st, fluxArray=fluxArray,
                                                element='e-', elementCount=chargeCount)

    DFMetIdNames = get_list_metabolite_ids_names(model)

    # get a list of metabolites that can be exchanged (ids and names)
    exchangeMetID = []
    exchangeName = []
    exchangeRxnID = []
    exchangeFLux = []
    for exRxn in exchRxn:
        exchangeRxnID.append(exRxn.id)
        exchangeFLux.append(exRxn.flux)

        for exMet in exRxn.metabolites:
            exchangeMetID.append(exMet.id)
            exchangeName.append(exMet.name)

    exchangeDict = {'Name': exchangeName, 'metabolite id': exchangeMetID, 'reaction id': exchangeRxnID,
                    'flux': exchangeFLux}
    DFexchange = pd.DataFrame(data=exchangeDict)

    if print2Excel:
        saveLocation = r'C:\Users\lucas\PycharmProjects\Alquimia\SBML screening\Excel analysis\{}'.format(saveName)
        with pd.ExcelWriter(saveLocation) as writer:
            StoiMatrixDF.to_excel(writer, sheet_name='Stoichiometric_matrix')
            fluxArray.to_excel(writer, sheet_name='Reaction_fluxes')
            DFexchange.to_excel(writer, sheet_name='Exchange_reactions')
            CarbonsDF.to_excel(writer, sheet_name='Carbons_Per_Metbolite')
            DFmetRnx.to_excel(writer, sheet_name='Missing_Formula_Reactions')
            DFMetIdNames.to_excel(writer, sheet_name='ID_2_name')
            DFcarbon.to_excel(writer, sheet_name='Carbon_Balance')
            DFoxygen.to_excel(writer, sheet_name='Oxygen_Balance')
            DFhydrogen.to_excel(writer, sheet_name='Hydrogen_Balance')
            DFcharge.to_excel(writer, sheet_name='Electron_Balance')
            inputDF.to_excel(writer, sheet_name='Carbon_Input')
            outputDF.to_excel(writer, sheet_name='Carbon_Output')


# functions to help fix the models
def balance_element(reaction, element):
    products = reaction.products
    reactants = reaction.reactants

    Cproducts = 0
    for prod in products:
        if element in prod.elements:
            nC = prod.elements[element]
        elif element == 'e-':
            nC = prod.charge
        else:
            nC = 0
        stoiFactor = reaction.metabolites[prod]
        Cproducts += nC * stoiFactor

    Creactants = 0
    for react in reactants:
        if element in react.elements:
            nC = react.elements[element]
        elif element == 'e-':
            nC = react.charge
        else:
            nC = 0
        stoiFactor = reaction.metabolites[react]
        Creactants += nC * stoiFactor
    return Creactants, Cproducts


def check_reaction(model, reactionID):
    """ Prints the reactions equations and look if the elements are balanced """
    # elements to look at
    elements = ['C', 'H', 'O', 'e-']
    reaction = model.reactions.get_by_id(reactionID)
    reactionEq = string_reactions(reaction, case='formulas')
    reactionEqNames = string_reactions(reaction=reaction, case='names')

    print(reactionEq)
    print(reactionEqNames)
    print(reaction)  # reaction with id codes
    for elm in elements:
        reac, prod = balance_element(reaction=reaction, element=elm)
        missing = prod + reac

        print('{} left of reaction = {}'.format(elm, reac))
        print('{} right of reaction = {}'.format(elm, prod))
        print('{} missing = {}'.format(elm, missing))
        print('')

def check_charge(model,reactionID):
    """ Checks which molecules are charged in the reactants and products to more easily identify missing charge """

    # get the reaction
    reaction = model.reactions.get_by_id(reactionID)

    # print reactions
    reactionEq = string_reactions(reaction, case='formulas')
    reactionEqNames = string_reactions(reaction=reaction, case='names')
    print(reactionEq)
    print(reactionEqNames)

    # get products and reactans
    products = reaction.products
    reactants = reaction.reactants

    # iniciate lists
    prodList = []
    chargeProdList = []
    stoiProdList = []
    totalChargeProdList = []
    for prod in products:
        # get the name
        prodList.append(prod.name)

        # get the charge
        charge = prod.charge
        chargeProdList.append(charge)

        # get the stoichiometric value
        stoiFactor = reaction.metabolites[prod]
        stoiProdList.append(reaction.metabolites[prod])

        # get the total charge
        totalCharge = charge * stoiFactor
        totalChargeProdList.append(totalCharge)

    # make the dataframe
    dfProducts = pd.DataFrame({'Product':prodList, 'Charge': chargeProdList, 'Stoichiometry':stoiProdList,
                                   'Total charge': totalChargeProdList})

    reacList = []
    chargeReacList = []
    stoiReacList = []
    totalChargeReacList = []
    for react in reactants:
        # get the name
        reacList.append(react.name)

        # get the charge
        charge = react.charge
        chargeReacList.append(charge)

        # get the stoichiometric value
        stoiFactor = reaction.metabolites[react]
        stoiReacList.append(reaction.metabolites[react])

        # get the total charge
        totalCharge = charge * stoiFactor
        totalChargeReacList.append(totalCharge)

    # make the dataframe
    dfReactants = pd.DataFrame({'Product': reacList, 'Charge': chargeReacList, 'Stoichiometry': stoiReacList ,
                                   'Total charge': totalChargeReacList})

    print('')
    print('the charge of the individual reactants are:')
    print(dfReactants)
    print('')
    print('total charge Reactants: ', sum(totalChargeReacList))

    print('')
    print('the charge of the individual products are:')
    print(dfProducts)
    print('')
    print('total charge Products: ', sum(totalChargeProdList) )
    print('')

    return sum(totalChargeProdList) + sum(totalChargeReacList)



def get_metabolites_whith_missing_formula(model):
    '''returns all metabolites without a formula'''
    metabolites = model.metabolites
    missing = []
    for met in metabolites:
        if not met.formula:
            missing.append(met)
    return missing


def estimate_biomass_formula_from_products(model, reactionID):
    """calculates formula for biomass based on the reactants when no explicit biomass metabolite is given
    (carbon hydrogen and oxygen)

    input:
        model (model)
        reaction (reactionID): reaction ID to estimnate the formula from

    returns:
        formula (str): estimate of the foprmula

    """
    rxn = model.reactions.get_by_id(reactionID)
    stoichiometry = rxn.metabolites
    rxnReactants = rxn.reactants

    elements = ['C', 'H', 'O']
    missingDict = {}
    for ele in elements:
        elementInProducts = count_element_in_list(reaction=rxn, reactionList=rxnReactants, element=ele)
        # also accounts for the stoichiometry
        missingDict.update({ele: elementInProducts})

    formulaEstimate = 'C{}H{}O{}'.format(missingDict['C'], missingDict['H'], missingDict['O'])
    MW = 12* missingDict['C'] + missingDict['H'] + 15.99 * missingDict['O']
    return formulaEstimate, MW


def estimate_formula(model, estimationDict):
    '''
    estimates the formula (CxHyOz) based of the stoichiometry of the given reaction in estimationDict
    inputs:
        model (CobraPy model)
        estimationDict (Dict): dictionary { metabolite you want to estimate the formula of : reaction to estimate from}
    '''
    formulaDict = {}
    for metID in estimationDict:
        reaction = model.reactions.get_by_id(estimationDict[metID])
        metabolite = model.metabolites.get_by_id(metID)
        # reaction = ProducingRct[0]  # just need one producing reaction so you can stop at the first one

        stoiCoef_rct = abs(reaction.metabolites[metabolite])
        reaction_reactants = reaction.reactants
        reaction_products = reaction.products
        # reaction_products.remove(metabolite)
        # remove the metabolite which we are looking that way you
        # can check with a reaction that does have the correct carbons if it works
        elements = ['C', 'H', 'O']
        missingDict = {}
        for ele in elements:
            elementsInReactants = count_element_in_list(reaction, reactionList=reaction_reactants, element=ele)
            elementsInProducts = count_element_in_list(reaction, reactionList=reaction_products, element=ele)
            elementMissingMet = (elementsInReactants - elementsInProducts) / stoiCoef_rct
            if elementMissingMet < 0:
                warnings.warn('the element {} is overproduced in reaction '
                              '{} plz check it out'.format(ele, reaction.id), category=UserWarning)
                elementMissingMet = 0
            if isinstance(elementMissingMet, float):
                elementMissingMet = round(elementMissingMet, 2)
            missingDict.update({ele: elementMissingMet})

        formulaEstimate = 'C{}H{}O{}'.format(round(missingDict['C']), round(missingDict['H']), round(missingDict['O']))
        formulaDict.update({metID: formulaEstimate})
    return formulaDict


def count_missing_formulas_in_rxn(rxnId, model):
    reaction = model.reactions.get_by_id(rxnId)
    # Get the list of metabolites involved in the reaction
    metabolites = reaction.metabolites
    # Initialize a counter for missing formulas
    countMissingFormulas = 0

    # Iterate over the metabolites
    for metabolite in metabolites:
        # Check if the metabolite formula is missing
        if not metabolite.formula:  # if the string is empty
            # If the formula is missing, increment the counter
            countMissingFormulas += 1

    # Return the count of metabolites with missing formulas
    return countMissingFormulas


def fix_missing_formulas(model, fixDict, maxIterations=10):
    '''
    fixes metaboliets that don't have a formula by looking at the reaction and counting the missing elements C H and O
    this is a very rough estimation!!

    Inputs:
    model (COBRA model): model
    fixDict (Dict): dictionary with metabolites as keys and rxn IDs as values (the reaction id you want to fix the metabolite with)

    output: model (COBRA model)
    '''
    estimateFormulas = {}
    stopCriteria = True
    iteration = 0
    while stopCriteria and iteration < maxIterations:
        toFix = {}
        notSolved = {}
        countCompleet = 0
        iteration += 1
        for metId in fixDict:
            rxnId = fixDict[metId]
            nMissing = count_missing_formulas_in_rxn(rxnId=rxnId, model=model)
            if nMissing == 1:
                toFix.update({metId: rxnId})

            elif nMissing == 0:
                countCompleet += 1

            else:
                notSolved.update({metId: rxnId})

        a = countCompleet
        b = len(fixDict)
        if countCompleet == len(fixDict):  # all the reactions have formulas now
            stopCriteria = False
        elif iteration == (maxIterations - 1):
            print(notSolved)
            raise Exception(
                'the above metabolite formulas could not be found, consider using a different reaction to find them')

        else:
            # estimate the formulas with reactions only missing 1 formula
            batchEstimastes = estimate_formula(model=model, estimationDict=toFix)
            estimateFormulas.update(batchEstimastes)
            # update the metabolites of the model
            for metId in batchEstimastes:
                met2change = model.metabolites.get_by_id(metId)
                formula = batchEstimastes[metId]
                met2change.formula = formula

    return model, estimateFormulas


def ATP_Biomass_Ratio(model, biomassRxnID, ATPmetID, modelName = 'NONE', printResults = True):
    """
    prints the flux of biomass and the total flux of produced ATP
    Params:
        model (COBRA model): the cobra model
        biomassRxnID (str):  the exchange reaction ID of biomass
        ATPmetID (str): the metabolite ID of ATP
        modelName (str): name of the model
    Returns:
        printed results of the biomass flux, ATP flux and the BM/ATP ratio
    """
    # get the flux solution
    solutions = model.optimize()
    fluxArray = solutions.fluxes

    # get the stoichiometric matrix
    StoiMatrixDF = cobra.util.create_stoichiometric_matrix(model, array_type='DataFrame')  # [:,0:posExchangeRxn]

    # find the ATP producing flux
    RowATP = StoiMatrixDF.loc[ATPmetID]
    fluxATP = RowATP * fluxArray

    fluxList = []
    for flux in fluxATP:
        if flux > 0:
            fluxList.append(flux)
    totalATPFlux = sum(fluxList)

    # find the biomass flux
    biomassExReaction = model.reactions.get_by_id(biomassRxnID)
    biomassFlux = biomassExReaction.flux

    # the ratio
    BM_ATP_ratio =  biomassFlux/ totalATPFlux

    if printResults:
        print('model:', modelName)
        print('the total flux (mmol/g/h) of ATP is: ', totalATPFlux)
        print('the total flux (mmol/g/h) of BM is: ', biomassFlux)
        print('the total ratio BM/ATP (mmolBM/mmolATP) is: ', BM_ATP_ratio)
        print('')

    return  BM_ATP_ratio


def find_yield(model, substrateExchangeRxnID, productExchangeRxnID, printResults = False, Biomass = False):
    """
    Finds the yields of a product derived from a given substrate

    Params:
        model (model): GEM model
        substrateID (str): id of the exchange reaction of the substrate
        productID (str): id of the exchange reaction of the product

    returs
        yield (float): returns the yield of the specified compounds in g/g

    """
    #model.optimize()
    cobra.flux_analysis.pfba(model)
    #solutionFVA = cobra.flux_analysis.flux_variability_analysis(model,processes= 1)
    #print(solutionFVA)
    # get the reactions:
    RxnSubstrate = model.reactions.get_by_id(substrateExchangeRxnID)
    RxnProduct = model.reactions.get_by_id(productExchangeRxnID)

    # get the metabolite objects
    metSubstrate =  RxnSubstrate.reactants[0]
    metProduct = RxnProduct.reactants[0]

    # get the molecular weight
    mwSubstrate = metSubstrate.formula_weight
    mwProduct = metProduct.formula_weight

    # get the fluxes
    fluxSubstrate = RxnSubstrate.flux
    fluxProduct = RxnProduct.flux
    print(metProduct.name)
    if Biomass: #'biomass' in metProduct.name:
        ratio = - (fluxProduct) / (fluxSubstrate * mwSubstrate * 0.001) # bio mass in g/g/h substrate in mmol/g/h
        MetabioliteName = 'Biomass'
    else:
        MetabioliteName = metProduct.name
        try:
            ratio = - (fluxProduct * mwProduct) / (fluxSubstrate * mwSubstrate)
        except: # if the division is by zero (a nonsense result) return ratio = 0
            ratio = 0

    if printResults:
        print('the yield (g/g) of {} is: {} \n'.format(MetabioliteName, ratio))

    return ratio

def is_protein_met(metabolite):
    """ checks if the given metaboltie is a protein base on the chemical formula"""
    elements = ['C', 'N', 'H', 'O']
    elementDict = {}

    for e in elements:
        nElement = count_atom_in_formula(metabolite= metabolite,atom=e)
        elementDict.update({e:nElement})

    if elementDict['C'] >= 2 and elementDict['N'] >= 1 and elementDict['H'] >= 4 and elementDict['O'] >= 2 :
        protein = True
    else:
        protein = False

    return protein