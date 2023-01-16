'''
This function import information of the model to an Excel file to check the missing carbon in the metabolic reactions
the percent of missing carbon can also  be attributed to the reactions to check on their importants of the reactions to
the outcome of the FBA
'''
import pandas as pd
import cobra
import numpy as np
from f_usefull_functions import get_location
from f_find_carbons import carbon_balance_in_out, countCarbonInFormula
#from cobra import Metabolite

def string_reactions(reaction, case = 'names'):
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
def print_SBML_info_2_excel(modelName):
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
    #print(StoiMatrixDF)

    FBA = model.optimize()
    fluxArray = FBA.fluxes #[0:posExchangeRxn] #.to_numpy()

    #fluxArray.drop(ExRxn, inplace=True)
    #rxnFluxes = np.dot(StoiMatrixDF, fluxArray)
    #metabolicFluxDF = pd.DataFrame(metabolicFlux, index=StoiMatrixDF.index)
    #print(metabolicFluxDF)

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
        carbon = countCarbonInFormula(metFormula= formula)
        carbonCount.append(carbon)
        if not formula:
            metID.append(met.id)
            metName.append(met.name)
            metNameWithSpaces.append(met.name)
            frozenRct = met.reactions
            # get the reaction that produces the metabolite
            allRxn = list(frozenRct)

            extraNames = [''] * (len(allRxn) -1)
            metNameWithSpaces = metNameWithSpaces + extraNames
            #metName = metName + extraNames

            # helpNames = metName + extraNames
            # metNameWithSpaces = metNameWithSpaces + helpNames
            # metName.append(extraNames)

            for rxn in allRxn:
                reactionStr = string_reactions(reaction= rxn)
                strRxns.append(reactionStr)

                fluxRxn.append(rxn.flux)
                stoiFactor = rxn.metabolites
                stoiMetMissingFormula.append(stoiFactor[met])


    DictCarbons ={'ID metabolite':list(StoiMatrixDF.index), '# Carbons' : carbonCount}
    CarbonsDF = pd.DataFrame(DictCarbons)

    DictMetabolites = {'ID': metID, 'Name' : metName}
    DFmetabolites = pd.DataFrame(DictMetabolites)
    #print(DFmetabolites)

    DictMetRnx = {'Name' : metNameWithSpaces, 'Reactions':strRxns, 'Flux' : fluxRxn, 'Stoichiometry' : stoiMetMissingFormula }
    DFmetRnx = pd.DataFrame(DictMetRnx)
    #print(DFmetRnx)

    # find the ingoing and outgoing fluxes
    objectiveMetID = 'S_biomass_ext'
    inputDF, outputDF  = carbon_balance_in_out(modelLocation=model, metIDsMissingCarbon=objectiveMetID, tol= 0.0001)

    # calculate the percentage of carbon at goes missing in each reaction (exculed the exchage reactions?)
    #exclude the transfer (exchange reactions) reactions
    exchRxn = model.exchanges
    # transform id's to a list
    keysListExRxn = [rx.id for rx in exchRxn]

    # to drop or not?
    fluxArray.drop(keysListExRxn, inplace=True)
    fluxArrayNp = np.array(fluxArray).reshape(-1,1)
    st = StoiMatrixDF.drop(keysListExRxn, axis= 'columns')# ,inplace=False)
    #st = StoiMatrixDF

    carbonArray = np.array(carbonCount)
    carbonArray = carbonArray.reshape(-1,1)
    StoiMatrixTranspose = np.transpose(st)

    carbonBalance = np.dot(StoiMatrixTranspose,carbonArray) # in mols of Carbon
    # multiply by 12 gC/mol
    carbonMassBalance = carbonBalance * 12
    carbonMissing = carbonMassBalance * fluxArrayNp
    rxnIDreduced = list(fluxArray.index)

    DictMissingCarbonALL = {'Reaction Id': rxnIDreduced, 'carbon (gCarbon)':carbonMissing.reshape(len(carbonMissing),) }
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
    DictMissingCarbon = {'Reaction Id': unbalancedId,'carbon (gCarbon)': unbalanced, 'reaction': rxnUnblanaced}
    DFUnbalancedCarbonReactions = pd.DataFrame(DictMissingCarbon)

    #UnbalancedCarbonReactions
    #carbonStoiMatrix
    #carbonStoiMatrix = np.dot(np.array(carbonCount),np.transpose(StoiMatrixDF))
    # numpy.array(carbonCount)

    with pd.ExcelWriter(saveName) as writer:
        StoiMatrixDF.to_excel(writer, sheet_name='StoichiometricMatrix')
        fluxArray.to_excel(writer, sheet_name= 'ReactionFluxes')
        CarbonsDF.to_excel(writer, sheet_name= 'CarbonsPerMetbolite')
        DFmetRnx.to_excel(writer, sheet_name='MissingFormulaReactions')
        DFUnbalancedCarbonReactions.to_excel(writer, sheet_name= 'Carbon Balance')
        inputDF.to_excel(writer, sheet_name='Carbon Input')
        outputDF.to_excel(writer, sheet_name='Carbon output')

if __name__ == '__main__':
    modelName = "P_sherm_model.xml"
    loc_sher = get_location(modelName)
    model = cobra.io.read_sbml_model(loc_sher)
    model.name = modelName

    # give biomass metabolite a name
    metbiomassId = 'S_biomass_ext'
    metbiomass = model.metabolites.get_by_id(metbiomassId)
    metbiomass.name = 'Biomass'
    metbiomass.formula = 'C15.12HNO'

    print_SBML_info_2_excel(modelName= model)