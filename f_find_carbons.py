# find the carbons of a certain metabolite by looking at the producing reaction(s)
import cobra.io
import os
import re
import numpy as np
import warnings
import pandas as pd

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

def countCarbonInFormula(metFormula):
    splitFormula = re.split('(\d+)', metFormula)
    nrOfCarbons = 0  # just in case something wierd happens
    if 'C' not in metFormula:  # if there is no carbon in the formula
        nrOfCarbons = 0
    else:
        for j, element in enumerate(splitFormula):
            if 'C' in element and len(element) == 1:
                nrOfCarbons = int(splitFormula[j + 1])
            elif 'C' in element and len(element) > 1:
                posCarbon = element.index('C')
                if element[posCarbon + 1].isupper():  # for case like CH4 there is no 1 next to the C
                    nrOfCarbons = 1  # only one carbon
                else:
                    continue  # for cases like Co (cobalt) just skip
            else:
                continue
    return nrOfCarbons

def countCarbonInList(reaction, reactionList):
    carbonNrAll = []
    for met in reactionList:
        stoiCoef_i = abs(reaction.metabolites[met])
        metFormula = met.formula
        splitFormula = re.split('(\d+)', metFormula)
        nrOfCarbons = 0  # just in case something wierd happens
        if 'C' not in metFormula:  # if there is no carbon in the formula
            nrOfCarbons = 0
        else:
            for j, element in enumerate(splitFormula):
                if 'C' in element and len(element) == 1:
                    nrOfCarbons = int(splitFormula[j + 1])
                elif 'C' in element and len(element) > 1:
                    posCarbon = element.index('C')
                    if element[posCarbon + 1].isupper():  # for case like CH4 there is no 1 next to the C
                        nrOfCarbons = 1  # only one carbon
                    else:
                        continue  # for cases like Co (cobalt) just skip
                else:
                    continue

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
            nCarbon = countCarbonInFormula(formula)
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

def carbon_balance_in_out(modelLocation, metIDsMissingCarbon=None, tol = 1e-4):
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

def carbonBalance(model, reactionDF, missingCarbonDict, tol=1e-4):
    metNamesAll = []
    carbonNrAll = []
    gramsCAll = []
    for i, metID in enumerate(reactionDF.metabolite):
        if abs(reactionDF.flux[i]) > tol:
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
                nrOfCarbons  = countCarbonInFormula(metFormula)

            carbonNrAll.append(nrOfCarbons)
            gramsC = nrOfCarbons * 12 * reactionDF.flux[i]
            gramsCAll.append(gramsC)


    dictReactions = {'Metabolite': metNamesAll,
                     '# of C': carbonNrAll,
                     'flux (gram-C/g-DW/h)': gramsCAll}

    dataFrameReactions = pd.DataFrame(dictReactions)
    return dataFrameReactions

if __name__ == '__main__':
    # check if I can find the # carbons in the cell wall
    loc = os.getcwd()
    loc_acidi = loc + r'\SBML models\PAC_4875_model.xml'
    model = cobra.io.read_sbml_model(loc_acidi)
    cellwallID = 'S_cpdnew27_c0'
    c = findCarbonsMissingMetabolite(model= model, metID= cellwallID)
    print(c)