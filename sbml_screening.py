import cobra.io
import os
import pandas as pd
from sbml_test import productFluxes, getNamesMediaMetabolites
import re
import numpy as np

def carbonBalance(model, case = 'secretion' ,tol = 1e-4):
    df = model.summary()

    if case == 'secretion':
        reactionDF = df.secretion_flux
    else:
        reactionDF = df.uptake_flux

    metNamesAll = []
    carbonNrAll = []
    gramsCAll = []
    for i,metID in enumerate(reactionDF.metabolite):
        if abs(reactionDF.flux[i]) > tol:
            met = model.metabolites.get_by_id(metID)
            metName = met.name
            metNamesAll.append(metName)
            metFormula = met.formula
            splitFormula = re.split('(\d+)', metFormula)
            nrOfCarbons = 0 # just in case something wierd happens
            if 'C' not in metFormula: # if there is no carbon in the formula
                nrOfCarbons = 0
            else:
                for j, element in enumerate(splitFormula):
                    if 'C' in element and len(element) == 1:
                        nrOfCarbons = int(splitFormula[j+1])
                    elif 'C' in element and len(element) > 1:
                        posCarbon = element.index('C')
                        if element[posCarbon+1].isupper(): # for case like CH4 there is no 1 next to the C
                            nrOfCarbons = 1 # only one carbon
                        else: continue  # for cases like Co (cobalt) just skip
                    else: continue

            carbonNrAll.append(nrOfCarbons)
            # TODO should also be multiplied by the stoichiometric coef.
            gramsC = nrOfCarbons * 12 * reactionDF.flux[i]
            gramsCAll.append(gramsC)

    dictSecretion = {'Metabolite': metNamesAll ,
            '# of C': carbonNrAll,
           'flux (gram-C/g-DW/h)': gramsCAll}
    dataFrameSecretion = pd.DataFrame(dictSecretion)

    return  dataFrameSecretion

def carbonBalanceReactions(model, metaboliteReactDict ,tol = 1e-4):
    metNamesAll = []
    carbonNrAll = []
    gramsCAll = []
    for i,metID in enumerate(metaboliteReactDict):
        rctID = metaboliteReactDict[metID]
        flux = model.reactions.get_by_id(rctID).flux
        met = model.metabolites.get_by_id(metID)
        stoiCoef = model.reactions.get_by_id(rctID).metabolites[met]
        #print(stoiCoef.keys())
        if abs(flux) > tol:
            met = model.metabolites.get_by_id(metID)
            metName = met.name
            metNamesAll.append(metName)
            metFormula = met.formula
            splitFormula = re.split('(\d+)', metFormula)
            nrOfCarbons = 0 # just in case something wierd happens
            if 'C' not in metFormula: # if there is no carbon in the formula
                nrOfCarbons = 0
            else:
                for j, element in enumerate(splitFormula):
                    if 'C' in element and len(element) == 1:
                        nrOfCarbons = int(splitFormula[j+1])
                    elif 'C' in element and len(element) > 1:
                        posCarbon = element.index('C')
                        if element[posCarbon+1].isupper(): # for case like CH4 there is no 1 next to the C
                            nrOfCarbons = 1 # only one carbon
                        else: continue  # for cases like Co (cobalt) just skip
                    else: continue

            carbonNrAll.append(nrOfCarbons)
            gramsC = nrOfCarbons * 12 * flux * stoiCoef
            gramsCAll.append(gramsC)

    dictSecretion = {'Metabolite': metNamesAll ,
            '# of C': carbonNrAll,
           'flux (gram-C/g-DW/h)': gramsCAll}
    dataFrameReturn = pd.DataFrame(dictSecretion)

    return  dataFrameReturn

def countCarbonInReactants(model,rctID, productID):
    #TODO if you want to figure out how much carbon is in a specific product, you also need to look at where the carbon is in the other products
    rct = model.reactions.get_by_id(rctID)
    flux = rct.flux
    listMets = rct.reactants
    productMet = model.metabolites.get_by_id(productID)
    productName = productMet.name
    stoiCoef = model.reactions.get_by_id(rctID).metabolites[productMet]
    carbonNrAll = []
    for met in listMets:
        stoiCoef_i = abs(model.reactions.get_by_id(rctID).metabolites[met])
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
    fluxCgrams = sum(carbonNrAll) * flux * stoiCoef
    totaalCarbonAtoms = sum(carbonNrAll)
    #gramsC = nrOfCarbons * 12 * flux * stoiCoef
    #gramsCAll.append(gramsC)
    return productName, totaalCarbonAtoms, fluxCgrams
#read model
loc = os.getcwd()
loc_acidi = loc + r'\SBML models\PAC_4875_model.xml'
model = cobra.io.read_sbml_model(loc_acidi)

###uncomment to activate again
# # df = model.summary()
# # print(df)
#
uptakeDF = carbonBalance(model, case = 'uptake')
secretionDF = carbonBalance(model)
print(secretionDF)
print(uptakeDF)
#
CgramsOut = sum(secretionDF['flux (gram-C/g-DW/h)'])
CgramsIn = sum(uptakeDF['flux (gram-C/g-DW/h)'])    # should be 10*6*12 = 720  # 10 mols of glu, 6 C atoms, 12 gC/mol
print("{} of the carbon is accounted for, let's see how much goes into bio mass".format(-CgramsOut/CgramsIn))

################################################
# 46 % of the carbon is accounted for, let's see how much goes into bio mass
# #(see colab script, checkObjectiveReaction as well)
model.metabolites.get_by_id('S_biomass_ext').summary()

mets = ['S_cpd11461_c0', # 'DNA', '
'S_cpd11463_c0' , # Protein'
'S_cpd11613_c0' , # , 'RNA'
'S_cpdnew18_c0' , # Propionibacterium lipids',
'S_cpdnew26_c0' , #'Small molecule pool',
'S_cpdnew27_c0'] # 'Propionibacterium cell wall'

searchNames = False
if searchNames:
    names = []
    for r in mets:
      names.append(model.metabolites.get_by_id(r).name)
    print(names)

################################################################################################################
# for metabolites like DNA in biomass cell wall count all carbons used to make the cell wall
rctWithMetId = {'S_cpd11461_c0': 'biomass_c0', # 'DNA', '
'S_cpd11463_c0': 'biomass_c0' , # Protein'
'S_cpd11613_c0':'biomass_c0'}  # , 'RNA'

df = carbonBalanceReactions(model= model, metaboliteReactDict= rctWithMetId)
#print(df)

################################################################################################################
# for metabolites like protein cell wall count all carbons used to make the cell wall, pool of molecules and lipids
smallPoolMet = countCarbonInReactants(model, rctID= 'rxnnew72_c0', productID= 'S_cpdnew26_c0')
print(smallPoolMet)

### cell wall
# sm = model.metabolites.get_by_id('S_cpdnew27_c0').summary()
#print(sm)
# polyschacarides
polyS = countCarbonInReactants(model, rctID= 'rxnnew64_c0',productID = 'S_cpdnew20_c0')
#print(polyS)
# Lipomannan
lipoM = countCarbonInReactants(model, rctID= 'rxnnew65_c0',productID = 'S_cpdnew21_c0')
#print(lipoM)
# peptidoglycan
pep = countCarbonInReactants(model, rctID= 'rxnnew66_c0',productID = 'S_cpdnew22_c0')
#print(pep)

carbonFluxInBiomass = polyS[2] + lipoM[2] + pep[2] + smallPoolMet[2] + sum(df['flux (gram-C/g-DW/h)'])


print('{} grams C is in Biomass representing {}% of the total carbon'.format(carbonFluxInBiomass,100*carbonFluxInBiomass/CgramsIn))
#smallPoolMet = countCarbonInReactants(model, rctID= 'rxnnew72_c0')
# #newRow ={'Metabolite':smallPoolMet[0],'# of C':smallPoolMet[1], 'flux (gram-C/g-DW/h)':smallPoolMet[2] }
#df2add = pd.DataFrame(newRow)
#dfnew = pd.concat([df,df2add],ignore_index=True)
#df.append(newRow, ignore_index=True)