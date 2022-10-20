import cobra.io
import os
import pandas as pd
from sbml_test import productFluxes, getNamesMediaMetabolites, mediaFluxes
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
            gramsC = nrOfCarbons * 12 * reactionDF.flux[i]
            gramsCAll.append(gramsC)

    dictSecretion = {'Metabolite': metNamesAll ,
            '# of C': carbonNrAll,
           'flux (gram-C/g-DW/h)': gramsCAll}
    dataFrameSecretion = pd.DataFrame(dictSecretion)

    return  dataFrameSecretion


#read model
loc = os.getcwd()
loc_acidi = loc + r'\SBML models\PAC_4875_model.xml'
model = cobra.io.read_sbml_model(loc_acidi)
# df = model.summary()
# print(df)
uptakeDF = carbonBalance(model, case = 'uptake')
secretionDF = carbonBalance(model)
print(secretionDF)

CgramsOut = sum(secretionDF['flux (gram-C/g-DW/h)'])
CgramsIn = sum(uptakeDF['flux (gram-C/g-DW/h)'])    # should be 10*6*12 = 780  # 10 mols of glu, 6 C atoms, 12 gC/mol
print(CgramsOut/CgramsIn)
print(uptakeDF)

# medium
# mediumNames = getNamesMediaMetabolites(model=model) # names in medium
# print(mediumNames)
# medium = model.medium
# print(medium)

# get secretion and uptake fluxes
# secretion = productFluxes(model,tol= 0.001)
# uptake = mediaFluxes()
#
# print(len(secretion))
# print(secretion)
#

