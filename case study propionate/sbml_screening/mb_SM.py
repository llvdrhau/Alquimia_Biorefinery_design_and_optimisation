
import cobra
import numpy as np
from f_find_carbons import countCarbonInFormula
import pandas as pd


def rxnExIDs_2_metExIDs(reactionList):
    listMetabolites = []
    reactionDict = {}
    for i in reactionList:
        #mediumReaction = model.reactions.get_by_id(i)
        reactant = i.reactants
        nameMetabolite = reactant[0].id
        listMetabolites.append(nameMetabolite)
        reactionDict.update({nameMetabolite:i})
    return listMetabolites, reactionDict
def find_metabolite_index(model, metID):
    mets = model.metabolites
    for i, met in enumerate(mets):
        if met.id == metID:
            break
    return i

def find_reaction_index(model, rxnID = None):
    rxns = model.metabolites
    if rxnID is None:
        for i, rxn in enumerate(rxns):
            if 'Ex' in rxn.id: # omit the exchange reactions (they should be all sequential)
                break
        return i
    else:
        for i, rxn in enumerate(rxns):
            if rxn.id == rxnID:
                break
        return i

loc_sher = r"C:\Users\lucas\PycharmProjects\Alquimia\SBML models\P_sherm_model.xml"
model = cobra.io.read_sbml_model(loc_sher)

exchangeRnxDictList = model.exchanges
ExRxn = list(exchangeRnxDictList.__dict__['_dict'])
ExMet = rxnExIDs_2_metExIDs(exchangeRnxDictList)[0]

stoichiometric_matrix = cobra.util.create_stoichiometric_matrix(model, array_type='DataFrame') # [:,0:posExchangeRxn]
stoichiometric_matrix.drop(ExRxn, axis=1, inplace=True) # get rid of the exchange reactions
stoichiometric_matrix.drop(ExMet, axis=0, inplace=True) # get rid of the exchange reactions

FBA = model.optimize()
fluxArray = FBA.fluxes #[0:posExchangeRxn] #.to_numpy()
fluxArray.drop(ExRxn, inplace=True)

metabolicFlux = np.dot(stoichiometric_matrix, fluxArray)
metabolicFluxDF = pd.DataFrame(metabolicFlux, index=stoichiometric_matrix.index)

metabolite_carbon_flux_array = []
positive = []
negative = []

for i,met in enumerate(model.metabolites):
    if met.id in ExMet:
        pass
    else:
        metID = met.id
        nCarbon = countCarbonInFormula(met.formula)
        carbonFlux =  metabolicFluxDF.loc[metID][0] * nCarbon * 12
        metabolite_carbon_flux_array.append(carbonFlux)  # in grams-carbon / h / gDW
        if carbonFlux > 0:
            positive.append(carbonFlux)
        elif carbonFlux < 0:
            negative.append(carbonFlux)

#glu_row = stoichiometric_matrix.loc['S_cpd00027_ext'] # id D-glucose is 'S_cpd00027_ext'
#flux_carbon_glu = np.dot(glu_row, fluxArray)
#print('the carbon flux of glucose is {}'.format(flux_carbon_glu))

print('the sum of all carbon fluxes in grams carbon / gDW/h is:')
print(sum(metabolite_carbon_flux_array))
metabolite_carbon_flux_DF = pd.DataFrame(metabolite_carbon_flux_array, index=stoichiometric_matrix.index )
print('')
print('the ratio in out is: {}'.format(abs(sum(positive)/sum(negative))))
print('')
print('the sum of carbon of that is created internally (grams carbon/gDW/h)  is: ')
print(sum(positive))