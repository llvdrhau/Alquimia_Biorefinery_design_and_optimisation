import pandas as pd
import cobra
from cobra import Metabolite
import numpy as np
from f_usefull_functions import get_location


modelName = "P_sherm_model.xml"
loc_sher = get_location(modelName)
model = cobra.io.read_sbml_model(loc_sher)

# give biomass metabolite a name
metbiomassId = 'S_biomass_ext'
metbiomass = model.metabolites.get_by_id(metbiomassId)
metbiomass.name = 'Biomass'


stoichiometric_matrix = cobra.util.create_stoichiometric_matrix(model, array_type='DataFrame')  # [:,0:posExchangeRxn]
print(stoichiometric_matrix)

FBA = model.optimize()
fluxArray = FBA.fluxes #[0:posExchangeRxn] #.to_numpy()
#fluxArray.drop(ExRxn, inplace=True)

metabolicFlux = np.dot(stoichiometric_matrix, fluxArray)
metabolicFluxDF = pd.DataFrame(metabolicFlux, index=stoichiometric_matrix.index)
print(metabolicFluxDF)

metabolites = model.metabolites
metID = []
metName = []
strRxns = []
for met in metabolites:
    formula = met.formula
    if not formula:
        metID.append(met.id)
        metName.append(met.name)
        frozenRct = met.reactions
        # get the reaction that produces the metabolite
        allRxn = list(frozenRct)

        # todo fix the append of filler names
        extraNames = [''] * len(allRxn)
        metName.append(extraNames)

        for rxn in allRxn:
            reactants = rxn.reactants
            reactantStr = ''
            for met in reactants:
                reactName = met.name
                reactantStr += reactName + ' + '

            productStr = ""
            products = rxn.products
            for met in products:
                prodName = met.name
                productStr += prodName + ' + '

            reactionStr = "{} = {}".format(reactantStr, productStr)
            strRxns.append(reactionStr)

DictMetabolites = {'ID': metID, 'Name' : metName}
DFmetabolites = pd.DataFrame(DictMetabolites)
print(DFmetabolites)

#with pd.ExcelWriter('P_sherm_analysis.xlsx') as writer:
#    stoichiometric_matrix.to_excel(writer, sheet_name='stoichiometric_matrix')

