
import cobra
import pandas as pd
from f_usefull_functions import get_location


modelName = "P_sherm_model.xml"
loc_sher = get_location(modelName)
model = cobra.io.read_sbml_model(loc_sher)
metabolites = model.metabolites
# blockedRxn = cobra.flux_analysis.find_blocked_reactions(model= model)
# print(len(blockedRxn))

metID = []
metName = []
for met in metabolites:
    formula = met.formula
    if not formula:
        metID.append(met.id)
        metName.append(met.name)
DictMetabolites = {'ID': metID, 'Name' : metName}
DFmetabolites = pd.DataFrame(DictMetabolites)
print(DFmetabolites)
DFmetabolites.to_excel('missing_metabolites2.xlsx')

