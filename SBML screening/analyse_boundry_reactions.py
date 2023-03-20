from f_screen_SBML import *

# read the model
modelName = "P_acnes_model.xml"
save = False

loc_model = get_location(modelName)
model = cobra.io.read_sbml_model(loc_model)
model.optimize()

print('----------- Glucose in Extracellular -----------  \n')
ExGlucoseId = 'S_cpd00027_ext'
ExGlucoseMet = model.metabolites.get_by_id(ExGlucoseId)
print_all_rxn_of_metabolite(ExGlucoseMet, case='names', printFlux=True)

print('----------- Glucose in Cytosol -----------  \n')
cytoGlucose = model.metabolites.get_by_id('S_cpd00027_c0')
print_all_rxn_of_metabolite(cytoGlucose, case='names', printFlux=True)


# look at currrent transport mechanism for acetate and propionate
print('----------- ACETATE in Cytosol -----------  \n')
acID = 'S_cpd00029_c0'
intraMet = model.metabolites.get_by_id(acID)
print_all_rxn_of_metabolite(intraMet, case='names', printFlux=True)

print('-----------  ACETATE in Extracellular -----------  \n')
CytoEthId = 'S_cpd00029_ext'
extraMet = model.metabolites.get_by_id(CytoEthId)
print_all_rxn_of_metabolite(extraMet, case='names', printFlux=True)

