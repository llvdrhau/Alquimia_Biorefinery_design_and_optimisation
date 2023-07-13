from f_screen_SBML import *

modelName = 'p-thermo.xml'
loc_sher = get_location(modelName)
model = cobra.io.read_sbml_model(loc_sher)
model.name = modelName

# # ----------------------------------------------------- look at anaerobic conditions' hello
# make anaerobic, o2 exchange reaction to zero
o2_exchRxn_ID = 'EX_o2_e'
o2_exchRxn = model.reactions.get_by_id(o2_exchRxn_ID).bounds = 0, 0
carbon_balance_in_out(model)
model.optimize()
#print_SBML_info_2_excel(modelName=model, saveName= 'analysis_p_thermo.xlsx', print2Excel= True)

print('Extracellular \n')
ExGlucoseId = 'glc__D_e'
ExGlucoseMet = model.metabolites.get_by_id(ExGlucoseId)
print_all_rxn_of_metabolite(ExGlucoseMet, case='id', printFlux=True)

print('Cytosol \n')
CytoGlucoseId = 'glc__D_c'
CytoGlucoseMet = model.metabolites.get_by_id(CytoGlucoseId)
print_all_rxn_of_metabolite(CytoGlucoseMet, case='id', printFlux=True)

print('ETHANOL in Cytosol \n')
CytoEthId = 'etoh_c'
CytoEthMet = model.metabolites.get_by_id(CytoEthId)
print_all_rxn_of_metabolite(CytoEthMet, case='id', printFlux=True)

print('ETHANOL in Extracellular \n')
CytoEthId = 'etoh_e'
CytoEthMet = model.metabolites.get_by_id(CytoEthId)
print_all_rxn_of_metabolite(CytoEthMet, case='id', printFlux=True)

# # ----------------------------------------------------- get the atp biomass ratio
# # metabolite and reaction id's
biomassExReactionID = 'EX_Biomass_e'
ATPMetaboliteID = 'atp_c'
# print the ratios
ATP_Biomass_Ratio(model=model, biomassRxnID=biomassExReactionID,
                  ATPmetID=ATPMetaboliteID, modelName=modelName)





