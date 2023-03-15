from f_screen_SBML import *

modelName = 'p-thermo.xml'
loc_sher = get_location(modelName)
model = cobra.io.read_sbml_model(loc_sher)
model.name = modelName
model.optimize()
# carbon_balance_in_out(model)


print('Extracellular \n')
ExGlucoseId = 'glc__D_e'
ExGlucoseMet = model.metabolites.get_by_id(ExGlucoseId)
print_all_rxn_of_metabolite(ExGlucoseMet, case='names', printFlux=True)

print('Cytosol \n')
CytoGlucoseId = 'glc__D_c'
CytoGlucoseMet = model.metabolites.get_by_id(CytoGlucoseId)
print_all_rxn_of_metabolite(CytoGlucoseMet, case='names', printFlux=True)
print('sw')

# # ----------------------------------------------------- get the atp biomass ratio
# # metabolite and reaction id's
# biomassExReactionID = 'EX_Biomass_e'
# ATPMetaboliteID = 'atp_c'
# # print the ratios
# ATP_Biomass_Ratio(model=model, biomassRxnID=biomassExReactionID,
#                   ATPmetID=ATPMetaboliteID, modelName=modelName)
#
# # ----------------------------------------------------- look at anaerobic conditions' hello
# # make anaerobic, o2 exchange reaction to zero
# o2_exchRxn_ID = 'EX_o2_e'
# o2_exchRxn = model.reactions.get_by_id(o2_exchRxn_ID).bounds = 0, 0
# carbon_balance_in_out(model)

# # print_SBML_info_2_excel(modelName=model, saveName= 'analysis_p_thermo.xlsx', print2Excel= False)

