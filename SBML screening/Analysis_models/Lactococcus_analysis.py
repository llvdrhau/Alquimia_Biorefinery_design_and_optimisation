from f_screen_SBML import *

modelName = 'iNF517.xml'
# Lactococcus lactis subsp. cremoris MG1363
# from BIGGMODELS: http://bigg.ucsd.edu/models/iNF517

loc_sher = get_location(modelName)
model = cobra.io.read_sbml_model(loc_sher)
model.name = modelName
model.optimize()
# print_SBML_info_2_excel(modelName=model, saveName= 'analysis_iNF517.xlsx', print2Excel= True)

# ----------------------------------------------------- look at anaerobic condition

print('Extracellular GLUCOSE \n')
ExGlucoseId = 'glc__D_e'
ExGlucoseMet = model.metabolites.get_by_id(ExGlucoseId)
print_all_rxn_of_metabolite(ExGlucoseMet, case='id', printFlux=True)

print('Cytosol GLUCOSE \n')
CytoGlucoseId = 'glc__D_c'
CytoGlucoseMet = model.metabolites.get_by_id(CytoGlucoseId)
print_all_rxn_of_metabolite(CytoGlucoseMet, case='id', printFlux=True)

print('Acetate in Cytosol \n')
acID = 'ac_c'
CytoEthMet = model.metabolites.get_by_id(acID)
print_all_rxn_of_metabolite(CytoEthMet, case='id', printFlux=True)

print('Acetate in Extracellular \n')
CytoEthId = 'ac_e'
CytoEthMet = model.metabolites.get_by_id(CytoEthId)
print_all_rxn_of_metabolite(CytoEthMet, case='id', printFlux=True)

# # ----------------------------------------------------- get the atp biomass ratio
# metabolite and reaction id's
biomassExReactionID = 'BIOMASS_LLA'
ATPMetaboliteID = 'atp_c'

bmFormula , MWbm = estimate_biomass_formula_from_products(model = model, reactionID=biomassExReactionID)
print(bmFormula)
# print the ratios
ratio = ATP_Biomass_Ratio(model=model, biomassRxnID=biomassExReactionID,
                  ATPmetID=ATPMetaboliteID, modelName=modelName)

ratio2 = ratio * MWbm
print(ratio2)


