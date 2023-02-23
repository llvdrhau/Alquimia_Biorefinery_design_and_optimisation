from old_scripts.functions import *

modelList = ["P_sherm_model.xml", 'P_propionicum_model.xml', 'PAC_4875_model.xml', 'P_acnes_model.xml', 'P_avidum_model.xml']

for modelName in modelList:
    loc_sher = get_location(modelName)
    model = cobra.io.read_sbml_model(loc_sher)
    model.name = modelName

    # change the biomass formula
    metbiomassId = 'S_biomass_ext'
    metbiomass = model.metabolites.get_by_id(metbiomassId)
    metbiomass.name = 'Biomass'
    #metbiomass.formula = 'C15.12HNO'

    print_SBML_info_2_excel(modelName=model, idMissingCarbon=None, print2Excel= False)