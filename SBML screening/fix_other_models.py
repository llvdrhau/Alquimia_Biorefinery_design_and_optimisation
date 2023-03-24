from old_scripts.functions import *
from cobra.io import write_sbml_model
# this is just a quick fix to roughly make the massbalances correct (not in detail like fix_model_sherm.py)

modelName = ["P_avidum_model.xml",
             "P_propionicum_model.xml",
             "PAC_4875_model.xml"]

saveName = ['P_avidum_V2.xml',
            "P_propionicum_V2.xml",
            "PAC_4875_V2.xml"]

for i,modelName in enumerate(modelName):
    loc_model = get_location(modelName)
    model = cobra.io.read_sbml_model(loc_model)
    model.name = saveName[i]

    # bioMass metabolite does not have a name
    metbiomassId = 'S_biomass_ext'
    metbiomass = model.metabolites.get_by_id(metbiomassId)
    metbiomass.name = 'Biomass'
    metbiomass.formula = 'C36H43O24' # not rounded off 'C36.74H43.76O24.42'
    print_SBML_info_2_excel(modelName=model, print2Excel=False, saveName='test666.xlsx')

    # ----------------- split into compartments by re-moving and adding new boundry reactions
    # Get the original (faulty) exchange reactions (even though there is only one compartment for the time being,
    # It does a good job at guessing what they are in this case)
    model._compartments = {'c': 'Cytoplasma', 'e': 'Extracellular'}
    listExchRxn = model.exchanges
    for original_reaction in listExchRxn:
        # Create a new exchange reaction with the same metabolites as the original reaction
        metabolite = original_reaction.reactants[0]  # there is only going to be one metabolite in the reaction reactants
        metabolite.compartment = 'e'

    # save to correct file location
    write_sbml_model(model, r"C:\Users\lucas\PycharmProjects\Alquimia\SBML models\{}".format(saveName[i]))
