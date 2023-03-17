from cobra import Metabolite

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




# ----------------- split into compartments by re-moving and adding new boundry reactions
# TODO make a loop over this code to change each boundry reaction
# see chatGTP as well
# Get the original (faulty) exchange reactions (even though there is only one compartment for the time being,
# It does a good job at guessing what they are in this case)
listExchRxn = model.exchanges

# Create a new exchange reaction with the same metabolites as the original reaction
new_reaction = cobra.Reaction(f'{original_reaction.id}_exchange')
new_reaction.name = f'{original_reaction.name} exchange'
new_reaction.add_metabolites(original_reaction.metabolites)
new_reaction.lower_bound = -1000
new_reaction.upper_bound = 1000

# Add the new exchange reaction to the model and remove the original reaction
model.add_reactions([new_reaction])
model.reactions.remove(original_reaction)


model.add_metabolites(Metabolite('nan', name='nan', compartment='e'))
# this is super hacky and probably not the best way to get a new compartment in the model

model.add_boundary(model.metabolites.get_by_id("nan"), type="exchange")
reactionNan = model.reactions.get_by_id('EX_nan')
#model.add_compartment('e', 'Extra-Cellular Space')
model._compartments = {'c': 'Cytoplasma', 'e': 'Extracellular'}


# define compartments to a reaction
glucoseExchID = 'Ex_S_cpd00027_ext'
reaction = model.reactions.get_by_id(glucoseExchID)
reaction.compartments = 'e'



# model.add_metabolites(Metabolite('co2_e', name='CO2', compartment='e'))
# model.add_boundary(model.metabolites.get_by_id("co2_e"), type="exchange")


# 'Ex_S_cpd00029_ext',    # propionate
# model.add_boundary(model.metabolites.get_by_id("co2_e"), type="exchange")
# model.add_metabolites(Metabolite('co2_e', name='CO2', compartment='e'))
# model.add_compartment('e', 'Extra-Cellular Space')
# model.add_compartment()
