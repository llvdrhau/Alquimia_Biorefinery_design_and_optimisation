from cobra import Metabolite

from f_screen_SBML import *

# read the model
modelName = "P_acnes_model.xml"
save = False

loc_model = get_location(modelName)
model = cobra.io.read_sbml_model(loc_model)
model.optimize()
# ----------------------------------------------------------- split into compartments
# define compartments a compartment

cytoGlucose = model.metabolites.get_by_id('S_cpd00027_c0')
print_all_rxn_of_metabolite(cytoGlucose, case='names', printFlux=True)

print(model.exchanges)
model.add_metabolites(Metabolite('co2_e', name='CO2', compartment='e'))
model.add_boundary(model.metabolites.get_by_id("co2_e"), type="exchange")
model._compartments = {'c': 'Cytoplasma', 'e': 'Extracellular'}
print(model.exchanges)

# 'Ex_S_cpd00029_ext',    # propionate
# model.add_boundary(model.metabolites.get_by_id("co2_e"), type="exchange")
# model.add_metabolites(Metabolite('co2_e', name='CO2', compartment='e'))
# model.add_compartment('e', 'Extra-Cellular Space')
# model.add_compartment()
