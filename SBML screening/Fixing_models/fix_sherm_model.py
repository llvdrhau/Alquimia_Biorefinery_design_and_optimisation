'''
This script intends to fix the SBML (GEMs) of "P_sherm_model.xml"
propionibacterium freudenreichii (also known as shermanii)

lucas.vanderhauwaert@usc.es
07/07/2023
'''

from f_screen_SBML import *

modelName = "P_sherm_model.xml"
loc_sher = get_location(modelName)
model = cobra.io.read_sbml_model(loc_sher)

# check the yields, are they realistic?
exchangeRxnId_Glucose = 'Ex_S_cpd00027_ext'
exchangeRxnId_Propionate = 'Ex_S_cpd00141_ext'
exchangeRxnId_Acetate = 'Ex_S_cpd00029_ext'
exchangeRxnId_Biomass = 'Ex_S_biomass_ext'

# yield of propionate
yProp = find_yield(model=model, substrateExchangeRxnID=exchangeRxnId_Glucose,
                   productExchangeRxnID=exchangeRxnId_Propionate)
print('The yield of propionate is: {} \n'.format(yProp))

# yield of acetate
yAce = find_yield(model=model, substrateExchangeRxnID=exchangeRxnId_Glucose,
                  productExchangeRxnID=exchangeRxnId_Acetate)
print('The yield of acetate is: {} \n'.format(yAce))

# yield of Biomass
yBm = find_yield(model=model, substrateExchangeRxnID=exchangeRxnId_Glucose,
                 productExchangeRxnID=exchangeRxnId_Biomass, Biomass=True)
print('The yield of biomass is: {} \n'.format(yBm))


# there is no flux going through the cell maintanace reaction. Tippically annotated by atpm (m= maintanace)
# The models are orriginally from ModelSeed, so you can look up atpm in ModelSeed and find the reaction id: rxn00062_c0
atpm = model.reactions.get_by_id('rxn00062_c0')
# let's set the bound of the flux to that of another gram negative bacteria: ecoli! iJO1366.xml 3.15 mmol/gDW/h
#atpm.bounds = (3.15, 1000)
atpm.bounds = (6, 1000)
#print(atpm)


print("#################### after changing bounds #######################")
yBm = find_yield(model=model, substrateExchangeRxnID=exchangeRxnId_Glucose,
                 productExchangeRxnID=exchangeRxnId_Biomass, Biomass=True)
print('The NEW yield of biomass is: {} \n'.format(yBm))

yProp = find_yield(model=model, substrateExchangeRxnID=exchangeRxnId_Glucose,
                   productExchangeRxnID=exchangeRxnId_Propionate)
print('The NEW yield of propionate is: {} \n'.format(yProp))




#
# print(model.summary())
# # copy atpm from ecoli, or other gram negative bacteria, tip see modelSeed
# # this looks like the cell maintenance reaction, but not really, is it?



# '''
# Bigg
# Kegg
# ModelSeed
# MetaCyc
# types of databases for annotation
#
# '''
# '''
# -1.0 ATP + -1.0 AMP +  = 2.0 ADP +
# the flux of reaction rxn00097_c0 is : 1.5640107976400617 mmol/gDW/h
#  the bounds are (-1000.0, 1000.0) mmol/gDW/h
# the compartment(s) of the reactions are {'c'}
# '''
#
# '''
# could consider changing the objective function eg
# minimise atp consumption
# maximise propionate production, setting bm to minimum flux (best option we think)
# minimise glucose consumption setting the biomass flux
#
# '''
# model.objective = 'Ex_S_cpd00141_ext'
# model.objective_direction = 'max'
# # set the bounds of biomass!!!!
#
# # for minimising glucose consumption you have to maximise the objective (the flux is negative!!)
# model.objective = 'Ex_S_cpd00027_ext'
# model.objective_direction = 'max'
# # set the bounds of biomass!!!!
#
#
# # ----------------- split into compartments by re-moving and adding new boundry reactions
# # Get the original (faulty) exchange reactions (even though there is only one compartment for the time being,
# # It does a good job at guessing what they are in this case)
# model._compartments = {'c': 'Cytoplasma', 'e': 'Extracellular'}
# listExchRxn = model.exchanges
# for original_reaction in listExchRxn:
#     # Create a new exchange reaction with the same metabolites as the original reaction
#     metabolite = original_reaction.reactants[0]  # there is only going to be one metabolite in the reaction reactants
#     metabolite.compartment = 'e'
#
# # save the model
# if save:
#     write_sbml_model(model, r"C:\Users\lucas\PycharmProjects\Alquimia\SBML models\P_sherm_test.xml")
#
