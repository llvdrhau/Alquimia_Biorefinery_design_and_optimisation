
from f_screen_SBML import find_yield, string_reactions
from f_usefull_functions import get_location
import cobra

# load the model
#modelName = "PAC_4875_model.xml"
#modelName = "PAC_4875_V2.xml"
modelName = "P_sherm_V2.xml"
loc_sher = get_location(modelName)
model = cobra.io.read_sbml_model(loc_sher)

# ---------------------------------------------------------------------------------------------------------------------
# split into compartments by re-moving and adding new boundry reactions
# Get the original (faulty) exchange reactions (even though there is only one compartment for the time being,
# It does a good job at guessing what they are in this case)
model._compartments = {'c': 'Cytoplasma', 'e': 'Extracellular'}
listExchRxn = model.exchanges
for original_reaction in listExchRxn:
    # Create a new exchange reaction with the same metabolites as the original reaction
    metabolite = original_reaction.reactants[0]  # there is only going to be one metabolite in the reaction reactants
    metabolite.compartment = 'e'

# print in green color that the model now has two compartments
print('\033[92mAn extra compartment has been added to the model\033[0m')

# there is no flux going through the cell maintenance reaction. Typically annotated by atpm (m= maintenance)
# The models are originally from ModelSeed, so you can look up atpm in ModelSeed and find the reaction id: rxn00062_c0
atpm = model.reactions.get_by_id('rxn00062_c0')
strATPm = string_reactions(atpm, printFlux=True)
print('')
print("The maintenance reaction is:")
print(strATPm[0])
print('the bounds of the reaction are: {}'.format(atpm.bounds))
print('The flux of the maintanance reaction is: {}'.format(strATPm[1]))

# let's set the bound of the flux to that of another gram negative bacteria: ecoli! iJO1366.xml 3.15 mmol/gDW/h
lbATPm = 9 # 3.15 # mmol/gDW/h
atpm.bounds = (lbATPm, 1000)

# glucose exchange reaction
exchangeRxnId_Glucose = 'Ex_S_cpd00027_ext'
# propionate exchange reaction
exchangeRxnId_Propionate = 'Ex_S_cpd00141_ext'

# set the objective function to biomass

# check if lactate can be consumed by the model
print(" --------------------------------------------- \n")
print("Check if the model can consume L-Lactate \n")

lactateExchangeReaction = model.reactions.get_by_id('Ex_S_cpd00159_ext')
model.reactions.get_by_id(exchangeRxnId_Glucose).bounds = (0, 1000)
lactateExchangeReaction.bounds = (-10, 1000)
print('the bounds of lactate are',lactateExchangeReaction.bounds)
print('the bounds of glucose are',model.reactions.get_by_id(exchangeRxnId_Glucose).bounds)

# find the yield of propionate from lactate
yProp = find_yield(model=model, substrateExchangeRxnID='Ex_S_cpd00159_ext', productExchangeRxnID=exchangeRxnId_Propionate,
                   printResults= True)
print('the yield of propionate from lactate is: {}'.format(yProp))

# reset the bounds of lactate
lactateExchangeReaction.bounds = (0, 1000)

# ---------------------------------------------------------------------------------------------------------------------
# check if the model can consume glycerol and what is the yield of propionate
print(" --------------------------------------------- \n")
print("Check if the model can consume Glycerol \n")
GlycerolExchangeReaction = model.reactions.get_by_id('Ex_S_cpd00100_ext')
GlycerolExchangeReaction.bounds = (-10, 1000)

# find the yield of propionate from glycerol
yProp = find_yield(model=model, substrateExchangeRxnID='Ex_S_cpd00100_ext',
                   productExchangeRxnID=exchangeRxnId_Propionate, printResults= True)
