'''
This script intends to fix the SBML (GEMs) of "PAC_4875_model.xml"
propionibacterium acidipropionici

lucas.vanderhauwaert@usc.es
12/07/2023
'''

from f_screen_SBML import find_yield, string_reactions
from f_usefull_functions import get_location
import cobra
from cobra.io import write_sbml_model

# save switch
saveModel = False
# load the model
modelName = "PAC_4875_model.xml"
loc_sher = get_location(modelName)
model = cobra.io.read_sbml_model(loc_sher)
# ---------------------------------------------------------------------------------------------------------------------
# give the metabolite biomass_c0 a name
model.metabolites.get_by_id('S_biomass_ext').name = 'Biomass'
# check the biomass reaction of the model
biomassRxn = model.reactions.get_by_id('biomass_c0')
formulaBiomass = string_reactions(biomassRxn, case='names')
print(formulaBiomass[0])
print('')
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

print(' --------------------------------------------- \n')

# check the yields, are they realistic?
exchangeRxnId_Glucose = 'Ex_S_cpd00027_ext'
exchangeRxnId_Propionate = 'Ex_S_cpd00141_ext'
exchangeRxnId_Acetate = 'Ex_S_cpd00029_ext'
exchangeRxnId_Biomass = 'Ex_S_biomass_ext'

# yield of propionate
yProp = find_yield(model=model, substrateExchangeRxnID=exchangeRxnId_Glucose,
                   productExchangeRxnID=exchangeRxnId_Propionate)
# yield of acetate
yAce = find_yield(model=model, substrateExchangeRxnID=exchangeRxnId_Glucose,
                  productExchangeRxnID=exchangeRxnId_Acetate)
# yield of Biomass
yBm = find_yield(model=model, substrateExchangeRxnID=exchangeRxnId_Glucose,
                 productExchangeRxnID=exchangeRxnId_Biomass, biomass=True)

# print the yields
print('The yield of propionate is: {}'.format(yProp))
print('The yield of acetate is: {}'.format(yAce))
print('The yield of biomass is: {}'.format(yBm))

# get the flux of biomass by finding the reaction id of the biomass exchange reaction
reactionBm = model.reactions.get_by_id(exchangeRxnId_Biomass)
fluxBm = reactionBm.flux
print('The growth rate of biomass is: {} gBM/gDW/h '.format(fluxBm))

# ---------------------------------------------------------------------------------------------------------------------
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


print(" --------------------------------------------- \n")
print("The minimum flux for the maintenance reaction has been set to {} mmol/gDW/h \n".format(lbATPm))

# find the yield of propionate and biomass
yProp = find_yield(model=model, substrateExchangeRxnID=exchangeRxnId_Glucose,
                   productExchangeRxnID=exchangeRxnId_Propionate)
yAce = find_yield(model=model, substrateExchangeRxnID=exchangeRxnId_Glucose,
                  productExchangeRxnID= exchangeRxnId_Acetate)
yBm = find_yield(model=model, substrateExchangeRxnID=exchangeRxnId_Glucose,
                 productExchangeRxnID=exchangeRxnId_Biomass, biomass=True)

#print the yields
print('The yield of propionate is: {}'.format(yProp))
print('The yield of acetate is: {}'.format(yAce))
print('The yield of biomass is: {} '.format(yBm))


# get the flux of biomass by finding the reaction id of the biomass exchange reaction
reactionBm = model.reactions.get_by_id(exchangeRxnId_Biomass)
fluxBm = reactionBm.flux
print('The growth rate of biomass is: {} gBM/gDW/h '.format(fluxBm))


# ---------------------------------------------------------------------------------------------------------------------
# the Biomass yield is still too high, now let's set the lower bound of the biomass exchange reaction to 95% of the flux
# found in the original model and set the objective to the propionate exchange reaction
print(" --------------------------------------------- \n")
print("The minimum flux for maintenance + objective set to propionate \n")

# the flux of the biomass reaction when the model is optimized for biomass
maxFluxBiomass = fluxBm * 0.99
model.reactions.get_by_id('Ex_S_biomass_ext').bounds = (maxFluxBiomass, 1000)
model.objective = exchangeRxnId_Propionate
model.objective_direction = 'max' # maximize the production of propionate


# find the yield of propionate
yProp = find_yield(model=model, substrateExchangeRxnID=exchangeRxnId_Glucose,
                   productExchangeRxnID=exchangeRxnId_Propionate)
# find the yield of acetate
yAce = find_yield(model=model, substrateExchangeRxnID=exchangeRxnId_Glucose,
                  productExchangeRxnID=exchangeRxnId_Acetate)
# find biomass yield
yBm = find_yield(model=model, substrateExchangeRxnID=exchangeRxnId_Glucose,
                 productExchangeRxnID=exchangeRxnId_Biomass, biomass=True)

#print the yields
print('The yield of propionate is: {}'.format(yProp))
print('The yield of acetate is: {}'.format(yAce))
print('The yield of biomass is: {} '.format(yBm))


# ---------------------------------------------------------------------------------------------------------------------
# let's still set the lower bound of the biomass exchange reaction to 95% of the flux
# found in the original model and set the objective to minimise the glucose exchange reaction!!
print(" --------------------------------------------- \n")
print("The minimum flux for maintenance + objective set to minimise the substrate uptake rate \n")

# the flux of the biomass reaction when the model is optimized for biomass
model.reactions.get_by_id('Ex_S_biomass_ext').bounds = (maxFluxBiomass, 1000)
model.objective = exchangeRxnId_Glucose
model.objective_direction = 'max' # maximize because we want to minimize the glucose uptake rate which is negative


# find the yield of propionate
yProp = find_yield(model=model, substrateExchangeRxnID=exchangeRxnId_Glucose,
                   productExchangeRxnID=exchangeRxnId_Propionate)
# find the yield of acetate
yAce = find_yield(model=model, substrateExchangeRxnID=exchangeRxnId_Glucose,
                  productExchangeRxnID=exchangeRxnId_Acetate)
# find biomass yield
yBm = find_yield(model=model, substrateExchangeRxnID=exchangeRxnId_Glucose,
                 productExchangeRxnID=exchangeRxnId_Biomass, biomass=True)

#print the yields
print('The yield of propionate is: {}'.format(yProp))
print('The yield of acetate is: {}'.format(yAce))
print('The yield of biomass is: {} '.format(yBm))
print("")
print("The glucose uptake rate is: {} mmol/gDW/h".format(model.reactions.get_by_id(exchangeRxnId_Glucose).flux))


#----------------------------------------------------------------------------------------------------------------------


# set the model whith the desired objective function and bounds and save it
# maximise the biomass production
model.objective = 'biomass_c0'
model.objective_direction = 'max'
# reset the bounds of the biomass exchange reaction
model.reactions.get_by_id('Ex_S_biomass_ext').bounds = (0, 1000)

# save the model
if saveModel:
    newModelName = r"C:\Users\lucas\PycharmProjects\Alquimia\SBML models\PAC_4875_V1.xml"
    write_sbml_model(model, newModelName)
    print("The model has been saved with biomass as the OF \n")


# ---------------------------------------------------------------------------------------------------------------------

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
yAce = find_yield(model=model, substrateExchangeRxnID='Ex_S_cpd00159_ext', productExchangeRxnID=exchangeRxnId_Acetate,
                  printResults= True)
yBm = find_yield(model=model, substrateExchangeRxnID='Ex_S_cpd00159_ext', productExchangeRxnID=exchangeRxnId_Biomass,
                 biomass=True, printResults= True)
#print('the yield of propionate from lactate is: {}'.format(yProp))

# reset the bounds of lactate
lactateExchangeReaction.bounds = (0, 1000)

# ---------------------------------------------------------------------------------------------------------------------
# check if the model can consume glycerol and what is the yield of propionate
print(" --------------------------------------------- \n")
print("Check if the model can consume Glycerol \n")
GlycerolExchangeReaction = model.reactions.get_by_id('Ex_S_cpd00100_ext')
GlycerolExchangeReaction.bounds = (-10, 1000)

print('Check bounds:')
print('the bounds of glycerol are',GlycerolExchangeReaction.bounds)
print('the bounds of glucose are',model.reactions.get_by_id(exchangeRxnId_Glucose).bounds)
print('the bounds of lactate are',lactateExchangeReaction.bounds)
print('')

# find the yield of propionate from glycerol
yProp = find_yield(model=model, substrateExchangeRxnID='Ex_S_cpd00100_ext',
                   productExchangeRxnID=exchangeRxnId_Propionate, printResults= True)
yAce = find_yield(model=model, substrateExchangeRxnID='Ex_S_cpd00100_ext',
                  productExchangeRxnID=exchangeRxnId_Acetate,printResults= True)
yBm = find_yield(model=model, substrateExchangeRxnID='Ex_S_cpd00100_ext', productExchangeRxnID=exchangeRxnId_Biomass,
                 biomass=True, printResults= True)