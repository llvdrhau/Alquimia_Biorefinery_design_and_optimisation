from f_screen_SBML import *

# read the model
modelName = 'P_acnes_lvdh.xml'
save = False

loc_model = get_location(modelName)
model = cobra.io.read_sbml_model(loc_model)

# ----------------- split into compartments by re-moving and adding new boundry reactions

# Get the original (faulty) exchange reactions (even though there is only one compartment for the time being,
# It does a good job at guessing what they are in this case)
model._compartments = {'c': 'Cytoplasma', 'e': 'Extracellular'}
listExchRxn = model.exchanges
for original_reaction in listExchRxn:
    # Create a new exchange reaction with the same metabolites as the original reaction
    metabolite = original_reaction.reactants[0]  # there is only going to be one metabolite in the reaction reactants
    metabolite.compartment = 'e'

# define compartments to a reaction
# print(model.exchanges) #check if the same error occurs

model.optimize()

print('----------- Glucose in Extracellular -----------  \n')
ExGlucoseId = 'S_cpd00027_ext'
ExGlucoseMet = model.metabolites.get_by_id(ExGlucoseId)
print_all_rxn_of_metabolite(ExGlucoseMet, case='names', printFlux=True)

print('----------- Glucose in Cytosol -----------  \n')
cytoGlucose = model.metabolites.get_by_id('S_cpd00027_c0')
print_all_rxn_of_metabolite(cytoGlucose, case='names', printFlux=True)

# look at currrent transport mechanism for acetate and propionate

# print('----------- ACETATE reactions in Cytosol -----------  \n')
# acID = 'S_cpd00029_c0'
# intraMet = model.metabolites.get_by_id(acID)
# print_all_rxn_of_metabolite(intraMet, case='names', printFlux=True)

print('-----------  ACETATE going to Extracellular -----------  \n')
CytoEthId = 'S_cpd00029_ext'
extraMet = model.metabolites.get_by_id(CytoEthId)
print_all_rxn_of_metabolite(extraMet, case='names', printFlux=True)

print('----------- Propionate going to Extracellular -----------  \n')
propID = 'S_cpd00141_ext'
extraMet = model.metabolites.get_by_id(propID)
print_all_rxn_of_metabolite(extraMet, case='names', printFlux=True)

# -------------------------------------------------------------------------------------------------------------------
# let's add an atp consumption to the reaction rxn05488_c0 (transport) to account for active transport
transportRxnId_Acetate = 'rxn05488_c0'
transportRxnId_Propionate = 'rxn05634_c0'

exchangeRxnId_Acetate = 'Ex_S_cpd00029_ext'
exchangeRxnId_Propionate = 'Ex_S_cpd00141_ext'
exchangeRxnId_Glucose = 'Ex_S_cpd00027_ext'

mwGlu = 180     # g/mol
mwAce = 59.044  # g/mol
mwProp = 73.07  # g/mol

# get the old yield
fluxGlu = model.reactions.get_by_id(exchangeRxnId_Glucose).flux
fluxAce = model.reactions.get_by_id(exchangeRxnId_Acetate).flux
flux_prop = model.reactions.get_by_id(exchangeRxnId_Propionate).flux

yield_ace = - (fluxAce * mwAce) / (mwGlu * fluxGlu)
yield_prop = (flux_prop * mwAce) / (mwGlu * fluxGlu)
print('the old yield (g/g) of acetate is: {} \n'
      'the old yield (g/g) of popionate is {}'.format(yield_ace, yield_prop))
print('')
ratioOld = ATP_Biomass_Ratio(model=model, biomassRxnID='Ex_S_biomass_ext', ATPmetID='S_cpd00002_c0',
                             modelName='---OLD---')

# change the transport reaction (add atp consumption)
# get the reaction
transportRxn_Acetate = model.reactions.get_by_id(transportRxnId_Acetate)
transportRxn_Propionate = model.reactions.get_by_id(transportRxnId_Propionate)

ATP_Id = 'S_cpd00002_c0'
ADP_Id = 'S_cpd00008_c0'
PP_id = 'S_cpd00009_c0'  # id of phosphate

metabolites2add = {
    ATP_Id: -1.0,  # consumed
    ADP_Id: 1.0,  # produced
    PP_id: 1.0  # produced
}

transportRxn_Acetate.add_metabolites(metabolites2add)
transportRxn_Propionate.add_metabolites(metabolites2add)

# check if adding metabolites to reactions worked
# eq = string_reactions(reaction= transportRxn_Acetate, case= 'names')
# print(eq)

# ---------------------------------------------- check the effect of adding ATP consumption to the transport reactions
model.optimize()  # re-optimise
# get the NEW yield
fluxGlu = model.reactions.get_by_id(exchangeRxnId_Glucose).flux
fluxAce = model.reactions.get_by_id(exchangeRxnId_Acetate).flux
flux_prop = model.reactions.get_by_id(exchangeRxnId_Propionate).flux

yield_ace = - (fluxAce * mwAce) / (mwGlu * fluxGlu)
yield_prop = (flux_prop * mwAce) / (mwGlu * fluxGlu)
print('the NEW yield (g/g) of acetate is: {} \n'
      'the NEW yield (g/g) of popionate is {}'.format(yield_ace, yield_prop))

print('')
ratioONEw = ATP_Biomass_Ratio(model=model, biomassRxnID='Ex_S_biomass_ext', ATPmetID='S_cpd00002_c0',
                              modelName='---NEW---')

print('')
print('----------- Propionate going to Extracellular NEW -----------  \n')
propID = 'S_cpd00141_ext'
extraMet = model.metabolites.get_by_id(propID)
print_all_rxn_of_metabolite(extraMet, case='names', printFlux=True)


# find MW of Biomass
metBM = model.metabolites.get_by_id('S_biomass_ext')
mwBM = metBM.formula_weight
print(mwBM)
print('')