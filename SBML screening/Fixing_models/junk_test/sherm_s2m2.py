"""
This script was quickly made to brainstrom some possible solutions/ improvements to the model during the s2m2 workshop
"""
from f_screen_SBML import *

modelName = "P_sherm_model.xml"
loc_sher = get_location(modelName)
model = cobra.io.read_sbml_model(loc_sher)

# check the yields, more realistic now?
exchangeRxnId_Glucose = 'Ex_S_cpd00027_ext'
exchangeRxnId_Propionate = 'Ex_S_cpd00141_ext'
exchangeRxnId_Acetate = 'Ex_S_cpd00029_ext'
exchangeRxnId_Biomass = 'Ex_S_biomass_ext'

yProp = find_yield(model= model, substrateExchangeRxnID= exchangeRxnId_Glucose,
                   productExchangeRxnID= exchangeRxnId_Propionate)
print('The yield of propionate is: {} \n'.format(yProp))
yAce = find_yield(model= model, substrateExchangeRxnID= exchangeRxnId_Glucose,
                  productExchangeRxnID= exchangeRxnId_Acetate)
print('The yield of acetate is: {} \n'.format(yAce))
yBm = find_yield(model= model, substrateExchangeRxnID= exchangeRxnId_Glucose,
                 productExchangeRxnID= exchangeRxnId_Biomass,Biomass = True)
print('The yield of biomass is: {} \n'.format(yBm))

atpMet = model.metabolites.get_by_id('S_cpd00002_c0')
#print_all_rxn_of_metabolite(metabolite=atpMet, case='names', printFlux=True)
#find_maintanace_reaction(metabolite=atpMet, case='names', printFlux=True)
atpm = model.reactions.get_by_id('rxn00062_c0')
atpm.bounds = (6,1000)
print(atpm)

yBm = find_yield(model= model, substrateExchangeRxnID= exchangeRxnId_Glucose,
                 productExchangeRxnID= exchangeRxnId_Biomass,Biomass = True)
print('The NEW yield of biomass is: {} \n'.format(yBm))

yProp = find_yield(model= model, substrateExchangeRxnID= exchangeRxnId_Glucose,
                   productExchangeRxnID= exchangeRxnId_Propionate)
print('The NEW yield of propionate is: {} \n'.format(yProp))

print(model.summary())
# copy atpm from ecoli, or other gram negative bacteria, tip see modelSeed
# this looks like the cell maintenance reaction, but not really, is it?
'''
Bigg
Kegg
ModelSeed
MetaCyc
types of databases for annotation 

'''
'''
-1.0 ATP + -1.0 AMP +  = 2.0 ADP + 
the flux of reaction rxn00097_c0 is : 1.5640107976400617 mmol/gDW/h 
 the bounds are (-1000.0, 1000.0) mmol/gDW/h 
the compartment(s) of the reactions are {'c'} 
'''

'''
could consider changing the objective function eg 
minimise atp consumption 
maximise propionate production, setting bm to minimum flux (best option we think) 
minimise glucose consumption setting the biomass flux 

'''
model.objective = 'Ex_S_cpd00141_ext' # propionate production
model.objective_direction = 'max'
# set the bounds of biomass!!!!

# for minimising glucose consumption you have to maximise the objective (the flux is negative!!)
model.objective = 'Ex_S_cpd00027_ext'
model.objective_direction = 'max'
# set the bounds of biomass!!!!
