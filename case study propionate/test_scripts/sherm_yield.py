
from f_usefull_functions import get_location
from f_make_surrogate_model import find_yield
import cobra

modelName  = 'P_propionicum_V2.xml'
#modelName  = 'P_sherm_test.xml'

loc = get_location(modelName)
model = cobra.io.read_sbml_model(loc)

# substrate exchange rate
exchangeRxnId_Glucose = 'Ex_S_cpd00027_ext'
substrateExchangeRxnID =  exchangeRxnId_Glucose
glucoseRnx = model.reactions.get_by_id(substrateExchangeRxnID)

print(glucoseRnx.bounds)


# check the yields, more realistic now?
exchangeRxnId_Propionate = 'Ex_S_cpd00141_ext'
exchangeRxnId_Acetate = 'Ex_S_cpd00029_ext'
exchangeRxnId_Biomass = 'Ex_S_biomass_ext'


print('-----------------using the function ---------------------')

yProp = find_yield(model= model, substrateExchangeRxnID= exchangeRxnId_Glucose, productExchangeRxnID= exchangeRxnId_Propionate)
print('The yield of propionate is: {} \n'.format(yProp))
yAce = find_yield(model= model, substrateExchangeRxnID= exchangeRxnId_Glucose, productExchangeRxnID= exchangeRxnId_Acetate)
print('The yield of acetate is: {} \n'.format(yAce))
yBm = find_yield(model= model, substrateExchangeRxnID= exchangeRxnId_Glucose, productExchangeRxnID= exchangeRxnId_Biomass)
print('The yield of biomass is: {} \n'.format(yBm))



# what happens in the function is that instead of finding the yields per substrate, we find the yields per product
# in other words more efficient to loop over the substrates in stead of products
# # let's check without using the function
# model.optimize()
# # get reactions
# acetateRxn = model.reactions.get_by_id(exchangeRxnId_Acetate)
# propionateRxn = model.reactions.get_by_id(exchangeRxnId_Propionate)
# biomassRxn = model.reactions.get_by_id(exchangeRxnId_Biomass)
#
# # get fluxes
# glucoseFlux = glucoseRnx.flux
# acetateFlux = acetateRxn.flux
# propionateFlux = propionateRxn.flux
# bmFlux = biomassRxn.flux
#
# # get MW
# # get the metabolite objects
# glucoseMetabolite = glucoseRnx.reactants[0]
# acetateMetabolite = acetateRxn.reactants[0]
# propionateMetabolite = propionateRxn.reactants[0]
# bmMetabolite = biomassRxn.reactants[0]
#
#
#
# glucoseMW = glucoseMetabolite.formula_weight
# acetateMW = acetateMetabolite.formula_weight
# propionateMW = propionateMetabolite.formula_weight
# bmMW = bmMetabolite.formula_weight
#
# # yields
# y_prop = (propionateFlux/glucoseFlux)* (propionateMW/glucoseMW)
# y_ace = (acetateFlux/glucoseFlux)* (acetateMW/glucoseMW)
# y_bm = (bmFlux/glucoseFlux)*(bmMW/glucoseMW)
#
#
# print('The yield of propionate is: {} \n'.format(y_prop))
# print('The yield of acetate is: {} \n'.format(y_ace))
# print('The yield of biomass is: {} \n'.format(y_bm))
#
# model.optimize()

