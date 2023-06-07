
from f_usefull_functions import get_location
from f_make_surrogate_model import find_yield
import cobra

modelName  = 'P_sherm_test.xml'
loc = get_location(modelName)
model = cobra.io.read_sbml_model(loc)

#
substrateExchangeRxnID = 'Ex_S_cpd00027_ext' # glucose
productExchangeRxnID = 'Ex_S_cpd00029_ext' # acetaqte

glucoseRnx = model.reactions.get_by_id(substrateExchangeRxnID)
print(glucoseRnx.bounds)
glucoseRnx.bounds = (-10, 1000)

# check the yields, more realistic now?
exchangeRxnId_Acetate = 'Ex_S_cpd00029_ext'
exchangeRxnId_Propionate = 'Ex_S_cpd00141_ext'
exchangeRxnId_Glucose = 'Ex_S_cpd00027_ext'
exchangeRxnId_Biomass = 'Ex_S_biomass_ext'

yProp = find_yield(model= model, substrateExchangeRxnID= exchangeRxnId_Glucose, productExchangeRxnID= exchangeRxnId_Propionate)
print('The yield of propionate is: {} \n'.format(yProp))

yAce = find_yield(model= model, substrateExchangeRxnID= exchangeRxnId_Glucose, productExchangeRxnID= exchangeRxnId_Acetate)
print('The yield of acetate is: {} \n'.format(yAce))

yBm = find_yield(model= model, substrateExchangeRxnID= exchangeRxnId_Glucose, productExchangeRxnID= exchangeRxnId_Biomass)
print('The yield of biomass is: {} \n'.format(yBm))


