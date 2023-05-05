from f_screen_SBML import *

modelName = "P_acnes_V2.xml"
save = False

loc_model = get_location(modelName)
model = cobra.io.read_sbml_model(loc_model)

GlucoseExchRnxId ='Ex_S_cpd00027_ext'
PropionateRnxID = 'Ex_S_cpd00141_ext'
AcetateRxnId  = 'Ex_S_cpd00029_ext'
BiomassRnxID = 'Ex_S_biomass_ext'

biomassRxn = model.reactions.get_by_id(BiomassRnxID)
PropRxn = model.reactions.get_by_id(PropionateRnxID)
PropRxn.bounds = 0, 1000
biomassRxn.bounds = 0,1000
print('the bounds of propioante are:', PropRxn.bounds)

uptakeRates = [1000]
for rates in uptakeRates:

    glucoseRxnExch = model.reactions.get_by_id(GlucoseExchRnxId)
    glucoseRxnExch.bounds = -rates, 1000 # make sure it's negative
    model.optimize()
    y_prop = find_yield(model=model, substrateExchangeRxnID= GlucoseExchRnxId, productExchangeRxnID= PropionateRnxID)
    y_bm = find_yield(model=model, substrateExchangeRxnID= GlucoseExchRnxId, productExchangeRxnID= BiomassRnxID)
    print('')
    print('the uptake rate of substrate is:', glucoseRxnExch.flux)
    print('the bm flux is:', biomassRxn.flux)
    print('the flux of product is:',PropRxn.flux )
    print('the yield (g/g) of propionate is:', y_prop)
    print('the yield (g/g) of biomass is:', y_bm)
    print('the yield (mol/mol) of biomass is:', biomassRxn.flux/glucoseRxnExch.flux)


# let's try and map the model
import json
import cobra.io.json
# Serialize the model to JSON format using cobra.io.json.to_json
model_json = cobra.io.json.to_json(model)

# Write the JSON-formatted string to a file
with open('acnesJ.json', 'w') as f:
    json.dump(model_json, f)
