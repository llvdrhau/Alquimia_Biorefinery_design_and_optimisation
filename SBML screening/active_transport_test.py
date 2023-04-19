
from f_screen_SBML import *

# read the model
modelNames = ['P_acnes_lvdh.xml', 'P_avidum_lvdh.xml', 'P_propionicum_lvdh.xml', 'P_sherm_lvdh.xml',
              'PAC_4875_lvdh.xml']

BM_atp_old = []
y_prop_old = []
y_ace_old = []
y_bm_old = []
dictOld = {}

BM_atp_new = []
y_prop_new = []
y_ace_new = []
y_bm_new = []
dictNew = {}
for name in modelNames:
    loc_model = get_location(name)
    model = cobra.io.read_sbml_model(loc_model)

    # ----------------- split into compartments by re-moving and adding new boundry reactions
    # Get the original (faulty) exchange reactions (even though there is only one compartment for the time being,
    # It does a good job at guessing what they are in this case)
    model._compartments = {'c': 'Cytoplasma', 'e': 'Extracellular'}
    listExchRxn = model.exchanges
    for original_reaction in listExchRxn:
        # Create a new exchange reaction with the same metabolites as the original reaction
        metabolite = original_reaction.reactants[
            0]  # there is only going to be one metabolite in the reaction reactants
        metabolite.compartment = 'e'

    # ----------------------------------------------- get the data of interest of the old model
    # get the molecular weight of biomass
    metBM = model.metabolites.get_by_id('S_biomass_ext')
    biomassMW = metBM.formula_weight
    # get the ratio
    bm_atp_old = ATP_Biomass_Ratio(model=model, biomassRxnID='Ex_S_biomass_ext', ATPmetID='S_cpd00002_c0',
                                   modelName=name, printResults=False)
    BM_atp_old.append(bm_atp_old * biomassMW)

    # get the yields of propionate and acetate
    exchangeRxnId_Acetate = 'Ex_S_cpd00029_ext'
    exchangeRxnId_Propionate = 'Ex_S_cpd00141_ext'
    exchangeRxnId_Glucose = 'Ex_S_cpd00027_ext'
    exchangeRxnId_Biomass = 'Ex_S_biomass_ext'

    yield_prop = find_yield(model=model, substrateExchangeRxnID=exchangeRxnId_Glucose,
                            productExchangeRxnID=exchangeRxnId_Propionate)
    yield_ace = find_yield(model=model, substrateExchangeRxnID=exchangeRxnId_Glucose,
                           productExchangeRxnID=exchangeRxnId_Acetate)
    yield_bm = find_yield(model=model, substrateExchangeRxnID=exchangeRxnId_Glucose,
                           productExchangeRxnID=exchangeRxnId_Biomass)

    y_prop_old.append(yield_prop)
    y_ace_old.append(yield_ace)
    y_bm_old.append(yield_bm)

    # ----------------------------------------------- update the model
    # change the transport reaction (add atp consumption)
    transportRxnId_Acetate = 'rxn05488_c0'
    transportRxnId_Propionate = 'rxn05634_c0'

    # get the reaction
    transportRxn_Acetate = model.reactions.get_by_id(transportRxnId_Acetate)
    transportRxn_Propionate = model.reactions.get_by_id(transportRxnId_Propionate)

    ATP_Id = 'S_cpd00002_c0'
    ADP_Id = 'S_cpd00008_c0'
    PP_id = 'S_cpd00009_c0'  # id of phosphate
    H_id_intra = 'S_cpd00067_c0' # id of intracellular H+
    H_id_ext = 'S_cpd00067_ext' # id of extracellular H+

    # Left side of chemical eq => extracellular space
    # Right side of chemical eq => intracellular space
    metabolites2add = {
        ATP_Id: 1.0,  # consume ATP intracelular (right side of the chemical eq)
        ADP_Id: -1.0,  # produce ADP intracelular (left side of the chemical eq)
        PP_id: -1.0,  # produced Phosphorus
        H_id_intra: 2, # intracellular, there is already one H in the reaction so 3 H+ in total react
        H_id_ext:-2 # extracellular, there is already one H in the reaction
    }

    transportRxn_Acetate.add_metabolites(metabolites2add)
    transportRxn_Propionate.add_metabolites(metabolites2add)

    # ----------------------------------------------- get the data of interest of the NEW model
    # get the ratio
    bm_atp_new = ATP_Biomass_Ratio(model=model, biomassRxnID='Ex_S_biomass_ext', ATPmetID='S_cpd00002_c0',
                                   modelName=name, printResults=False)
    BM_atp_new.append(bm_atp_new * biomassMW)

    yield_prop = find_yield(model=model, substrateExchangeRxnID=exchangeRxnId_Glucose,
                            productExchangeRxnID=exchangeRxnId_Propionate)

    yield_ace = find_yield(model=model, substrateExchangeRxnID=exchangeRxnId_Glucose,
                           productExchangeRxnID=exchangeRxnId_Acetate)

    yield_bm =  find_yield(model=model, substrateExchangeRxnID=exchangeRxnId_Glucose,
                           productExchangeRxnID=exchangeRxnId_Biomass)

    y_prop_new.append(yield_prop)
    y_ace_new.append(yield_ace)
    y_bm_new.append(yield_bm)

DICT_old = {'BM_ATP (g/mol)': BM_atp_old, 'yield biomass': y_bm_old ,'yield acetate (g/g)': y_ace_old ,
            'yield propionate (g/g)': y_prop_old}

DICT_new = {'BM_ATP (g/mol)': BM_atp_new, 'yield biomass (g/g)': y_bm_new, 'yield acetate (g/g)': y_ace_new ,
            'yield propionate (g/g)': y_prop_new}

DF_old = pd.DataFrame(DICT_old, index=modelNames)
DF_new = pd.DataFrame(DICT_new, index=modelNames)

print(DF_old)
print('')
print(DF_new)