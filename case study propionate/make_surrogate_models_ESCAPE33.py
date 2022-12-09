
from f_make_surrogate_model import *

# open fermentation
excelFile =  'Glucose_PH_Data.xlsx' # with a 5th polynomial it is acctually quite nice....
polynomial = 5
out = regression_2_json(excelFile, normalise= False ,save=True, saveName='open_fermentation_polynomial_ESCAPE33.json',
                         showPLot= False, polynomial= {'pH':polynomial}, case= 'Ridge')

# SBML surrogates
substrates = ['Ex_S_cpd00027_ext', # glucose
              'Ex_S_cpd00082_ext'] # fructose
products = ['Ex_S_cpd00141_ext',   # acetate
            'Ex_S_cpd00029_ext']   # propionate

loc_acidi = 'PAC_4875_model.xml'
loc_acnes = 'P_acnes_model.xml'
loc_prop = 'P_propionicum_model.xml'
loc_avidum = 'P_avidum_model.xml'
loc_sher = 'P_sherm_model.xml'
microorganisms = [loc_acidi, loc_acnes, loc_prop, loc_avidum, loc_sher]  # all microorganisms
for i in microorganisms:
    sbml = SBML_2_json(modelName= i, substrate_exchange_rnx= substrates, product_exchange_rnx= products,
                       save = True, checkCarbon= False, printEq= True)