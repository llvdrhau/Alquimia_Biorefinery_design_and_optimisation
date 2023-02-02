from f_make_surrogate_model import *

loc_acidi = 'PAC_4875_lvdh.xml'
loc_acnes = 'P_acnes_lvdh.xml'
loc_prop = 'P_propionicum_lvdh.xml'
loc_avidum = 'P_avidum_lvdh.xml'
loc_sher = 'P_sherm_lvdh.xml'

# SBML surrogates
substrates = ['Ex_S_cpd00027_ext', # glucose
              'Ex_S_cpd00082_ext', # fructose
               'Ex_S_cpd00100_ext' # glycerol
              ]

products = ['Ex_S_cpd00141_ext',   # acetate
            'Ex_S_cpd00029_ext',    # propionate
            'Ex_S_biomass_ext']     # biomass

microorganisms = [loc_acidi, loc_acnes, loc_prop, loc_avidum]  # all microorganisms
saveNames = ['PAC.json', 'acnes.json', 'propionicum.json' ,'avidum.json']
for i, organism in enumerate(microorganisms):
    sbml = SBML_2_json(modelName= organism, substrate_exchange_rnx= substrates, product_exchange_rnx= products,
                       saveName= saveNames[i], save = True, checkCarbon= False, printEq= True)

# the model P_sherm_model can not consume glycerol so only screen it for glucose and fructose 'sherm'

organism = loc_sher
saveName = 'sherm.json'

substrates = ['Ex_S_cpd00027_ext', # glucose
              'Ex_S_cpd00082_ext', # fructose
              ]

products = ['Ex_S_cpd00141_ext',   # acetate
            'Ex_S_cpd00029_ext',    # propionate
            'Ex_S_biomass_ext']     # biomass

sbmlSherm = SBML_2_json(modelName=organism, substrate_exchange_rnx=substrates, product_exchange_rnx=products,
                       saveName=saveName, save=True, checkCarbon=False, printEq=True)

# # open fermentation
# excelFile =  'Glucose_PH_Data.xlsx' # with a 5th polynomial it is acctually quite nice....
# polynomial = 5
# out = regression_2_json(excelFile, normalise= False ,save=False, saveName='open_fermentation_polynomial_ESCAPE33.json',
#                          showPLot= False, polynomial= {'pH':polynomial}, case= 'Ridge')
