from f_make_surrogate_model import *

# ----- script specifications

saveSwitch = True
printSwitch = False
carbonCheckSwitch = False

# -------------------------- SBML models
loc_acidi = 'PAC_4875_lvdh.xml'
loc_acnes = 'P_acnes_lvdh.xml'
loc_prop = 'P_propionicum_lvdh.xml'
loc_avidum = 'P_avidum_lvdh.xml'
loc_sher = 'P_sherm_lvdh.xml'

# SBML surrogates
substrates = ['Ex_S_cpd00027_ext', # glucose
              'Ex_S_cpd00082_ext', # fructose
              ]

products = ['Ex_S_cpd00141_ext',   # acetate
            'Ex_S_cpd00029_ext',    # propionate
            'Ex_S_biomass_ext']     # biomass

microorganisms = [loc_acidi, loc_acnes, loc_prop, loc_avidum]  # all microorganisms
saveNames = ['v1_PAC.json', 'v1_acnes.json', 'v1_propionicum.json' ,'v1_avidum.json'] # make save names

# define maximum allowed concentrations
maxConcentration = {'Propionate' : 0.018} # max concentration of propionate is 0.018 kg/L

for i, organism in enumerate(microorganisms):
    sbml = SBML_2_json(modelName= organism, substrate_exchange_rnx= substrates, product_exchange_rnx= products,
                       case= 'mass_yield', saveName= saveNames[i], save = saveSwitch,
                       printEq= printSwitch, maxConcentration= maxConcentration)


# --------------------
# the model P_sherm_model can not consume glycerol so only screen it for glucose and fructose 'sherm'

organism = loc_sher
saveName = 'v1_sherm.json'

substrates = ['Ex_S_cpd00027_ext', # glucose
              'Ex_S_cpd00082_ext', # fructose
              ]

products = ['Ex_S_cpd00141_ext',   # acetate
            'Ex_S_cpd00029_ext',    # propionate
            'Ex_S_biomass_ext']     # biomass

maxProductConcentration = {'Propionate': 0.030} # in kg/L
sbmlSherm = SBML_2_json(modelName=organism, substrate_exchange_rnx=substrates, product_exchange_rnx=products,
                       case= 'mass_yield', saveName=saveName, save = saveSwitch,
                       printEq= printSwitch, maxConcentration= maxProductConcentration)