from f_make_surrogate_model import *

# ----- script specifications

saveSwitch = True
printSwitch = False
carbonCheckSwitch = False

# -------------------------- SBML models
loc_acidi = 'PAC_4875_V2.xml'
loc_sher = 'P_sherm_V2.xml'

loc_acnes = 'P_acnes_lvdh.xml'
loc_prop = 'P_propionicum_lvdh.xml'
loc_avidum = 'P_avidum_lvdh.xml'


# SBML surrogates
substrates = ['Ex_S_cpd00027_ext', # glucose
              'Ex_S_cpd00082_ext', # fructose
              ]

products = ['Ex_S_cpd00141_ext',   # acetate
            'Ex_S_cpd00029_ext',    # propionate
            'Ex_S_biomass_ext']     # biomass

#microorganisms = [loc_acidi, loc_acnes, loc_prop, loc_avidum, loc_sher]  # all microorganisms
#saveNames = ['v1_PAC.json', 'v1_acnes.json', 'v1_propionicum.json' ,'v1_avidum.json', 'v1_sherm.json'] # make save names

microorganisms = [loc_acidi, loc_sher]  # all microorganisms
saveNames = ['v1_PAC.json', 'v1_sherm.json'] # make save names


# define maximum allowed concentrations
maxConcentration = {'Propionate': 42.37e-3}  # max concentration of propionate is 0.042 kg/L or 42 g/L

# todo make sure this script works with SBML_2_json_V2!!! (it does not work yet)
for i, organism in enumerate(microorganisms):
    sbml = SBML_2_json(modelName= organism, substrate_exchange_rnx= substrates, product_exchange_rnx= products,
                       case= 'mass_yield', saveName= saveNames[i], save = saveSwitch,
                       printEq= printSwitch, maxConcentration= maxConcentration)


