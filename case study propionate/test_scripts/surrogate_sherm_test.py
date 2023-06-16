from f_make_surrogate_model import *

# ----- script specifications

saveSwitch = True
printSwitch = True
carbonCheckSwitch = True

# -------------------------- SBML models

modelName = 'P_sherm_v2.xml'
microorganisms = [modelName]  # all microorganisms
saveNames = ['sherm_test.json'] # make save names

# SBML surrogates
substrates = ['Ex_S_cpd00027_ext', # glucose
              'Ex_S_cpd00082_ext', # fructose
              'Ex_S_cpd00159_ext'  # lactate
              ]

products = ['Ex_S_cpd00141_ext',   # propionate
            'Ex_S_cpd00029_ext',    # acetate
            'Ex_S_biomass_ext']     # biomass



# define maximum allowed concentrations
maxConcentration = {'Propionate': 42.37e-3}  # max concentration of propionate is 0.042 kg/L or 42 g/L

for i, organism in enumerate(microorganisms):

    sbml = SBML_2_json_v2(modelName=organism, substrate_exchange_rnx=substrates, product_exchange_rnx=products,
                        saveName=saveNames[i], save=saveSwitch,
                        maxConcentration=maxConcentration)


