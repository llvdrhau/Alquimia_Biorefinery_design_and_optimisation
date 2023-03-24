from f_make_surrogate_model import *
"""
This script creates the surrogate models of the GEMs. Json files are made from which the superstructure can then read 
the appropriate mathematical expression 

In this script you will be ask to include/exclude substrates. Based on the history of your choice the same substrates 
will be included or excluded during the loop to minimise time selecting substrates manually  
"""
# ----- script specifications

saveSwitch = True
printSwitch = False
carbonCheckSwitch = False

# -------------------------- SBML models
loc_acidi = 'PAC_4875_V2.xml'
loc_acnes = 'P_acnes_V2.xml'
loc_prop = 'P_propionicum_V2.xml'
loc_avidum = 'P_avidum_V2.xml'
loc_sher = 'P_sherm_V2.xml'

# SBML surrogate options
substrates = 'select'

products = ['Ex_S_cpd00141_ext',  # propionate
            'Ex_S_cpd00029_ext',  # acetate
            'Ex_S_biomass_ext']  # biomass

# propionate needs a minimum yield of 0.1
tolerance = {'Ex_S_cpd00141_ext': 0.15}
# define maximum allowed concentrations
maxConcentration = {'Propionate': 0.018}  # max concentration of propionate is 0.018 kg/L or 18 g/L
# substrates = ['Ex_S_cpd00027_ext']  # glucose
microorganisms = [loc_acnes, loc_acidi, loc_prop, loc_avidum, loc_sher]  # all microorganisms
saveNames = ['acnes_v2.json', 'PAC_v2.json', 'propionicum_v2.json', 'avidum_v2.json',
             'sherm_v2.json']  # make save names
alreadyConsidered = ['D-Glucose', 'D-Fructose', 'Glycerol', 'L-Lactate', 'Maltose']
toIgnore = ['D-Lyxitol', 'alpha-Methyl-D-glucoside', 'CELB', 'TRHL', 'GLCN', 'Succinate', 'Glycerol-3-phosphate',
            'beta-Methylglucoside', 'Glycerone', 'Pyruvate', 'LACT', 'D-Mannitol', 'ACTN', 'L-Inositol', 'Ribitol',
            'Thyminose', 'D-Ribose', 'Glucuronate', 'Fumarate', '2-Oxobutyrate', 'L-Malate', 'DTYL', '2-Oxoglutarate',
            'D-Arabinose', 'D-Mannose', 'Salicin']

allConsideredSubstrates = []
allIgnoredSubstrate = []
for i, organism in enumerate(microorganisms):
    objec, considered, ignored = SBML_2_json_v2(modelName=organism, substrate_exchange_rnx='select',
                                                product_exchange_rnx=products,
                                                yieldTol=tolerance, saveName=saveNames[i], save=True,
                                                toIgnore=toIgnore, alreadyConsidered=alreadyConsidered)

    allConsideredSubstrates += considered
    allIgnoredSubstrate += ignored

    alreadyConsidered = list(set(allConsideredSubstrates))
    toIgnore = list(set(allIgnoredSubstrate))

# print the original
print(list(set(allConsideredSubstrates)))
print(list(set(allIgnoredSubstrate)))
