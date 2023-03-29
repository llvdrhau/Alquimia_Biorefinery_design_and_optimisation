from f_make_surrogate_model import *
"""
This script creates the surrogate models of the GEMs. Json files are made from which the superstructure can then read 
the appropriate mathematical expression 

In this script you will be ask to include/exclude substrates. Based on the history of your choice the same substrates 
will be included or excluded during the loop to minimise time selecting substrates manually  
"""
# ----- script specifications
saveSwitch = True

# -------------------------- SBML models
loc_acidi = 'PAC_4875_V2.xml'
loc_acnes = 'P_acnes_V2.xml'
loc_prop = 'P_propionicum_V2.xml'
loc_avidum = 'P_avidum_V2.xml'
loc_sher = 'P_sherm_V2.xml'

# select substrates qnd products
substrates = 'select'  # substrates = ['Ex_S_cpd00027_ext']  # glucose
products = ['Ex_S_cpd00141_ext',  # propionate
            'Ex_S_cpd00029_ext',  # acetate
            'Ex_S_biomass_ext']  # biomass

# propionate needs a minimum yield of 0.1 for the possible substrate selection procedure
tolerance = {'Ex_S_cpd00141_ext': 0.15}

# define maximum allowed concentrations in reactor
maxConcentration = {'Propionate': 0.071}  # max concentration of propionate is 0.018 kg/L or 18 g/L

# give the names of the organisms and the save names of the json files
microorganisms = [loc_acnes, loc_acidi, loc_prop, loc_avidum, loc_sher]  # all microorganisms
saveNames = ['acnes_v2.json', 'PAC_v2.json', 'propionicum_v2.json', 'avidum_v2.json','sherm_v2.json']

# give list of components to already ignore or include
alreadyConsidered = ['D-Glucose', 'D-Fructose', 'Glycerol', 'L-Lactate', 'Maltose']
toIgnore = ['D-Lyxitol', 'alpha-Methyl-D-glucoside', 'CELB', 'TRHL', 'GLCN', 'Succinate', 'Glycerol-3-phosphate',
            'beta-Methylglucoside', 'Glycerone', 'Pyruvate', 'LACT', 'D-Mannitol', 'ACTN', 'L-Inositol', 'Ribitol',
            'Thyminose', 'D-Ribose', 'Glucuronate', 'Fumarate', '2-Oxobutyrate', 'L-Malate', 'DTYL', '2-Oxoglutarate',
            'D-Arabinose', 'D-Mannose', 'Salicin']

# create the json files
allConsideredSubstrates = []
allIgnoredSubstrate = []
for i, organism in enumerate(microorganisms):
    objec, considered, ignored = SBML_2_json_v2(modelName=organism, substrate_exchange_rnx='select',
                                                product_exchange_rnx=products, maxConcentration= maxConcentration,
                                                yieldTol=tolerance, saveName=saveNames[i], save=saveSwitch,
                                                toIgnore=toIgnore, alreadyConsidered=alreadyConsidered)

    allConsideredSubstrates += considered
    allIgnoredSubstrate += ignored

    alreadyConsidered = list(set(allConsideredSubstrates))
    toIgnore = list(set(allIgnoredSubstrate))

# print all the considered and ignored substrates to the terminal
print(alreadyConsidered)
print(toIgnore)

# save the input clusters to a json file as well
inputCluster = {'inputs': alreadyConsidered}
save_2_json(saveName='input_cluster.json', saveObject=inputCluster)
