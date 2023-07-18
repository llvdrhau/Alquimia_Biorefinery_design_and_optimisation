"""
This script creates the surrogate models of the GEMs. Json files are made from which the superstructure can then read 
the appropriate mathematical expression 

In this script you will be asked to include/exclude substrates. Based on the history of your choice the same substrates
will be included or excluded during the loop to minimise time selecting substrates manually  
"""

from f_make_surrogate_model import *

# ----- script specifications
saveSwitch = True

# -------------------------- SBML models
loc_acidi = 'PAC_4875_V2.xml'
loc_acnes = 'P_acnes_V2.xml'
loc_prop = 'P_propionicum_V2.xml'
loc_avidum = 'P_avidum_V2.xml'
loc_sher = 'P_sherm_V2.xml'

# select substrates qnd products
substrates = 'select'               # substrates = ['Ex_S_cpd00027_ext']  # glucose
products = ['Ex_S_cpd00141_ext',    # propionate
            'Ex_S_cpd00029_ext',    # acetate
            'Ex_S_biomass_ext']     # biomass

# propionate needs a minimum yield of 0.1 for the possible substrate selection procedure
tolerance = {'Ex_S_cpd00141_ext': 0.15,   # tolerance for the yield of propionate
             'Ex_S_cpd00029_ext': 0, # tolerance for the yield of acetate -1e-10
             'Ex_S_biomass_ext': 0.05 }  # tolerance for the yield of biomass

# define maximum allowed concentrations in reactor
maxConcentration = {'Propionate': 42.37e-3}  # max concentration of propionate is 0.035 kg/L or 35 g/L

# give the names of the organisms and the save names of the json files
microorganisms = [loc_acidi, loc_acnes, loc_prop, loc_avidum, loc_sher]  # all microorganisms
saveNames = ['v2_PAC.json', 'v2_acnes.json', 'v2_propionicum.json', 'v2_avidum.json', 'v2_sherm.json']

# give list of components to already ignore or include
alreadyConsidered = ['D-Glucose', 'D-Fructose', 'L-Lactate', 'Glycerol', 'Maltose', 'Sucrose', 'Xylose', 'Dextrin']
toIgnore = ['D-Lyxitol', 'alpha-Methyl-D-glucoside', 'CELB', 'TRHL', 'GLCN', 'Succinate', 'Glycerol-3-phosphate',
            'beta-Methylglucoside', 'Glycerone', 'Pyruvate', 'LACT', 'D-Mannitol', 'ACTN', 'L-Inositol', 'Ribitol',
            'Thyminose', 'D-Ribose', 'Glucuronate', 'Fumarate', '2-Oxobutyrate', 'L-Malate', 'DTYL', '2-Oxoglutarate',
            'D-Arabinose', 'D-Mannose', 'Salicin', 'L-Lyxitol', 'Xylitol', 'D-GLUCOSE-6-PHOSPHATE', 'GLUCOSE-1-PHOSPHATE']

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
