"""
Octobre 2022
This script is to check the carbon balance of the Propioni bacteria
"""

#import modules
import os
from f_find_carbons import *

# get the location of the models
loc = os.getcwd()
loc_acidi = loc + r'\SBML models\PAC_4875_model.xml'
loc_acnes = loc + r'\SBML models\P_acnes_model.xml'
loc_prop = loc + r'\SBML models\P_propionicum_model.xml'
loc_avidum = loc + r'\SBML models\P_avidum_model.xml'
loc_sher = loc + r'\SBML models\P_sherm_model.xml'

microorganisms = [loc_acidi, loc_acnes, loc_prop,loc_avidum,loc_sher]
objectiveMetID = 'S_biomass_ext'

for i in microorganisms:
    carbonBalanceInOut(modelLocation=i, metIDsMissingCarbon=objectiveMetID, tol=1e-4)

print('')