"""
Octobre 2022
This script is to check the carbon balance of the Propioni bacteria
"""

# import modules
from f_find_carbons import *

# get the location of the models
loc = os.getcwd()
loc_acidi = loc + r'\SBML models\PAC_4875_model.xml'
loc_acnes = loc + r'\SBML models\P_acnes_model.xml'
loc_prop = loc + r'\SBML models\P_propionicum_model.xml'
loc_avidum = loc + r'\SBML models\P_avidum_model.xml'
loc_sher = loc + r'\SBML models\P_sherm_model.xml'

microorganisms = [loc_acidi, loc_acnes, loc_prop, loc_avidum, loc_sher]
microorganisms = []
objectiveMetID = 'S_biomass_ext'

for i in microorganisms:
    carbonBalanceInOut(modelLocation=i, metIDsMissingCarbon=objectiveMetID, tol=1e-4)

#####################################################################################
#test other micro organisms
microViviane = loc + r'\SBML models\p-thermo.xml'  #iNF517.xml' # latcate producing organism
validat = cobra.io.sbml.validate_sbml_model(microViviane)
print(validat[0])
print('')
print(validat[1]['SBML_FATAL']) #
#model = cobra.io.read_sbml_model(microViviane)
carbonBalanceInOut(modelLocation=microViviane, metIDsMissingCarbon=[], tol=1e-4)

"""
ok so we can see the following, for the P. species have a carbon balance which is 60% to 70% compleet
so compare to other species  P. thermoglucosidasius (Viviane et al 2020) and lactococus (from BIGG) 70 and 80% compleet 
respectively so where does the mas sin those models go?? Do internal metabolites consitute the rest of the metabolites?
wat is the deal man?? 
"""