from f_find_carbons import *

# load model acidi bacteria
loc = os.getcwd()
loc_acidi = loc + r'\SBML models\PAC_4875_model.xml'
model = cobra.io.read_sbml_model(loc_acidi)

# metaboliets that make the biomass model.
# model.metabolites.get_by_id('S_biomass_ext').summary()
mets = ['S_cpd11461_c0', # 'DNA', '
'S_cpd11463_c0' , # Protein'
'S_cpd11613_c0' , # , 'RNA'
'S_cpdnew26_c0' , #'Small molecule pool',
'S_cpdnew18_c0' , # Propionibacterium lipids',
'S_cpdnew27_c0'] # 'Propionibacterium cell wall'

objectiveReactionID = 'biomass_c0' #
c = findCarbonsOfReaction(model= model ,reactionID= objectiveReactionID)

print(c)