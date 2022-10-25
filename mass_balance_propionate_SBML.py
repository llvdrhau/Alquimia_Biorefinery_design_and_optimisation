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
'S_cpdnew18_c0' , # Propionibacterium lipids',
'S_cpdnew26_c0' , #'Small molecule pool',
'S_cpdnew27_c0'] # 'Propionibacterium cell wall'

for metID in mets:
    carbonOfEachMolecule = []
    if model.metabolites.get_by_id(metID).formula:   # if it has formula count the carbons
        formula = model.metabolites.get_by_id(metID).formula
        nCarbon = countCarbonInMetabolite(formula)
        carbonOfEachMolecule.append(nCarbon)
    else:
        nCarbon = findCarbonsMissingMetabolite(model=model, metID=metID)
        if nCarbon == 0:
            # get missing metabolites and run
            subMetabolites, coef = getIDlistOfProducingMetabolties(model, metID)
            carbonSubReactions= []
            for subMetID in subMetabolites:
                nSubCarbon = findCarbonsMissingMetabolite(model=model, metID=subMetID)
                carbonSubReactions.append(nSubCarbon)
                ##okye so now i have to account for the coeficients because i'm not carring though across the functions if we go a layer deeper
                # sum up all the carbons of the metabolites from subreactions and take into account original stoichiometry


            subNCarbon = 1

# Small molecule pool
smp = mets[-2]
print(smp)
c = findCarbonsMissingMetabolite(model=model, metID=smp)
print(c)

#
# cellwallID = mets[-1]
# c = findCarbonsMissingMetabolite(model=model, metID=cellwallID)
# print(c)