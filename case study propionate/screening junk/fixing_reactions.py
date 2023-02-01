
from functions import *

'''
fixing the reactions where a lot of carbon goes missing look at the Excel files generated by print_model_2_excel
to see which reactions are leaking carbon
'''

def balance(reaction):
    products = reaction.products
    reactants = reaction.reactants

    Cproducts = 0
    for prod in products:
        if 'C' in prod.elements:
            nC = prod.elements['C']
        else:
            nC = 0
        stoiFactor = reaction.metabolites[prod]
        Cproducts += nC * stoiFactor

    Creactants = 0
    for react in reactants:
        if 'C' in react.elements:
            nC = react.elements['C']
        else:
            nC = 0
        stoiFactor = reaction.metabolites[react]
        Creactants += nC * stoiFactor
    return Creactants, Cproducts

def check_reaction(model, reactionID):
    reaction = model.reactions.get_by_id(reactionID)
    reactionEq = string_reactions(reaction, case='formulas')
    reactionEqNames = string_reactions(reaction=reaction, case='names')
    print(reactionEq)
    print(reactionEqNames)
    reac, prod = balance(reaction=reaction)
    print(reac, prod)
    missing = prod + reac
    print(missing)
    return missing

# p_sherm_model
# rxnnew66_c0

modelName = "P_sherm_model.xml"
loc_sher = get_location(modelName)
model = cobra.io.read_sbml_model(loc_sher)
model.name = "P_sherm_model"
rxnID = ['rxnnew73_c0','rxnnew66_c0' ]
for i in rxnID:
    Cmissing = check_reaction(model= model, reactionID= i)

#'S_cpdnew22_c0'
#'Propionibacterium peptidoglycan'

idProt = 'S_cpd11463_c0'
#'Protein'

metProtein = model.metabolites.get_by_id(idProt)
print(metProtein.formula)

# TODO find the missing carbons of the specified metabolites in the reaction
#  with unbalanced carbon and ajust the formula with replace, make it so that O , H can also be balanced