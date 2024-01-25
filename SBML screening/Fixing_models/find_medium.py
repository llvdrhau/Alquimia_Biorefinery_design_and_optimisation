'''
This script is used to find the medium of the models.
The medium is defined as the set of exchange reactions that are not blocked and have a negative boundry (i.e., uptake).

'''

from f_screen_SBML import find_yield, string_reactions
from f_usefull_functions import get_location
import cobra


# list of models to be analyzed
ListModels = [ 'P_acnes_V2.xml', 'PAC_4875_V2.xml', 'P_propionicum_V2.xml','P_avidum_V2.xml',  'P_sherm_V2.xml']

def find_medium(modelName):
    '''
    This function finds the medium of a model
    :param model: cobra model
    :return: list of exchange reactions that are not blocked and have a negative boundry (i.e., uptake).
    '''

    # read the model
    locationModel = get_location(modelName)
    model = cobra.io.read_sbml_model(locationModel)

    # find the exchange reactions
    listExchRxn = model.exchanges
    # find the medium
    medium = []
    lb = [] # lower bound
    ub = [] # upper bound

    for rxn in listExchRxn:
        if rxn.bounds[0] < 0:
            substrateMetabolite = rxn.reactants[0]
            medium.append(substrateMetabolite.name)
            lb.append(rxn.bounds[0])
            ub.append(rxn.bounds[1])
            # make a dictionary of the medium and the bounds
            mediumDict = dict(zip(medium, zip(lb, ub)))

    # if mediumDict does not exist, then return error
    if not 'mediumDict' in locals():
        raise Exception('No medium found. Check the bounds of the exchange reactions.')

    # print the dictionary to the console as a table

    print('The medium of the model {} is:'.format(modelName))
    print('---------------------------------------------')
    print('Metabolite \t\t\t lower bound \t\t upper bound')
    for i in mediumDict:
        print('{:<20s} {:<20s} {:<20s}'.format(str(i), str(mediumDict[i][0]), str(mediumDict[i][1])))

    return mediumDict

# run the function for all models
DictMediumDict = {}
for modelName in ListModels:
    mediumDict = find_medium(modelName=modelName)
    DictMediumDict.update({modelName: mediumDict})
    print('--------------------------------------------- \n')
    print('')

#I now want to make one table with all the unique media (i.e., uniquie metabolites) and their bounds of all the models
#I will use pandas to do this
import pandas as pd
#make a dataframe of the dictionary
df = pd.DataFrame.from_dict(DictMediumDict)
savePath = r'C:\Users\Lucas\PycharmProjects\Alquimia\SBML screening\Fixing_models\medium.xlsx'
df.to_excel(savePath)




