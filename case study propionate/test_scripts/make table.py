'''
Makes the table of the propionate yields to excel
'''

from f_usefull_functions import get_location
import json
import pandas as pd
import numpy as np

modelNames = ['v2_acnes.json',
              'v2_PAC.json',
              'v2_propionicum.json',
              'v2_avidum.json',
              'v2_sherm.json']

organismNames = ['P. acnes',
                 'P. acidipropionici',
                 'P. propionicum',
                 'P. avidum',
                 'P. freudenreichii']

substrates = ["Sucrose",
              "D-Glucose",
              "Maltose",
              "Dextrin",
              "L-Lactate",
              "Glycerol",
              "Xylose",
              "D-Fructose"]

products = ['Propionate', 'Acetate', 'Biomass']

dictionary = {}
for prod in products:
    matrixYields = np.zeros((1,8))
    for i, model in enumerate(modelNames):
        listYields = []
        modelLocation = get_location(model)
        with open(modelLocation) as file:
            mdl = json.load(file)

        for sub in substrates:
            try:
                y = mdl['coef'][prod][sub]
                y = round(y, 2)
            except:
                y = 'n.a.'
            listYields.append(y)
        arrayYields = np.array(listYields)
        arrayYields = arrayYields.reshape(1,-1)
        matrixYields = np.vstack((matrixYields,arrayYields))
    #print(matrixYields)
    matrixYields = matrixYields[1:].T

    dfYield = pd.DataFrame(matrixYields, columns= organismNames , index=substrates)
    dictionary[prod] = dfYield


excelFiles = ['output_prp.xlsx', 'output_ace.xlsx', 'output_bm.xlsx']
for i, exl in enumerate(excelFiles):
    df = dictionary[products[i]]
    df.to_excel(exl, sheet_name=products[i], index=True)

# # Create a Pandas Excel writer using openpyxl engine
# writer = pd.ExcelWriter('xxxxxxx.xlsx')
#
# # Iterate over each key in the dictionary
# for key, values in dictionary.items():
#     # Create a DataFrame using the values and substrates list
#     df = values
#     # Write the DataFrame to a separate sheet in the Excel file
#     df.to_excel(writer, sheet_name=key)
#
