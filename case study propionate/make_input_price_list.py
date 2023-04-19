from f_usefull_functions import get_location, save_2_json
import json
import pandas as pd

jsonFile = 'input_cluster.json'
jsonLoc = get_location(file=jsonFile)
with open(jsonLoc) as file:
    inputs = json.load(file)


# Prices searched on sigma ALdich (see Excel file, sheet price_list)
Excelfile = 'propionate_case_study_v2.xlsx'
locationExcel = get_location(Excelfile)
priceDF = pd.read_excel(locationExcel, sheet_name= 'price_list', index_col= 'inputs')


inputList = inputs['inputs']
dictInputs = {}
for i in inputList:
    price = round(priceDF.price[i],2)
    dictInputs.update({i: price})

save_2_json(saveName='inputs_v2.json', saveObject= dictInputs)
