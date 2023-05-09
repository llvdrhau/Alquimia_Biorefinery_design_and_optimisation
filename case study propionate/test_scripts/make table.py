'''
Makes the table of the propionate yields to excel
'''

from f_usefull_functions import get_location
import json

modelNames = ['v2_acnes.json',
              'v2_PAC.json',
              'v2_propionicum.json',
              'v2_avidum.json',
              'v2_sherm.json']

lables = [ "Sucrose",
    "D-Glucose",
    "Maltose",
    "Dextrin",
    "L-Lactate",
    "Glycerol",
    "Xylose",
    "D-Fructose"]

for model in modelNames:
    modelLocation = get_location(model)
    with open(modelLocation) as file:
        mdl = json.load(file)
    print(mdl)
    yieldPropMdl = mdl['coef']['Propionate']