"""
This script is intended to make surrogate models for the open fermentation model of Alberte
"""
import numpy as np
import matplotlib.pyplot as plt
from f_makeMLmodel import make_elasticNet_model, regression

# open fermentation only using glucose formulation 1
#excelFile =  'Glucose_open_fermentation.xlsx'
#model = ridge_regression(excelFile, save=True, saveName='Glucose_open_fermentation.json', showPLot= False)
# output = makeElasticNetModel(excelFile)

##################################################################################################
# open fermentation only using glucose formulation 2
#excelFile =  'Glucose_PH_200_Data_Points.xlsx'
excelFile = 'Glucose_pH_200_max_range_8_5.xlsx'
polynomial = 6
out = regression(excelFile, normalise= False ,save=False, saveName='Glucose_pH_open_fermentation.json',
                         showPLot= True, polynomial= {'pH':polynomial}, case= 'Ridge')

model = out[1]
modelProp = model['CV_Propionate'][0]
xpHDataProp =  model['CV_Propionate'][1]
ydataProp = model['CV_Propionate'][2]

points = 50
pH = np.zeros(shape=(points,polynomial))
ones = np.ones(shape=(points,1))
pH1 = np.linspace(4, 8.5, 50)
for i in range(polynomial):
    ph_help = pH1**(1+i)
    pH[:,i] = ph_help
ph =  np.hstack(tup=(ones,pH))

pred = modelProp.predict(ph)

plt.figure(2)
plt.plot(ph[:,1], pred)
plt.plot(xpHDataProp['pH_1'], ydataProp,'r*')
plt.show()


print('df')
# output = makeElasticNetModel(excelFile)

##################################################################################################
switchVar = False
if switchVar:
    output = make_elasticNet_model('GelatineData_elastic_net.xlsx') #, modelName= 'SAVE TEST')
    #output = makeElasticNetModel('Gelatine_Glucose_Data.xlsx')
    # print equations
    eq = output[0]
    for i in eq:
        print(i)
    # show plot of
    plt = output[1]
    plt.show()
