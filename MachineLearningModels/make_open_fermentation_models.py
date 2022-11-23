"""
This script is intended to make surrogate models for the open fermentation model of Alberte
"""

from f_makeMLmodel import make_elasticNet_model, ridge_regression

# open fermentation only using glucose formulation 1
#excelFile =  'Glucose_open_fermentation.xlsx'
#model = ridge_regression(excelFile, save=True, saveName='Glucose_open_fermentation.json', showPLot= False)
# output = makeElasticNetModel(excelFile)

##################################################################################################
# open fermentation only using glucose formulation 2
excelFile =  'Glucose_PH_Data.xlsx'
model = ridge_regression(excelFile, normalise= False ,save=False, saveName='Glucose_pH_open_fermentation.json', showPLot= True)
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
