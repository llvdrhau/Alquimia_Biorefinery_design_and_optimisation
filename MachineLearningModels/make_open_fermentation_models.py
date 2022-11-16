"""
This script is intended to make surrogate models for the open fermentation model of Alberte
"""

from f_makeMLmodel import make_elasticNet_model, ridge_regression

# open fermentation only using glucose
excelFile =  'Glucose_open_fermentation.xlsx'
model = ridge_regression(excelFile, save=True,saveName='Glucose_open_fermentation.json' ,showPLot= False)
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
