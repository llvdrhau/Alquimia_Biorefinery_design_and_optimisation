import pandas as pd

def manipulate_excel(ExcelName,feature,threshold, newExcelName):
    X = pd.read_excel(ExcelName, sheet_name='inputs')
    Y = pd.read_excel(ExcelName, sheet_name='outputs')
    inputNames = list(X.keys())
    outputNames = list(Y.keys())

    values2inspect = X[feature]
    loc = values2inspect < threshold

    helpDictIn = {}
    for i in X:
        key = i
        values = X[i][loc]
        helpDictIn.update({key:values})

    helpDictOut = {}
    for i in Y:
        key = i
        values = Y[i][loc]
        helpDictOut.update({key:values})


    inDF = pd.DataFrame(helpDictIn)
    outDF = pd.DataFrame(helpDictOut)
    with pd.ExcelWriter("{}.xlsx".format(newExcelName)) as writer:
        inDF.to_excel(writer, sheet_name='inputs',  index=False)
        outDF.to_excel(writer, sheet_name='outputs', index=False)



if __name__ == '__main__':
    ExcelName = 'Glucose_PH_200_Data_Points.xlsx'
    manipulate_excel(ExcelName =  'Glucose_PH_200_Data_Points.xlsx', feature= 'pH',threshold= 8.5,
                     newExcelName='Glucose_pH_max_range_8_5')
