from functions import *

excelFile = r'\propionate_case_study.xlsx'

excelDict = read_excel_sheets4_superstructure(excelName=excelFile )

objectsOutputDict1  = make_output_intervals(excelDict)
objectsOutputDict2  = make_output_intervals(excelDict)

if objectsOutputDict1 == objectsOutputDict2:
    print('coooool')
else:
    print('fuuccccck')