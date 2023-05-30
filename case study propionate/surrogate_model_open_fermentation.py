'''
This script create the open fermentation surrogate model for the superstructure optimisation. the model is saved as a
Json file and later read by the function that make the superstructure.
'''

from f_make_surrogate_model import *

save2json = False

# Glucose_PH_Data.xlsx # without biomass as variable
# excelFile = 'Glucose_PH_Data.xlsx'
excelFile = 'Glucose_PH_60_Data_Points.xlsx'  # with a 5th polynomial it is acctually quite nice....

# ---- plot the data
dataLocation = get_location(file=excelFile, case='ML')
x_pH = pd.read_excel(dataLocation, sheet_name='inputs')
y_outputs = pd.read_excel(dataLocation, sheet_name='outputs')
# plot_data_subplots(y_data=y_outputs, x_data=x_pH)
# deleet the row where the pH is larger than 8.49
# (otherwise the fit is not going to be great, see the propionate plot)

# get the indexes of rows with values greater than 8.48
x_pH = x_pH[x_pH < 8.5]
x_pH = x_pH.dropna()  # drop the nan's
indexes = x_pH.index
y_outputs = y_outputs.loc[x_pH.index]

# plot the data to see the patterns
# plot_data_subplots(y_data=y_outputs, x_data=x_pH)

# ---- fit poly data
polynomial = 6
model = regression_open_fermentation(xdata=x_pH, ydata=y_outputs, polynomialDegree=polynomial,
                                     case='Ridge', plot=True)
# compare to other scenarios
# model = regression_open_fermentation(xdata=x_pH, ydata=y_outputs, polynomialDegree=5, case='Linear')
# model = regression_open_fermentation(xdata=x_pH, ydata=y_outputs, polynomialDegree=6, case='Linear')

# output names are the column names (could also change)
outputNames = list(y_outputs.columns)

# the inputs are the polynomials of the feature pH
featureNames = []
for i in range(polynomial):
    featureNames.append('pH**{}'.format(i))

# CV_Propionate is the name of propionate in the model!! carefull to get this right, or you'll get an error
maxProductConcentration = {'CV_Propionate': 25e-3}  # 25 g/L or 0.025 kg/L
regression_2_json_v2(outputNames, featureNames, model, maxConcentration=maxProductConcentration,
                     saveName='open_fermentation_polynomial_case_study.json', save=save2json)
