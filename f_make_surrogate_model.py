"""
functions to find surrogate models with machinelearning models and SBML models
the models are transformed into a JSON  file so it can be read by the superstructure functions
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from f_find_carbons import carbon_balance_in_out
import cobra.io
import re

from sklearn.model_selection import train_test_split
from sklearn.linear_model import ElasticNetCV
from sklearn.linear_model import RidgeCV
from sklearn.linear_model import LassoCV
from sklearn.linear_model import Lasso
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures

import json
import math

def find_carbons_in_formula(formula):
    metFormula = formula
    splitFormula = re.split('(\d+)', metFormula)
    nrOfCarbons = 0  # just in case something wierd happens
    if 'C' not in metFormula:  # if there is no carbon in the formula
        nrOfCarbons = 0
    else:
        for j, element in enumerate(splitFormula):
            if 'C' in element and len(element) == 1:
                nrOfCarbons = int(splitFormula[j + 1])
            elif 'C' in element and len(element) > 1:
                posCarbon = element.index('C')
                if element[posCarbon + 1].isupper():  # for case like CH4 there is no 1 next to the C
                    nrOfCarbons = 1  # only one carbon
                else:
                    continue  # for cases like Co (cobalt) just skip
            else:
                continue
    return nrOfCarbons


class SurrogateModel:
    def __init__(self,name, outputs, coef, lable, intercept= None):
        self.name = name
        self.outputs = outputs
        self.coef = coef
        if intercept is None:
            interceptDict = {}
            for i in outputs:
                interceptDict.update({i:''})
            self.intercept = interceptDict
        else:
            self.intercept = intercept
        self.lable = lable


def regression_2_json(ExcelName, showPLot = True, save = False, saveName = 'data.json', normalise = False,
                     case = 'Ridge',polynomial=None):
    loc = os.getcwd()
    posAlquimia = loc.find('Alquimia')
    loc = loc[0:posAlquimia + 8]
    # modelLocations = loc + r'\excel files\' + modelName
    dataLocation = loc + r"\machine learning models\excel_data\{}".format(ExcelName)

    if polynomial is None:
        polynomial = {}

    # features
    X = pd.read_excel(dataLocation, sheet_name='inputs')
    inputNames = list(X.keys())

    for name in inputNames:
        polynomialKeys = list(polynomial.keys())
        if len(polynomialKeys) > 1:
            raise Exception('so something to fix, you can not have more then two variables makes it messy fix the code '
                            'if you have time, if not put the features als polynomials in the excel file xoxox Lucas of the past')
        if name in polynomialKeys:
            x = X[name].to_numpy()
            nPolynoms = polynomial[name]
            X_new = PolynomialFeatures(nPolynoms).fit_transform(x.reshape((len(x), 1)))
            #X_new = X_new[:,1:len(X_new)]
            dict2Pandas = {}
            for nr, col in enumerate(X_new.T): # don't forget to transpose the matrix to loop over the cols
                key = name + '**{}'.format(nr)
                dict2Pandas.update({key:col})
            X = pd.DataFrame(dict2Pandas)
            #X = X_new # todo shiiiit what if more then one variable!!


    # target values (the reactor outputs)
    Y = pd.read_excel(dataLocation, sheet_name='outputs')
    outputNames = list(Y.keys())

    if normalise:
        # normalise
        for input in X:
            meanInput = X[input].mean()
            stdInput = X[input].std()
            inputNormalised = (X[input] - meanInput) / stdInput
            X[input] = inputNormalised

    if polynomial and normalise:
        inputNames = list(X.keys())
        dropName = inputNames[0] # the bais is what you want to make into ones
        X[dropName] = np.zeros(shape=(len(X),1)) # bais should be ones not nan if normalised


        # normalise
        # for output in Y:
        #     meanInput = Y[output].mean()
        #     stdInput = Y[output].std()
        #     outputNormalised = (Y[output] - meanInput) / stdInput
        #     Y[output] = outputNormalised

    # plot variables
    rows = 2
    cols = math.ceil(Y.shape[1] / 2)

    # ridge regression for each output
    outputVariables = []
    coefficients = {}
    intercepts = {}
    coefMatrix = np.zeros(shape=(len(Y.keys()), len(X.keys()))) # columns are the inputs (= features) rows the output variables
    modelDictionary = {}
    jsonDict = {}

    for i,outName in enumerate(Y):
        # select output
        output = Y[outName]
        # split training data/ test data
        X_train, X_test, y_train, y_test = train_test_split(X, output, test_size=0.2,random_state = 2) # random_state = 2, so consistent results are obtained (2 being the seed)
        # define model
        alfas_ridge = (1e-2, 1e-1, 1.0, 10.0, 100.0)
        alfas_lasso = (1e-12, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 500, 800,1000)
        if case == 'Ridge':
            model = RidgeCV(alphas= alfas_ridge)
        elif case == 'Lasso':
            model = LassoCV(alphas= alfas_lasso, max_iter=10000)
            #model = Lasso(alpha= 1, max_iter= 4000)
        elif case == 'Linear':
            model = LinearRegression()
        else:
            raise Exception("The string variable _case_ can only be 'Linear, 'Lasso' or 'Ridge'")
        #model = Ridge(alpha=1)
        # fit model
        model.fit(X_train, y_train)

        ##### check if the fit is OK
        # Make predictions using the testing set
        y_predicted =  model.predict(X)
        y_observed = output

        # subplot to evaluate goodness of fit
        ax = plt.subplot(rows, cols, i + 1)
        ax.plot(y_observed, y_predicted, 'k*')
        minimum = min([min(y_observed),min(y_predicted)])
        maximum = max([max(y_observed),max(y_predicted)])
        ax.plot([minimum, maximum], [minimum, maximum], 'r') # plot diagonal
        ax.set_title(outName)
        ax.set_xlabel("real")
        ax.set_ylabel("predicted")
        # add equation to the vector of strings => 'y = ax + b'
        yName = outName
        outputVariables.append(yName)
        eq = yName + ' == '
        for j, xname in enumerate(X):
            eq = eq + ' + ' + xname + ' * {0}'.format(model.coef_[j])
            coefficients.update()
            coefMatrix[i,j] = model.coef_[j]
        eq = eq + ' + {}'.format(model.intercept_)

        intercepts.update({outName: model.intercept_})

        if case == 'Lasso'  or case == 'Ridge':
            print(model.alpha_)
        print('the model coef are {}'.format(model.coef_))
        print('the model intercept is {}'.format(model.intercept_))
        normFactor = 1/(max(y_predicted)- min(y_predicted))
        NMSE = math.sqrt(sum((y_observed - y_predicted) ** 2) / len(y_observed)) * normFactor
        print('the NMSE is: {}'.format(NMSE))
        print(eq)
        modelDictionary.update({outName:(model,X,output)})  # X are the inputs

    if showPLot:
        plt.show()

    coefDict = {}
    for i, out in enumerate(outputNames):
        featureCoefDict = {}
        for j, feature in enumerate(X):
            coefOfFeature = coefMatrix[i,j]
            featureCoefDict.update({feature:coefOfFeature}) # can't put a np array in a json file
        coefDict.update({out:featureCoefDict})
    #coefficients.update({'coefficients':coefDict})
    jsonDict.update({'lable': case,
                     'Name': saveName,
                     'CV_Equations': {'coef':coefDict, 'intercept': model.intercept_}})

    surrogateModel = SurrogateModel(name=saveName ,outputs= outputVariables, coef=coefDict, intercept=intercepts,
                                    lable = 'other')
    if save:
        loc = os.getcwd()
        posAlquimia = loc.find('Alquimia')
        loc = loc[0:posAlquimia + 8]
        loc = loc + r'\json models' + r'\{}'.format(saveName)
        with open(loc, 'w+', encoding='utf-8') as f:
            json.dump(surrogateModel.__dict__, f, ensure_ascii=False, indent=4)

        #with open("/path/to/file.json", "w+") as f:
        #    json.dump(object_to_write, f)


    return surrogateModel,modelDictionary

def SBML_2_json(modelName, substrate_exchange_rnx, product_exchange_rnx, newObjectiveReaction = None, saveName = None,
                substrate2zero= 'Ex_S_cpd00027_ext', missingCarbonId = None, printEq = False, checkCarbon = True, save = False):
    loc = os.getcwd()
    posAlquimia = loc.find('Alquimia')
    loc = loc[0:posAlquimia + 8]
    # modelLocations = loc + r'\excel files\' + modelName
    modelLocation = loc + r"\SBML models\{}".format(modelName)

    if saveName is None:
        saveName = modelName.replace('.xml','.json')

    allYields_FBA =[]
    allEquations = []

    model = cobra.io.read_sbml_model(modelLocation)
    # make sure the right objective is set
    if newObjectiveReaction:
        model.objective = newObjectiveReaction

    # change the standard exchange reaction to zero
    exchange_rnx_2_zero = substrate2zero
    model.reactions.get_by_id(exchange_rnx_2_zero).bounds = 0, 0

    coefDict = {}
    outputVariables = []
    for product in product_exchange_rnx:
        productMet = model.reactions.get_by_id(product).reactants[0]
        productName = productMet.name
        productFormula = productMet.formula
        Cprod = find_carbons_in_formula(productFormula)
        strEq = '{} == '.format(productName)
        outputVariables.append(productName)
        substrateCoefDict = {}
        for substrate in substrate_exchange_rnx:
            # set all substrates to zero
            for exchRnx in substrate_exchange_rnx:
                model.reactions.get_by_id(exchRnx).bounds = 0, 0
            # change bound of new desired substrate to -10 mol/h/gDW
            model.reactions.get_by_id(substrate).bounds = -10, 0

            # do regular FBA
            solution = model.optimize()
            FBA_substrate_flux = solution.fluxes[substrate]
            substrateMet = model.reactions.get_by_id(substrate).reactants[0]
            substrateName = substrateMet.name
            substrateFormula = substrateMet.formula
            Csub = find_carbons_in_formula(substrateFormula)

            FBA_product_flux = solution.fluxes[product]
            FBA_yield = abs((FBA_product_flux / FBA_substrate_flux) * (Cprod * 12) / (Csub * 12))  # in gramsC / grams C: 12 gCarbon/mol
            allYields_FBA.append(FBA_yield)
            strEq += ' + {} * {}'.format(FBA_yield, substrateName)
            substrateCoefDict.update({substrateName:FBA_yield})
        coefDict.update({productName:substrateCoefDict})
        allEquations.append(strEq)
        if printEq:
            print(strEq)

    if printEq:
        print(modelName)

    if checkCarbon:
        carbon_balance_in_out(modelLocation=model, metIDsMissingCarbon=missingCarbonId, tol=1e-4)

    surrogateModel = SurrogateModel(name=modelName, outputs=outputVariables, coef=coefDict,
                                    lable='SBML')

    if save:
        loc = os.getcwd()
        posAlquimia = loc.find('Alquimia')
        loc = loc[0:posAlquimia + 8]
        loc = loc + r'\json models' + r'\{}'.format(saveName)
        with open(loc, 'w+', encoding='utf-8') as f:
            json.dump(surrogateModel.__dict__, f, ensure_ascii=False, indent=4)

    return allEquations, allYields_FBA