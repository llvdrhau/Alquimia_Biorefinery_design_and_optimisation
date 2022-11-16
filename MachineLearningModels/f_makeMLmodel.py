# code to automatically make an elastic net regression model
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
from sklearn.model_selection import train_test_split
from sklearn.linear_model import ElasticNetCV
from sklearn.linear_model import RidgeCV
import pickle
import json

class SurrogateModel:
    def __init__(self,name, outputs, coef, intercept):
        self.name = name
        self.outputs = outputs
        self.coef = coef
        self.intercept = intercept


# info at https://scikit-learn.org/stable/modules/linear_model.html#elastic-net
# https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.ElasticNetCV.html#sklearn.linear_model.ElasticNetCV
def make_elasticNet_model(ExcelName, iterations = 1000, modelName = '', plotSwitch = 0, case = 'ridge'):
    X = pd.read_excel(ExcelName,sheet_name='inputs')
    inputnames = X.keys()
    y = pd.read_excel(ExcelName,sheet_name='outputs')
    rows = 2
    cols = math.ceil(y.shape[1]/2)
    stringEquationVector= []
    for id,i in enumerate(y):
        # might need to normalise data....
        # mean standard normilsation if inputs are not in the same decimal rank?
        # select output
        Y = y[i]
        # split training data/ test data
        X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.2,random_state = 2) # random_state = 2, so consistent results are obtained (2 being the seed)

        ######## find correct hyper parameters: alfa (regularisation parameter)
        # and l1_ratio => how is the alfa parameter devided over the L1 and L2 norm
        if case == 'elastic':
            ratios = np.arange(0, 1.1, 0.1)
            # ratios = [0, .1, .5, .7, .9, .95, .99, 1]
        elif case == 'lasso':
            ratios = 1
        else: # if case is ridge regression
            ratios = 0

        alphas = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.0, 1.0, 10.0, 100.0]
        alphas = [1e-1, 0.0, 1.0, 10.0, 100.0]
        model = ElasticNetCV(l1_ratio=ratios, alphas=alphas, normalize = True , max_iter = iterations,n_jobs=-1)

        # fit model
        model.fit(X_train, y_train)
        #a = model.alpha_  #check() out hyperparameters
        #b = model.l1_ratio_ #check() out hyperparameters
        # summarize chosen configuration
        print('alpha: %f' % model.alpha_)
        print('l1_ratio_: %f' % model.l1_ratio_)

        # Make predictions using the testing set
        y_pred_en = model.predict(X_test)

        # subplot to evaluate goodness of fit
        ax = plt.subplot(rows, cols, id+1)
        ax.plot(y_test,y_pred_en, 'k*')
        ax.set_title(i)
        ax.set_xlabel("")

        # add equation to the vector of strings => 'y = ax +b'
        yName = i
        eq = yName + ' == '
        for j,xname in enumerate(inputnames):
            eq = eq + '+' + xname + '*{0}'.format(model.coef_[j])

        stringEquationVector.append(eq)

    if plotSwitch:
        plt.show()

    if modelName:
        with open("stringEquationVector.bin", "wb") as modelName:  # "wb" because we want to write in binary mode
            pickle.dump(stringEquationVector, modelName)

        # if we want to relaod in other scripts:
        # with open("state.bin", "rb") as f:  # "rb" because we want to read in binary mode
        #     state = pickle.load(f)

    return stringEquationVector,plt

def ridge_regression(ExcelName, showPLot = True, save = False, saveName = 'data.json'):

    # features
    X = pd.read_excel(ExcelName, sheet_name='inputs')
    # normalise
    for input in X:
        meanInput = X[input].mean()
        stdInput = X[input].std()
        inputNormalised = (X[input] - meanInput) / stdInput
        X[input] = inputNormalised
    inputNames = list(X.keys())
    # target values (the reactor outputs)
    Y = pd.read_excel(ExcelName, sheet_name='outputs')
    outputNames =  list(Y.keys())
    # normalise
    for output in Y:
        meanInput = Y[output].mean()
        stdInput = Y[output].std()
        outputNormalised = (Y[output] - meanInput) / stdInput
        Y[output] = outputNormalised

    # plot variables
    rows = 2
    cols = math.ceil(Y.shape[1] / 2)

    # ridge regression for each output
    outputVariables = []
    coefficients = {}
    intercepts = {}
    coefArry = np.zeros(shape=(len(Y.keys()), len(X.keys()))) # columns are the inputs (= features) rows the output variables
    for i,outName in enumerate(Y):
        # select output
        output = Y[outName]
        # split training data/ test data
        X_train, X_test, y_train, y_test = train_test_split(X, output, test_size=0.2,random_state = 2) # random_state = 2, so consistent results are obtained (2 being the seed)
        # define model
        model = RidgeCV(alphas=(1e-2, 1e-1, 1.0, 10.0, 100.0))
        # fit model
        model.fit(X_train, y_train)

        ##### check if the fit is OK
        # Make predictions using the testing set
        y_pred_en = model.predict(X_test)
        # subplot to evaluate goodness of fit
        ax = plt.subplot(rows, cols, i + 1)
        ax.plot(y_test, y_pred_en, 'k*')
        xlims = list(ax.get_xlim())
        ylims = list(ax.get_ylim())
        ax.plot(xlims, ylims, 'r')
        ax.set_title(outName)
        ax.set_xlabel("real")
        ax.set_ylabel("predicted")
        # add equation to the vector of strings => 'y = ax + b'
        yName = outName
        outputVariables.append(yName)
        eq = yName + ' == '
        for j, xname in enumerate(inputNames):
            eq = eq + ' + ' + xname + ' * {0}'.format(model.coef_[j])
            coefficients.update()
            coefArry[i,j] = model.coef_[j]
        eq = eq + ' + {}'.format(model.intercept_)
        intercepts.update({outName: model.intercept_})
        print(eq)

    if showPLot:
        plt.show()

    coefDict = {}
    for i, out in enumerate(outputNames):
        featureCoefDict = {}
        for j, feature in enumerate(inputNames):
            coefOfFeature = coefArry[i,j]
            featureCoefDict.update({feature:coefOfFeature}) # can't put a np array in a json file
        coefDict.update({out:featureCoefDict})
    #coefficients.update({'coefficients':coefDict})

    surrogateModel = SurrogateModel(name=saveName ,outputs= outputVariables, coef=coefDict, intercept=intercepts)
    if save:
        loc = os.getcwd()
        posAlquimia = loc.find('Alquimia')
        loc = loc[0:posAlquimia + 8]
        loc = loc + r'\json models' + r'\{}'.format(saveName)
        with open(loc, 'w+', encoding='utf-8') as f:
            json.dump(surrogateModel.__dict__, f, ensure_ascii=False, indent=4)

        #with open("/path/to/file.json", "w+") as f:
        #    json.dump(object_to_write, f)

    return surrogateModel





