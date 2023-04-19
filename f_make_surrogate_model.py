import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import cobra.io

import seaborn as sns

from sklearn.model_selection import train_test_split
from sklearn.linear_model import RidgeCV, LassoCV, Ridge, Lasso, LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.metrics import r2_score, mean_squared_error

import json
import math
from f_usefull_functions import get_location, save_2_json
from f_screen_SBML import count_atom_in_formula, carbon_balance_in_out, find_yield, is_protein_met

"""
functions to find surrogate models with machinelearning models and SBML models
the models are transformed into a JSON  file so it can be read by the superstructure functions
"""


# --------------------------------------------------------------------------------------
# Surogate model class
# ------------------------------------------- -------------------------------------------
# Class to make all surrogate models uniform!
class SurrogateModel:
    def __init__(self, name, inputs, outputs, coef, lable,lightKey= None, maxConcentration=None, intercept=None):
        self.name = name
        self.inputs = inputs
        self.outputs = outputs
        self.coef = coef
        if intercept is None:
            interceptDict = {}
            for i in outputs:
                interceptDict.update({i: ''})
            self.intercept = interceptDict
        else:
            self.intercept = intercept
        self.lable = lable

        if maxConcentration is None:
            maxConcentration = {}

        if lightKey is not None:
            self.lightKey = lightKey


        if maxConcentration:
            key = list(maxConcentration.keys())[0]
            val = list(maxConcentration.values())[0]
            self.waterEq = 'water == {} / {}'.format(key, val)

            # for key, val in maxConcentration.items():


# # ------------------------------------------- -------------------------------------------
# Surogate model functions open_fermentation unit
# # ------------------------------------------- -------------------------------------------
def regression_open_fermentation(xdata, ydata, polynomialDegree, case='Lasso', plot=True):
    # make the polynomial data
    poly = PolynomialFeatures(degree=polynomialDegree, include_bias=True)
    X_poly = poly.fit_transform(xdata)

    # Split data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X_poly, ydata, test_size=0.35, random_state=42)

    # Fit linear regression model to training data
    if case == 'Ridge':
        reg = Ridge().fit(X_train, y_train)
    elif case == 'Lasso':
        reg = Lasso(alpha=0.002).fit(X_train, y_train)
        # model = Lasso(alpha= 1, max_iter= 4000)
    elif case == 'Linear':
        reg = LinearRegression().fit(X_train, y_train)
    else:
        raise Exception("The string variable _case_ can only be 'Linear, 'Lasso' or 'Ridge'")

    # Predict Qtot for test data
    y_pred = reg.predict(X_test)
    r2_scores = r2_score(y_test, y_pred, multioutput='raw_values')
    MSE = mean_squared_error(y_test, y_pred, multioutput='raw_values')
    coefs = reg.coef_

    # print info
    print("R2 score for each output variable:", r2_scores)
    print("mean_squared_error for each output variable:", MSE)
    print('')
    print("the ceof are:", coefs)

    # only dependent on 1 variable, so we can plot the outcome of the model vs the pH
    # create data to plot the model
    colName = xdata.columns[0]
    x_data_model = np.linspace(min(xdata[colName]), max(xdata[colName]), num=100)
    x_data_model = x_data_model.reshape((len(x_data_model), 1))  # reshape
    x_poly_data_model = poly.fit_transform(x_data_model)
    y_data_model = reg.predict(x_poly_data_model)

    # check out the  plots
    if plot:
        plot_parity_plots(yPred=y_pred, yObv=y_test)
        plot_model_vs_data(x_data=X_poly[:, 1], y_data=ydata, x_data_model=x_data_model.squeeze(),
                           y_data_model=y_data_model)
    return reg


def plot_data_subplots(x_data, y_data):
    """ plots the data of the model we want to regress """

    num_cols_x = x_data.shape[1]
    if num_cols_x > 1:
        raise Exception('The number of inputs should not be larger then 1')

    num_cols = y_data.shape[1]  # number of columns in y_data
    num_rows = (num_cols - 1) // 2 + 1  # calculate number of rows for subplot layout
    fig, axes = plt.subplots(nrows=num_rows, ncols=2, figsize=(10, 5 * num_rows))  # create subplots
    for i, ax in enumerate(axes.flatten()):  # iterate over subplots
        if i < num_cols:  # plot data if there are still columns left
            sns.set_style('darkgrid')
            sns.scatterplot(x=x_data.squeeze(), y=y_data.iloc[:, i],
                            ax=ax)  # plot i-th column of y_data against x_data using seaborn
            ax.set_title(y_data.columns[i])  # set title to column name
        else:  # remove unused subplots
            ax.remove()
    fig.tight_layout()  # adjust subplot spacing
    plt.show()  # display plot


def plot_parity_plots(yPred, yObv):
    num_cols_yPred = yPred.shape[1]
    num_cols_yObv = yObv.shape[1]
    try:
        subplotTitles = list(yObv.columns)
    except:
        subplotTitles = ['unspecified'] * len(yObv)
    if num_cols_yPred != num_cols_yObv:
        raise Exception('yPred and yObv should be the same size')

    num_cols = yPred.shape[1]  # number of columns in y_data
    num_rows = (num_cols - 1) // 2 + 1  # calculate number of rows for subplot layout
    fig, axes = plt.subplots(nrows=num_rows, ncols=2, figsize=(10, 5 * num_rows))  # create subplots
    sns.set_style('darkgrid')
    for i, ax in enumerate(axes.flatten()):  # iterate over subplots
        if i < num_cols:  # plot data if there are still columns left
            observed = yObv[subplotTitles[i]].to_numpy()
            predicted = yPred[:, i]
            sns.scatterplot(x=observed, y=predicted, ax=ax)  # plot i-th column of y_data against x_data using seaborn
            sns.lineplot(x=predicted, y=predicted, color='red', ax=ax)
            ax.set_xlabel('Observed Qtot')
            ax.set_ylabel('Predicted Qtot')
            ax.set_title(f'Parity plot for {subplotTitles[i]}')  # set title to column name
        else:  # remove unused subplots
            fig.delaxes(ax)
    fig.tight_layout()  # adjust subplot spacing
    plt.show()  # display plot


def plot_model_vs_data(x_data, y_data, x_data_model, y_data_model):
    """
    plots the data of the regression model and the data it was trained on
    """
    try:
        y_data = y_data.to_numpy() # change to numpy array if it is a dataframe
    except:
        pass

    num_cols = y_data.shape[1]  # number of columns in y_data
    num_rows = (num_cols - 1) // 2 + 1  # calculate number of rows for subplot layout
    fig, axes = plt.subplots(nrows=num_rows, ncols=2, figsize=(10, 5 * num_rows))  # create subplots
    for i, ax in enumerate(axes.flatten()):  # iterate over subplots
        if i < num_cols:  # plot data if there are still columns left
            sns.set_style('darkgrid')
            sns.scatterplot(x=x_data, y=y_data[:, i], ax=ax)  # plot i-th column of y_data against x_data using seaborn
            sns.lineplot(x=x_data_model, y=y_data_model[:, i], ax=ax)
            #ax.set_title(y_data.columns[i])  # set title to column name
        else:  # remove unused subplots
            ax.remove()
    fig.tight_layout()  # adjust subplot spacing
    plt.show()  # display plot


def regression_2_json(data, showPLot=True, save=False, saveName='data.json', normalise=False,
                      case='Ridge', polynomial=None):
    if isinstance(data, str) and ':xlsx' in data:
        pass
    dataLocation = get_location(file=data, case='ML')

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
            # X_new = X_new[:,1:len(X_new)]
            dict2Pandas = {}
            for nr, col in enumerate(X_new.T):  # don't forget to transpose the matrix to loop over the cols
                key = name + '**{}'.format(nr)
                dict2Pandas.update({key: col})
            X = pd.DataFrame(dict2Pandas)

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
        dropName = inputNames[0]  # the bais is what you want to make into ones
        X[dropName] = np.zeros(shape=(len(X), 1))  # bais should be ones not nan if normalised

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
    coefMatrix = np.zeros(
        shape=(len(Y.keys()), len(X.keys())))  # columns are the inputs (= features) rows the output variables
    modelDictionary = {}
    jsonDict = {}

    for i, outName in enumerate(Y):
        # select output
        output = Y[outName]
        # split training data/ test data
        X_train, X_test, y_train, y_test = train_test_split(X, output, test_size=0.2,
                                                            random_state=2)  # random_state = 2, so consistent results are obtained (2 being the seed)
        # define model
        alfas_ridge = (1e-2, 1e-1, 1.0, 10.0, 100.0)
        alfas_lasso = (1e-12, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 500, 800, 1000)
        if case == 'Ridge':
            model = RidgeCV(alphas=alfas_ridge)
        elif case == 'Lasso':
            model = LassoCV(alphas=alfas_lasso, max_iter=10000)
            # model = Lasso(alpha= 1, max_iter= 4000)
        elif case == 'Linear':
            model = LinearRegression()
        else:
            raise Exception("The string variable _case_ can only be 'Linear, 'Lasso' or 'Ridge'")
        # model = Ridge(alpha=1)
        # fit model
        model.fit(X_train, y_train)

        ##### check if the fit is OK
        # Make predictions using the testing set
        y_predicted = model.predict(X)
        y_observed = output

        # subplot to evaluate goodness of fit
        ax = plt.subplot(rows, cols, i + 1)
        ax.plot(y_observed, y_predicted, 'k*')
        minimum = min([min(y_observed), min(y_predicted)])
        maximum = max([max(y_observed), max(y_predicted)])
        ax.plot([minimum, maximum], [minimum, maximum], 'r')  # plot diagonal
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
            coefMatrix[i, j] = model.coef_[j]
        eq = eq + ' + {}'.format(model.intercept_)

        intercepts.update({outName: model.intercept_})

        if case == 'Lasso' or case == 'Ridge':
            print(model.alpha_)
        print('the model coef are {}'.format(model.coef_))
        print('the model intercept is {}'.format(model.intercept_))
        normFactor = 1 / (max(y_predicted) - min(y_predicted))
        NMSE = math.sqrt(sum((y_observed - y_predicted) ** 2) / len(y_observed)) * normFactor
        print('the NMSE is: {}'.format(NMSE))
        print(eq)
        modelDictionary.update({outName: (model, X, output)})  # X are the inputs

    if showPLot:
        plt.show()

    coefDict = {}
    for i, out in enumerate(outputNames):
        featureCoefDict = {}
        for j, feature in enumerate(X):
            coefOfFeature = coefMatrix[i, j]
            featureCoefDict.update({feature: coefOfFeature})  # can't put a np array in a json file
        coefDict.update({out: featureCoefDict})
    # coefficients.update({'coefficients':coefDict})
    jsonDict.update({'lable': case,
                     'Name': saveName,
                     'CV_Equations': {'coef': coefDict, 'intercept': model.intercept_}})

    surrogateModel = SurrogateModel(name=saveName, outputs=outputVariables, coef=coefDict, intercept=intercepts,
                                    lable='yield_equation')
    if save:
        loc = os.getcwd()
        posAlquimia = loc.find('Alquimia')
        loc = loc[0:posAlquimia + 8]
        loc = loc + r'\json models' + r'\{}'.format(saveName)
        with open(loc, 'w+', encoding='utf-8') as f:
            json.dump(surrogateModel.__dict__, f, ensure_ascii=False, indent=4)

        # with open("/path/to/file.json", "w+") as f:
        #    json.dump(object_to_write, f)

    return surrogateModel, modelDictionary


def regression_2_json_v2(outputNames, featureNames, model, saveName, inputNames = None, save=True,
                         maxConcentration=None, lightKey = None, lable = 'yield_equation'):
    """ Saves the regresion model as a readable json file for the superstructure
    Params

    Returns

    """
    if isinstance(outputNames, str):
        outputNames = [outputNames] # make sure it's a list so you interrate correctly

    if inputNames is None:
        inputNames = []
    if isinstance(inputNames, str):
        inputNames = [inputNames] # make sure it's a list

    coef_ = model.coef_
    interpect_ = model.intercept_
    coefDict = {}
    interpectDict = {}
    for i, out in enumerate(outputNames):
        featureCoefDict = {}
        interpectDict.update({out: interpect_[i]})
        for j, feature in enumerate(featureNames):
            coefOfFeature = coef_[i, j]
            featureCoefDict.update({feature: coefOfFeature})  # can't put a np array in a json file
        coefDict.update({out: featureCoefDict})

    surrogateModel = SurrogateModel(name=saveName, inputs = inputNames, outputs=outputNames, coef=coefDict, intercept=interpectDict,
                                    lable=lable, maxConcentration=maxConcentration, lightKey=lightKey)
    if save:
        save_2_json(saveName=saveName, saveObject=surrogateModel)


# # ------------------------------------------- -------------------------------------------
# Surogate model functions for GEMS
# # ------------------------------------------- -------------------------------------------

def SBML_2_json(modelName, substrate_exchange_rnx, product_exchange_rnx, case='carbon_yield', maxConcentration=None,
                newObjectiveReaction=None, saveName=None, substrate2zero='Ex_S_cpd00027_ext',
                printEq=False, save=False):
    """ Starting from the SBML model a json file is created so that the equations can be quickly constructed in pyomo
    Params:
        * modelName(str): the name of the model
        * substrate_exchange_rnx (list): list of substrate id's that are of interest
        * product_exchange_rnx (list): list of product id's that are of interest
        * case (str): what the yield schould be based on nl carbon yield or mass yield
        * maxConcentration (dict): the maximum concentration that a reactor can have for a certain product
        this parameter determines how much water needs to be added to the reactor (!!! kg/L !!!)
        * newObjectiveReaction (str): ID of the reaction you maximise (default is alway biomass)
        * missingCarbonId (str): ID of a metabolite that has a missing formula and you want to estimate it with the
        check carbon function
        * yieldTol (array): yield tolerances to accept an exchange metabolite as a potential substrate

    returns:
        a json file save in 'json models'
        allEquations
        allYields_FBA
    """

    #  read in the SBML model
    modelLocation = get_location(modelName)
    model = cobra.io.read_sbml_model(modelLocation)

    # make sure the right objective is set
    if newObjectiveReaction:
        model.objective = newObjectiveReaction

    # change the standard exchange reaction to zero
    exchange_rnx_2_zero = substrate2zero
    model.reactions.get_by_id(exchange_rnx_2_zero).bounds = 0, 0

    # make the save name
    if saveName is None:
        saveName = modelName.replace('.xml', '.json')

    # preallocate lists and dictionaries
    allYields_FBA = []
    allEquations = []
    coefDict = {}
    outputVariables = []

    # start the loop over the products
    if printEq:
        print(modelName)
    for product in product_exchange_rnx:
        productMet = model.reactions.get_by_id(product).reactants[0]
        productName = productMet.name
        Cprod = count_atom_in_formula(metabolite=productMet, atom='C')
        MW_prod = productMet.formula_weight

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

            # if the model does not consume the substrate then there is a problem
            if FBA_substrate_flux >= 0:  # consuming reactions are negative, that's why >= is used
                raise Exception(
                    'The model {} does not consume the substrate {} so a yield can not be obtained'.format(modelName,
                                                                                                           substrate))

            substrateMet = model.reactions.get_by_id(substrate).reactants[0]
            substrateName = substrateMet.name
            Csub = count_atom_in_formula(metabolite=substrateMet, atom='C')
            MW_sub = substrateMet.formula_weight

            # get the flux solutions
            FBA_product_flux = solution.fluxes[product]

            # get the right type of yield base on the given case
            if case == 'carbon_yield':
                FBA_yield = abs((FBA_product_flux / FBA_substrate_flux) * (Cprod * 12) / (
                        Csub * 12))  # in gramsC / grams C: 12 gCarbon/mol

            elif case == 'mass_yield':
                FBA_yield = abs(
                    (FBA_product_flux / FBA_substrate_flux) * MW_prod / MW_sub)  # in gramsC / grams C: 12 gCarbon/mol

            else:
                raise Exception(
                    "The variable 'case' (= {}) can only take the string 'carbon_yield' or 'mass_yield' ".format(case))

            allYields_FBA.append(FBA_yield)
            strEq += ' + {} * {}'.format(FBA_yield, substrateName)
            substrateCoefDict.update({substrateName: FBA_yield})
        coefDict.update({productName: substrateCoefDict})
        allEquations.append(strEq)
        if printEq:
            print(strEq)

    if printEq:
        print('')  # extra space to make it more readable

    surrogateModel = SurrogateModel(name=modelName, outputs=outputVariables, coef=coefDict,
                                    lable='SBML', maxConcentration=maxConcentration)

    if save:
        loc = os.getcwd()
        posAlquimia = loc.find('Alquimia')
        loc = loc[0:posAlquimia + 8]
        loc = loc + r'\json models' + r'\{}'.format(saveName)
        with open(loc, 'w+', encoding='utf-8') as f:
            json.dump(surrogateModel.__dict__, f, ensure_ascii=False, indent=4)

    return allEquations, allYields_FBA


def SBML_2_json_v2(modelName, substrate_exchange_rnx, product_exchange_rnx, maxConcentration=None,
                   newObjectiveReaction=None, saveName=None, exchRnx2zero='Ex_S_cpd00027_ext', yieldTol=None,
                   save=False, toIgnore=None, alreadyConsidered=None):
    """ Starting from the SBML model a json file is created so that the equations can be quickly constructed in pyomo
    Params:
        * modelName(str): the name of the model
        * substrate_exchange_rnx (list): list of substrate id's that are of interest
        * product_exchange_rnx (list): list of product id's that are of interest
        * case (str): what the yield schould be based on nl carbon yield or mass yield
        * maxConcentration (dict): the maximum concentration that a reactor can have for a certain product
        this parameter determines how much water needs to be added to the reactor (!!! kg/L !!!)
        * newObjectiveReaction (str): ID of the reaction you maximise (default is alway biomass)
        * missingCarbonId (str): ID of a metabolite that has a missing formula, and you want to estimate it with the
        check carbon function
        * yieldTol (array): yield tolerances to accept an exchange metabolite as a potential substrate
        * toIgnore (list): names of potential substrates that can be ignored
        * alreadyConsidered (list): names of potential substrates that can be automatically considered
    returns:
        a json file save in 'json models'
        allEquations
        allYields_FBA
    """

    #  read in the SBML model
    modelLocation = get_location(modelName)
    model = cobra.io.read_sbml_model(modelLocation)

    # make sure the right objective is set
    if newObjectiveReaction:
        model.objective = newObjectiveReaction

    # change the standard exchange reaction to zero
    # exchange_rnx_2_zero = substrate2zero
    # model.reactions.get_by_id(exchange_rnx_2_zero).bounds = 0, 0

    # make the save name
    if saveName is None:
        saveName = modelName.replace('.xml', '.json')

    # determine how to make the coefficient dictionary base on the input
    coefDict, considered, ignored = get_coef_all_substrates_SBML(modelName=model,
                                                                 substrateExchRxnIDs=substrate_exchange_rnx,
                                                                 productExchRxnIDs=product_exchange_rnx,
                                                                 yieldTol=yieldTol, exchRnx2zero=exchRnx2zero,
                                                                 ignore=toIgnore, include=alreadyConsidered)
    outputNames = list(coefDict.keys())
    surrogateModel = SurrogateModel(name=modelName, outputs=outputNames, coef=coefDict,
                                    lable='SBML', maxConcentration=maxConcentration)

    if save:
        loc = os.getcwd()
        posAlquimia = loc.find('Alquimia')
        loc = loc[0:posAlquimia + 8]
        loc = loc + r'\json models' + r'\{}'.format(saveName)
        with open(loc, 'w+', encoding='utf-8') as f:
            json.dump(surrogateModel.__dict__, f, ensure_ascii=False, indent=4)

    return surrogateModel, considered, ignored  # allEquations, allYields_FBA


def get_coef_all_substrates_SBML(modelName, substrateExchRxnIDs, productExchRxnIDs, yieldTol,
                                 exchRnx2zero='Ex_S_cpd00027_ext', ignore=None, include=None):
    """ Get the list of possible substrates from a model: Substrates have at least 3 carbons and produce a yield which is
    at least as big as the yield tolerance and is not a protein

    Params:
        model (str, model): can be a string of the model name or the model its self
        productExchRxnIDs (list): list of id strings of the desired products
        substrateExchRxnIDs (list or str): list of id strings of the desired substrates or str: 'select'
        exchRnx2zero (str): string of the original substrate exchange reaction that needs to be set to zero. Glusoe in
                            the case of the propioni bacteria
        yieldTol (array): yield tolarances to accept an exchange metabolite as a potential substrate
        ignore (list): list of substrate names to ignore
        include (list): list of substrate names to automatically include

    Returns:
        substrateList (list): list of possible substrates

    """

    # read in the model if necessary
    if isinstance(modelName, str):
        modelLocation = get_location(modelName)
        modelStrName = modelName
        #  read in the SBML model
        model = cobra.io.read_sbml_model(modelLocation)
    else:
        model = modelName
        modelStrName = model.name

    # change the original substrate to zero
    model.reactions.get_by_id(exchRnx2zero).bounds = 0, 1000

    # get the exchange reactions based on the input of substrateExchRxnIDs
    if isinstance(substrateExchRxnIDs, str) and substrateExchRxnIDs == 'select':
        allExchRxn = model.exchanges
        selectSwitch = True
    elif isinstance(substrateExchRxnIDs, list):
        allExchRxn = substrateExchRxnIDs
        selectSwitch = False
    else:
        raise Exception("the input variable 'substrateExchRxnIDs' must be a list or the string: 'select' ")

    # preallocate variables
    substrateLists = []
    yieldDict = {}

    # loop over all the products and possible substrates
    for i, prodId in enumerate(productExchRxnIDs):
        productRxn = model.reactions.get_by_id(prodId)
        productName = productRxn.reactants[0].name
        substrateCoefDict = {}
        yields = []
        substrateNames = []

        try:
            tolerance = yieldTol[prodId]
        except:
            tolerance = 0

        for rxnExch in allExchRxn:
            if isinstance(rxnExch, str):
                rxnExch = model.reactions.get_by_id(rxnExch)  # get rxn object from the model if a list of strings
            substrateMetabolite = rxnExch.reactants[0]
            substrateName = substrateMetabolite.name
            nCarbon = count_atom_in_formula(substrateMetabolite, atom='C')
            proteinCheck = is_protein_met(metabolite=substrateMetabolite)

            # filter out the substrates that aren't carbon based or are proteins
            if nCarbon >= 2 and not proteinCheck:

                # save the original bounds
                originalReactionBounds = rxnExch.bounds

                # change bound of new desired substrate to -10 mol/h/gDW
                rxnExch.bounds = -10, 1000
                FBA_yield = find_yield(model, substrateExchangeRxnID=rxnExch.id, productExchangeRxnID=prodId,
                                       printResults=False)

                # check if the yield is above the given tolerance
                if FBA_yield >= tolerance and FBA_yield < 1:  # bigger than the tolerance and smaller then 1
                    yields.append(FBA_yield)  # already in g/g
                    substrateNames.append(substrateName)
                    substrateCoefDict.update({substrateName: FBA_yield})

                # reset the bounds to the original bounds again
                rxnExch.bounds = originalReactionBounds

        yieldDict.update({productName: substrateCoefDict})
        substrateLists.append(substrateNames)

    # now we're going to select the same substrates for each product as the list of substrates with the least amount of
    # substrate names
    shortestList = min(substrateLists, key=len)

    # selectSwitch = False # comment this line, just temporaily so don't have to go through this each fucking time
    if selectSwitch:
        finalSubstrateList, ignored = substrate_run_through(substrateNames=shortestList, modelName=modelStrName,
                                                            ignore=ignore
                                                            , include=include)
    else:
        finalSubstrateList = shortestList
        ignored = []
        included = []

    finalDict = {}
    for key in yieldDict:
        substrateDict_help = yieldDict[key]
        listOfSubstrateKeys = substrateDict_help.keys()
        list2deleet = list(set(listOfSubstrateKeys) ^ set(finalSubstrateList))
        for key2deleet in list2deleet:
            del substrateDict_help[key2deleet]
        finalDict.update({key: substrateDict_help})

    return finalDict, finalSubstrateList, ignored


def substrate_run_through(substrateNames, modelName, ignore, include):
    """ this function loops over the names of the substrates and asks in the terminal if a certain element need not be
    considered as a substrate for the reactor model

    Params:
        substrateNames (list): list of all the names of the substrates
        substrateIds (list): list of all the id's of the given substrates

    Returns:
         deleetList (list): list of substrate names that need to be deleeted

    """

    if include is None:
        include = []
    if ignore is None:
        ignore = []

    # initiate lists (substrates you want to keep and ask permission to delete/keep)
    keepList = []
    excludeList = []
    alreadyIncluded = []

    # ----------------------------------- ask to include (include) list
    print("Here is a preview of the remaining substrates of model {}: ".format(modelName))
    print(substrateNames)
    print('')

    print("In the previous model(s) the previous substrates where considered: {} \n"
          "do you wish to include them aswell (y/n) ?".format(include))

    answerInclude = input()

    switchReAsk = True
    if answerInclude == 'y' or answerInclude == 'n':
        switchReAsk = False

    while switchReAsk:
        print("you can only type 'y' or 'n' ")
        answerInclude = input()
        if answerInclude == 'y' or answerInclude == 'n':
            switchReAsk = False

    if answerInclude == 'y':
        toInclude = list(set(substrateNames) & set(include))
        keepList += toInclude
        alreadyIncluded = toInclude

    # ----------------------------------- ask to exclude previous (ignore) list
    print("In the previous model(s) the previous substrates where excluded: {} \n"
          "do you wish to exclude them aswell (y/n) ?".format(ignore))

    answerExclude = input()

    switchReAsk = True
    if answerExclude == 'y' or answerExclude == 'n':
        switchReAsk = False

    while switchReAsk:
        print("you can only type 'y' or 'n' ")
        answerExclude = input()
        if answerExclude == 'y' or answerExclude == 'n':
            switchReAsk = False

    if answerExclude == 'y':
        excludeList += list(set(substrateNames) & set(ignore))

    # -------------------------------------------------------- make the final list to ask for
    askList = set(substrateNames) - set(excludeList) - set(alreadyIncluded)

    # if not askList: # so if it is empty
    #     askList = substrateNames

    # -------------------------------------------------------- loop over remaining substrates
    nList = len(askList)
    for i, substrate in enumerate(askList):
        print("Do you want to consider the following metabolite as a substrate: {}\n"
              "type 'y' (yes) or 'n' (no)".format(substrate.upper()))
        answer = input()
        percentage = (i + 1) / nList * 100
        print("{} % of the list completed".format(percentage))

        switchReAsk = True
        if answer == 'y' or answer == 'n':
            switchReAsk = False

        while switchReAsk:
            print("you can only type 'y' or 'n' ")
            answer = input()
            if answer == 'y' or answer == 'n':
                switchReAsk = False

        if answer == 'y':
            keepList.append(substrate)

    ignored = list(set(keepList) ^ set(substrateNames))
    return keepList, ignored


# # ------------------------------------------- -------------------------------------------
# Surrogate model functions distillation units
# # ------------------------------------------- -------------------------------------------
def simulate_distilation(x_D, x_B, F, x_F, alfa_f,  # for mass balances
                         Hvap_LK, Hvap_HK,  # for condenser duty
                         T_F, T_D, T_B, Cp_LK, Cp_HK,  # for reboiler duty
                         printResults=False):
    """
        Calculates the flow of mass, reflux ratio, and energy requirements for a distillation column.

        Parameters:

            F (float): flow rate of the incoming stream (kg/hr)
            x_F (float): composition of the LK in the feed component (mass %)
            x_D (float): desired composition of the distillate component (mass% of the LK)
            x_B (float): desired composition of the bottom component (mass% of the LK)
            alfa_f (float): vapor pressure the relative volatility is given by the ratio of vapor pressures,
                             and thus is a function only of temperature. (-)
            Hvap_LK (float): Evaporation enthalpy of the light key (kJ/kg)
            Hvap_HK (float): Evaporation enthalpy of the heavy key (kJ/kg)
            T_F (float): Temperature of the Feed (ªC or K)
            T_D (float): Temperature of the Distillate stream (ªC or K)
            T_B (float): Temperature of the Bottom stream (ªC or K)
            Cp_LK (float): Heat capacity of the light key (kJ/K/kg)
            Cp_HK (float): Heat capacity of the heavy key (kJ/K/kg)


        Returns:
            D (float): flow of distillate (kg/h)
            B (float): flow of bottom (kg/h)
            Q (float): energy requirements of the re-boiler (J/hr)
        """

    # flow of mass
    D = F * (x_F - x_B) / (x_D - x_B)  # in kg/h
    B = F - D  # in kg/h

    # reflux ratio assumed at 1.3
    L = (F * ((D * x_D) / (F * x_F) - alfa_f * D * (1 - x_D) / (F * (1 - x_F))) / (alfa_f - 1)) * 1.3  # in kg/h
    V = L + D  # in kg/h

    # condenser
    Hvap = x_D * Hvap_LK + (1 - x_D) * Hvap_HK  # in kJ/kg
    Qc = Hvap * V  # in kJ/kg * kg/h = kJ/h

    # re-boiler
    hF = (x_F * Cp_LK + (1 - x_F) * Cp_HK) * (T_F - T_D)  # in kJ/kg
    hB = (x_B * Cp_LK + (1 - x_B) * Cp_HK) * (T_B - T_D)  # in kJ/kg
    Qr = B * hB + Qc - F * hF  # kJ/h

    # total energy requirements the same as that of the re-boiler
    # Qtot = Qr #- Qc

    # transform the Qr to kwh power consumption per kg
    kw = Qr / 3600  # kJ/h to kJ/s = kW
    # assume 1 hour of operation?
    powerConsumption = kw / F  # kwh/kg

    seperationCoefDist = (F * x_F) / (D * x_D)
    seperationCoefBtm = (F * x_F) / (B * x_B)
    seperationCoef = [seperationCoefDist, seperationCoefBtm]

    # print statments
    if printResults:
        print('')
        print('the flow of the LK in the feed (kg/h) is {} mols\n'.format(F * x_F))
        print('the flow of distilate leaving (kg/h): {} where \n'
              'the LK has {} kg'.format(D, D * x_D))
        print('')
        print('the flow of bottom (kg/h): {} where \n'
              'the LK has {} kg\n'.format(B, B * x_B))
        print('the sum of the LK in the bottom and distilate is: {} \n'.format(D * x_D + B * x_B))
        print('Hvap (kJ/kg): {}'.format(Hvap))
        print('Qc (kJ/h): {}'.format(Qc))
        print('Qr (kJ/h): {}'.format(Qr))
        print('the sum of the dutys: {}'.format(Qr - Qc))
        print('the power consumption in kWh/kg is', powerConsumption)

    return powerConsumption, seperationCoef


def make_surrogate_model_distillation(xdata, ydata, polynomialDegree, case='Linear', plot=True, alfa = 0.1):
    """ Create the surrogate model (linear or lasso regresion ) for distillation units, which is only dependent on the
    variable x_F which is the fraction of the LK element in the feed

     Params:
        powerConsumption(list/array): the response variable, power consumption in kwh/kg as calculated by the
        function simulate_distillation

        x_F_vals(list/array): the variable of the regression, the fraction of the light key in the feed (mass%)

     Returns:
        regression model: lasso, ridge or linear model
     """

    # reshape the data
    xdata = np.array(xdata).reshape(-1, 1)
    ydata = np.array(ydata).reshape(-1, 1)

    # make the polynomial data
    poly = PolynomialFeatures(degree=polynomialDegree, include_bias=True)
    X_poly = poly.fit_transform(xdata)

    # Split data into training and testing sets

    X_train, X_test, y_train, y_test = train_test_split(X_poly, ydata, test_size=0.35, random_state=42)

    # Fit linear regression model to training data
    if case == 'Ridge':
        reg = Ridge().fit(X_train, y_train)
    elif case == 'Lasso':
        reg = Lasso(alpha=alfa).fit(X_train, y_train)
        # model = Lasso(alpha= 1, max_iter= 4000)
    elif case == 'Linear':
        reg = LinearRegression().fit(X_train, y_train)
    else:
        raise Exception("The string variable _case_ can only be 'Linear, 'Lasso' or 'Ridge'")

    # Predict Qtot for test data
    y_pred = reg.predict(X_test)
    r2_scores = r2_score(y_test, y_pred, multioutput='raw_values')
    MSE = mean_squared_error(y_test, y_pred, multioutput='raw_values')
    coefs = reg.coef_

    # print info
    print("R2 score for each output variable:", r2_scores)
    print("mean_squared_error for each output variable:", MSE)
    print('')
    print("the ceof are:", coefs)

    # only dependent on 1 variable, so we can plot the outcome of the model vs the pH
    # create data to plot the model

    x_data_model = np.linspace(min(xdata), max(xdata), num=100)
    x_data_model = x_data_model.reshape((len(x_data_model), 1))  # reshape
    x_poly_data_model = poly.fit_transform(x_data_model)
    y_data_model = reg.predict(x_poly_data_model)

    # check out the  plots
    if plot:
        #plot_parity_plots(yPred=y_pred, yObv=y_test)
        if len(y_data_model.shape) < 2:
            y_data_model = y_data_model.reshape((len(y_data_model), 1))
        plot_model_vs_data(x_data=X_poly[:, 1], y_data=ydata, x_data_model=x_data_model.squeeze(), y_data_model=y_data_model)
    return reg
