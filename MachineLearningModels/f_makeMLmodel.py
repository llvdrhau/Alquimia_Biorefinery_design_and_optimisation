# code to automatically make an elastic net regression model
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
from sklearn.model_selection import train_test_split
from sklearn.linear_model import ElasticNetCV
import pickle


# info at https://scikit-learn.org/stable/modules/linear_model.html#elastic-net
# https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.ElasticNetCV.html#sklearn.linear_model.ElasticNetCV
def makeElasticNetModel (ExcelName, iterations = 1000, modelName = '', plotSwitch = 0, case = 'ridge'):
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


if __name__ == '__main__':
    #output = makeElasticNetModel('GelatineData_elastic_net.xlsx') #, modelName= 'SAVE TEST')
    output = makeElasticNetModel('Gelatine_Glucose_Data.xlsx')
    # print equations
    eq = output[0]
    for i in eq:
        print(i)
    # show plot of
    plt = output[1]
    plt.show()



