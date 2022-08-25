# code to automaticlly make an elastic net regression model
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import time
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split
from sklearn import linear_model
from sklearn.preprocessing import PolynomialFeatures

def makeElasticNetModel (ExcelName, hyperparameter=1):
    x = pd.read_excel(ExcelName,sheet_name='inputs')
    inputnames = x.keys()
    y = pd.read_excel(ExcelName,sheet_name='outputs')
    outputNames = y.keys()
    stringEquationVector= []
    for i in y:
        # might need to normalise data....
        # mean standard normilsation if inputs are not in the same decimal rank?
        # select output
        Y = y[i]
        # split training data/ test data
        X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.2)
        # Create elastic net regression objects
        en_regr = linear_model.ElasticNet(alpha= hyperparameter)
        # Train the model using the training sets
        en_regr.fit(X_train, y_train)
        # Make predictions using the testing set
        y_pred_en = en_regr.predict(X_test)
        # plot to evaluate goodness of fit

        plt.figure()
        plt.xlabel('observered')
        plt.ylabel('predicted')
        plt.plot( y_test,y_pred_en,'.')
        plt.plot([min(y_test),min(y_pred_en)],[max(y_test),max(y_pred_en)],'-')
        plt.title('Plot to evaluate fit of the model')
        plt.legend([ "elastic Net", "digonal"])
        plt.show()

        # add equation to the vector of strings = 'y = ax +b'
        yName = outputNames[i]
        eq = yName + '= '
        for j,xname in enumerate(inputnames):
            txt1 = "My name is {fname}, I'm {age}".format(fname="John", age=36)
            eq = eq + xname + '*{0}'.format(en_regr.coef_[j])







x = np.linspace(0, 2*np.pi, 50)
y = 5*np.sin(x) + 0.1*np.random.randn(50)

X = PolynomialFeatures(3).fit_transform(x.reshape((50,1)))  # polynomial features

print(X.shape)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

# Create linear regression object
regr = linear_model.LinearRegression()

# Train the model using the training sets
regr.fit(X_train, y_train)

# Make predictions using the testing set
y_pred = regr.predict(X_test)

plt.figure()
plt.plot(x, regr.predict(X), '.')
plt.plot(x, y, 'k')
plt.title("Linear Regression");
