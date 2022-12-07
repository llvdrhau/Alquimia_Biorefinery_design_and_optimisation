import numpy as np
from sklearn.model_selection import train_test_split

from sklearn import linear_model

import matplotlib.pyplot as plt


x = np.linspace(0, 2*np.pi, 50)
y = 5*np.sin(x) + 0.1*np.random.randn(50)

from sklearn.preprocessing import PolynomialFeatures

X = PolynomialFeatures(3).fit_transform(x.reshape((50,1)))  # polynomial features

print(X.shape)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

# Create linear regression object
regr = linear_model.RidgeCV()

# Train the model using the training sets
regr.fit(X_train, y_train)

# Make predictions using the testing set
y_pred = regr.predict(X_test)

plt.figure()
plt.plot(x, regr.predict(X), '.')
plt.plot(x, y, 'k')
plt.title("Linear Regression")
plt.show()