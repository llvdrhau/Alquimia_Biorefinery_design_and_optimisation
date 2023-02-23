import pandas as pd
import numpy as np
from standerdization import meanCenterStandardization, plotSeveralOutputs
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
import keras
from keras.models import Sequential
from keras.layers import Dense


# import tensorflow as tf
# from tensorflow import keras
# from tensorflow.keras import layers
# from tensorflow.keras.layers.experimental import preprocessing

# import data
df_features = pd.read_excel('Gelatinedata.xlsx', sheet_name = 'features')
df_outputs = pd.read_excel('Gelatinedata.xlsx', sheet_name = 'outputs')
df_mw =  pd.read_excel('Gelatinedata.xlsx', sheet_name = 'MW')
# see structure of the data
print(df_features.head())
print(df_outputs.head())
print(df_outputs.keys())
# standerdise the features
namesFeatures = df_features.keys()
arrayFeatures = df_features.to_numpy()
standerdizedFeatures = meanCenterStandardization(arrayFeatures,namesFeatures)


namesOutputs = df_outputs.keys()
arrayOutputs = df_outputs.to_numpy()
standerdizedOutputs = meanCenterStandardization(arrayOutputs,namesOutputs)

# pair wise plots to see the distribitions and potential correlations
plotParwise = False
if plotParwise:
    plt.figure()
    sns.pairplot(standerdizedFeatures)
    plt.show()

# split the data
X_train, X_test, y_train, y_test = train_test_split(standerdizedFeatures, standerdizedOutputs, test_size=0.2)

# act = "relu"
# act = "sigmoid"
act = "tanh"

model = Sequential()
model.add(Dense(3, input_dim=len(df_features.keys()), activation=act))
model.add(Dense(30, activation=act)) # 12
model.add(Dense(18, activation=act)) # 24 waren goei aantal nodes
model.add(Dense(len(df_outputs.keys())))
model.compile(loss='mean_absolute_error', optimizer='adam', metrics=['accuracy']) # Model compilation for training
a = model.summary()
print(a)

#train model
history = model.fit(X_train, y_train, epochs=100, batch_size=5, validation_data=(X_test.to_numpy(), y_test.to_numpy()))
# validate model
y_pred = model.predict(X_test.to_numpy())
#yPred = y_pred.to_numpy()
yObv = y_test.to_numpy()
t = y_test.keys()
plotSeveralOutputs(y_pred,yObv,(3,2),titleVec= t)

model.evaluate(X_test.to_numpy(),y_test.to_numpy())

# save model

model.save('NN_gelatine')