# -*- coding: utf-8 -*-
"""
Created on Mon May 16 16:51:51 2022

@author: lucas
"""

 # coregionalized regression modeL 
# Multiple output regresion ; TODO fix it 
 
import pandas as pd 
#pylab inline
import matplotlib.pyplot as plt
import numpy as np 
#pylab.ion()
import GPy

#%% GET data 
trainingDataX = pd.read_excel('TrainingData.xlsx', sheet_name = 'input')
trainingDataX.columns = [c.replace(' ', '_') for c in trainingDataX.columns]
trainingDataY = pd.read_excel('TrainingData.xlsx', sheet_name = 'output')
trainingDataY.columns = [c.replace(' ', '_') for c in trainingDataY.columns]

dimX = len(trainingDataX.columns)
dimY  =len(trainingDataY.columns)

# get input data
X = trainingDataX.to_numpy() 

# get output data
Y = trainingDataY.to_numpy()

xList = []
for i in range(dimY):
    xList.append(X)

yList = []               
for i in trainingDataY: 
    yHelp = trainingDataY[i].to_numpy()
    yHelp = np.reshape(yHelp, (105,1))
    yList.append(yHelp)  


#%% make model
 
kern = GPy.kern.Matern32(dimX)
icm = GPy.util.multioutput.ICM(input_dim=dimX, num_outputs=dimY, kernel= kern)
print(icm) 
    
m = GPy.models.GPCoregionalizedRegression(xList,yList, kernel=icm)
m['.*Mat32.var'].constrain_fixed(1.) #For this kernel, B.kappa encodes the variance now.
m.optimize()
print(m)

#%% 
import pickle
pickle.dump( m, open( "save.coModel", "wb" ) )

# =============================================================================
# with open('save.pkl', 'wb') as file:
#     pickle.dump(m, file)
# with open('save.pkl', 'rb') as file:
# loaded_model = pickle.load(file)
# https://stackoverflow.com/questions/64557786/how-to-save-load-optimized-gpy-regression-model
# =============================================================================



#%% 

    
