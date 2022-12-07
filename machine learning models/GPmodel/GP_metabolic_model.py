# -*- coding: utf-8 -*-
"""
Created on Thu May 12 16:58:00 2022

@author: Reboots
"""

import pandas as pd 
#pylab inline
import matplotlib.pyplot as plt
import numpy as np 
#pylab.ion()
import GPy
#from example_multi_ouput_GP_model import plot_2outputs


trainingDataX = pd.read_excel('TrainingData.xlsx', sheet_name = 'input')
trainingDataX.columns = [c.replace(' ', '_') for c in trainingDataX.columns]
trainingDataY = pd.read_excel('TrainingData.xlsx', sheet_name = 'output')
trainingDataY.columns = [c.replace(' ', '_') for c in trainingDataY.columns]

dimX = len(trainingDataX.columns)
dimY  =len(trainingDataY.columns)

# get input data
X1 = trainingDataX.pH.to_numpy() 
X2 = trainingDataX.percentage_protein.to_numpy()  
# reshape 
X1 = np.reshape(X1, (105,1))
X2 = np.reshape(X2, (105,1))

# get output data
Y = trainingDataY.acetate.to_numpy()
Y = np.reshape(Y, (105,1))

plotOn = True 
subPLotOn = False
if plotOn and subPLotOn: 
    plt.figure(figsize=(15, 12), dpi=100)
    for i,colNames in enumerate(trainingDataX):
        data = trainingDataX[colNames]
        plt.subplot(3,7,i+1)
        plt.hist(data, bins = 13)
        plt.xlabel(colNames, fontsize = 12)
        
elif plotOn and not subPLotOn :
    plt.figure()
    for i,colNames in enumerate(trainingDataX):
        data = trainingDataX[colNames]
        plt.hist(data, bins = 13)
        plt.xlabel(colNames, fontsize = 12)
        plt.show()


kernel = GPy.kern.RBF(input_dim=1, variance=1., lengthscale=1.)
m = GPy.models.GPRegression(X2,Y,kernel)
print(m)
m.optimize(messages=True)
print(m)

# plot
plotOn2 = False 
if plotOn2: 
    fig = m.plot()
    GPy.plotting.show(fig, filename='basic_gp_regression_density_notebook_optimized')
#%%  n dimentional model 

names = ['pH','percentage_protein','percentage_Carbohydrates']
names = trainingDataX.columns
dimIn = len(names)

Xmulti = trainingDataX[names].to_numpy()
# define kernel
ker = GPy.kern.Matern52(dimIn,ARD=True) + GPy.kern.White(dimIn)
# create simple GP modelb
m = GPy.models.GPRegression(Xmulti,Y,ker)
print(m)
# optimize and plot
if dimIn == 2: 
    m.optimize(messages=True,max_f_eval = 1000)
    fig = m.plot()
    print(GPy.plotting.show(fig, filename='basic_gp_regression_notebook_2d'))
    print(m)
    
#%% validation test 

validationDataX = pd.read_excel('ValidationData.xlsx', sheet_name = 'input')
validationDataX.columns = [c.replace(' ', '_') for c in trainingDataX.columns]
validationDataY = pd.read_excel('ValidationData.xlsx', sheet_name = 'output')
validationDataY.columns = [c.replace(' ', '_') for c in trainingDataY.columns]

# get output data
Yobv = validationDataY.acetate.to_numpy()
Yobv = np.reshape(Yobv, (20,1))

newX =validationDataX.to_numpy()
Ynew = m.predict(newX)
Ynew = Ynew[0]

diagX = range(20)

plt.figure(1)
plt.plot(Yobv, Ynew, 'r+')
plt.plot(diagX,diagX)



 #%% Multiple output regresion ; TODO fix it 
 
# =============================================================================
# kern = GPy.kern.Matern32(1)
# icm = GPy.util.multioutput.ICM(input_dim=dimX, num_outputs=dimY, kernel= kern)
# print(icm) 
# 
# X1 =trainingDataX.to_numpy()
# X = trainingDataX.to_numpy().transpose()
# Y = trainingDataY.to_numpy().transpose() 
# 
# xList = []
# for i in range(len(X)):
#     xList.append(X[i].transpose())
# 
# yList = []               
# for i in range(len(Y)): 
#     yList.append(Y[i].transpose())    
# 
# xTest = []
# for i in range(len(Y)): 
#     xTest.append(X.transpose())
#     
# m = GPy.models.GPCoregionalizedRegression(X1,yList,kernel=icm)
# m['.*Mat32.var'].constrain_fixed(1.) #For this kernel, B.kappa encodes the variance now.
# m.optimize()
# print (m)
# =============================================================================

