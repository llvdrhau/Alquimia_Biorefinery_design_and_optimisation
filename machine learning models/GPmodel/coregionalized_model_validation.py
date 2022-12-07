# -*- coding: utf-8 -*-
"""
Created on Wed May 18 09:41:38 2022

@author: Reboots
"""
import numpy as np 
import pickle 
import matplotlib.pyplot as plt
import pandas as pd
import math



#%% 
validationDataX = pd.read_excel('ValidationData.xlsx', sheet_name = 'input')
validationDataX.columns = [c.replace(' ', '_') for c in validationDataX.columns]

validationDataY = pd.read_excel('ValidationData.xlsx', sheet_name = 'output')
validationDataY.columns = [c.replace(' ', '_') for c in validationDataY.columns]
validationDataY.columns = [c.replace('_gCOD/L', '') for c in validationDataY.columns]


#%%

with open('save.coModel', 'rb') as file:
    m = pickle.load(file)
    
#%% get output data
newXx = validationDataX.to_numpy() #[:,None]
row, col = validationDataY.shape

yAll = [] 
for i in range(col):
    toStack = np.ones((row,1)) * i 
    newX = np.hstack([newXx,toStack])
    noise_dict = {'output_index':newX[:,-1:].astype(int)}
    y = m.predict(newX,Y_metadata=noise_dict)
    yAll.append(y)



#%% plots 
plt.figure(1) 
plt.figure(figsize=(15, 12), dpi=100)
for i,colName in enumerate(validationDataY):
    plt.subplot(2,4,i+1)
    y = validationDataY[colName].to_numpy()[:,None]
    y_predicted = yAll[i][0]
    y_err = yAll[i][1]
    y_err = np.reshape(y_err,(20,))
    #y_err = y_err.tolist()
    #plt.plot(y,y_predicted,'rx')
    plt.errorbar(y, y_predicted, yerr = y_err , fmt = 'x',color = 'orange',)
    diagDim = max(y) + max(y)*0.1  #round(max(y) + max(y)*0.1)  
    diagX = np.linspace(0, diagDim)
    plt.plot(diagX,diagX)
    plt.xlabel('y_modelAlberte')
    plt.ylabel('Y_predicted')
    plt.title(colName)
plt.show()    


#%% 
        
plt.figure()
for i,colNames in enumerate(validationDataX):
    data = validationDataX[colNames]
    plt.hist(data, bins = 13, color = 'red')
    plt.xlabel(colNames, fontsize = 12)
    plt.show()



