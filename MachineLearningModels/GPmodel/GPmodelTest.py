# -*- coding: utf-8 -*-
"""
Created on Thu May 12 10:25:36 2022

@author: Reboots
"""

import GPy 
import numpy as np
import chart_studio.plotly as plt
#from matplotlib import pyplot as plt

X = np.random.uniform(-3.,3.,(20,1))
Y = np.sin(X) + np.random.randn(20,1)*0.05

kernel = GPy.kern.RBF(input_dim=1, variance=1., lengthscale=1.)
                      
#type GPy.kern.<tab> here:
#GPy.kern.RBF?

m = GPy.models.GPRegression(X,Y,kernel)
# to find the values of the parameters that maximize the likelihood of the data 
m.optimize(messages=True) 
#m.optimize_restarts(num_restarts = 20)

# display model parameters 
print(m)

# plot 
fig = m.plot()
GPy.plotting.show(fig, filename='basic_gp_regression_density_notebook_optimized')

#%% 
import GPy 
import numpy as np
import chart_studio.plotly as plt
# sample inputs and outputs
X = np.random.uniform(-3.,3.,(50,2))
Y = np.sin(X[:,0:1]) * np.sin(X[:,1:2])+np.random.randn(50,1)*0.05

# define kernel
ker = GPy.kern.Matern52(2,ARD=True) + GPy.kern.White(2)

# create simple GP model
m = GPy.models.GPRegression(X,Y,ker)

# optimize and plot
m.optimize(messages=True,max_f_eval = 1000)
fig = m.plot()
print(GPy.plotting.show(fig, filename='basic_gp_regression_notebook_2d'))
print(m)


