# -*- coding: utf-8 -*-
"""
Created on Fri May 13 10:20:55 2022

@author: Reboots
"""
import numpy as np
import matplotlib.pyplot as pb
import GPy

def plot_2outputs(m,xlim,ylim):
    fig = pb.figure(figsize=(18,10))
    #Output 1
    ax1 = fig.add_subplot(211)
    ax1.set_xlim(xlim)
    ax1.set_title('Output 1')
    m.plot(plot_limits=xlim,fixed_inputs=[(1,0)],which_data_rows=slice(0,100),ax=ax1)
    ax1.plot(Xt1[:,:1],Yt1,'rx',mew=1.5)
    #Output 2
    ax2 = fig.add_subplot(212)
    ax2.set_xlim(xlim)
    ax2.set_title('Output 2')
    m.plot(plot_limits=xlim,fixed_inputs=[(1,1)],which_data_rows=slice(100,200),ax=ax2)
    ax2.plot(Xt2[:,:1],Yt2,'rx',mew=1.5)
    

if __name__ == '__main__': 
  
    #This functions generate data corresponding to two outputs
    f_output1 = lambda x: 4. * np.cos(x/5.) - .4*x - 35. + np.random.rand(x.size)[:,None] * 2.
    f_output2 = lambda x: 6. * np.cos(x/5.) + .2*x + 35. + np.random.rand(x.size)[:,None] * 8.
    
    
    #{X,Y} training set for each output
    X1 = np.random.rand(100)[:,None]; X1=X1*75
    X2 = np.random.rand(100)[:,None]; X2=X2*70 + 30
    Y1 = f_output1(X1)
    Y2 = f_output2(X2)
    #{X,Y} test set for each output
    Xt1 = np.random.rand(100)[:,None]*100
    Xt2 = np.random.rand(100)[:,None]*100
    Yt1 = f_output1(Xt1)
    Yt2 = f_output2(Xt2)
    
    a = [X1,X2]
    b = [Y1,Y2]
    
    K = GPy.kern.Matern32(1)
    icm = GPy.util.multioutput.ICM(input_dim=1,num_outputs=2,kernel=K)
    
    m = GPy.models.GPCoregionalizedRegression([X1,X2],[Y1,Y2],kernel=icm)
    m['.*Mat32.var'].constrain_fixed(1.) #For this kernel, B.kappa encodes the variance now.
    m.optimize()
    print(m)
    plot_2outputs(m,xlim=(0,100),ylim=(-20,60))