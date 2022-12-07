import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def meanCenterStandardization(data,colNames):
    # Substract mean, normalized by std.dev
    mu = np.mean(data, axis=0)  # mean vector
    #print('Mean: {}'.format(mu))
    Xmu = data - mu  # mean-centered
    sd = np.std(data, axis=0)
    Xmustd = Xmu / sd  # mean centered and normalized by std. dev.
    Xmustd = pd.DataFrame(Xmustd,columns = colNames)
    return Xmustd

def plotSeveralOutputs(dataPrediction, dataObserved,shapePlots,titleVec):
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.4, wspace=0.4)
    for i in range(dataPrediction.shape[1]):
        plt.subplot(shapePlots[0], shapePlots[1], i+1)
        x = dataObserved[:, i]
        y = dataPrediction[:, i]
        diag = np.linspace(min(y), max(y))
        plt.plot(x, y, 'o')
        plt.plot(diag, diag)
        plt.title(titleVec[i])
    plt.show()
    # fig, ax = plt.subplots(shapePlots[0], shapePlots[1])
    # for i in range(dataPrediction.shape[1]) : #go over columns
    #     x = dataObserved[:,i]
    #     y = dataPrediction[:,i]
    #     diag = np.linspace(min(y), max(y))
    #     ax[i] = plt.plot(x,y,'o')
    #     ax[i] = plt.plot(diag, diag)
    #     plt.show()
    # plt.show()