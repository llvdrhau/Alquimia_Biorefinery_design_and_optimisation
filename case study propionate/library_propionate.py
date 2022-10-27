'''
Created on Tue sep 20 2022

library of reactors for the case study of propionate production

@author: Lucas Van der hauwaert
email: lucas.vanderhauwaert@usc.es
'''

from f_makeIntervalObjects import *

# Input
objInputDict = makeInputIntervals(r'\data_propionibacteria.xlsx')
objReqctorDict = makeReactorInterval(r'\data_propionibacteria.xlsx')

#estObje = makeReactorInterval(r'\data_propionibacteria.xlsx')
# loop over the dictionary to put the names of the intervals in the script
for i in objInputDict:
    locals()[i] = objInputDict[i]

for i in objReqctorDict:
    locals()[i] = objReqctorDict[i]

#print(Glucose.compositionDict)

# glucose = inputCharaterisation()
# # Reactors
# P_acidi_batch = makeReactor()
# P_freu_batch = makeReactor()
# P_avi_batch = makeReactor()
# P_acn_batch = makeReactor()
# P_prop_batch = makeReactor()
# mix_culture = makeReactor()
