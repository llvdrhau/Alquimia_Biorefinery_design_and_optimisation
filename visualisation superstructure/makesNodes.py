# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 22:31:13 2022
This code reads an excel file with the data of th super structure
@author: lucas
"""

# import sys 
# print( sys.executable )

import pandas as pd
#import numpy as np
#import matplotlib as plt
import networkx as nx

strm = pd.read_excel('superstructureTest.xlsx', sheet_name = 'streams', index_col=0)  
conectionMatrix = pd.read_excel('superstructureTest.xlsx', sheet_name = 'connections', index_col=0)  

allNodes = pd.read_excel('superstructureTest.xlsx', sheet_name = 'NODES')
inputList = allNodes.Inputs.dropna()
outputList = allNodes.Products.dropna()
reactorList = allNodes.Reactors.dropna()


# print(allNodes[:Inputs])

# node data 
nodeDF = pd.read_excel('superstructureTest.xlsx', sheet_name= 'nodeMatrix', index_col=0)
nodeNames = nodeDF.index.values
nodeMatrix = nodeDF.to_numpy()


#any('in1' in inputList)


# creat reprisentation of superstructure
G = nx.Graph()
#G.add_nodes_from(nodeNames)

# conected nodes and position nodes correctly
posInpunts = 1
posOutputs = 1
posReactor = 1
for i in nodeNames:
    a = i == inputList
    b = i == outputList
    c = i == reactorList
    if any(a):
       G.add_node(i,pos=(1,posInpunts))
       posInpunts +=1
    if any(b):
        G.add_node(i,pos=(4,posOutputs))
        posOutputs += 1 
    if any(c):
        checkConectionRow = nodeDF[i:].to_numpy()
        nInputs = len(inputList)
        if any(checkConectionRow[0,0:nInputs]) :
            G.add_node(i, pos =(2,posReactor))
        else: 
            G.add_node(i, pos =(3,posReactor))
        posReactor += 1
        
        
        
# create conections 
n = len(nodeNames)
for i in range(n):  
    for j in range(n):
        if nodeMatrix[i,j] == 1 :  #could change this to is string, so edge cqn be labled
            G.add_edge(nodeNames[i],nodeNames[j])

#pos = nx.planar_layout(G)
pos=nx.get_node_attributes(G,'pos')

nx.draw(G,pos,with_labels=True)

#nx.draw_networkx_labels(G)
# G = nx.petersen_graph()
# subax1 = plt.subplot(121)
# nx.draw(G, with_labels=True, font_weight='bold')
# subax2 = plt.subplot(122)
# nx.draw_shell(G, nlist=[range(5, 10), range(5)], with_labels=True, font_weight='bold')


# import networkx as nx

# G=nx.Graph()

 
# G.add_node(1,pos=(1,1))

 
# G.add_node(2,pos=(2,2))

 
# G.add_edge(1,2)

 
# pos=nx.get_node_attributes(G,'pos')

# nx.draw(G,pos)


