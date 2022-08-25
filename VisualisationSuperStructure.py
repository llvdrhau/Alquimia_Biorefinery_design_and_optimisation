# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 09:47:07 2022

Visualisation of the superstructure 

@author: Lucas Van der hauwaert 
"""

import pandas as pd
import networkx as nx
import importlib
import numpy as np 
import warnings 

def visualiseSuperStructure(excelFile,reactorLibrary,positions): 
     
    components = pd.read_excel(excelFile, sheet_name = 'components')  
    layers = components.layer.to_numpy()
    compounds = components.components.tolist()
    conectionMatrix = pd.read_excel(excelFile, sheet_name = 'connectionMatrix')  
    reactorNames = conectionMatrix.Reactor.tolist()
    
    
    layerDict = {compounds[i]:layers[i] for i in range(len(layers))} # helping dict to make the one here under 
    layerAndElementDict = {} # prealllocate 
    for i in range(max(layers)+1):
        helpList =[] # reset for each new layer 
        for j in layerDict.keys():
            if layerDict[j] == i :
                helpList.append(j)    
        layerAndElementDict.update({i:helpList})     
            
      
    #layerDict2  = {layers[i]:}
    
    inputNodeList = []
    outputNodeList = []
    streamList = []
    
    for i  in compounds:
        if layerDict[i] == 1 :
            inputNodeList.append(i) 
        elif layerDict[i] == max(layers): 
            outputNodeList.append(i)
        else: 
            streamList.append(i)
    
    if len(inputNodeList) ==  len(positions) :
        startPositionsDict = {inputNodeList[i]:positions[i] for i in range(len(positions))}
    else: 
        warnings.warn('the starting positions and amount of inputs nodes do not coinside')
        
    
    # creat object G for information on the superstructure
    G = nx.Graph()
    
 # input layer ################################################################
    # create the input nodes and assign them a position 
    postionDict = {}
    for i in inputNodeList: 
        x = 1
        y = startPositionsDict[i]
        G.add_node(i,pos=(x,y))
        postionDict.update({i: [x, y]})
   
        
# characterise edges as a dictionary  ##########################################
        
    connectionDict = {}    
    for i in conectionMatrix:
        a = conectionMatrix[i]
        locateEntry = a == 1 
        locateBoolean =  a == 2  
        locateSplit = a == 3 
        
        if sum(locateBoolean) > 1 : 
            typeConnection = 'boolean'
            pos = np.where(locateBoolean == True)[0]
            reactorList = [reactorNames[x] for x in pos]   
            reactorList.append(typeConnection)
            connectionDict.update({i:reactorList})
            
        if sum(locateSplit) > 1 :
            typeConnection = 'split'
            pos = np.where(locateSplit == True)[0]
            reactorList = [reactorNames[x] for x in pos] 
            reactorList.append(typeConnection)
            connectionDict.update({i:reactorList})    
            
        if sum(locateEntry) >= 1 :
            typeConnection = 'singleEnrty'
            pos = np.where(locateEntry == True)[0]
            reactorList = [reactorNames[x] for x in pos] 
            reactorList.append(typeConnection)
            connectionDict.update({i:reactorList}) 
            
            
# Next layer after input  #####################################################
#edges   
    
    conectionMatrix.pop('Reactor') 
    connectionList = list(connectionDict.keys()) 
    
    nLayers = max(layers)
   
    layerAndNodeDict= {}
    for i in range(1,nLayers+1) : 
        for j in layerAndElementDict[i]: # elements in each layer 
            reactorListHelp = []
            for r in reactorNames:
                get_reactor = getattr(importlib.import_module(reactorLibrary), r)
                outReactor = get_reactor.outputs
                if j in outReactor :
                    reactorListHelp.append(r)
            layerAndNodeDict.update({i:r})
                     
   # TODO  okye so thjis part just trying to get a dictionary with the nodes on the right layer 
    # but i think this should be specifide in the  excel file I think                  
        
        
            
    
    
    
    
    for i in reactorNames:
        get_reactor = getattr(importlib.import_module(reactorLibrary), i)
        inReactor = get_reactor.inputs
        outReactor = get_reactor.outputs
        # TODO focous here on joining the first and subsequent layer 

        for j in inReactor :
            if j in inputNodeList: # connection to the first layer 
                x = layerDict[j]
                
                G.add_node(i,pos=(x,y))
        
        ##################################################### #####################################################
        #####################################################
        #####################################################
        for j in outReactor: 
            if j in connectionList :
                x = layerDict[j]

                G.add_node(i,pos=(x,y))
                edge1 = i # the reactor one node behind 
                reactor2conect = connectionDict[j]
                case = reactor2conect[-1]
                reactor2conect.pop()
                for x in reactor2conect:
                    edge2 = x 
                    if case == 'boolean':
                        G.add_edge(edge1,edge2, color='r')
                        #nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
                    else: 
                        G.add_edge(edge1,edge2, color='g')
                
# =============================================================================
#                 for x in connectionDict[j]
#                     G.add_node(i,pos=(x,y))
#                     G.add_edge(i,j, color='r')
#                 
#          
# =============================================================================
        
        
        
        
# =============================================================================
#         
#         
#         for j in inReactor :
#             if j in inReactor and j in inputList :
#                 x = 2 #layerDict[2]
#                 y = postionDict[j][1]
#                 G.add_node(i,pos=(x,y))
#                 G.add_edge(j,i)
#                 postionDict.update({i: [x, y]})
#                 
#             
#                 
#             if j in outsReactor and j in streamList :
#                 a= 1
#                 
#                 
#         
#    
#         
#         
# =============================================================================


# =============================================================================
#     
#     
#     G.add_edge(nodeNames[i],nodeNames[j])
#     
#     pos=nx.get_node_attributes(G,'pos')
#     nx.draw(G,pos,with_labels=True)
#     
#     
#             
#         if nrInputs == 1 :
#             y = yMax/2 
#         else : 
#             y = 2
#             
# =============================================================================
    
    
    
# =============================================================================
#     inputList = allNodes.Inputs.dropna()
#     outputList = allNodes.Products.dropna()
#     reactorList = allNodes.Reactors.dropna()
# =============================================================================
    
    
# =============================================================================
#     # print(allNodes[:Inputs])
#     
#     # node data 
#     nodeDF = pd.read_excel('superstructureTest.xlsx', sheet_name= 'nodeMatrix', index_col=0)
#     nodeNames = nodeDF.index.values
#     nodeMatrix = nodeDF.to_numpy()
#     
#     
#     #any('in1' in inputList)
#     
#     
#     # creat reprisentation of superstructure
#     G = nx.Graph()
#     #G.add_nodes_from(nodeNames)
#     
#     # conected nodes and position nodes correctly
#     posInpunts = 1
#     posOutputs = 1
#     posReactor = 1
#     for i in nodeNames:
#         a = i == inputList
#         b = i == outputList
#         c = i == reactorList
#         if any(a):
#            G.add_node(i,pos=(1,posInpunts))
#            posInpunts +=1
#         if any(b):
#             G.add_node(i,pos=(5,posOutputs))
#             posOutputs += 1 
#         if any(c):
#             checkConectionRow = nodeDF[i:].to_numpy()
#             nInputs = len(inputList)
#             if any(checkConectionRow[0,0:nInputs]) :
#                 G.add_node(i, pos =(2,posReactor))
#             else: 
#                 G.add_node(i, pos =(3,posReactor))
#             posReactor += 1
#             
#             
#             
#     # create conections 
#     n = len(nodeNames)
#     for i in range(n):  
#         for j in range(n):
#             if nodeMatrix[i,j] == 1 :  #could change this to is string, so edge cqn be labled
#                 G.add_edge(nodeNames[i],nodeNames[j])
#     
#     #pos = nx.planar_layout(G)
#     pos=nx.get_node_attributes(G,'pos')
#     
#     nx.draw(G,pos,with_labels=True)
#     
#     #nx.draw_networkx_labels(G)
#     # G = nx.petersen_graph()
#     # subax1 = plt.subplot(121)
#     # nx.draw(G, with_labels=True, font_weight='bold')
#     # subax2 = plt.subplot(122)
#     # nx.draw_shell(G, nlist=[range(5, 10), range(5)], with_labels=True, font_weight='bold')
#     
#     
#     # import networkx as nx
#     
#     # G=nx.Graph()
#     
#      
#     # G.add_node(1,pos=(1,1))
#     
#      
#     # G.add_node(2,pos=(2,2))
#     
#      
#     # G.add_edge(1,2)
#     
#      
#     # pos=nx.get_node_attributes(G,'pos')
#     
#     # nx.draw(G,pos)
#     
# =============================================================================

if __name__ == '__main__' :
    visualiseSuperStructure('testModel2.xlsx','testLibraryClass', [2.5,7.5])
    
    