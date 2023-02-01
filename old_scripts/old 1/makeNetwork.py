# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 13:33:00 2022

@author: lucas
"""

import pandas as pd
import numpy as np
import matplotlib as plt
from schemdraw import flow

strm = pd.read_excel('superstructureTest.xlsx', sheet_name = 'streams', index_col=0)  
conectionMatrix = pd.read_excel('superstructureTest.xlsx', sheet_name = 'connections', index_col=0)  

# node data 
nodeDF = pd.read_excel('superstructureTest.xlsx', sheet_name= 'nodeMatrix', index_col=0)
nodeNames = nodeDF.index.values
nodeMatrix = nodeDF.to_numpy()