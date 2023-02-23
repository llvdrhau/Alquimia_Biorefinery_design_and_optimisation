# -*- coding: utf-8 -*-
"""
Created on Tue May 17 18:33:41 2022

@author: lucas
"""
import pickle

m = pickle.load(open("save.p", "rb"))
print(m)