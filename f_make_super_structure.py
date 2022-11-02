# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 10:43:34 2022

Builds further upone the script makeModelWithClassObjects
The goal is to make a superstructure where the streams have different
compositions plus intergrating the generic process interval
@author: Lucas Van der Hauwaert
"""
# imports
import pandas as pd
import numpy as np
import pyomo.environ as pe
import pyomo.opt as po
import importlib
import os
from f_makeIntervalObjects import make_reactor_intervals, make_input_intervals


def make_super_structure():
    pass