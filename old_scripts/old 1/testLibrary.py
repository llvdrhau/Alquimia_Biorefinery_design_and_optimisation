# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 11:50:12 2022
test library of models 

@author: lucas van der hauwaert
"""

# for neural networks I'll need to sum over the vector 
# https://stackoverflow.com/questions/41833664/pyomo-summation-of-a-product-of-a-matrix-by-a-vector
# FROm NN you need to make string expresions of the various outputs!!!! so you get a list of expresions



def Randy(): 
    strExpr = ['str1 == in1 *0.8 ']   
               
    ins = ['in1']
    outs = ['str1']
    
    in_outDict = { 
            'inputs': ins,
            'outputs':outs 
            }
    
    return  in_outDict , strExpr 

def R1(): 
    strExpr = ['str1 == in1 *0.8 ',
               'str2 == in1*0.7' ]   
               
    ins = ['in1']
    outs = ['str1','str2']
    
    in_outDict = { 
            'inputs': ins,
            'outputs':outs 
            }
    
    return  in_outDict , strExpr 

def R2(): 
    strExpr = 'p2 == str2 *0.8 ' 
    ins = ['str2']
    outs = ['p2']
    
    in_outDict = { 
            'inputs': ins,
            'outputs':outs 
            }
    return  in_outDict , strExpr

def R3(): 
    strExpr = ['p3a == str2*0.8 + in2 *0.2',
                   'p3b == str2 *0.2 + in2 *0.1' ]
    ins = ['in2','str2']
    outs = ['p3a','p3b']
    
    utility = 0.1  # euros/kg electricity consumption 
    
    
    in_outDict = { 
            'inputs': ins,
            'outputs':outs 
            }
    return  in_outDict , strExpr, utility 

def R4(): 
    strExpr = ['str4 == str1 * 0.9 ']
                
    ins = ['str1']
    outs = ['str4']
    
    in_outDict = { 
            'inputs': ins,
            'outputs':outs 
            }
    return  in_outDict , strExpr

def R5(): 
    strExpr = ['p5 == str1 *0.3+5']
                
    ins = ['str1']
    outs = ['p5']
    
    in_outDict = { 
            'inputs': ins,
            'outputs':outs 
            }
    return  in_outDict , strExpr

def Rw(): 
    strExpr = ['pw == in2 *6 ']
                
    ins = ['in2']
    outs = ['pw']
    
    in_outDict = { 
            'inputs': ins,
            'outputs':outs 
            }
    return  in_outDict , strExpr

def Rq(): 
    strExpr = ['pq == str1 * 0.3  ']
                
    ins = ['str1']
    outs = ['pq']
    
    in_outDict = { 
            'inputs': ins,
            'outputs':outs 
            }
    return  in_outDict , strExpr

def Rx(): 
    strExpr = ['px == str4 * 0.1 ']
                
    ins = ['str4']
    outs = ['px']
    
    in_outDict = { 
            'inputs': ins,
            'outputs':outs 
            }
    return  in_outDict , strExpr

def Ry(): 
    strExpr = ['py == str4 *0.95']
                
    ins = ['str4']
    outs = ['py']
    
    in_outDict = { 
            'inputs': ins,
            'outputs':outs 
            }
    return  in_outDict , strExpr





'''
COULD BE VERY INTERESTING TO IMPLEMENT THE REACTORS WITH IN CLASSES 

# =============================================================================
# class makeReactor:
#     def __init__(self, inputs, outputs, eq):
#         self.inputs  = inputs 
#         self.outputs = outputs
#         
#         if isinstance(self.eq,str): # fail safe if you forget to code the string reactor expression as a list
#             strExpr_help = [self.eq]   # make it a list 
#             self.eq = strExpr_help
#         else: 
#             self.eq = eq 
#         
#     def makeDict(self):
#         in_outDict = { 
#                 'inputs': self.inputs,
#                 'outputs':self.outputs 
#                 }
#         return in_outDict
#    
#       
#       
# 
# exprStr = ['p3a == str2*0.8 + in2 *1.1+ 8.8',
#                'p3b == str2 *0.2 + in2 *0.9' ]
# inputs = ['in2','str2']
# outputs = ['p3a','p3b']
# 
# R1 = makeReactor(inputs, outputs, exprStr)  
# 
# =============================================================================

'''



