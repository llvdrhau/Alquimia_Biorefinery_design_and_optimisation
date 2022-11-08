# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 17:08:25 2022
test how I could use classes instead of seperate funtions 
@author: Reboots
"""

class makeReactor:
    def __init__(self, inputs, outputs, eq, utilities = []):
        self.inputs  = inputs 
        self.outputs = outputs
        self.utilities = utilities
        
        if isinstance(eq, str): # fail safe if you forget to code the string reactor expression as a list 
            self.eq = [eq]   # make it a list
        else: 
            self.eq = eq 
        
    def makeDict(self):
        in_outDict = { 
                'inputs': self.inputs,
                'outputs':self.outputs 
                }
        return in_outDict
    
    # alternative class method 
    @classmethod
    def fromNeuralNetwork(cls,neuralNetwork):
        # so here comes the alternative constructor for neural networks wich will make the string equations from a machine learning model 
        # https://www.youtube.com/watch?v=rq8cL2XMM5M&ab_channel=CoreySchafer 
        pass 
   

# think of sizing example   
############################################################ Randy 

strExpr = ['str1 == in1 *0.8 ']   
         
ins = ['in1']
outs = ['str1']
  
Randy = makeReactor(ins, outs, strExpr)  


############################################################ R1

strExpr = ['str1 == in1 *0.8 ',
           'str2 == in1*0.7' ]             
ins = ['in1']
outs = ['str1','str2']
utilities  =0.0014 

R1 = makeReactor(ins, outs, strExpr, utilities)  



############################################################ R2
strExpr = 'p2 == str2 *0.8 ' 
ins = ['str2']
outs = ['p2']

R2 = makeReactor(ins, outs, strExpr)  

############################################################ R3
strExpr = ['p3a == str2*0.8 + in2 *0.2',
               'p3b == str2 *0.2 + in2 *0.1' ]
ins = ['in2','str2']
outs = ['p3a','p3b']
utilities  =0.0014 
R3 = makeReactor(ins, outs, strExpr,utilities)  

############################################################ R4 
strExpr = ['str4 == str1 * 0.9 ']
            
ins = ['str1']
outs = ['str4']
R4 = makeReactor(ins, outs, strExpr)  


############################################################ R5
strExpr = ['p5 == str1 *0.3+5']
            
ins = ['str1']
outs = ['p5']

R5 = makeReactor(ins, outs, strExpr)  

############################################################ Rw

strExpr = ['pw == in2 * 0.6 ']
            
ins = ['in2']
outs = ['pw']
Rw = makeReactor(ins, outs, strExpr)  
   

############################################################ Rq

strExpr = ['pq == str1 * 0.3  ']
ins = ['str1']
outs = ['pq']
Rq = makeReactor(ins, outs, strExpr)  


############################################################ Rq1

strExpr = ['strq == str1 * 0.3  ']
            
ins = ['str1']
outs = ['strq']

RQ1 = makeReactor(ins, outs, strExpr)  


############################################################ Rq2 

strExpr = ['pq == strq * 0.3 +5']
            
ins = ['strq']
outs = ['pq']

RQ2 = makeReactor(ins, outs, strExpr)  


############################################################ Rx
strExpr = ['px == str4 * 0.1 ']
            
ins = ['str4']
outs = ['px']

Rx = makeReactor(ins, outs, strExpr)  

############################################################ Ry               
strExpr = ['py == str4 *0.95']
            
ins = ['str4']
outs = ['py']
Ry = makeReactor(ins, outs, strExpr)


if __name__ == '__main__' :
    
    test = getattr(Rx, 'inputs')
    a = Rx.eq