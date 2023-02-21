import mpmath
import numpy as np
import os
import json
import seaborn as sns
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Lasso
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import PolynomialFeatures
import matplotlib.pyplot as plt

def distillation_check(x_D, x_B, F, x_F, alfa_f,    # for mass balances
                      Hvap_LK, Hvap_HK,             # for condenser duty
                      T_F, T_D, T_B, Cp_LK, Cp_HK,  # for reboiler duty
                       printResults = True):
    """
        Calculates the flow of mass, number of stages, reflux ratio, and energy requirements for a distillation column.

        Parameters:

            F (float): flow rate of the incoming stream (kg/hr)
            x_F (float): composition of the LK in the feed component (mass %)
            x_D (float): desired composition of the distillate component (mass% of the LK)
            x_B (float): desired composition of the bottom component (mass% of the LK)
            alfa_f (float): vapor pressure the relative volatility is given by the ratio of vapor pressures,
                             and thus is a function only of temperature. (-)
            Hvap_LK (float): Evaporation enthalpy of the light key (kJ/kg)
            Hvap_HK (float): Evaporation enthalpy of the heavy key (kJ/kg)
            T_F (float): Temperature of the Feed (ªC or K)
            T_D (float): Temperature of the Distillate stream (ªC or K)
            T_B (float): Temperature of the Bottom stream (ªC or K)
            Cp_LK (float): Heat capacity of the light key (kJ/K/kg)
            Cp_HK (float): Heat capacity of the heavy key (kJ/K/kg)


        Returns:
            D (float): flow of distilaate (kg/h)
            B (float): flow of bottom (kg/h)
            Q (float): energy requirements (J/hr)
        """

    # flow of mass
    D = F*(x_F - x_B)/(x_D - x_B)
    B = F - D

    # reflux ratio
    L = (F*( (D*x_D)/(F*x_F) - alfa_f * D*(1-x_D)/ ( F*(1-x_F)) ) / (alfa_f -1)) * 1.3
    V = L + D

    # condenser
    Hvap = x_D * Hvap_LK + (1-x_D) * Hvap_HK
    Qc = Hvap * V

    # reboiler
    hF = (x_F* Cp_LK + (1-x_F) * Cp_HK) *(T_F - T_D)
    hB = (x_B* Cp_LK + (1-x_B) * Cp_HK) *(T_B - T_D)
    Qr = B* hB + Qc - F*hF

    print('')
    Qtot = Qr - Qc

    # print statments
    if printResults:
        print('the flow of the LK in the feed is {} mols\n'.format(F * x_F))
        print('the flow of distilate leaving: {} where \n'
              'the LK has {} mols'.format(D, D * x_D))
        print('')
        print('the flow of bottom: {} where \n'
              'the LK has {} mols\n'.format(B, B * x_B))
        print('the sum of the LK in the bottom and distilate is: {} \n'.format(D * x_D + B * x_B))
        print('Hvap: {}'.format(Hvap))
        print('Qc: {}'.format(Qc))
        print('Qr: {}'.format(Qr))
        print('the sum of the dutys: {}'.format(Qtot))

    return Qtot

def make_distillation_equations(T_F, T_D, T_B,                       # temperatures
                               VaporPressureLK, VaporPressureHK,    # vapor presures
                               Cp_LK, Cp_HK, Hvap_LK, Hvap_HK,      # constants
                               intervalName):

    """
    abbv: Light Key (LK), Heavy Key (HK)

    TO gather from pyomo:
    1) FEED_FlOW (Kg/h): "Feed"
    2) FEED_FLOW compostion (% mass of the LK): "x_F"
    TO be given to the model
    3) Distillate composition (% mass of the LK): "x_D"
    4) Bottoms composition (% mass of the LK): "x_B"

    # Parameters 2 pass on
    1) alfa_f: relative volativity = VP_LK / VP_HK (at the temperature of the reboiler more or less)
    VP_LK, VP_HK (Pa)
    Variable
    2) T_F : temperature of the feed (K)
    3) T_D : temperature of the distilate (K)
    4) T_B : temperature of the bottom (K)
    Constant
    5) Cp_LK: heat capacity LK (kJ/kg/K)
    6) Cp_HK: heat capacity HK (kJ/kg/K)
    7) Hvap_LK: vapor enthalpy LK (kJ/Kg)
    8) Hvap_HK: vapor enthalpy HK (kJ/Kg)


    """

    # determin realtive volatility
    alfa_f = VaporPressureLK / VaporPressureHK

    # preallocate all variables necessary
    pyomoVariables = []
    pyomoEquations = []

    # -----------  flow distillate and bottom equations
    # variables
    DistilateVar = 'Distillate_{}'.format(intervalName)
    BottomVar = 'Bottom_{}'.format(intervalName)

    # equations
    DistilateEq = "model.var['{}'] == Feed * (x_F - x_B) / (x_D - x_B)".format(DistilateVar)
    # Feed, x_F, x_B, x_D need to be replaced by parameters form the Excel file
    BottomsEq = "model.var['{}'] == Feed - model.var['{}'] ".format(BottomVar,DistilateVar)

    # add to variable and equaition lists
    pyomoVariables += [DistilateVar,BottomVar]
    pyomoEquations += [DistilateEq,BottomsEq]

    # ----------- Internal flows
    # variables
    Lvar = "L_{}".format(intervalName)
    Vvar = "V_{}".format(intervalName)

    # equations
    Internal_L_Eq = "model.var['{0}'] == (Feed * ((model.var['{1}'] * x_D) / (Feed * x_F) - " \
                    "{2} * model.var['{1}'] * (1 - x_D) / (Feed * (1 - x_F))) / ({2} - 1)) * 1.3".format(Lvar,DistilateVar, alfa_f)
    Internal_V_Eq =  "model.var['{}'] == model.var['{}'] + model.var['{}'] ".format(Vvar, Lvar, DistilateVar)

    # add to variable and equaition lists
    pyomoVariables += [Lvar, Vvar]
    pyomoEquations += [Internal_L_Eq, Internal_V_Eq]

    # ----------- condenser duty
    QcVar = "Qc_{}".format(intervalName) # condenser variable in kJ/h
    #Hvap = (x_D * {} + (1 - x_D) * {})
    CondenserDutyEq = "model.var['{}'] == (x_D * {} + (1 - x_D) * {}) * model.var['{}']".format(QcVar, Hvap_LK, Hvap_HK, Vvar)

    # add to variable and equaition lists
    pyomoVariables += [QcVar]
    pyomoEquations += [CondenserDutyEq]

    # ----------- reboiler duty
    EntalpyFeedVar = "hF_{}".format(intervalName)
    EntalpyBottomVar = "hB_{}".format(intervalName)
    QrVar = "Qr_{}".format(intervalName) # reboiler variable in kJ/h

    EnthalpyFeedEq = "model.var['{}'] == (x_F * {} + (1 - x_F) * {}) * ({} - {})".format(EntalpyFeedVar, Cp_LK, Cp_HK, T_F, T_D)
    EnthalpyBottomEq = "model.var['{}'] == (x_B * {} + (1 - x_B) * {}) * ({} - {})".format(EntalpyBottomVar, Cp_LK, Cp_HK, T_B, T_D)
    DutyReboilerEq = "model.var['{}'] == model.var['{}'] * model.var['{}'] + model.var['{}'] - Feed * model.var['{}'] ".format(QrVar, BottomVar, EntalpyBottomVar,QcVar, EntalpyFeedVar  )

    # add to variable and equaition lists
    pyomoVariables += [EntalpyFeedVar, EntalpyBottomVar, QrVar]
    pyomoEquations += [EnthalpyFeedEq, EnthalpyBottomEq, DutyReboilerEq]


    # ----------- total duty
    QtotVar = "Q_tot_{}".format(intervalName)
    QtotEq = "model.var['{}'] == (model.var['{}'] - model.var['{}']) / 3600".format(QtotVar, QrVar, QcVar)
    # /3600 to get the units into kWh (assuming 1 hour of operation)

    # add to variable and equaition lists
    pyomoVariables += [QtotVar]
    pyomoEquations += [QtotEq]

    return pyomoVariables, pyomoEquations

# check witht the example given in the text book!
distillation_check(x_D=0.90, x_B=0.05, F= 100, x_F= 0.6, alfa_f=3.1,
                  Hvap_LK= 33800, Hvap_HK= 38000,
                  T_F= 90, T_D= 82 , T_B= 108, Cp_LK= 133, Cp_HK= 157 )

# make preliminary object for distilation2
var, eq = make_distillation_equations(T_F= 80, T_D= 90, T_B= 100,               # temperatures (K)
                                VaporPressureLK = 9, VaporPressureHK= 3.5,    # vapor presures (Pa)
                                Cp_LK= 130, Cp_HK = 157,                       # heat capacity (KJ/Kg/K)
                                Hvap_LK = 33000 , Hvap_HK = 40000,              # constants (KJ/Kg)
                                intervalName= 'DIST1')

jsonSwitch = True
if jsonSwitch:
    jsonDict = {"variables": var,
                "equations":eq,
                "light_key": 'water',
                "heavy_key": 'ace,prop',
                "lable": "distillation",
                "name": 'DIST1'
                }

    saveName = "testDistilation.json"
    loc = os.getcwd()
    posAlquimia = loc.find('Alquimia')
    loc = loc[0:posAlquimia + 8]
    loc = loc + r'\json models' + r'\{}'.format(saveName)
    with open(loc, 'w+', encoding='utf-8') as f:
        json.dump(jsonDict, f, ensure_ascii=False, indent=4)
