import mpmath
import numpy as np


def distilation_check(x_D, x_B, F, x_F, alfa_f,  # for mass balances
                      Hvap_a, Hvap_b,       # FOR condenser duty
                      T_F, T_D, T_B, Cp_a, Cp_b): # for reboiler duty
    """
        Calculates the flow of mass, number of stages, reflux ratio, and energy requirements for a distillation column.

        Parameters:

            F (float): flow rate of the incoming stream (kg/hr)
            x_F (float): composition of the LK in the feed component (mass %)

            x_D (float): desired composition of the distillate component (mass% of the LK)
            x_B (float): desired composition of the bottom component (mass% of the LK)

            alfa_f (float): vapor presure the relative volatility is given by the ratio of vapor pressures,
                             and thus is a function only of temperature.
            nex


        Returns:
            D (float): flow of distilaate (kg/h)
            B (float): flow of bottom (kg/h)
            Q (float): energy requirements (J/hr)
        """

    # flow of mass
    D = F*(x_F - x_B)/(x_D - x_B)
    B = F - D
    print('the flow of the LK in the feed is {} mols\n'.format(F*x_F))
    print('the flow of distilate leaving: {} where \n'
          'the LK has {} mols'.format(D, D*x_D))
    print('')
    print('the flow of bottom: {} where \n'
          'the LK has {} mols\n'.format(B, B*x_B))
    print('the sum of the LK in the bottom and distilate is: {} \n'.format(D*x_D + B*x_B))

    # reflux ratio
    L = (F*( (D*x_D)/(F*x_F) - alfa_f * D*(1-x_D)/ (F*(1-x_F)) ) / (alfa_f -1)) * 1.3
    V = L + D
    print(L)
    print(V)

    # condenser
    Hvap = x_D * Hvap_a + (1-x_D) * Hvap_b
    Qc = Hvap * V
    print('Hvap: {}'.format(Hvap))
    print('Qc: {}'.format(Qc))


    # reboiler
    hF = (x_F* Cp_a + (1-x_F) * Cp_b) *(T_F - T_D)
    hB = (x_B* Cp_a + (1-x_B) * Cp_b) *(T_B - T_D)
    Qr = B* hB + Qc - F*hF
    print('Qr: {}'.format(Qr))

    print('')
    print('the sum of the dutys: {}'.format(Qr - Qc))



    return D, B, L

def make_distilation_equations(T_F, T_D, T_B,                       # temperatures
                               VaporPressureLK, VaporPressureHK,    # vapor presures
                               Cp_LK, Cp_HK, Hvap_LK, Hvap_HK,      # constants
                               intervalName):

    """
    abbv: Light Key (LK), Heavy Key (HK)

    TO gather from pyomo:
    1) FEED_FlOW (Kg/h): "Feed"
    2) FEED_FLOW compostion (% mass of the LK): "x_F"
    3) Distillate composition (% mass of the LK): "x_D"
    4) Bottoms composition (% mass of the LK): "x_B"

    # Parameters 2 pass on
    1) alfa_f: relative volativity = VP_LK / VP_HK (at the temperature of the reboiler more or less)
    VP_LK, VP_HK
    Variable
    2) T_F : temperature of the feed
    3) T_D : temperature of the distilate
    4) T_B : temperature of the bottom
    Constant
    5) Cp_LK: heat capacity LK
    6) Cp_HK: heat capacity HK
    7) Hvap_LK: vapor enthalpy LK
    8) Hvap_HK: vapor enthalpy HK


    """

    # determin realtive volatility
    alfa_f = VaporPressureLK / VaporPressureHK

    # preallocate all variables necessary
    pyomoVariables = []
    pyomoEquations = []

    # -----------  flow distillate and bottom equations
    # variables
    DistilateVar = 'Distilate_{}'.format(intervalName)
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
    QcVar = "Qc_{}".format(intervalName)
    #Hvap = (x_D * {} + (1 - x_D) * {})
    CondenserDutyEq = "model.var['{}'] == (x_D * {} + (1 - x_D) * {}) * model.var['{}']".format(QcVar, Hvap_LK, Hvap_HK, Vvar)

    # add to variable and equaition lists
    pyomoVariables += [QcVar]
    pyomoEquations += [CondenserDutyEq]

    # ----------- reboiler duty
    EntalpyFeedVar = "hF_{}".format(intervalName)
    EntalpyBottomVar = "hB_{}".format(intervalName)
    DutyReboilerVar = "Qr_{}".format(intervalName)

    EnthalpyFeedEq = "model.var['{}'] == (x_F * {} + (1 - x_F) * {}) * ({} - {})".format(EntalpyFeedVar, Cp_LK, Cp_HK, T_F, T_D)
    EnthalpyBottomEq = "model.var['{}'] == (x_B * {} + (1 - x_B) * {}) * ({} - {})".format(EntalpyBottomVar, Cp_LK, Cp_HK, T_B, T_D)
    DutyReboilerEq = "model.var['{}'] == model.var['{}'] * model.var['{}'] + model.var['{}'] - Feed * model.var['{}'] ".format(DutyReboilerVar, BottomVar, EntalpyBottomVar,QcVar, EntalpyFeedVar  )

    # add to variable and equaition lists
    pyomoVariables += [EntalpyFeedVar, EntalpyBottomVar, DutyReboilerVar]
    pyomoEquations += [EnthalpyFeedEq, EnthalpyBottomEq, DutyReboilerEq]

    return pyomoVariables, pyomoEquations


distilation_check(x_D=0.90, x_B=0.05, F= 100, x_F= 0.6, alfa_f=3.1,
                  Hvap_a= 33800, Hvap_b= 38000,
                  T_F= 90, T_D= 82 , T_B= 108, Cp_a= 133, Cp_b= 157 )


var, eq = make_distilation_equations(T_F= 80, T_D= 90, T_B= 100,                       # temperatures
                               VaporPressureLK = 9, VaporPressureHK= 3.5,    # vapor presures
                               Cp_LK= 130, Cp_HK = 157, Hvap_LK = 33000 , Hvap_HK = 40000,      # constants
                               intervalName= 'DIST1')

for i in eq:
    print(i)
    print('')