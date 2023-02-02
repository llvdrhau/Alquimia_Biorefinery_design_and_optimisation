import mpmath
import numpy as np


def distilation_check(x_D, x_B, F, x_F, alfa_f,  # for mass balances
                      Hvap_a, Hvap_b,       # FOR condenser duty
                      T_F, T_D, T_B, Cp_a, Cp_b): # for reboiler duty
    """
        Calculates the flow of mass, number of stages, reflux ratio, and energy requirements for a distillation column.

        Parameters:
            x_D (float): desired composition of the distillate component (mass%)
            F (float): flow rate of the incoming stream (kg/hr)
            x_F (float): composition of the feed component (mass %)
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
    print(D)
    print(B)

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
    Qr = B* hB + Qc -F*hF
    print('Qr: {}'.format(Qr))



    return D, B, L


distilation_check(x_D=0.95, x_B=0.05, F= 450, x_F= 0.6, alfa_f=3.1,
                  Hvap_a= 33800, Hvap_b= 38000,
                  T_F= 90, T_D= 82 , T_B= 108, Cp_a= 133, Cp_b= 157 )