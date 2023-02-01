import numpy as np


def distillation_design(x_D, q, Cp, T_F, T_D, T_B):
    """
    Calculates the number of stages, reflux ratio, and energy requirements for a distillation column.

    Parameters:
        x_D (float): desired composition of the distillate component
        q (float): flow rate of the incoming stream (kg/hr)
        Cp (float): heat capacity of the components (J/kg-K)
        T_F (float): temperature of the incoming stream (K)

        T_D (float): desired temperature of the distillate (K)
        T_B (float): desired temperature of the bottoms (K)

    Returns:
        N (int): number of stages
        R (float): reflux ratio
        Q (float): energy requirements (J/hr)
    """

    # calculate the flows of the Bottom and Distilate



    # Calculate the number of stages using De Fenske equation
    N = (np.log((1 - x_D) / x_D)) / (np.log(q * (Cp * (T_F - T_D)) / (q * Cp * (T_B - T_D))))

    # Calculate the reflux ratio using Underwood equation
    R = (1 - x_D) / (x_D - (1 - x_D) * np.exp(-N))

    # Calculate the energy requirements
    Q = q * Cp * (T_F - T_D) * N * R / (1 + R)

    return int(round(N)), R, Q


def distilation_check(x_D, F, x_F, alfa_f):
    """
        Calculates the flow of mass, number of stages, reflux ratio, and energy requirements for a distillation column.

        Parameters:
            x_D (float): desired composition of the distillate component (mass%)
            F (float): flow rate of the incoming stream (kg/hr)
            x_F (float): composition of the feed component (mass %)
            alfa_f (float): vapor presure the relative volatility is given by the ratio of vapor pressures, a1;2 ¼ Ps
                             and thus is a function only of temperature.
            nex


        Returns:
            D (float): flow of distilaate (kg/h)
            B (float): flow of bottom (kg/h)
            Q (float): energy requirements (J/hr)
        """

    # flow of mass
    x_B = 1-x_D
    D = F*(x_F - x_B)/(x_D - x_B)
    B = F - D
    print(D)
    print(B)

    # reflux ratio
    Lmin = F*( (D*x_D)/(F*x_F) - alfa_f * D*(1-x_D)/ (F*(1-x_F)) ) / (alfa_f -1)
    L = 1.3 * Lmin
    print(L)


    return D, B, L

distilation_check(x_D=0.95, F= 450, x_F= 0.6, alfa_f=3.1)