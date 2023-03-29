
from f_make_surrogate_model import simulate_distilation
import numpy as np
import math
# -------------------------------------------------
#           Distillation unit 2
# -------------------------------------------------
# VFAs - water separation
# light-key: water
# heavy-key: VFA (grouped together)

# data on the feed (variable)
x_F = np.linspace(0.01, 0.99) # in mass %
# F =  np.linspace(100, 500)    # kg/h
#x_F = 0.3
F = 10     # kg/h

# desired separation outcome
x_D = 0.9  # VFA's mass % in distillate
x_B = 0.05 # VFA's mass % in bottom

# vapour pressure LK
# https://www.engineeringtoolbox.com/water-vapor-saturation-pressure-d_599.html
VP_LK = 1.96 * 1.0135 # bar at 120 °C

# vapour pressure HK
# https://webbook.nist.gov/cgi/cbook.cgi?ID=C79094&Mask=4&Type=ANTOINE&Plot=on#ANTOINE
T = 120 + 273 # temperature in K near distillate temperature
VP_HK_log10 = 4.74558 - (1679.869 / (T -59.832))
VP_HK = 10**VP_HK_log10 # in bar.

# vapor presure ratio
alfa = VP_LK / VP_HK # vapour Pressure ratio

# vapour enthalpies heavy-key light-key
Hvap_LK = 2202 # kJ/kg # https://www.engineeringtoolbox.com/water-properties-d_1573.html?vA=120&units=C#
Hvap_HK = 755 # kJ/kg  # https://webbook.nist.gov/cgi/cbook.cgi?ID=C79094&Mask=4

# temperatures
T_F = 25       # feed temperature (C)
T_D = 95       # distillate temperature (C) # close to the boiling point of the light key
T_B = 120      # bottom temperature (C)     # close to the boiling point of the light key (propionate)

# heat capacities
Cp_LK = 4.184 # (kJ/K/kg)
Cp_HK = 2.334 # (kJ/K/kg)

Q = simulate_distilation(x_D= x_D, x_B= x_B, F= F, x_F= x_F, alfa_f= alfa,          # for mass balances
                         Hvap_LK= Hvap_LK, Hvap_HK= Hvap_HK,                        # for condenser duty
                         T_F= T_F, T_D=T_D, T_B=T_B, Cp_LK= Cp_LK, Cp_HK= Cp_LK, printResults=True)    # for reboiler duty