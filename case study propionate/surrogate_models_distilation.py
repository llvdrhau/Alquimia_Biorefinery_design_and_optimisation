
from f_make_surrogate_model import simulate_distilation, make_surrogate_model_distillation, regression_2_json_v2
import numpy as np

# -------------------------------------------------
#           Distillation unit 2
# -------------------------------------------------
# VFAs - water separation
# light-key: water
# heavy-key: VFA (grouped together)

# data on the feed (variable)
n_samples = 50
np.random.seed(5)
x_F = np.random.uniform(0.2, 0.8, size=n_samples)

# give input feed
F = 1000   # kg/h

# desired separation outcome
x_D = 0.9  # VFA's mass % in distillate
x_B = 0.05 # VFA's mass % in bottom

# temperatures
T_F = 25       # feed temperature (C)
T_D = 95       # distillate temperature (C) # close to the boiling point of the light key (water 100 C)
T_B = 120      # bottom temperature (C)     # close to the boiling point of the heavy key (propionic acid, 141 C)

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
Hvap_LK = 46/72*1e3 # kJ/kg # https://www.engineeringtoolbox.com/water-properties-d_1573.html?vA=120&units=C#
Hvap_HK = 3774 # kJ/kg  # https://webbook.nist.gov/cgi/cbook.cgi?ID=C79094&Mask=4

# heat capacities
Cp_LK = 4.184 # (kJ/K/kg)
Cp_HK = 2.334 # (kJ/K/kg)

powerConsumption =[]
for xf in x_F:
    Q = simulate_distilation(x_D= x_D, x_B= x_B, F= F, x_F= xf, alfa_f= alfa,          # for mass balances
                         Hvap_LK= Hvap_LK, Hvap_HK= Hvap_HK,                        # for condenser duty
                         T_F= T_F, T_D=T_D, T_B=T_B, Cp_LK= Cp_LK, Cp_HK= Cp_LK, printResults=False)    # for reboiler duty
    powerConsumption.append(Q)


# fit model
n = 4 # ploynomial degree
reg = make_surrogate_model_distillation(xdata= x_F, ydata=powerConsumption, polynomialDegree=n, case='Linear', alfa= 1,
                                        plot=True)


# save the model
featureNames = []
for i in range(n+1):
    featureNames.append('x_F**{}'.format(i))

regression_2_json_v2(outputNames= 'energy_consumption', inputNames = 'x_F' ,featureNames=featureNames, model= reg,
                     saveName= 'Distillation_2.json', lable= 'Distillation_Regresion', lightKey = 'water')



# -------------------------------------------------
#           Distillation unit 3
# -------------------------------------------------
# acetate - propionate separation
# light-key: acetate
# heavy-key: propionate

# data on the feed (variable)
n_samples = 50
np.random.seed(42)
x_F = np.random.uniform(0.1, 0.9, size=n_samples)

# give input feed
F = 1000   # kg/h

# desired separation outcome
x_D = 0.95  # LK acetic acid mass % in distillate
x_B = 0.05 # LK acetic acid mass % in bottom

# operating temperatures
T_F = 100       # feed temperature (C)
T_D = 115      # distillate temperature (C) # close to the boiling point of the light key (118 °C acetic acid)
T_B = 138      # bottom temperature (C)     # close to the boiling point of the heavy key (141.2 °C propionic acid)

# vapour pressure LK (acetate)
# http://www.ddbst.com/en/EED/PCP/VAP_C84.php
VP_LK =  57.06 * 1e-2 # bar at 100 °C

# vapour pressure HK (propionate)
# https://webbook.nist.gov/cgi/cbook.cgi?ID=C79094&Mask=4&Type=ANTOINE&Plot=on#ANTOINE
T = 100 + 273 # temperature in K near distillate temperature
VP_HK_log10 = 4.74558 - (1679.869 / (T -59.832))
VP_HK = 10**VP_HK_log10 # in bar.

# vapor pressure ratio
alfa = VP_LK / VP_HK # vapour Pressure ratio

# vapour enthalpies heavy-key light-key
Hvap_LK = 2202 # kJ/kg; acetic acid  # https://www.engineeringtoolbox.com/water-properties-d_1573.html?vA=120&units=C#
Hvap_HK = 3774 # kJ/kg; propionic acid # https://webbook.nist.gov/cgi/cbook.cgi?ID=C79094&Mask=4

# heat capacities
Cp_LK = 2.050 # (kJ/K/kg) acetic acid    https://rapidn.jrc.ec.europa.eu/substance/acetic-acid
Cp_HK = 2.334 # (kJ/K/kg) propionic acid

powerConsumption =[]
for xf in x_F:
    Q = simulate_distilation(x_D= x_D, x_B= x_B, F= F, x_F= xf, alfa_f= alfa,          # for mass balances
                         Hvap_LK= Hvap_LK, Hvap_HK= Hvap_HK,                        # for condenser duty
                         T_F= T_F, T_D=T_D, T_B=T_B, Cp_LK= Cp_LK, Cp_HK= Cp_LK, printResults=False)    # for reboiler duty
    powerConsumption.append(Q)


# fit model
n = 4 # ploynomial degree
reg = make_surrogate_model_distillation(xdata= x_F, ydata=powerConsumption, polynomialDegree=n, case='Linear', alfa= 1,
                                        plot=True)


# save the model
featureNames = []
for i in range(n+1):
    featureNames.append('x_F**{}'.format(i))

regression_2_json_v2(outputNames= 'energy_consumption', inputNames = 'x_F' ,featureNames=featureNames, model= reg,
                     saveName= 'Distillation_3.json', lable= 'Distillation_Regresion', lightKey = 'ace')