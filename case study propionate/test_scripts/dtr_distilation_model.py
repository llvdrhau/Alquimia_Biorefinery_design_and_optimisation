from f_make_surrogate_model import simulate_distilation
import numpy as np
import matplotlib.pyplot as plt
from sklearn.tree import DecisionTreeRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from scipy.spatial import Delaunay, ConvexHull


# script options
saveSwitch = False
nDist2 = 2 #polynomial
nDist3 = 2
# -------------------------------------------------
#           Distillation unit 2
# -------------------------------------------------
# VFAs - water separation
# light-key: water
# heavy-key: VFA (grouped together)

# data on the feed (variable)
n_samples = 50
np.random.seed(5)
x_F = np.random.uniform(0.70, 0.95, size=n_samples)

# give input feed
F = 1000   # kg/h

# desired separation outcome (i.e. compositions for the light key in bottom and distillate stream)
x_D = 0.98  # water mass % in distillate
x_B = 0.05  # water mass % in bottom

# temperatures
T_F = 25       # feed temperature (C)
T_D = 95       # distillate temperature (C) # close to the boiling point of the light key (water 100 C)
T_B = 120      # bottom temperature (C)     # close to the boiling point of the heavy key (propionic acid, 141 C)

# vapour pressure LK water
# https://www.engineeringtoolbox.com/water-vapor-saturation-pressure-d_599.html
VP_LK = 1.96 * 1.0135 # bar at 120 °C

# vapour pressure HK propionate (+ acetate)
# https://webbook.nist.gov/cgi/cbook.cgi?ID=C79094&Mask=4&Type=ANTOINE&Plot=on#ANTOINE
T = 120 + 273 # temperature in K near distillate temperature 120 °C
VP_HK_log10 = 4.74558 - (1679.869 / (T -59.832))
VP_HK = 10**VP_HK_log10 # in bar.

# vapor presure ratio
alfa = VP_LK / VP_HK # vapour Pressure ratio

# vapour enthalpies heavy-key light-key
Hvap_LK = 46/72*1e3 # kJ/kg # https://www.engineeringtoolbox.com/water-properties-d_1573.html?vA=120&units=C#
Hvap_HK = 3774 # kJ/kg  # https://webbook.nist.gov/cgi/cbook.cgi?ID=C79094&Mask=4

# heat capacities
Cp_LK = 4.184 # (kJ/K/kg) # water
Cp_HK = 2.334 # (kJ/K/kg) # propionate

powerConsumption =[]
seperationCoef_LK_Dist = []
seperationCoef_LK_Btm = []
for xf in x_F:
    Q, sepeationCoef = simulate_distilation(x_D= x_D, x_B= x_B, F= F, x_F= xf, alfa_f= alfa,          # for mass balances
                         Hvap_LK= Hvap_LK, Hvap_HK= Hvap_HK,                        # for condenser duty
                         T_F= T_F, T_D=T_D, T_B=T_B, Cp_LK= Cp_LK, Cp_HK= Cp_LK, printResults=False)    # for reboiler duty
    powerConsumption.append(Q)
    seperationCoef_LK_Dist.append(sepeationCoef[0])
    seperationCoef_LK_Btm.append(sepeationCoef[1])

# plot the % of water seperated to the distilation stream (compared to the feed)
# import matplotlib.pyplot as plt
# # create scatter plot
# # create figure and axes objects
# fig, ax1 = plt.subplots()
#
# # plot data on primary axis
# ax1.scatter(x_F, seperationCoef_LK_Dist, color='tab:red')
# ax1.set_xlabel('X Data')
# ax1.set_ylabel('water in dist %', color='tab:red')
# ax1.tick_params(axis='y', labelcolor='tab:red')
#
# # create secondary axis and plot data on it
# ax2 = ax1.twinx()
# ax2.scatter(x_F, seperationCoef_LK_Btm, color='tab:blue')
# ax2.set_ylabel('water in btm %', color='tab:blue')
# ax2.tick_params(axis='y', labelcolor='tab:blue')
#
# plt.title('Scatter Plot of X and Y Data')
# # show plot
# plt.show()


# fit model



def triangulation(dictData):
    data = dictData["X"]
    prediction = dictData["Y"]
    xy = np.hstack((data,prediction))

    # triangulate input data with Delaunay function, obtain simplices and convex hull
    tridata = Delaunay(xy)
    dictData["triangulation_0"] = tridata

    # assign vertices to dictData
    dictData["vindices_0"] = list(range(len(data)))
    dictData["vertices_0"] = data
    dictData["yvertices_0"] = xy

    # assign simplices to dictData
    dictData["sindices_0"] = tridata.simplices

    # predicion values
    dictData["prediction_0"] = prediction

    return dictData

# reshape the data
X = x_F.reshape(-1,1)
Y = np.array(powerConsumption)
Y = Y.reshape(-1,1)

dictData = {"X": X, "Y": Y}
dataOut = triangulation(dictData=dictData)

print(dataOut)

# # Split data into training and testing sets
# X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.5, random_state=42)
#
# # Create DTR model with default parameters
# dtr = DecisionTreeRegressor()
#
# # Fit DTR model to training data
# dtr.fit(X_train, y_train)
#
# # Make predictions on testing data
# y_pred = dtr.predict(X_test)
#
# # Evaluate the model's performance on testing data
# mse = mean_squared_error(y_test, y_pred)
# print("Mean Squared Error:", mse)
#

# # Plot data and predictions
# plt.scatter(X_test, y_test, color='black')
# plt.scatter(X_test, y_pred, color='blue', linewidth=3)
# plt.xlabel('X')
# plt.ylabel('Y')
# plt.title('Decision Tree Regressor')
# plt.show()