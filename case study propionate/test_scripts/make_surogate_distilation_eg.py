import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.metrics import r2_score
from make_distilation_model import distillation_check

# Define the function distilation_check()

# Create data
np.random.seed(123)
n_samples = 1000
F_vals = np.random.uniform(100, 500, size=n_samples)
x_F_vals = np.random.uniform(0.2, 0.8, size=n_samples)

# Calculate Q for each set of F and x_F values
Q_vals = []
for F, x_F in zip(F_vals, x_F_vals):
    Q = distillation_check(F=F, x_F=x_F, x_D=0.9, x_B=0.1, alfa_f=1.5,
                           Hvap_LK=100, Hvap_HK=50,
                           T_F=80, T_D=150, T_B=180,
                           Cp_LK=2, Cp_HK=1)
    Q_vals.append(Q)

# Convert the lists to NumPy arrays
F_vals = np.array(F_vals).reshape(-1, 1)
x_F_vals = np.array(x_F_vals).reshape(-1, 1)
Q_vals = np.array(Q_vals).reshape(-1, 1)

# Split the data into training and testing sets
split_ratio = 0.7
split_idx = int(n_samples * split_ratio)
F_train = F_vals[:split_idx]
F_test = F_vals[split_idx:]
x_F_train = x_F_vals[:split_idx]
x_F_test = x_F_vals[split_idx:]
Q_train = Q_vals[:split_idx]
Q_test = Q_vals[split_idx:]

# Fit a linear regression model
lin_reg = LinearRegression()
lin_reg.fit(np.concatenate((F_train, x_F_train), axis=1), Q_train)
Q_linreg_pred_train = lin_reg.predict(np.concatenate((F_train, x_F_train), axis=1))
Q_linreg_pred_test = lin_reg.predict(np.concatenate((F_test, x_F_test), axis=1))

# Fit a polynomial regression model
poly_reg = PolynomialFeatures(degree=2)
X_poly_train = poly_reg.fit_transform(np.concatenate((F_train, x_F_train), axis=1))
X_poly_test = poly_reg.fit_transform(np.concatenate((F_test, x_F_test), axis=1))
poly_reg.fit(X_poly_train, Q_train)
lin_reg_2 = LinearRegression()
lin_reg_2.fit(X_poly_train, Q_train)
Q_polyreg_pred_train = lin_reg_2.predict(X_poly_train)
Q_polyreg_pred_test = lin_reg_2.predict(X_poly_test)

# Create the parity plots
sns.set_style("whitegrid")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

# Plot the linear regression parity plot
sns.scatterplot(x=Q_test.squeeze(), y=Q_linreg_pred_test.squeeze(), ax=ax1)
sns.lineplot(x=Q_test.squeeze(), y=Q_test.squeeze(), color='r', ax=ax1)
ax1.set_xlabel("Observed Qtot")
ax1.set_ylabel("Predicted Qtot")
ax1.set_title("Linear regression")

# Plot the polynomial regression parity plot
sns.scatterplot(x=Q_test.squeeze(), y=Q_polyreg_pred_test.squeeze(), ax=ax2)
sns.lineplot(x=Q_test.squeeze(), y=Q_test.squeeze(), color='r', ax=ax2)
ax2.set_xlabel("Observed Qtot")
ax2.set_ylabel("Predicted Qtot")
ax2.set_title("Polynomial regression")
plt.show()
