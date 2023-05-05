import numpy as np
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator
import matplotlib.pyplot as plt

# Generate hypothetical data
x = np.random.rand(50) * 10
y = np.sin(x) + np.random.randn(50) * 0.1

# Perform Delaunay triangulation
tri = Delaunay(np.column_stack((x, y)), qhull_options='QJ10')

# Create interpolation function
interp = LinearNDInterpolator(tri, y)

# Generate predictions
x_pred = np.linspace(0, 10, 100)
y_pred = interp(x_pred, np.sin(x_pred))

# Plot data and predictions
plt.scatter(x, y, color='blue', label='data')
plt.plot(x_pred, y_pred, color='red', label='predictions')
plt.legend()
plt.show()
