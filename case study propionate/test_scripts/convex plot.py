import numpy as np
import matplotlib.pyplot as plt

# Define the convex function: f(x) = x^2
def convex_function(x):
    return x**2

# Define the non-convex function: f(x) = sin(x) + 0.1x^2
def nonconvex_function(x):
    return np.sin(x) + 0.1 * x**2

# Create an array of x-values to plot the functions
x_vals = np.linspace(-10, 10, 200)

# Plot the convex function
plt.plot(x_vals, convex_function(x_vals), label='Convex function')

# Plot the non-convex function
plt.plot(x_vals, nonconvex_function(x_vals), label='Non-convex function')

# Add titles and labels
plt.title('Convex and non-convex function graphs')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend()

# Display the plot
plt.show()
