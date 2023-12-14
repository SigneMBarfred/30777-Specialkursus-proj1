#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

# Green points
x_green = np.array([-3.02, -3.2, -3.5, -3.75, -3.95, -4.2, -4.4, -4.65, -4.9, -5.1, -5.35, -5.55, -5.85, -6.05, -6.28])
y_green = np.array([1.3, 1.75, 2.2, 2.55, 3.2, 3.55, 3.95, 4.2, 4.45, 4.6, 4.75, 4.9, 5.05, 5.2, 5.4])
# Remove duplicate
x_green, indices = np.unique(x_green, return_index=True)
y_green = y_green[indices]

# Polynomial fit
degree = 10  # You can adjust the degree of the polynomial
coefficients = np.polyfit(x_green, y_green, degree)
poly_fit = np.poly1d(coefficients)

# Plot the polynomial fit
x_fit = np.linspace(min(x_green), max(x_green), 1000)
y_fit = poly_fit(x_fit)

plt.figure()
plt.plot(x_green, y_green, 'go-', label='FR06 = Franceschini et al. 2006', zorder=1)
plt.plot(x_fit, y_fit, label=f'Polynomial Fit (Degree {degree})', color='red', zorder=2)
plt.xlabel(r'$\log(F_{3.6}\, \mu\mathrm{m}\, [\mathrm{Jy}])$')
plt.ylabel(r'$\log(\#(>F_{3.6}\, \mu\mathrm{m}))\, [\mathrm{deg}^{-2}]$')
plt.legend()
plt.show()

# Load the interpolated flux values file
input_file_interpolated = 'C:/Users/signe/Documents/DTU/Specialkursus_radiogal/flux_values_interpolated.txt'
