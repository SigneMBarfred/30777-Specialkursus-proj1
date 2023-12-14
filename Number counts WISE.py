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
# spline fit
spline = UnivariateSpline(x_green, y_green, s=0)

x_spline = np.linspace(min(x_green), max(x_green), 1000)
y_spline = spline(x_spline)

plt.figure()
plt.plot(x_green, y_green, 'go-', label='FR06 = Franceschini et al. 2006', zorder=1)
plt.plot(x_spline, y_spline, label='Spline Fit', color='red', zorder=2)
plt.xlabel(r'$\log(F_{3.6}\, \mu\mathrm{m}\, [\mathrm{Jy}])$')
plt.ylabel(r'$\log(\#(>F_{3.6}\, \mu\mathrm{m}))\, [\mathrm{deg}^{-2}]$')
plt.legend()
plt.show()

# Read the flux values from the input file
input_file = 'C:/Users/olive/Dropbox/5. semester/Special kursus/flux_values.txt'

data = np.genfromtxt(input_file, skip_header=1, dtype=float)

# Extract flux values and take the log10
flux_values_data = np.log10(data[:, 1:5])  # Assuming columns 1 to 4 are flux_W1 to flux_W4 and taking the log10
print(flux_values_data[:,3])
 
# Extract unique obj_id values
obj_ids = data[:, 0]

# Use the spline to interpolate the y value
y_value_interpolated = 10**spline(flux_values_data)

# Save y_value_interpolated to a text file
output_file_y_values = 'ns.txt'
header_line_y_values = 'obj_id ' + ' '.join([f'Flux_W{i+1}' for i in range(y_value_interpolated.shape[1])])
np.savetxt(output_file_y_values, np.column_stack([obj_ids] + [y_value_interpolated[:, i] for i in range(y_value_interpolated.shape[1])]), fmt='%d ' + ' '.join(['%.8f' for _ in range(y_value_interpolated.shape[1])]), header=header_line_y_values, comments='')


# Find the index of the minimum y_value_interpolated along each row
row_indices, column_indices = np.unravel_index(np.argmin(y_value_interpolated, axis=None), y_value_interpolated.shape)

# Extract the minimum value, obj_id, and column index
# Find the minimum y_value_interpolated along each column
nlim = np.min(y_value_interpolated, axis=0)
obj_id_nlim = obj_ids[np.argmin(y_value_interpolated, axis=0)]  # Get obj_id corresponding to the minimum value in each column
column_index_nlim = np.arange(1, y_value_interpolated.shape[1] + 1)  # Column indices starting from 1

print(f'Lowest y values: {nlim}')
print(f'obj_ids: {obj_id_nlim}')
print(f'Column indices: {column_index_nlim}')

#print(f'For x={flux_values_data}, interpolated y={y_value_interpolated}')

#plt.figure()
#plt.plot(flux_values_data, y_value_interpolated, 'bo', label='Interpolated Points', zorder=1)
#colors = ['red', 'blue', 'black', 'orange']
#for i in range(4):
#    indices = np.arange(i, len(flux_values_data), 4)
#    plt.plot(flux_values_data[indices], y_value_interpolated[indices], 'o', color=colors[i], zorder=i+2)
#plt.plot(x_green, y_green, 'go-', label='FR06 = Franceschini et al. 2006', zorder=6)
#plt.plot(x_spline, y_spline, label='Spline Fit', color='red', zorder=7)
#legend_entries = [
#    plt.Line2D([0], [0], marker='o', color='w', label='Flux_W1', markerfacecolor='red', markersize=10),
#    plt.Line2D([0], [0], marker='o', color='w', label='Flux_W2', markerfacecolor='blue', markersize=10),
#    plt.Line2D([0], [0], marker='o', color='w', label='Flux_W3', markerfacecolor='black', markersize=10),
#    plt.Line2D([0], [0], marker='o', color='w', label='Flux_W4', markerfacecolor='orange', markersize=10),
#    plt.Line2D([0], [0], marker='o', color='w', label='FR06 = Franceschini et al. 2006', markerfacecolor='green', markersize=10),
#    plt.Line2D([0], [0], marker='', linestyle='-', color='red', label='Spline Fit', markersize=10),
#]
#plt.xlabel(r'$\log(F_{3.6}\, \mu\mathrm{m}\, [\mathrm{Jy}])')
#plt.ylabel(r'$\log(\#(>F_{3.6}\, \mu\mathrm{m}))\, [\mathrm{deg}^{-2}]$')
#plt.legend(handles=legend_entries)
#plt.show()

# Load angular distances from the text file
input_file_radius = 'C:/Users/olive/Dropbox/5. semester/Special kursus/angular_distances_without_16.2.txt'
data_radius = np.genfromtxt(input_file_radius, skip_header=1, dtype=float, delimiter='\t', missing_values='', filling_values=np.nan, invalid_raise=False)

# Filter out rows with unexpected number of columns
data_radius = data_radius[~np.isnan(data_radius).any(axis=1)]

# Extract object IDs and angular distances
obj_ids_radius = data_radius[:, 0].astype(int)
angular_distances_str = data_radius[:, 1]

# Ensure the same number of rows for both arrays
min_rows = min(len(angular_distances_str), len(y_value_interpolated))
angular_distances_str = angular_distances_str[:min_rows]
y_value_interpolated = y_value_interpolated[:min_rows, :]

# Calculate critical probability (constant for all columns)
r_s = 15.0 / 3600.0
p_crit = np.pi * r_s**2 * nlim

# Initialize arrays to store mu_r, Poisson probabilities, and additional calculations for each column
mu_r_values = np.zeros_like(y_value_interpolated)
poisson_probabilities = np.zeros_like(y_value_interpolated)
mu_cor_values = np.zeros_like(y_value_interpolated)
final_prob_values = np.zeros_like(y_value_interpolated)

# Calculate mu_r, Poisson probabilities, and additional calculations for each column
for i in range(y_value_interpolated.shape[1]):
    mu_r_values[:, i] = np.pi * angular_distances_str**2 * y_value_interpolated[:, i]

    poisson_probabilities[:, i] = 1 - np.exp(-mu_r_values[:, i])
    
    # Calculate corrected mu
    mu_cor_values[:, i] = poisson_probabilities[:, i] * (1 + np.log(p_crit[i] / poisson_probabilities[:, i]))
    
    # Calculate final probability
    final_prob_values[:, i] = 1 - np.exp(-mu_cor_values[:, i])

# Save obj_id, mu_r, Poisson probabilities, and additional calculations to a text file
output_file = 'additional_calculations.txt'
header_line = 'obj_id ' + ' '.join([f'mu_r_{i+1} poisson_probability_{i+1} mu_cor_{i+1} final_prob_{i+1}' for i in range(y_value_interpolated.shape[1])])
np.savetxt(output_file, np.column_stack([obj_ids_radius] + [mu_r_values[:, i] for i in range(y_value_interpolated.shape[1])] + [poisson_probabilities[:, i] for i in range(y_value_interpolated.shape[1])] + [mu_cor_values[:, i] for i in range(y_value_interpolated.shape[1])] + [final_prob_values[:, i] for i in range(y_value_interpolated.shape[1])]), fmt='%d ' + ' '.join(['%.8f %.8f %.8f %.8f' for _ in range(y_value_interpolated.shape[1])]), header=header_line, comments='')



# values for 16.2
#Flux
#16: Source 2 
Ch1 = 4.095592869919121*10**(-5)
Ch2 = 3.1427833660837705*10**(-5)
Ch3 = 0.0002056637519693072
CH4 = 0.001568780190748341

data_16_2 = np.array(np.log10([Ch1, Ch2, Ch3, CH4]))
y_value_interpolated_16_2 = 10**spline(data_16_2)
print(y_value_interpolated_16_2)

r=0.002518286671670176
mu_r_16_2 = np.pi * r**2 * y_value_interpolated_16_2
print(mu_r_16_2)
P_16_2 = 1 - np.exp(-mu_r_16_2)
print(P_16_2)
mu_cor_16_2 = P_16_2 * (1 + np.log(p_crit / P_16_2))
print(mu_cor_16_2)
final_prob_16_2 = 1 - np.exp(-mu_cor_16_2)
print(final_prob_16_2)


# %%
