import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.signal import savgol_filter

# Load the data from the file
filename = "npt.dump.gofr.dat"
data = pd.read_csv(filename, delim_whitespace=True)

# Extract relevant columns (dist, Fe-Fe RDF, Int(Fe-Fe))
dist = data['dist'].values  # Convert to numpy array
fe_fe_rdf = data['Fe-Fe'].values  # Convert to numpy array
coordination_number = data['Int(Fe-Fe)'].values  # Assuming 'Int(Fe-Fe)' stores coordination number

# Apply a Savitzky-Golay filter to smooth the data (window size must be odd)
smoothed_rdf = savgol_filter(fe_fe_rdf, window_length=11, polyorder=3)

# Find the global maximum and its value
max_index = np.argmax(smoothed_rdf)
max_value = smoothed_rdf[max_index]
max_dist_value = dist[max_index]

# Get the two nearest values around the maximum (left and right neighbors)
left_idx = max(0, max_index - 1)
right_idx = min(len(smoothed_rdf) - 1, max_index + 1)
neighboring_max_values = smoothed_rdf[[left_idx, max_index, right_idx]]
neighboring_max_distances = dist[[left_idx, max_index, right_idx]]

# Find the first minimum within 5 Å after the maximum
min_after_max_index = None
min_value = None

# Limit the search range to 5 Å after the global maximum
max_dist_value_plus_5 = max_dist_value + 2
if max_dist_value_plus_5 < dist[-1]:  # Ensure that the maximum + 5 Å is within the range
    # Search for the minimum within the range [max_dist_value, max_dist_value + 5]
    valid_range = (dist >= max_dist_value) & (dist <= max_dist_value_plus_5)
    min_after_max_index = np.argmin(smoothed_rdf[valid_range]) + np.where(valid_range)[0][0]
    min_value = smoothed_rdf[min_after_max_index]

# Get the two nearest values around the minimum
left_idx_min = max(0, min_after_max_index - 1)
right_idx_min = min(len(smoothed_rdf) - 1, min_after_max_index + 1)
neighboring_min_values = smoothed_rdf[[left_idx_min, min_after_max_index, right_idx_min]]
neighboring_min_distances = dist[[left_idx_min, min_after_max_index, right_idx_min]]

# Extract the corresponding coordination numbers for the minimum values
coordination_at_min = coordination_number[[left_idx_min, min_after_max_index, right_idx_min]]

# Calculate average and standard deviation of the coordination numbers
avg_coordination = np.mean(coordination_at_min)
std_coordination = np.std(coordination_at_min)

# Print the results for the maximum and minimum values
print("Maximum value and its two nearest values:")
for d, rdf in zip(neighboring_max_distances, neighboring_max_values):
    print(f"Dist: {d}, Fe-Fe RDF: {rdf}")

print("\nMinimum value within 2 Å after the maximum and its two nearest values:")
if min_after_max_index is not None:
    for d, rdf, coord in zip(neighboring_min_distances, neighboring_min_values, coordination_at_min):
        print(f"Dist: {d}, Fe-Fe RDF: {rdf}, Coordination Number: {coord}")
else:
    print("No minimum found within 5 Å after the maximum.")

# Print the average and standard deviation of the coordination number
print("\nCoordination Number Statistics:")
print(f"Average Coordination Number: {avg_coordination:.2f}")
print(f"Standard Deviation of Coordination Number: {std_coordination:.2f}")

# Save the results to a file
output_filename = "rdf_max_min_values_with_coordination.txt"
with open(output_filename, 'w') as f:
    f.write("Maximum Values (Dist, Fe-Fe RDF):\n")
    for d, rdf in zip(neighboring_max_distances, neighboring_max_values):
        f.write(f"{d}\t{rdf}\n")
    f.write("\nMinimum Values (Dist, Fe-Fe RDF, Coordination Number):\n")
    if min_after_max_index is not None:
        for d, rdf, coord in zip(neighboring_min_distances, neighboring_min_values, coordination_at_min):
            f.write(f"{d}\t{rdf}\t{coord}\n")
    f.write("\nCoordination Number Statistics:\n")
    f.write(f"Average Coordination Number: {avg_coordination:.2f}\n")
    f.write(f"Standard Deviation of Coordination Number: {std_coordination:.2f}\n")

# Plot the RDF with markers at the maximum and minimum
plt.plot(dist, fe_fe_rdf, label="Original Fe-Fe RDF", color='black', alpha=0.5)
plt.plot(dist, smoothed_rdf, label="Smoothed Fe-Fe RDF", color='red', linewidth=2)
plt.scatter(neighboring_max_distances, neighboring_max_values, color='blue', label="Maximum", zorder=5)
if min_after_max_index is not None:
    plt.scatter(neighboring_min_distances, neighboring_min_values, color='green', label="Minimum within 2 Å", zorder=5)

# Add vertical lines at the identified maximum and minimum
plt.axvline(x=neighboring_max_distances[1], color='blue', linestyle='--', linewidth=1)
if min_after_max_index is not None:
    plt.axvline(x=neighboring_min_distances[1], color='green', linestyle='--', linewidth=1)

# Add labels and legend
plt.xlabel('Distance (Å)')
plt.ylabel('Fe-Fe RDF')
plt.legend()

# Save the plot as a PNG file
plot_filename = "fe_fe_rdf_with_markers_and_coordination.pdf"
plt.tight_layout()
plt.savefig(plot_filename, dpi=450)
plt.show()

