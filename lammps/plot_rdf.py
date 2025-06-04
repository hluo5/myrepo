import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Set global font to Arial
mpl.rcParams['font.family'] = 'Arial'

# Read the data file
file_path = 'npt.dump.gofr.dat'  # Replace with your file path

# Load the data into a structured array, skipping the header row
data = np.loadtxt(file_path, skiprows=1)

# Extract distance (first column) and RDF data (columns 2, 4, ..., 34 for different pairs)
distances = data[:, 0]
rdf_values = data[:, 1::2]  # All RDF values (columns 1, 3, 5, ..., 33)

# Extract coordination numbers (the Int values)
coordination_numbers = data[:, 2::2]  # Columns 2, 4, 6, ..., 34 for coordination numbers

# Create figure for RDFs
#plt.figure(figsize=(10, 6))
#for i, label in enumerate([
#    "Mg-Mg", "Mg-Si", "Mg-O", "Mg-Fe", "Si-Mg", "Si-Si", "Si-O", "Si-Fe", 
#    "O-Mg", "O-Si", "O-O", "O-Fe", "Fe-Mg", "Fe-Si", "Fe-O", "Fe-Fe"]):
#    plt.plot(distances[distances <= 8.0], rdf_values[distances <= 8.0, i], label=label)

# List of all pairs based on your column pattern
# For example, if your columns follow a known pattern like "Mg-Mg", "Mg-Si", "Mg-O", etc.
# You could define the pairs dynamically based on your column structure.
all_pairs = [
    "Mg-Mg", "Mg-Si", "Mg-O", "Mg-Fe", "Si-Mg", "Si-Si", "Si-O", "Si-Fe",
    "O-Mg", "O-Si", "O-O", "O-Fe", "Fe-Mg", "Fe-Si", "Fe-O", "Fe-Fe"
]

# Function to get pair indices
def get_pair_indices(pairs_to_plot):
    if pairs_to_plot == 'all':
        return range(len(all_pairs))  # Return all indices if 'all' is specified
    else:
        pair_indices = []
        for pair in pairs_to_plot:
            if pair in all_pairs:
                pair_indices.append(all_pairs.index(pair))
            else:
                print(f"Warning: '{pair}' not found in the list of available pairs.")
        return pair_indices

# Input your pairs to plot here, or use 'all' to plot everything
#pairs_to_plot = 'all'  # Modify this to 'all' or specify pairs like ["Mg-O", "Si-O"]
pairs_to_plot = ["Mg-O", "Si-O", "O-O", "Fe-O", "Fe-Fe"]
pair_indices = get_pair_indices(pairs_to_plot)

# Create figure for RDFs within 8 Å for selected pairs
plt.figure(figsize=(10, 6))
for idx in pair_indices:
    plt.plot(distances[distances <= 10.0], rdf_values[distances <= 10.0, idx], label=all_pairs[idx])

plt.xlabel('r (Å)', fontsize=14)
plt.ylabel('g(r)', fontsize=14)
plt.title('Partial Radial Distribution Functions', fontsize=16)
plt.legend(loc='upper right', fontsize=14)
plt.xticks(fontsize=12)  # Increase x-axis tick number size
plt.yticks(fontsize=12)  # Increase y-axis tick number size
plt.grid(True)
plt.tight_layout()
plt.savefig("RDF.pdf", format="pdf", dpi=450)

# Create figure for CNs
#plt.figure(figsize=(10, 6))
#for i, label in enumerate([
#    "Mg-Mg", "Mg-Si", "Mg-O", "Mg-Fe", "Si-Mg", "Si-Si", "Si-O", "Si-Fe", 
#    "O-Mg", "O-Si", "O-O", "O-Fe", "Fe-Mg", "Fe-Si", "Fe-O", "Fe-Fe"]):
#    plt.plot(distances[distances <= 4.0], coordination_numbers[distances <= 4.0, i], label=label)

# Create figure for CNs within 8 Å for selected pairs
plt.figure(figsize=(10, 6))
for idx in pair_indices:
    plt.plot(distances[distances <= 4.0], coordination_numbers[distances <= 4.0, idx], label=all_pairs[idx])

plt.xlabel('r (Å)', fontsize=14)
plt.ylabel('Coordination Number', fontsize=14)
plt.title('Coordination Numbers', fontsize=16)
plt.legend(loc='upper right', fontsize=14)
plt.xticks(fontsize=12)  # Increase x-axis tick number size
plt.yticks(fontsize=12)  # Increase y-axis tick number size
plt.grid(True)
plt.tight_layout()
plt.savefig("CN.pdf", format="pdf", dpi=450)

# Show the plots
plt.show()

