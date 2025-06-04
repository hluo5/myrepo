import numpy as np
import pandas as pd

# Boltzmann constant in J/mol/K (scaled for convenience in mixing entropy calculation)
k_B = 1.380649e-23 * 6.02214076e23  # J/mol/K
r_cut = 8.0

# Load the RDF data from the file
filename = "npt.dump.gofr.dat"
data = pd.read_csv(filename, delim_whitespace=True)

# Extract the RDFs and the distance values
dist = data['dist'].values  # Distance (r) values
rdf_dict = {
    'Mg-Mg': data['Mg-Mg'].values,
    'Mg-Si': data['Mg-Si'].values,
    'Mg-O': data['Mg-O'].values,
    'Mg-Fe': data['Mg-Fe'].values,
    'Si-Si': data['Si-Si'].values,
    'Si-O': data['Si-O'].values,
    'Si-Fe': data['Si-Fe'].values,
    'O-O': data['O-O'].values,
    'O-Fe': data['O-Fe'].values,
    'Fe-Fe': data['Fe-Fe'].values
}

# Define mole fractions for each species (these need to be provided or calculated)
# Example mole fractions (you should replace them with the actual mole fractions for your system)
x_Mg = 0.1764706
x_Si = 0.1617647
x_O = 0.5
x_Fe = 0.1617647

# Function to calculate mixing entropy
def calculate_mixing_entropy(rdf_dict, dist, x_Mg, x_Si, x_O, x_Fe):
    # Initialize mixing entropy
    S_mix = 0.0
    
    # List of species pairs and their corresponding mole fractions
    species_pairs = [
        ('Mg-Mg', x_Mg, x_Mg), ('Mg-Si', x_Mg, x_Si), ('Mg-O', x_Mg, x_O), ('Mg-Fe', x_Mg, x_Fe),
        ('Si-Si', x_Si, x_Si), ('Si-O', x_Si, x_O), ('Si-Fe', x_Si, x_Fe),
        ('O-O', x_O, x_O), ('O-Fe', x_O, x_Fe), ('Fe-Fe', x_Fe, x_Fe)
    ]
    
    # Iterate over all species pairs to calculate their contribution to the mixing entropy
    for pair, x_i, x_j in species_pairs:
        g_ij = rdf_dict[pair]  # Partial RDF for pair (i, j)
        
        # Add a small epsilon to avoid log(0) and prevent errors
        epsilon = 1e-10
        g_ij_safe = g_ij + epsilon
        
        # Calculate the integral of g_ij(r) * ln(g_ij(r)) over r
        mask  = dist <= r_cut
        integrand = g_ij_safe[mask] * np.log(g_ij_safe[mask])
        integral = np.trapz(integrand, dist[mask])  # Trapezoidal integration over the distance r
        # Add the contribution to the mixing entropy
        S_mix += x_i * x_j * integral
    
    # Multiply by -k_B to get the entropy in J/mol/K
    S_mix *= -k_B
    return S_mix

# Calculate the mixing entropy for the system
S_mix = calculate_mixing_entropy(rdf_dict, dist, x_Mg, x_Si, x_O, x_Fe)

# Print the mixing entropy
print(f"Mixing Entropy: {S_mix:.4e} J/mol/K")

