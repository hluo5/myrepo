import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# Define the files and their labels
files = {
    "Mg_Vol.txt": ("Mg", 0.0),
    "Si_Vol.txt": ("Si", 0.0),
    "O_Vol.txt": ("O", 0.0),
    "Fe_Vol.txt": ("Fe", 0.0)
}

# Create a figure with 2x2 subplots
fig, axs = plt.subplots(2, 2, figsize=(12, 10))

# Flatten axes for easy iteration
axs = axs.flatten()

for i, (filename, (label, subtract_value)) in enumerate(files.items()):
    # Load the data
    data = np.loadtxt(filename)
    
    # Subtract (reverse direction)
    adjusted_data = data - subtract_value
    
    # Use Gaussian KDE for smooth curve
    kde = gaussian_kde(adjusted_data)
    
    # Create a range for x values
    x_vals = np.linspace(adjusted_data.min(), adjusted_data.max(), 500)
    
    # Plot on the corresponding subplot
    axs[i].plot(x_vals, kde(x_vals), color='black')
    axs[i].set_title(f'{label}')
    axs[i].set_xlabel('Bader volume')
    axs[i].set_ylabel('Abundance')
    axs[i].grid(True)

# Adjust layout with more vertical spacing
plt.tight_layout(h_pad=2.5)

# Save the whole figure as a single PDF
plt.savefig("Volume_Distributions.pdf", dpi=450, format='pdf')

# Show the figure
plt.show()

