import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# Define the files and their labels
files = {
    "Mg_Val.txt": ("Mg", 8.0),
    "Si_Val.txt": ("Si", 4.0),
    "O_Val.txt": ("O", 6.0),
    "Fe_Val.txt": ("Fe", 14.0)
}

# Create a figure with 2x2 subplots
fig, axs = plt.subplots(2, 2, figsize=(12, 10))

# Flatten axes for easy iteration
axs = axs.flatten()

for i, (filename, (label, subtract_value)) in enumerate(files.items()):
    # Load the data
    data = np.loadtxt(filename)
    
    # Subtract (reverse direction)
    adjusted_data = subtract_value - data
    
    # Use Gaussian KDE for smooth curve
    kde = gaussian_kde(adjusted_data)
    
    # Create a range for x values
    x_vals = np.linspace(adjusted_data.min(), adjusted_data.max(), 500)
    
    # Plot on the corresponding subplot
    axs[i].plot(x_vals, kde(x_vals), color='black')
    axs[i].set_title(f'{label}')
    axs[i].set_xlabel('Bader Charge')
    axs[i].set_ylabel('Abundance')
    axs[i].grid(True)

# Adjust layout with more vertical spacing
plt.tight_layout(h_pad=2.5)

# Save the whole figure as a single PDF
plt.savefig("Charge_Distributions.pdf", dpi=450, format='pdf')

# Show the figure
plt.show()

