import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# Define elements and their subtract values
elements = {
    "Mg_Val.txt": ("Mg", 8.0),
    "Si_Val.txt": ("Si", 4.0),
    "O_Val.txt": ("O", 6.0),
    "Fe_Val.txt": ("Fe", 14.0)
}

# Define temperatures and their corresponding folders
temperature_folders = {
    "3500K": "/scratch/gpfs/JIEDENG/Haiyang/projects/elem_part/superEarth/3500K/DFT_MD/dry_544atoms/544atoms_run3/bader",
    "6500K": "/scratch/gpfs/JIEDENG/Haiyang/projects/elem_part/superEarth/6500K/DFT_MD/dry_544atoms/544atoms_run3/bader",
    "9000K": "/scratch/gpfs/JIEDENG/Haiyang/projects/elem_part/superEarth/9000K/DFT_MD/dry_544atoms/544atoms_run3/bader",
    "13000K": "/scratch/gpfs/JIEDENG/Haiyang/projects/elem_part/superEarth/13000K/DFT_MD/dry_544atoms/544atoms_run2/bader",
    "15000K": "/scratch/gpfs/JIEDENG/Haiyang/projects/elem_part/superEarth/15000K/DFT_MD/dry_544atoms/544atoms_run2/bader"
}

# Set up the figure
fig, axs = plt.subplots(2, 2, figsize=(14, 12))
axs = axs.flatten()

# Colors for temperatures (you can adjust them)
colors = ['black', 'red', 'blue', 'green', 'orange']

for idx, (filename, (label, subtract_value)) in enumerate(elements.items()):
    ax = axs[idx]
    
    for (temp, folder), color in zip(temperature_folders.items(), colors):
        file_path = f"{folder}/{filename}"
        
        try:
            data = np.loadtxt(file_path)
            
            # Subtract (reverse direction)
            adjusted_data = subtract_value - data
            
            # KDE
            kde = gaussian_kde(adjusted_data)
            x_vals = np.linspace(adjusted_data.min(), adjusted_data.max(), 500)
            
            # Plot
            ax.plot(x_vals, kde(x_vals), label=temp, color=color)
        
        except Exception as e:
            print(f"Failed to load {file_path}: {e}")
    
    ax.set_title(label, fontsize=14)
    ax.set_xlabel('Bader charge', fontsize=12)
    ax.set_ylabel('Abundance', fontsize=12)
    ax.grid(True)
    ax.legend(fontsize=10)

# Adjust spacing
plt.tight_layout(h_pad=3.0)

# Save the whole figure
plt.savefig("Charge_Distributions_with_Temperatures.pdf", dpi=450, format='pdf')

# Show plot
plt.show()

