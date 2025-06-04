import matplotlib.pyplot as plt
import numpy as np

# Number of atoms in the system
num_atoms = 544 # Number of atoms in the system
T = 6500 # Temperature in K
vol = 3286.065 # System volume in Å³
eV_A3_2_GPa  = 160.21766208 # 1 eV/Å3 = 160.2176621 GPa


# Define file names and corresponding labels
files = {
    f"{T}K_result.e.out": ("DFT Energy (eV)", "DP Energy (eV)", "eV", "eV/atom"),
    f"{T}K_result.f.out": ("DFT Forces (eV/Å)", "DP Forces (eV/Å)", "eV/Å"),
    f"{T}K_result.v.out": ("DFT Stress (GPa)", "DP Stress (GPa)", "GPa")
}

# Initialize a figure for saving
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Loop through each file and process the data
for i, (file_name, labels) in enumerate(files.items()):
    true_values = []
    pred_values = []
    unit = labels[2]  # Extract unit for RMSE

    # Read the file
    with open(file_name, "r") as file:
        lines = file.readlines()

    # Process lines (skip the first line if it's a header)
    for line in lines[1:]:  # Assuming the first line is a header
        values = list(map(float, line.split()))
        num_values = len(values) // 2  # Number of true/predicted pairs

        true_values.extend(values[:num_values])  # First half: true values
        pred_values.extend(values[num_values:])  # Second half: predicted values

    # Convert to NumPy arrays
    true_values = np.array(true_values)
    pred_values = np.array(pred_values)

    # Compute RMSE
    rmse = np.sqrt(np.mean((true_values - pred_values) ** 2))

    # Adjust RMSE for energy (convert to eV/atom)
    if file_name == f"{T}K_result.e.out":
        rmse /= num_atoms
        unit = labels[3]  # Use eV/atom for energy

    # Adjust RMSE for stress (convert to GPa)
    if file_name == f"{T}K_result.v.out":
        true_values = (true_values / vol) * eV_A3_2_GPa
        pred_values = (pred_values / vol) * eV_A3_2_GPa
        rmse = np.sqrt(np.mean((true_values - pred_values) ** 2))
        unit = labels[2]  # Use GPa for virial

    # Define min and max for a consistent y=x line across the entire range
    min_val = min(true_values.min(), pred_values.min())
    max_val = max(true_values.max(), pred_values.max())

    # Scatter plot
    ax = axes[i]
    ax.scatter(true_values, pred_values, alpha=0.5, label=f"Data vs Prediction (RMSE = {rmse:.3e} {unit})")
    ax.plot([min_val, max_val], [min_val, max_val], 'k--', label="y = x")  # Ensures full-range y=x line

    # Set axis limits to match full data range
    ax.set_xlim(min_val, max_val)
    ax.set_ylim(min_val, max_val)

    # Labels and title
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    #ax.set_title(f"Comparison: {file_name}")
    ax.legend()
    #ax.grid(True)

# Save figure as PDF
plt.tight_layout()
plt.savefig("comparison_plots.pdf", dpi=600)
plt.show()

