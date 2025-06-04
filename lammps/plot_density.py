import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl

# Set global font to Arial
mpl.rcParams['font.family'] = 'Arial'

# Constants
N_A = 6.022e23  # Avogadro's number (atoms/mol)
angstrom_to_cm = 1e-8  # Conversion factor from Angstroms to cm

# Load the LAMMPS dump file with explicit topology format
u = mda.Universe('npt.dump', format='LAMMPSDUMP')

# Load the LAMMPS dump file again manually to read BOX BOUNDS
dump_file = 'npt.dump'

# Define atom type mapping (based on your LAMMPS atom types)
atom_type_mapping = {1: 'Mg', 2: 'Si', 3: 'O', 4: 'Fe'}  # Adjust based on your file's atom type IDs
molar_masses = {'Mg': 24.305, 'Si': 28.085, 'O': 16.000, 'Fe': 55.847}

# Extract positions (coordinates) and atom types
atom_types = u.atoms.types  # Get atom types as numeric IDs (1, 2, 3, 4, etc.)
atom_names = [atom_type_mapping[int(atom_id)] for atom_id in atom_types]

# Exclude initial steps (e.g., 100 steps)
start_step = 440
end_step = u.trajectory.n_frames  # You can adjust this depending on the number of frames
num_bins = 5  # You can adjust this based on your needs

# Create a dictionary to store density for each atom type
density = {'Mg': np.zeros((num_bins, num_bins, num_bins)),
           'Si': np.zeros((num_bins, num_bins, num_bins)),
           'O': np.zeros((num_bins, num_bins, num_bins)),
           'Fe': np.zeros((num_bins, num_bins, num_bins))}

# Loop through the frames and calculate the density for each atom type
n_frames = 0
# --- read the box bounds separately ---
with open(dump_file, 'r') as f:
    lines = f.readlines()

frame_idx = -1
line_idx = 0
while line_idx < len(lines):
    if lines[line_idx].startswith('ITEM: TIMESTEP'):
        timestep = int(lines[line_idx + 1].strip())
        frame_idx += 1
        if frame_idx < start_step:
            # Skip lines for this frame
            line_idx += 9 + int(lines[line_idx + 3].strip())
            continue

        # Get box bounds
        box_bounds_line = line_idx + 5
        xlo, xhi = map(float, lines[box_bounds_line + 0].strip().split())
        ylo, yhi = map(float, lines[box_bounds_line + 1].strip().split())
        zlo, zhi = map(float, lines[box_bounds_line + 2].strip().split())
        print(zlo, zhi)

        # Now, MDAnalysis will be synchronized at this frame
        ts = u.trajectory[frame_idx]
        positions = u.atoms.positions  # Get the positions of all atoms
        atom_types = atom_names  # Get atom types (Mg, Si, O, Fe)

        # Create dynamic bin edges for this frame
        x_bins = np.linspace(xlo, xhi, num_bins + 1)
        y_bins = np.linspace(ylo, yhi, num_bins + 1)
        z_bins = np.linspace(zlo, zhi, num_bins + 1)

        # Bin indices
        x_bin_indices = np.digitize(positions[:, 0], x_bins) - 1
        y_bin_indices = np.digitize(positions[:, 1], y_bins) - 1
        z_bin_indices = np.digitize(positions[:, 2], z_bins) - 1

        # Update density by incrementing counts in the appropriate bins for each atom type
        for i, j, k, atom in zip(x_bin_indices, y_bin_indices, z_bin_indices, atom_types):
            if 0 <= i < num_bins and 0 <= j < num_bins and 0 <= k < num_bins:
                if atom in density:
                    density[atom][i, j, k] += 1
        n_frames += 1

        # Skip to next frame
        line_idx += 9 + int(lines[line_idx + 3].strip())
    else:
        line_idx += 1
print(f"Processed {n_frames} frames.")

density = {atom: density[atom] / n_frames for atom in density}

bin_volume = ((xhi - xlo) * (yhi - ylo) * (zhi - zlo)) * (angstrom_to_cm**3) / (num_bins**3)

# Correct calculation: sum the mass contributions from all atom types in each bin
average_mass_density = np.zeros((num_bins, num_bins, num_bins))  # To store total mass density in g/cm³
for atom in density:
    # Convert molar mass to atomic mass (grams per atom)
    atomic_mass = molar_masses[atom] / N_A  # g/atom
    # Sum the mass densities for each atom type in each bin
    average_mass_density += density[atom] * atomic_mass / bin_volume  # Add the contribution for each atom type

# Plot the 3D density profile for the average density
fig = plt.figure(figsize=(12, 10))

# Create a meshgrid for the bins
x_grid, y_grid, z_grid = np.meshgrid(x_bins[:-1] + (x_bins[1] - x_bins[0]) / 2,
                                     y_bins[:-1] + (y_bins[1] - y_bins[0]) / 2,
                                     z_bins[:-1] + (z_bins[1] - z_bins[0]) / 2)

# Plot a 3D scatter plot for average density (points will be plotted at the center of each bin)
ax = fig.add_subplot(111, projection='3d')
sc = ax.scatter(x_grid.flatten(), y_grid.flatten(), z_grid.flatten(),
                c=average_mass_density.flatten(), cmap='viridis', s=50)

# Labels and title
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Average Density Profile (g/cm³)')

# Add a color bar
plt.colorbar(sc)

plt.tight_layout()
plt.savefig("Density grid.pdf", format="pdf", dpi=450)
plt.show()

mid_x = num_bins // 2
plt.imshow(average_mass_density[mid_x, :, :], origin='lower', cmap='viridis',
           extent=[y_bins[0], y_bins[-1], z_bins[0], z_bins[-1]])
plt.colorbar(label='Density (g/cm³)')
plt.xlabel('Y (Å)')
plt.ylabel('Z (Å)')
plt.title('Density slice at middle X-plane')
plt.tight_layout()
plt.savefig("Density slice.pdf", format="pdf", dpi=450)
plt.show()
