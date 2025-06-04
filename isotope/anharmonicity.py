import matplotlib.pyplot as plt
import numpy as np

# Define the Lennard-Jones potential function
def lennard_jones(r, epsilon, sigma):
    return 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)

# Define the range for r
r = np.linspace(0.8, 3, 400)

# Parameters for low pressure (higher anharmonicity)
epsilon_low = 1.0
sigma_low = 1.0

# Parameters for high pressure (lower anharmonicity)
epsilon_high = 1.5
sigma_high = 0.9

# Calculate Lennard-Jones potential for low and high pressure
lj_low = lennard_jones(r, epsilon_low, sigma_low)
lj_high = lennard_jones(r, epsilon_high, sigma_high)

# Create the plot
plt.figure(figsize=(10, 6))

# Plot low pressure Lennard-Jones potential curve
plt.plot(r, lj_low, label='Low Pressure (Higher Anharmonicity)', color='blue')

# Plot high pressure Lennard-Jones potential curve
plt.plot(r, lj_high, label='High Pressure (Lower Anharmonicity)', color='red')

# Adding labels and legend
plt.xlabel('Interatomic Distance (r)')
plt.ylabel('Potential Energy (V)')
plt.title('Lennard-Jones Potential and Effect of Pressure on Anharmonicity')
plt.legend()

# Show plot
#plt.grid(True)
plt.ylim(-2, 2)
plt.xlim(0.8, 3)

# Save the figure with 600 dpi
plt.savefig("Lennard_Jones_Potential_Pressure_Effect.png", dpi=600)

# Display the plot
plt.show()
