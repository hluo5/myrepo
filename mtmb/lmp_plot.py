import lammps_logfile
import matplotlib.pyplot as plt

# Load LAMMPS log file
log = lammps_logfile.File("log.lammps")
N = 1  # Assuming no specific run number is used

# Extract data
step = log.get("Step", run_num=N)
pot_eng = log.get("PotEng", run_num=N)
press = log.get("Press", run_num=N)
temp = log.get("Temp", run_num=N)
volume = log.get("Volume", run_num=N)

# Create a 2x2 grid of subplots
fig, axs = plt.subplots(2, 2, figsize=(10, 8))

# Plot PotEng vs Step
axs[0, 0].plot(step, pot_eng, label="PotEng", color="blue")
axs[0, 0].set_title("Potential Energy (PotEng)")
axs[0, 0].set_xlabel("Step")
axs[0, 0].set_ylabel("PotEng (eV)")

# Plot Press vs Step
axs[0, 1].plot(step, press, label="Press", color="red")
axs[0, 1].set_title("Pressure (Press)")
axs[0, 1].set_xlabel("Step")
axs[0, 1].set_ylabel("Press (bar)")

# Plot Volume vs Step
axs[1, 0].plot(step, volume, label="Volume", color="green")
axs[1, 0].set_title("Volume")
axs[1, 0].set_xlabel("Step")
axs[1, 0].set_ylabel("Volume (A³)")

# Plot Temp vs Step for reference
axs[1, 1].plot(step, temp, label="Temp", color="purple")
axs[1, 1].set_title("Temperature")
axs[1, 1].set_xlabel("Step")
axs[1, 1].set_ylabel("Temp (K)")

# Adjust layout to prevent overlap
plt.tight_layout()

# Save the figure as a PDF with 300 dpi
plt.savefig("lammps_plots.pdf", dpi=300)

# Show the plots
plt.show()

