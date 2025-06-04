#!/usr/bin/env python3
"""
vasp_md_traces.py
-----------------
Plot TOTEN (from OSZICAR), ionic temperature (OSZICAR), and total pressure
(OUTCAR) versus time.  Writes vasp_md_traces.pdf.

Usage
-----
python plot_vasp_md_tp.py OSZICAR OUTCAR POTIM_fs
                                ^^^^^^^^^^^^^^^^
                                MD time step in femtoseconds
"""
import re, sys, numpy as np, matplotlib.pyplot as plt
from pathlib import Path

# ---------------- command-line ----------------
if len(sys.argv) != 4:
    sys.exit("Usage:  python vasp_md_traces.py OSZICAR OUTCAR POTIM_fs")

osz, outc = Path(sys.argv[1]), Path(sys.argv[2])
dt_ps     = float(sys.argv[3]) * 1e-3          # fs → ps

for f in (osz, outc):
    if not f.is_file():
        sys.exit(f"File not found: {f}")

# ---------------- regex patterns -------------
# OSZICAR line:   25 F= -8.56E+02 E0= ... T= 4217.53  P= ...
num = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][+-]?\d+)?"

re_osz = re.compile(
    rf"\s*\d+\s+T=\s*({num}).*?F=\s*({num})",
    re.I,
)
# OUTCAR line: total pressure  =  1169.58 kB
re_press = re.compile(r"total\s+pressure\s*=\s*([-+]?\d+\.\d+)\s*kB", re.I)

# ---------------- read OSZICAR ----------------
toten, temp = [], []
with osz.open() as f:
    for ln in f:
        m = re_osz.match(ln)
        if m:
            toten.append(float(m.group(1)))
            temp.append(float(m.group(2)))

if not toten:
    sys.exit("No MD-style lines found in OSZICAR – check the file.")

# ---------------- read OUTCAR -----------------
pres = []
with outc.open() as f:
    for ln in f:
        m = re_press.search(ln)
        if m:
            pres.append(float(m.group(1)))      # still in kbar

if not pres:
    sys.exit("No 'total pressure =' lines found in OUTCAR.")

# ---------------- synchronise lengths ----------
n = min(len(toten), len(temp), len(pres))
toten, temp, pres = toten[:n], temp[:n], pres[:n]
time = np.arange(n) * dt_ps          # ps

print(f"Found {n} MD steps  (Δt = {dt_ps:.4f} ps)")

# convert kbar → GPa if you prefer
# pres = [p/10.0 for p in pres]      # uncomment for GPa
# ylabel_P = "P (GPa)"
ylabel_P = "P (kbar)"

# ---------------- plotting ---------------------
fig, axs = plt.subplots(3, 1, figsize=(6.2, 9), sharex=True)

axs[0].plot(time, toten, lw=0.9)
axs[0].set_ylabel("T  (K)")
axs[0].set_title("Ionic temperature")

axs[1].plot(time, temp, lw=0.9)
axs[1].set_ylabel("TOTEN  (eV)")
axs[1].set_title("Electronic free energy per MD step")

axs[2].plot(time, pres, lw=0.9)
axs[2].set_ylabel(ylabel_P)
axs[2].set_xlabel("time  (ps)")
axs[2].set_title("Total pressure")

for ax in axs:
    ax.grid(alpha=0.3)

fig.tight_layout()
fig.savefig("vasp_md_ETP.pdf", dpi=450)
fig.show()

