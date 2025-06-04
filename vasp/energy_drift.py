#!/usr/bin/env python3
"""
check_drift.py  – monitor energy drift in VASP Nosé NVT MD

Usage
-----
python check_drift.py OUTCAR POTIM_fs [--out drift.pdf]

Arguments
---------
OUTCAR       path to VASP OUTCAR
POTIM_fs     integration time step in femtoseconds
--out FILE   (optional) write figure to FILE instead of / in addition to
             displaying on screen; extension decides format (png, pdf …)
"""

import re, sys, argparse, numpy as np, matplotlib.pyplot as plt

# ---------- CLI ----------
parser = argparse.ArgumentParser()
parser.add_argument("outcar")
parser.add_argument("potim", type=float, help="POTIM in fs")
parser.add_argument("--out", help="filename for saved figure")
args = parser.parse_args()
dt_ps = args.potim * 1e-3     # fs → ps

# ---------- scan OUTCAR ----------
nions, etotal = None, []
with open(args.outcar) as f:
    for ln in f:
        if nions is None:
            m = re.search(r"NIONS\s*=\s*(\d+)", ln)
            if m:
                nions = int(m.group(1))
        m = re.search(r"total energy\s+ETOTAL\s*=\s*(-?\d+\.\d+)", ln)
        if m:
            etotal.append(float(m.group(1)))

if nions is None or not etotal:
    sys.exit("Could not locate NIONS or ETOTAL in the file – check regex.")

et = np.asarray(etotal)                    # eV
drift = (et - et[0]) * 1000.0 / nions      # meV / atom
time  = np.arange(len(et)) * dt_ps         # ps

# ---------- plot ----------
plt.figure(figsize=(6,4))
plt.plot(time, drift, lw=1)
plt.xlabel("time (ps)")
plt.ylabel("drift (meV / atom)")
plt.title("Energy drift from ETOTAL (Nosé NVT)")
plt.tight_layout()
plt.savefig("drift.pdf", dpi=450)   # always write PDF
print("Figure saved to drift.pdf")
plt.show()

# ---------- linear drift ----------
slope, _ = np.polyfit(time, drift, 1)      # meV atom⁻¹ ps⁻¹
print(f"Linear drift ≈ {slope:.3f} meV atom⁻¹ ps⁻¹ over {time[-1]:.1f} ps")

