#!/usr/bin/env python3
"""
make_supercell.py   – build an N×N×N super‑cell from a VASP POSCAR.

Fixes:
  • wraps Direct coordinates back into [0,1) so the tiling is correct
  • groups atoms by species to satisfy the POSCAR convention
  • keeps ‘Selective dynamics’ lines if present
Usage:
  python make_supercell.py  POSCAR_in  N  POSCAR_out
"""
import sys
import numpy as np
from collections import Counter
from ase.io import read

def main() -> None:
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)

    fn_in, n_rep, fn_out = sys.argv[1], int(sys.argv[2]), sys.argv[3]

    # ------------------------------------------------------------------ read
    atoms      = read(fn_in, format="vasp")
    supercell  = atoms.repeat((n_rep, n_rep, n_rep))

    # ---------------------------------------------------------------- header
    with open(fn_in) as f:
        hdr = f.read().splitlines()

    comment = hdr[0].rstrip()

    ptr_species = 5                                         # default index
    if hdr[6].strip().lower().startswith(("s", "selective")):
        selective = True
        ptr_coord = 7
    else:
        selective = False
        ptr_coord = 6

    species_order = hdr[ptr_species].split()
    coord_mode_in = hdr[ptr_coord].strip().lower()
    use_direct    = coord_mode_in.startswith("d")

    # ---------------------------------------------------------------- counts
    counts = Counter(supercell.get_chemical_symbols())
    counts_line = " ".join(str(counts[s]) for s in species_order)

    # ---------------------------------------------------------------- sort by species
    order_index = {sp: i for i, sp in enumerate(species_order)}
    sorted_atoms = supercell[np.argsort([order_index[a.symbol] for a in supercell])]

    # ---------------------------------------------------------------- positions
    if use_direct:
        pos = sorted_atoms.get_scaled_positions(wrap=False) % 1.0  # wrap!
        coord_tag = "Direct"
    else:
        pos = sorted_atoms.get_positions()
        coord_tag = "Cartesian"

    # ---------------------------------------------------------------- write POSCAR
    with open(fn_out, "w") as f:
        f.write(comment + "\n")
        f.write("1.0\n")
        for v in sorted_atoms.get_cell():
            f.write(f"  {v[0]:.16f}  {v[1]:.16f}  {v[2]:.16f}\n")

        f.write("  " + "  ".join(species_order) + "\n")
        f.write("  " + counts_line + "\n")

        if selective:
            f.write("Selective dynamics\n")
        f.write(coord_tag + "\n")

        for x, y, z in pos:
            f.write(f"  {x:.10f}  {y:.10f}  {z:.10f}\n")

    print(f"✅  {fn_out} written ({len(sorted_atoms)} atoms, {n_rep}× super‑cell)")

if __name__ == "__main__":
    main()

