#!/usr/bin/env python3
"""
usage:  python shuffle_trim_poscar.py  POSCAR_in  N_to_delete

Example:
    python shuffle_trim_poscar.py POSCAR 500
"""

import sys, random
from collections import OrderedDict
from ase.io import read, write

# ------------- input parsing ------------------------------------------------
if len(sys.argv) != 3:
    sys.exit("usage:  shuffle_trim_poscar.py  POSCAR_in  N_to_delete")

fname      = sys.argv[1]
n_to_drop  = int(sys.argv[2])
rng_seed   = None            # set to an int for reproducible shuffles
if rng_seed is not None:
    random.seed(rng_seed)

# ------------- load structure ----------------------------------------------
atoms = read(fname, format="vasp")
nat   = len(atoms)
if n_to_drop >= nat:
    sys.exit(f"N_to_delete must be < total atoms ({nat})")

# ------------- 1) random deletion ------------------------------------------
drop_indices = set(random.sample(range(nat), n_to_drop))
kept_indices = [i for i in range(nat) if i not in drop_indices]
trimmed      = atoms[kept_indices]            # an ASE Atoms object

# ------------- 2) shuffle *within* each element type -----------------------
# keep the original element order as it first appears
species_order = list(OrderedDict.fromkeys(trimmed.get_chemical_symbols()))

shuffled_idx = []
for el in species_order:
    idx_list = [i for i, a in enumerate(trimmed) if a.symbol == el]
    random.shuffle(idx_list)                  # independent shuffle per element
    shuffled_idx.extend(idx_list)

shuffled_atoms = trimmed[shuffled_idx]

# ------------- 3) write out -------------------------------------------------
write("POSCAR_trimmed", shuffled_atoms, vasp5=True, direct=True)
print(f"Deleted {n_to_drop} atoms and shuffled each species → POSCAR_trimmed")

