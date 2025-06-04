#!/usr/bin/env python3
"""
merge_blocks.py

Merge repeated species‐blocks in a POSCAR into single blocks.

Usage:
    python merge_blocks.py POSCAR_in POSCAR_out
"""
import sys

if len(sys.argv) != 3:
    print(__doc__)
    sys.exit(1)

inp, outp = sys.argv[1], sys.argv[2]

# --- read entire file ---
with open(inp) as f:
    lines = [l.rstrip() for l in f]

# comment, scale, lattice
comment    = lines[0]
scale      = lines[1]
lattices   = lines[2:5]

# original species and counts (e.g. 8 entries each)
species_all = lines[5].split()
counts_all  = list(map(int, lines[6].split()))
coord_type  = lines[7]

# data lines
coords_all = lines[8:]

# detect single‐block length
# find first repetition of species_all[0] after index 0
blk_size = next(i for i in range(1, len(species_all)) if species_all[i]==species_all[0])
n_blocks = len(species_all)//blk_size

# unique species in order:
species = species_all[:blk_size]

# split counts into blocks and sum them
counts = [ sum(counts_all[i::blk_size]) for i in range(blk_size) ]

# convert coords_all into a list of lists, one per block
blocks = []
idx = 0
for b in range(n_blocks):
    block = coords_all[idx: idx+ sum(counts_all[b*blk_size:(b+1)*blk_size]) ]
    blocks.append(block)
    idx += len(block)

# now for each species, collect its coords across all blocks
merged_coords = []
for i, sp in enumerate(species):
    # in each block, the coords for species i occupy a contiguous slice
    # find start/end in that block
    for b in range(n_blocks):
        # compute offset of species i in block
        offs = sum(counts_all[b*blk_size + j] for j in range(i))
        n   = counts_all[b*blk_size + i]
        merged_coords.extend(blocks[b][offs: offs+n])

# --- write new POSCAR ---
with open(outp, 'w') as f:
    f.write(f"{comment}\n")
    f.write(f"{scale}\n")
    for L in lattices:
        f.write(f"{L}\n")
    f.write("   " + "   ".join(species) + "\n")
    f.write("   " + "   ".join(map(str, counts)) + "\n")
    f.write(f"{coord_type}\n")
    for line in merged_coords:
        f.write(f"{line}\n")

print(f"Written {outp}: species={species}, counts={counts}")

