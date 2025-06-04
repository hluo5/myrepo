#!/usr/bin/env python3
"""
Fast in-place reordering of each frame in a LAMMPS dump so that atoms
are grouped by TYPE (e.g. 111111… 222222… 333333…).

Usage
-----
    python sort_dump_by_type.py  bigdump.lammpstrj  bigdump_sorted.lammpstrj

If you need to sort by *id* instead, pass  --key id
"""

import sys
import argparse
import numpy as np

# ----------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Sort every frame in a LAMMPS dump by TYPE (default) or ID."
    )
    p.add_argument("input",  help="input dump file")
    p.add_argument("output", help="output dump file (sorted)")
    p.add_argument("--key",  choices=("type", "id"), default="type",
                   help="column to sort by (default: type)")
    return p.parse_args()

# ----------------------------------------------------------------------
def sort_dump_by_key(infile, outfile, key="type"):
    """
    Stream over `infile`, write a new dump to `outfile` where, in every frame,
    atoms are ordered by the chosen `key` (id or type).
    """
    with open(infile,  "rb") as fin, \
         open(outfile, "wb", buffering=2**20) as fout:

        while True:
            first = fin.readline()          # 'ITEM: TIMESTEP'
            if not first:                   # EOF
                break

            # ---- collect header lines until we hit 'ITEM: ATOMS ...' ----
            header_lines = [first]
            n_atoms = None
            while True:
                line = fin.readline()
                if not line:
                    raise EOFError("Unexpected EOF while reading header")
                header_lines.append(line)

                if line.startswith(b"ITEM: NUMBER OF ATOMS"):
                    n_atoms = int(fin.readline())
                    header_lines.append(str(n_atoms).encode() + b"\n")
                    continue                # keep going: we still need the rest

                if line.startswith(b"ITEM: ATOMS"):
                    atom_header = line      # save for later parsing
                    break                   # header complete

            # ---- know the atom field order ------------------------------
            fields = atom_header.split()[2:]        # after 'ITEM: ATOMS'
            try:
                key_idx = fields.index(key.encode())
            except ValueError:
                raise ValueError(
                    f"Column '{key}' not present in ATOMS line: {fields!r}"
                )

            # ---- read the atom block for this frame ---------------------
            atom_lines = [fin.readline() for _ in range(n_atoms)]
            if any(l == b"" for l in atom_lines):
                raise EOFError("Unexpected EOF while reading atom block")

            # ---- extract chosen key quickly with NumPy ------------------
            # Join lines -> one big ascii string, parse all numbers at once
            block_str = b"".join(atom_lines).decode("ascii", "replace")
            data = np.fromstring(block_str, sep=" ")   # 1-D float64
            ncols = len(atom_lines[0].split())         # tokens per line
            data = data.reshape(-1, ncols)             # (n_atoms, ncols)

            key_values = data[:, key_idx]
            order = np.argsort(key_values, kind="stable")

            # ---- write out: header + sorted atom lines ------------------
            fout.writelines(header_lines)
            #fout.writelines([atom_lines[i] for i in order])
            # ---- rewrite atom lines with new IDs (1 to N) ----
            id_idx = fields.index(b"id")  # get index of 'id' column
            new_atom_lines = []
            for new_id, i in enumerate(order, start=1):
                tokens = atom_lines[i].decode().strip().split()
                tokens[id_idx] = str(new_id)  # replace old ID
                new_line = " ".join(tokens) + "\n"
                new_atom_lines.append(new_line.encode())

            fout.writelines(new_atom_lines)

    print(f"Finished.  Sorted dump written to:  {outfile}")

# ----------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    sort_dump_by_key(args.input, args.output, key=args.key)

