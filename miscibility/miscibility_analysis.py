#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Integrated Pcond pipeline: dump parsing -> Fe-prob histogram/peaks -> phase windows ->
per-frame composition counts -> frame-averaged atomic fractions -> (optional) last-frame k-label dump.

What this script does
---------------------
1) Reads a LAMMPS dump (npt.dump by default) using your MolecularModTools `read_lammpstrj_dump`.
2) Computes conditional Fe-probabilities (Pcond) for every atom in each processed frame
   using neighborhood composition:
      - For a neighborhood with (nbA, nbB) counts, we estimate P(Fe | nbA, nbB) = (#central Fe)/(#atoms)
   Then for Fe atoms we store that Pcond directly; for non-Fe we also store the Fe-probability
   (1 - P(non-Fe | nbA, nbB)).
3) Writes two CSVs in OutputFiles/Pcond_Values/:
      - pFe.csv            (rows: Fe atoms, cols: frame indices)
      - 1 - pNonFe.csv     (rows: Non-Fe atoms, cols: frame indices)  -> Fe-probability of non-Fe atoms
   It also saves the atom type array and Fe/Non-Fe index arrays for later reconstruction.
4) From the last `ratio` fraction of frames (default 0.5), collects Fe-atom Pcond values and makes
   a 1-column file Graphs/all_pconds_raw.txt (header 'pcond').
5) Builds a histogram (density) with fixed bin size and smooths it (Gaussian). Finds the two tallest
   peaks and their distance (immiscibility measure). Saves plots and binned data in Graphs/.
6) If immiscible (distance >= immiscibility_cutoff), defines silicate and metal windows around the two peaks:
      - peak-centered windows of width ±(threshold_scale * peak_distance)
   Using *all atoms'* Fe-probabilities, classifies each atom per frame into silicate / interface / metal,
   tallies elemental composition for each phase, and aggregates across selected frames.
   Writes:
      - OutputFiles/composition_counts_per_frame.csv
      - OutputFiles/element_fractions_per_frame.csv
      - OutputFiles/composition_counts_aggregated.csv (summed over frames)
      - OutputFiles/element_fraction_summary_average.csv (fractions from aggregated counts)
   Also writes a last-frame visualization dump with a 'k' column (0=silicate, 2=interface, 1=metal):
      - OutputFiles/immisc_region_viz.dump
7) If miscible (distance < immiscibility_cutoff), reports it and still produces the histogram plot,
   but skips phase composition outputs.

You may run the parser once (expensive), then comment out `process_dump_file()` in main and re-run
only the analysis section as you tweak parameters (ratio, bin size, smoothing, thresholds, etc.).

NOTE
----
- Atom type mapping is assumed to be: 1=Mg, 2=Si, 3=O, 4=Fe. Adjust `ATOM_MAP` if yours differs.
- This script needs SciPy for best peak detection/smoothing but includes NumPy fallbacks.
- The `read_lammpstrj_dump` import path is set below; change it to your environment if needed.
"""

import os
default_n_threads = 8
os.environ['OPENBLAS_NUM_THREADS'] = f"{default_n_threads}"
os.environ['MKL_NUM_THREADS']    = f"{default_n_threads}"
os.environ['OMP_NUM_THREADS']    = f"{default_n_threads}"

import sys
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.spatial import KDTree
import matplotlib.pyplot as plt

# Set global font to Arial
#plt.rcParams['font.family'] = 'Arial'

# --- Optional SciPy-based smoothing/peaks; graceful fallbacks ---
try:
    from scipy.ndimage import gaussian_filter1d
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False

try:
    from scipy.signal import find_peaks
    HAVE_FIND_PEAKS = True
except Exception:
    HAVE_FIND_PEAKS = False

# ---- MolecularModTools dump reader (edit this path if needed) ----
sys.path.append('/home/hl0876/bin/miscib_analysis_final/MolecularModTools/OutputInfo/')
from ReadBox import read_lammpstrj_dump  # Uses .dump reader from MolecularModTools

# ====================== CONFIGURATION ======================
BASE_DIR = "./"              # Base working folder
DUMP_FILE = os.path.join(BASE_DIR, "npt.dump")
DUMP_INTERVAL = 2           # Read every Nth frame from dump
NNB = 100                    # Number of nearest neighbors to consider
DESIRED_FE_TYPE = 4          # Atom type that represents Fe

# Histogram / smoothing / peak detection
BIN_SIZE = 0.01             # Pcond bin width (units of probability)
SIGMA_BINS = 0.4             # Gaussian smoothing strength (in bin units)
AS_DENSITY = True            # Use density (integrates to 1); if False, uses counts + normalization
REF_TOTAL = 10               # Only used when AS_DENSITY=False

# Frame selection for analysis / averaging
RATIO_LAST_FRAMES = 0.1      # Use last 50% of frames for histogram/averaging

# Immiscibility threshold and window scales
IMMISCIBILITY_CUTOFF = 0.005  # If peak distance < this, treat as miscible
THRESHOLD_SCALE = 0.01       # Phase core window half-width = THRESHOLD_SCALE * peak_distance
VIZ_EXPAND = 1.0            # Visualization windows are ±(VIZ_EXPAND * THRESHOLD_SCALE * peak_distance)

# Directories
OUTPUT_DIR = os.path.join(BASE_DIR, "OutputFiles")
PCOND_DIR  = os.path.join(OUTPUT_DIR, "Pcond_Values")
GRAPHS_DIR = os.path.join(BASE_DIR, "Graphs")
os.makedirs(PCOND_DIR, exist_ok=True)
os.makedirs(GRAPHS_DIR, exist_ok=True)

# Atom type map (adjust to your system if different)
ATOM_MAP = {1: "Mg", 2: "Si", 3: "O", 4: "Fe"}

# ===========================================================


def smooth_gaussian(y: np.ndarray, sigma_bins: float) -> np.ndarray:
    """Gaussian smoothing in 'bin units'. Uses SciPy if available; falls back to NumPy convolution."""
    y = np.asarray(y, dtype=float)
    if sigma_bins <= 0:
        return y
    if HAVE_SCIPY:
        return gaussian_filter1d(y, sigma=sigma_bins, mode="reflect")
    # NumPy fallback
    half = int(np.ceil(4 * sigma_bins))
    kx = np.arange(-half, half + 1, dtype=float)
    kernel = np.exp(-0.5 * (kx / sigma_bins) ** 2)
    kernel /= kernel.sum()
    return np.convolve(y, kernel, mode="same")

def find_two_tallest_peaks_including_edges(
    x: np.ndarray,
    y: np.ndarray,
    snap_edge_to_01: bool = True,
    min_prominence: float | None = None,
):
    """
    Find two tallest peaks, treating array edges as valid peaks.
    Returns (sorted_peak_positions, distance) where positions are in x-units.
    If fewer than 2 peaks are found, returns ([], None).

    Args:
        x: bin centers (monotonic)
        y: values (density or counts)
        snap_edge_to_01: if True, when the top peak is at the first/last bin,
                         report its x as 0.0 or 1.0 (helpful for Pcond).
        min_prominence: optional prominence to suppress noise peaks.
    """
    # 1) standard interior peaks
    kwargs = {}
    if min_prominence is not None:
        kwargs["prominence"] = min_prominence
    try:
        from scipy.signal import find_peaks
        idx, _ = find_peaks(y, **kwargs)
    except Exception:
        idx = np.array([], dtype=int)

    # 2) explicit edge candidates
    candidates = [(i, y[i]) for i in idx]

    # left edge: y[0] >= y[1] and > 0 (so it's a “local” max at boundary)
    if y.size >= 2 and y[0] >= y[1] and y[0] > 0:
        candidates.append((0, y[0]))
    # right edge: y[-1] >= y[-2] and > 0
    if y.size >= 2 and y[-1] >= y[-2] and y[-1] > 0:
        candidates.append((len(y) - 1, y[-1]))

    # handle tiny arrays
    if y.size == 1 and y[0] > 0:
        candidates.append((0, y[0]))

    if not candidates:
        return [], None

    # 3) take the two tallest by height
    candidates.sort(key=lambda t: t[1])  # sort by y-height
    top = candidates[-2:] if len(candidates) >= 2 else candidates[-1:]

    # 4) map indices -> x positions; optionally snap edges to 0 or 1
    pos = []
    n = len(y)
    for (i, _) in top:
        xi = x[i]
        if snap_edge_to_01:
            if i == 0:
                xi = 0.0
            elif i == n - 1:
                xi = 1.0
        pos.append(xi)

    pos = np.sort(pos)
    dist = float(pos[1] - pos[0]) if len(pos) == 2 else None
    return pos.tolist(), dist


def process_dump_file() -> None:
    """
    Parse dump and compute Fe-probabilities (Pcond) per atom per processed frame.
    Saves:
        - OutputFiles/Pcond_Values/pFe.csv             (rows=Fe atoms, cols=frame indices)
        - OutputFiles/Pcond_Values/1 - pNonFe.csv      (rows=Non-Fe atoms, cols=frame indices)  [Fe-prob of non-Fe]
        - OutputFiles/npt_lastframe.dump               (raw last frame for visualization k-labeling)
        - OutputFiles/atom_types.npy                   (atom type per atom index; 1-based index matches dump id)
        - OutputFiles/fe_indices.npy                   (0-based indices of Fe atoms)
        - OutputFiles/nonfe_indices.npy                (0-based indices of non-Fe atoms)
    """
    print("\n=== Processing dump for Pcond ===")
    frames, atom_type_global = read_lammpstrj_dump(DUMP_FILE, interval=DUMP_INTERVAL)
    Type = np.asarray(atom_type_global, dtype=int)
    Natoms = Type.size

    fe_mask = (Type == DESIRED_FE_TYPE)
    nonfe_mask = ~fe_mask
    fe_indices = np.where(fe_mask)[0]
    nonfe_indices = np.where(nonfe_mask)[0]

    # Persist type/index info for later reconstruction
    np.save(os.path.join(OUTPUT_DIR, "atom_types.npy"), Type)
    np.save(os.path.join(OUTPUT_DIR, "fe_indices.npy"), fe_indices)
    np.save(os.path.join(OUTPUT_DIR, "nonfe_indices.npy"), nonfe_indices)

    time_list = []
    fe_pconds_by_frame = []       # each: array of length N_Fe (ordered by atom index among Fe)
    nonfe_pconds_by_frame = []    # each: array of length N_nonFe (Fe-probabilities for non-Fe atoms)

    for index, (step, coordinates, boxvecs) in enumerate(frames):
        print(f"  frame {index} (timestep {step})")
        xyz = np.asarray(coordinates, dtype=float)
        kdtree = KDTree(xyz)
        dist, points = kdtree.query(xyz, NNB)

        # Label central & neighbor types: Fe -> 1, non-Fe -> 2
        Type_label = np.where(Type == DESIRED_FE_TYPE, 1, 2)

        # (Natoms x NNB) type labels for neighbors (first column is the central atom itself)
        df_type = pd.DataFrame(Type_label[points], dtype=int)
        xi = pd.DataFrame(df_type[0])              # central atom type
        df_nb = df_type.drop(columns=[0], axis=1)  # neighbors only

        # Count neighbors by type
        xi['nbA'] = (df_nb == 1).sum(axis=1).astype(int)
        xi['nbB'] = (df_nb == 2).sum(axis=1).astype(int)

        # Group by neighbor counts, for central type Fe (1) and non-Fe (2)
        xi_groups = xi.groupby(['nbA', 'nbB']).size()
        xi_A = xi[xi[0] == 1].groupby(['nbA', 'nbB']).size()
        xi_B = xi[xi[0] == 2].groupby(['nbA', 'nbB']).size()

        # Conditional fractions
        group_A_frac = (xi_A / xi_groups).to_dict()  # P(central is Fe | nbA, nbB)
        group_B_frac = (xi_B / xi_groups).to_dict()  # P(central is non-Fe | nbA, nbB)

        # Build per-atom arrays (full length Natoms) then extract Fe / Non-Fe slices
        cond_A_full = np.full(Natoms, np.nan, dtype=float)
        cond_B_full = np.full(Natoms, np.nan, dtype=float)

        for i, row in xi.iterrows():
            nbA, nbB = int(row['nbA']), int(row['nbB'])
            cond_A_full[i] = group_A_frac.get((nbA, nbB), np.nan)
            cond_B_full[i] = group_B_frac.get((nbA, nbB), np.nan)

        # Fe-probability for *every* atom: Fe atoms take cond_A_full; non-Fe take 1 - cond_B_full
        feprob_full = np.where(fe_mask, cond_A_full, 1.0 - cond_B_full)

        # Store slice-ordered by atom index for persistence
        fe_pconds_by_frame.append(feprob_full[fe_indices])
        nonfe_pconds_by_frame.append(feprob_full[nonfe_indices])

        time_list.append(index)  # simple index; you can store 'step' if preferred

        # Save the raw last frame (for later k-label visualization dump)
        if index == len(frames) - 1:
            out_path = os.path.join(OUTPUT_DIR, "npt_lastframe.dump")
            with open(DUMP_FILE, 'r') as f:
                lines = f.readlines()
            frame_indices = [i for i, line in enumerate(lines) if "ITEM: TIMESTEP" in line]
            last_frame_start = frame_indices[-1]
            with open(out_path, 'w') as out:
                out.writelines(lines[last_frame_start:])

    # Persist Pcond matrices: rows=atoms subset, cols=frames
    Pcond_A = pd.DataFrame(np.column_stack(fe_pconds_by_frame), columns=time_list, index=fe_indices)
    Pcond_B = pd.DataFrame(np.column_stack(nonfe_pconds_by_frame), columns=time_list, index=nonfe_indices)

    Pcond_A.to_csv(os.path.join(PCOND_DIR, "pFe.csv"), index=False)
    Pcond_B.to_csv(os.path.join(PCOND_DIR, "1 - pNonFe.csv"), index=False)

    print(f"Saved pFe.csv and 1 - pNonFe.csv to {PCOND_DIR}")
    print("Saved atom_types.npy / fe_indices.npy / nonfe_indices.npy and npt_lastframe.dump to OutputFiles/")


def build_histogram_and_peaks(infile: str, out_prefix: str,
                              bin_size: float, sigma_bins: float,
                              as_density: bool = True, ref_total: float = 10):
    """
    Load 1-column 'pcond' file, make histogram and smoothed curve, find two tallest peaks.
    Saves binned data, smoothed curve, plot. Returns (centers, y_out, peak_positions, peak_distance).
    """
    data = np.loadtxt(infile, skiprows=1)  # single column after header
    p = np.asarray(data, dtype=float)
    p = p[np.isfinite(p)]
    if p.size == 0:
        raise RuntimeError("No valid Pcond values loaded from file.")

    # Edges/centers over observed range clamped to [0,1]
    pmin = max(0.0, np.floor(p.min() / bin_size) * bin_size)
    pmax = min(1.0, np.ceil(p.max() / bin_size) * bin_size)
    edges = np.arange(pmin, pmax + bin_size, bin_size)
    centers = 0.5 * (edges[:-1] + edges[1:])
    hist, _ = np.histogram(p, bins=edges)

    if as_density:
        total = hist.sum()
        y_base = hist.astype(float) / (total * bin_size) if total > 0 else np.zeros_like(hist, dtype=float)
        y_label = "Probability Density"
    else:
        y_base = hist.astype(float)
        y_label = "Counts"

    y_smooth = smooth_gaussian(y_base, sigma_bins)
    y_smooth = np.clip(y_smooth, 0, None)

    # Optional normalization for counts mode
    if not as_density:
        total_counts = y_base.sum()
        scale = (ref_total / total_counts) if total_counts > 0 else 1.0
        y_out_base = y_base * scale
        y_out = y_smooth * scale
        y_label_out = f"Normalized Counts (total={ref_total})"
    else:
        y_out_base = y_base
        y_out = y_smooth
        y_label_out = "Probability Density"

    # Peaks
    peak_positions, peak_distance = find_two_tallest_peaks_including_edges(centers, y_out, snap_edge_to_01=True, min_prominence=None)
    if peak_positions:
        print(f"Two tallest peaks at: {peak_positions}, distance = {peak_distance:.6f}")
    else:
        print("Less than two peaks found for the smoothed curve.")

    # Save outputs
    Path(out_prefix).parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(f"{out_prefix}_binned.txt",
               np.column_stack([centers, y_out_base]),
               fmt="%.10g",
               header=f"x_center  {'density' if as_density else f'normalized_counts(total={ref_total})'}",
               comments="")
    np.savetxt(f"{out_prefix}_curve.txt",
               np.column_stack([centers, y_out]),
               fmt="%.10g",
               header=f"x_center  {'smoothed_density' if as_density else f'normalized_smoothed_counts(total={ref_total})'}  (sigma={sigma_bins} bins)",
               comments="")

    # Plot
    plt.figure(figsize=(10, 6))
    plt.step(centers, y_out_base, where="mid", alpha=0.4, label="Histogram")
    plt.plot(centers, y_out, color="blue", linewidth=1.5, alpha=0.9, label=f"Smoothed)")
    for px in peak_positions:
        plt.axvline(px, color="red", linestyle="--")
    title = f"Pcond histogram (Δ={bin_size})"
    if not as_density:
        title += f", normalized to {ref_total}"
    if peak_distance is not None:
        title += f"\nTwo-peak distance = {peak_distance:.6f}"
    plt.title(title)
    plt.xlabel("Fe Conditional Probability", fontsize=14)
    plt.ylabel(y_label_out, fontsize=14)
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_plot.png", dpi=300)
    plt.close()

    return centers, y_out_base, y_out, peak_positions, peak_distance


def analyze_region_and_miscibility() -> None:
    """
    1) Loads pFe.csv and (1 - pNonFe).csv
    2) Selects the last RATIO_LAST_FRAMES fraction of frames
    3) Builds Graphs/all_pconds_raw.txt from Fe atoms only (per your request)
    4) Histogram + smoothing + peak detection (on Fe-only values)
    5) If immiscible: define phase windows, classify *all atoms* per frame, compute per-frame and aggregated
       composition counts and elemental fractions. Also write a last-frame k-label dump.
    """
    print("\n=== Analyzing regions and miscibility ===")
    file_A = os.path.join(PCOND_DIR, "pFe.csv")         # rows: Fe atoms, cols: frame indices
    file_B = os.path.join(PCOND_DIR, "1 - pNonFe.csv")  # rows: Non-Fe atoms, cols: frame indices (already Fe-prob)
    df_A = pd.read_csv(file_A)
    df_B = pd.read_csv(file_B)

    # Determine selected frame columns
    n_frames = df_A.shape[1]
    if n_frames == 0:
        raise RuntimeError("pFe.csv has no frames/columns.")
    k = max(1, int(n_frames * RATIO_LAST_FRAMES))
    selected_cols = list(df_A.columns[-k:])

    # 3) Build ALL-atoms Pcond file (Fe atoms + non-Fe atoms' Fe-prob) for histogram
    #    Use the same selected frame columns for both A and B.
    vals_A = df_A[selected_cols].to_numpy().ravel()      # Pcond for Fe atoms
    vals_B = df_B[selected_cols].to_numpy().ravel()      # Fe-prob for non-Fe atoms
    all_vals = np.concatenate([vals_A, vals_B])
    all_vals = all_vals[np.isfinite(all_vals)]

    all_pconds_path = os.path.join(GRAPHS_DIR, "all_pconds_raw.txt")
    np.savetxt(all_pconds_path, all_vals.reshape(-1, 1), fmt="%.10g", header="pcond", comments="")
    print(f"Wrote ALL-atoms Pcond values (Fe + non-Fe Fe-prob) to {all_pconds_path}")

    # 4) Histogram + peaks
    centers, y_base, y_smooth, peak_positions, peak_distance = build_histogram_and_peaks(
        infile=all_pconds_path,
        out_prefix=os.path.join(GRAPHS_DIR, "pcond_hist"),
        bin_size=BIN_SIZE,
        sigma_bins=SIGMA_BINS,
        as_density=AS_DENSITY,
        ref_total=REF_TOTAL
    )

    if (not peak_positions) or (peak_distance is None) or (peak_distance < IMMISCIBILITY_CUTOFF):
        print(f"System appears MISCIBLE based on two-peak distance: {peak_distance}")
        print(f"(cutoff = {IMMISCIBILITY_CUTOFF})")
        print(f"Saved histogram outputs in {GRAPHS_DIR}.")
        return

    # Immiscible case: define windows based on peaks
    silicate_peak = min(peak_positions)
    metal_peak    = max(peak_positions)
    immiscibility = peak_distance
    print(f"System appears IMMISCIBLE. distance = {immiscibility:.6f}")
    thresh = THRESHOLD_SCALE * immiscibility       # core windows (narrow) for composition stats
    viz_thresh = VIZ_EXPAND * thresh               # widened for visualization k-labels

    sil_lo, sil_hi = (silicate_peak - thresh), (silicate_peak + thresh)
    met_lo, met_hi = (metal_peak    - thresh), (metal_peak    + thresh)
    sil_hi_viz     = (silicate_peak + viz_thresh)
    met_lo_viz     = (metal_peak    - viz_thresh)

    # Load type/index maps to reconstruct full-atom arrays per frame
    Type = np.load(os.path.join(OUTPUT_DIR, "atom_types.npy"))
    fe_indices = np.load(os.path.join(OUTPUT_DIR, "fe_indices.npy"))
    nonfe_indices = np.load(os.path.join(OUTPUT_DIR, "nonfe_indices.npy"))
    Natoms = Type.size

    # Prepare composition tallies
    elements = ["O", "Mg", "Si", "Fe"]
    def fresh_counts():
        return {elem: 0 for elem in elements}

    # Per-frame outputs
    per_frame_counts = []     # dict per frame: phase -> element -> count
    per_frame_fracs  = []     # dict per frame: phase -> element -> fraction

    # Iterate over selected frames
    for col in selected_cols:
        # Reconstruct full Fe-probability array in atom order
        pcond_full = np.empty(Natoms, dtype=float)
        pcond_full[fe_indices]    = df_A[col].to_numpy()
        pcond_full[nonfe_indices] = df_B[col].to_numpy()  # already Fe-prob for non-Fe

        # Tally compositions using *core* phase windows
        metal_comp    = fresh_counts()
        silicate_comp = fresh_counts()
        interface_comp = fresh_counts()

        for i in range(Natoms):
            atom_type = int(Type[i])
            species = ATOM_MAP.get(atom_type, f"T{atom_type}")
            pv = pcond_full[i]

            if np.isnan(pv):
                continue

            if (pv >= sil_lo) and (pv <= sil_hi):
                silicate_comp[species] += 1
            elif (pv > sil_hi) and (pv < met_lo):
                interface_comp[species] += 1
            elif (pv >= met_lo) and (pv <= met_hi):
                metal_comp[species] += 1
            # else: outside both core windows -> ignore for composition stats

        # Convert to fractions per phase (normalized by that phase's total atoms)
        def to_fractions(counts: dict):
            total = sum(counts.values())
            if total == 0:
                return {k: 0.0 for k in counts}
            return {k: counts[k] / total for k in counts}

        fracs_metal    = to_fractions(metal_comp)
        fracs_silicate = to_fractions(silicate_comp)
        fracs_interface = to_fractions(interface_comp)

        per_frame_counts.append({
            "frame": col,
            "Metal": dict(metal_comp),
            "Silicate": dict(silicate_comp),
            "Interface": dict(interface_comp),
        })
        per_frame_fracs.append({
            "frame": col,
            "Metal": dict(fracs_metal),
            "Silicate": dict(fracs_silicate),
            "Interface": dict(fracs_interface),
        })

    # Aggregate across frames (sum counts, then derive fractions from sums)
    agg_counts = {
        "Metal": fresh_counts(),
        "Silicate": fresh_counts(),
        "Interface": fresh_counts(),
    }
    for rec in per_frame_counts:
        for phase in ["Metal", "Silicate", "Interface"]:
            for elem in elements:
                agg_counts[phase][elem] += rec[phase][elem]

    def frac_from_counts(counts: dict):
        total = sum(counts.values())
        if total == 0:
            return {k: 0.0 for k in counts}
        return {k: counts[k] / total for k in counts}

    agg_fracs = {
        "Metal": frac_from_counts(agg_counts["Metal"]),
        "Silicate": frac_from_counts(agg_counts["Silicate"]),
        "Interface": frac_from_counts(agg_counts["Interface"]),
    }

    # Save per-frame counts
    # Wide table: rows=frames, columns=phase.element count
    def dicts_to_wide(dict_list, value_kind="count"):
        rows = []
        for rec in dict_list:
            row = {"frame": rec["frame"]}
            for phase in ["Metal", "Silicate", "Interface"]:
                for elem in elements:
                    row[f"{phase}_{elem}"] = rec[phase][elem]
            rows.append(row)
        return pd.DataFrame(rows).sort_values("frame")

    df_counts_per_frame = dicts_to_wide(per_frame_counts, "count")
    df_counts_per_frame.to_csv(os.path.join(OUTPUT_DIR, "composition_counts_per_frame.csv"), index=False)

    # Per-frame fractions
    def dicts_to_wide_frac(dict_list):
        rows = []
        for rec in dict_list:
            row = {"frame": rec["frame"]}
            for phase in ["Metal", "Silicate", "Interface"]:
                for elem in elements:
                    row[f"{phase}_{elem}"] = rec[phase][elem]
            rows.append(row)
        return pd.DataFrame(rows).sort_values("frame")

    df_fracs_per_frame = dicts_to_wide_frac(per_frame_fracs)
    df_fracs_per_frame.to_csv(os.path.join(OUTPUT_DIR, "element_fractions_per_frame.csv"), index=False)

    # Aggregated counts & average (from aggregated counts)
    df_agg_counts = pd.DataFrame(agg_counts).T[elements]
    df_agg_counts.index.name = "Phase"
    df_agg_counts.to_csv(os.path.join(OUTPUT_DIR, "composition_counts_aggregated.csv"))

    df_agg_fracs = pd.DataFrame(agg_fracs).T[elements]
    df_agg_fracs.index.name = "Phase"
    df_agg_fracs.to_csv(os.path.join(OUTPUT_DIR, "element_fraction_summary_average.csv"))
    
    #compute 1 SD per column excluding exact zeros and append to element_fraction_summary_average.csv
    phases = ["Metal", "Silicate", "Interface"]
    cols_all = [f"{p}_{e}" for p in phases for e in elements]
    cols = [c for c in cols_all if c in df_fracs_per_frame.columns]
    X = df_fracs_per_frame[cols].replace(0, np.nan)
    sds = X.std(skipna=True, ddof=1)
    sd_rows = []
    for p in phases:
        row = {"Phase": p}
        for e in elements:
            col = f"{p}_{e}"
            row[f"{e}_sd"] = float(sds[col]) if col in sds.index else np.nan
        sd_rows.append(row)
    sd_df = pd.DataFrame(sd_rows, columns=["Phase"] + [f"{e}_sd" for e in elements])
    avg_path = os.path.join(OUTPUT_DIR, "element_fraction_summary_average.csv")
    avg_df = pd.read_csv(avg_path)
    out_df = avg_df.merge(sd_df, on="Phase", how="left")
    avg_cols_in_file = [e for e in elements if e in avg_df.columns]
    sd_cols_in_file  = [f"{e}_sd" for e in avg_cols_in_file if f"{e}_sd" in out_df.columns]
    desired_cols = ["Phase"] + avg_cols_in_file + sd_cols_in_file
    out_df = out_df[[c for c in desired_cols if c in out_df.columns]]
    out_df.to_csv(avg_path, index=False)
    print(f"Updated {avg_path} with columns: {', '.join([c for c in desired_cols if c in out_df.columns])}")

    print(f"Saved per-frame & aggregated composition summaries to {OUTPUT_DIR}")

    # Last-frame visualization dump with 'k' column for Ovito coloring
    last_frame_path = os.path.join(OUTPUT_DIR, "npt_lastframe.dump")
    k_dump_out = os.path.join(OUTPUT_DIR, "immisc_region_viz.dump")
    if os.path.exists(last_frame_path):
        # Use the *last selected frame*'s pcond_full for k-labels
        last_col = selected_cols[-1]
        pcond_full = np.empty(Natoms, dtype=float)
        pcond_full[fe_indices]    = df_A[last_col].to_numpy()
        pcond_full[nonfe_indices] = df_B[last_col].to_numpy()

        with open(last_frame_path, 'r') as f_in, open(k_dump_out, 'w') as f_out:
            in_atoms = False
            for line in f_in:
                if line.startswith("ITEM: ATOMS"):
                    f_out.write(line.strip() + " k\n")
                    in_atoms = True
                    continue
                elif in_atoms:
                    parts = line.strip().split()
                    atom_id = int(parts[0])     # LAMMPS ids are 1-based
                    pv = pcond_full[atom_id - 1]

                    # Visualization k using *widened* windows
                    if pv <= sil_hi_viz:
                        kval = 0  # Silicate
                    elif pv >= met_lo_viz:
                        kval = 1  # Metal
                    else:
                        kval = 2  # Interface
                    f_out.write(f"{line.strip()} {kval}\n")
                else:
                    f_out.write(line)
        print(f"Wrote {k_dump_out} for Ovito (k: 0=silicate, 2=interface, 1=metal)")
    else:
        print("npt_lastframe.dump not found; skipping immisc_region_viz.dump creation.")

    # Optional: overlay windows on histogram plot for clarity
    centers2 = centers
    y_base2  = y_base
    ys2      = y_smooth

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.step(centers2, y_base2, where="mid", alpha=0.4, label="Histogram")
    ax.plot(centers2, ys2, linewidth=1.5, label=f"Smoothed")

    # overlay “viz windows”
    #ax.axvspan(0, sil_hi_viz, alpha=0.2, label="Silicate (viz)")
    #ax.axvspan(met_lo_viz, 1, alpha=0.2, label="Metal (viz)")

    # peaks and core bounds
    ax.axvline(silicate_peak, color="blue", linestyle="--", label=f"Silicate peak: {silicate_peak:.4f}")
    ax.axvline(sil_lo,        color="blue", alpha=0.5, linestyle=":",  label="Silicate core bounds")
    ax.axvline(sil_hi,        color="blue", alpha=0.5, linestyle=":")
    ax.axvline(metal_peak,    color="red", linestyle="--", label=f"Metal peak: {metal_peak:.4f}")
    ax.axvline(met_lo,        color="red", alpha=0.5, linestyle=":",  label="Metal core bounds")
    ax.axvline(met_hi,        color="red", alpha=0.5, linestyle=":")

    ax.set_xlabel("Fe Conditional Probability", fontsize=14)
    ax.set_ylabel("Probability Density" if AS_DENSITY else "Counts", fontsize=14)
    ax.legend(loc="upper right", ncol=1, fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(GRAPHS_DIR, "pcond_hist_plot_with_overlay.png"), dpi=300)
    plt.close(fig)

    


def main():
    # 1) Run once to generate Pcond CSVs and lastframe dump. Comment out after first run if desired.
    #process_dump_file()

    # 2) Analyze histogram/peaks and compute per-frame & averaged compositions
    analyze_region_and_miscibility()

    print("\nAll done. Outputs written to:")
    print(f"  - {PCOND_DIR}")
    print(f"  - {GRAPHS_DIR}")
    print(f"  - {OUTPUT_DIR}")


if __name__ == "__main__":
    main()

