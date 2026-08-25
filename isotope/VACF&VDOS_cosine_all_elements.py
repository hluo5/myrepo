#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
General VACF -> cosine-transform VDOS -> mean force constant analysis
for ALL elements in a VASP XDATCAR trajectory.

Required user inputs (edit only these in the USER INPUT section):
    1) XDATCAR_FILE
    2) DT_FS
    3) TEMPERATURE_K
    4) ELEMENT_MASSES_AMU

The element names are normally read automatically from XDATCAR.
For VASP4-style XDATCAR files without element symbols, set ELEMENTS_OVERRIDE.

Outputs:
    - VACF for every element
    - cosine-transform VDOS for every element
    - force-constant convergence F(nu_max) for every element
    - final mean force constant for every element
    - CSV data files
    - PNG and PDF figures

Theory:
    C_vv(t) = <v_i(t0) . v_i(t0+t)>_(t0, atoms)
    g(nu)   ~ integral_0^infinity C_vv(t) cos(2*pi*nu*t) dt

    <nu^2> = integral nu^2 g(nu) dnu / integral g(nu) dnu
    (the plotted/integrated axis is cm^-1; internally it is converted to Hz)
    F       = m (2*pi)^2 <nu^2>

Notes:
    - The VACF is calculated with multi-time-origin averaging.
    - A mass-weighted total-system center-of-mass drift is removed at each time.
    - The VACF cosine transform uses a half-cosine taper and a DCT-I,
      which implements the trapezoidal endpoint weights correctly.
    - No empirical high-frequency baseline subtraction is used.
    - Negative finite-trajectory oscillations in the cosine spectrum are NOT clipped.
    - FINAL_CUTOFF_CM1 = 3500.0 cm^-1 is a default numerical choice, not a universal
      physical constant. Always inspect the F-vs-cutoff convergence plot.
"""

# =============================================================================
# USER INPUT: normally only edit this block
# =============================================================================

XDATCAR_FILE = "XDATCAR"

DT_FS = 1.0
TEMPERATURE_K = 3500.0

ELEMENT_MASSES_AMU = {
    "Mg": 24.305,
    "Si": 28.085,
    "O": 15.999,
    "Fe": 55.845,
}

# Only needed for an old VASP4-style XDATCAR without element symbols.
# Example: ELEMENTS_OVERRIDE = ["Mg", "Si", "O", "Fe"]
ELEMENTS_OVERRIDE = None


# =============================================================================
# OPTIONAL NUMERICAL SETTINGS
# These defaults can normally be left unchanged.
# =============================================================================

# Maximum correlation time used for the VACF.
# If the trajectory is too short, the code automatically shortens it.
MAX_CORR_FS = 500.0

# Default upper integration wavenumber used for the single "final" F value.
# This is intentionally separate from the convergence scan.
FINAL_CUTOFF_CM1 = 3500.0

# Force-constant convergence scan.
CUTOFF_MIN_CM1 = 300.0
CUTOFF_STEP_CM1 = 50.0

# Plot range for VDOS. None -> use at least FINAL_CUTOFF_CM1 and up to 150 cm^-1.
VDOS_PLOT_MAX_CM1 = None

# FFT chunk size for VACF calculation.
# Smaller value = lower memory, larger value = often faster.
VACF_CHUNK_ATOMS = 128

# Save both PNG and PDF.
SAVE_PDF = True

# Float32 substantially reduces memory for large trajectories.
POSITION_DTYPE = "float32"


# =============================================================================
# IMPORTS
# =============================================================================

import os
import csv
import math
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pathlib import Path
from scipy.fft import rfft, irfft, dct, next_fast_len


AMU = 1.66053906660e-27  # kg
C_CM_S = 2.99792458e10    # speed of light, cm/s


# =============================================================================
# XDATCAR READER
# =============================================================================

def _nonempty_line(lines, idx):
    """Advance over blank lines and return (line, new_index)."""
    n = len(lines)
    while idx < n and not lines[idx].strip():
        idx += 1
    if idx >= n:
        return None, idx
    return lines[idx].strip(), idx


def read_xdatcar(filename, elements_override=None, dtype=np.float32):
    """
    Read a standard VASP XDATCAR with a fixed simulation cell.

    Returns
    -------
    frac_all : ndarray, shape (n_frames, n_atoms, 3)
        Wrapped fractional coordinates.
    cell : ndarray, shape (3, 3)
        Cell vectors in Angstrom, stored by row.
    species : ndarray, shape (n_atoms,)
        Element symbol of each atom.
    """
    print("=" * 78)
    print(f"Reading XDATCAR: {filename}")

    with open(filename, "r") as f:
        lines = f.readlines()

    if len(lines) < 8:
        raise ValueError("XDATCAR is too short or malformed.")

    comment = lines[0].strip()
    scale_raw = float(lines[1].split()[0])

    cell0 = np.array(
        [list(map(float, lines[i].split()[:3])) for i in range(2, 5)],
        dtype=np.float64,
    )

    # VASP scale handling: positive = scale factor; negative = target volume.
    if scale_raw > 0:
        scale_factor = scale_raw
    elif scale_raw < 0:
        target_volume = abs(scale_raw)
        raw_volume = abs(np.linalg.det(cell0))
        scale_factor = (target_volume / raw_volume) ** (1.0 / 3.0)
    else:
        raise ValueError("XDATCAR scale factor cannot be zero.")

    cell = cell0 * scale_factor

    line5 = lines[5].split()

    # VASP5/6: line 6 has symbols, line 7 has counts.
    # VASP4: line 6 has counts and symbols are absent.
    is_vasp4 = all(tok.lstrip("+-").isdigit() for tok in line5)

    if is_vasp4:
        counts = list(map(int, line5))
        if elements_override is None:
            raise ValueError(
                "This looks like a VASP4-style XDATCAR without element symbols.\n"
                "Please set ELEMENTS_OVERRIDE, e.g. ['Mg','Si','O','Fe']."
            )
        elements = list(elements_override)
        if len(elements) != len(counts):
            raise ValueError(
                "ELEMENTS_OVERRIDE length does not match the number of atom-count groups."
            )
        ptr = 6
    else:
        elements = line5
        counts = list(map(int, lines[6].split()))
        if len(elements) != len(counts):
            raise ValueError("Element-name and atom-count lines have different lengths.")
        ptr = 7

    species = []
    for e, n in zip(elements, counts):
        species.extend([e] * n)
    species = np.asarray(species)
    natoms = len(species)

    frames = []
    coord_mode = None

    while ptr < len(lines):
        line, ptr0 = _nonempty_line(lines, ptr)
        if line is None:
            break

        low = line.lower()

        if low.startswith("direct"):
            coord_mode = "direct"
            ptr = ptr0 + 1
        elif low.startswith("cart"):
            coord_mode = "cartesian"
            ptr = ptr0 + 1
        else:
            # Some files may omit a repeated configuration marker.
            # If a coordinate mode has already been established, try reading coordinates.
            if coord_mode is None:
                ptr = ptr0 + 1
                continue
            ptr = ptr0

        if ptr + natoms > len(lines):
            break

        frame = np.empty((natoms, 3), dtype=dtype)
        ok = True

        for i in range(natoms):
            vals = np.fromstring(lines[ptr + i], sep=" ")
            if vals.size < 3:
                ok = False
                break
            frame[i] = vals[:3]

        if not ok:
            break

        if coord_mode == "cartesian":
            # Convert Cartesian Angstrom coordinates to fractional coordinates.
            # This code assumes a fixed cell.
            frame = (frame.astype(np.float64) @ np.linalg.inv(cell)).astype(dtype)

        frames.append(frame)
        ptr += natoms

    if len(frames) < 3:
        raise ValueError("Need at least 3 XDATCAR frames to compute central-difference velocities.")

    frac_all = np.asarray(frames, dtype=dtype)

    unique, counts_unique = np.unique(species, return_counts=True)

    print(f"Comment      : {comment}")
    print(f"Frames       : {frac_all.shape[0]}")
    print(f"Atoms        : {natoms}")
    print("Composition  :", dict(zip(unique, counts_unique)))
    print(f"Cell volume  : {abs(np.linalg.det(cell)):.6f} Angstrom^3")

    return frac_all, cell, species


# =============================================================================
# VELOCITY
# =============================================================================

def fractional_positions_to_velocity(frac_all, cell, dt_fs):
    """
    Compute Cartesian velocities using a central difference with PBC.

    Instead of explicitly unwrapping the entire trajectory, use the sum of
    two consecutive minimum-image displacements:

        v(t) = [dr(t-1 -> t) + dr(t -> t+1)] / (2*dt)

    Returns
    -------
    vel : ndarray, shape (n_frames-2, n_atoms, 3), Angstrom/fs
    """
    d_prev = frac_all[1:-1] - frac_all[:-2]
    d_next = frac_all[2:] - frac_all[1:-1]

    d_prev -= np.round(d_prev)
    d_next -= np.round(d_next)

    d_two = d_prev + d_next
    vel = np.einsum("tai,ij->taj", d_two, cell, optimize=True)
    vel /= (2.0 * dt_fs)

    return vel.astype(np.float32, copy=False)


def subtract_total_com_drift(vel, species, mass_dict):
    """
    Remove the instantaneous mass-weighted center-of-mass velocity
    of the entire system.
    """
    missing = sorted(set(species) - set(mass_dict))
    if missing:
        raise KeyError(
            "Missing atomic masses for element(s): " + ", ".join(missing)
        )

    masses = np.array([mass_dict[s] for s in species], dtype=np.float64)
    total_mass = masses.sum()

    # shape: (n_time, 3)
    v_com = np.einsum(
        "tai,a->ti",
        vel.astype(np.float64, copy=False),
        masses,
        optimize=True,
    ) / total_mass

    vel -= v_com[:, None, :].astype(vel.dtype)
    return vel


# =============================================================================
# MULTI-ORIGIN VACF
# =============================================================================

def compute_species_vacf_fft(vel, atom_indices, max_lag, chunk_atoms=128):
    """
    Compute a species-resolved multi-time-origin VACF efficiently by FFT.

    For each lag:
        VACF_raw(lag) =
            mean_{time origins, atoms} [ v_i(t0) dot v_i(t0+lag) ]

    The unbiased normalization uses all available time origins:
        number of origins = n_time - lag

    The final VACF is normalized so VACF(0) = 1.
    """
    n_time = vel.shape[0]
    n_atoms_species = len(atom_indices)

    if n_atoms_species == 0:
        raise ValueError("No atoms selected.")

    max_lag = min(int(max_lag), n_time - 1)

    # Zero-padding avoids circular correlation.
    n_fft = next_fast_len(2 * n_time - 1)

    corr_sum = np.zeros(n_time, dtype=np.float64)

    for start in range(0, n_atoms_species, chunk_atoms):
        idx = atom_indices[start:start + chunk_atoms]

        # Shape (time, 3*Nchunk)
        x = vel[:, idx, :].reshape(n_time, -1).astype(np.float64, copy=False)

        X = rfft(x, n=n_fft, axis=0)

        # Sum power over all selected atoms and Cartesian components BEFORE iFFT.
        power_sum = np.sum(
            X.real * X.real + X.imag * X.imag,
            axis=1,
            dtype=np.float64,
        )

        corr_sum += irfft(power_sum, n=n_fft)[:n_time]

    lags = np.arange(n_time, dtype=np.float64)
    n_origins = n_time - lags

    # Divide by time origins and number of atoms.
    # Cartesian x/y/z are summed, reproducing v(t0) dot v(t0+t).
    vacf_raw = corr_sum / (n_origins * n_atoms_species)

    vacf_raw = vacf_raw[:max_lag + 1]
    vacf = vacf_raw / vacf_raw[0]

    return vacf, vacf_raw


# =============================================================================
# VACF -> COSINE-TRANSFORM VDOS
# =============================================================================

def vacf_to_vdos_cosine(vacf, dt_fs):
    """
    Convert a one-sided VACF to a cosine-transform spectrum.

    A half-cosine taper is applied:
        w(t) = 0.5 * [1 + cos(pi*t/t_max)]

    The DCT-I has endpoint weights exactly matching a trapezoidal cosine
    integral. Therefore:

        integral C(t) cos(2*pi*nu*t) dt
        ~= 0.5 * dt * DCT-I[C(t)]

    after tapering.

    No empirical baseline correction and no negative-value clipping are used.
    """
    n = len(vacf)
    if n < 2:
        raise ValueError("VACF needs at least 2 points.")

    j = np.arange(n, dtype=np.float64)
    window = 0.5 * (1.0 + np.cos(np.pi * j / (n - 1)))

    x = vacf * window

    # dt_fs is included so the spectrum approximates the cosine integral.
    spectrum = 0.5 * dt_fs * dct(x, type=1)

    # DCT-I frequency grid:
    # nu_k = k / [2*(n-1)*dt]
    freq_hz = np.arange(n, dtype=np.float64) / (
        2.0 * (n - 1) * dt_fs * 1e-15
    )

    # Convert ordinary frequency (Hz) to spectroscopic wavenumber (cm^-1):
    # wavenumber = frequency / c
    wavenumber_cm1 = freq_hz / C_CM_S

    return wavenumber_cm1, spectrum, window


# =============================================================================
# FORCE CONSTANT
# =============================================================================

def trapz_compat(y, x):
    if hasattr(np, "trapezoid"):
        return np.trapezoid(y, x)
    return np.trapz(y, x)


def force_constant_from_vdos(wavenumber_cm1, vdos, cutoff_cm1, mass_amu):
    """
    Mean force constant from the second frequency moment:

        F = m (2*pi)^2 <nu^2>

        <nu^2> =
            integral_0^cutoff nu^2 g(nu) dnu
            / integral_0^cutoff g(nu) dnu

    Returns NaN if there are too few points or the integrated spectral
    weight is non-positive.
    """
    mask = (wavenumber_cm1 >= 0.0) & (wavenumber_cm1 <= cutoff_cm1)

    if np.count_nonzero(mask) < 3:
        return np.nan

    freq_hz = wavenumber_cm1[mask] * C_CM_S
    g = vdos[mask]

    area = trapz_compat(g, freq_hz)

    if not np.isfinite(area) or area <= 0:
        return np.nan

    second_moment = trapz_compat(freq_hz**2 * g, freq_hz) / area

    if not np.isfinite(second_moment) or second_moment < 0:
        return np.nan

    mass_kg = mass_amu * AMU
    F = mass_kg * (2.0 * np.pi)**2 * second_moment

    return F


# =============================================================================
# OUTPUT HELPERS
# =============================================================================

def save_figure(fig, outbase):
    fig.savefig(str(outbase) + ".png", dpi=300, bbox_inches="tight")
    if SAVE_PDF:
        fig.savefig(str(outbase) + ".pdf", bbox_inches="tight")
    plt.close(fig)


def write_column_csv(filename, first_name, first_values, data_dict):
    keys = list(data_dict.keys())
    arrays = [np.asarray(first_values)] + [np.asarray(data_dict[k]) for k in keys]
    n = min(len(a) for a in arrays)

    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([first_name] + keys)
        for i in range(n):
            writer.writerow([arrays[0][i]] + [a[i] for a in arrays[1:]])


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def main():
    xdatcar_path = Path(XDATCAR_FILE)

    if not xdatcar_path.exists():
        raise FileNotFoundError(f"Cannot find XDATCAR file: {xdatcar_path}")

    frac_all, cell, species = read_xdatcar(
        xdatcar_path,
        elements_override=ELEMENTS_OVERRIDE,
        dtype=np.float32,
    )

    # Preserve the order in which elements appear in the XDATCAR.
    elements = []
    for s in species:
        if s not in elements:
            elements.append(s)

    missing = [e for e in elements if e not in ELEMENT_MASSES_AMU]
    if missing:
        raise KeyError(
            "Please add atomic masses for: " + ", ".join(missing)
        )

    print("\nComputing central-difference velocities ...")
    vel = fractional_positions_to_velocity(frac_all, cell, DT_FS)

    # Free the large coordinate array before the VACF FFTs.
    del frac_all

    print("Removing instantaneous total-system COM drift ...")
    vel = subtract_total_com_drift(
        vel,
        species,
        ELEMENT_MASSES_AMU,
    )

    n_time = vel.shape[0]
    velocity_duration_fs = (n_time - 1) * DT_FS

    requested_lag = int(round(MAX_CORR_FS / DT_FS))

    # Avoid using a correlation window longer than half the available velocity trajectory.
    max_lag_allowed = max(2, (n_time - 1) // 2)
    max_lag = min(requested_lag, max_lag_allowed)
    actual_max_corr_fs = max_lag * DT_FS

    if max_lag < requested_lag:
        print(
            f"WARNING: requested MAX_CORR_FS={MAX_CORR_FS:g} fs is too long "
            f"for this trajectory; using {actual_max_corr_fs:g} fs instead."
        )

    print(f"Velocity frames    : {n_time}")
    print(f"Velocity duration  : {velocity_duration_fs:.3f} fs")
    print(f"VACF max lag       : {actual_max_corr_fs:.3f} fs")
    print(f"Temperature        : {TEMPERATURE_K:g} K")

    outdir = Path(
        f"{xdatcar_path.stem}_{TEMPERATURE_K:g}K_VACF_cosine_cm-1"
    )
    outdir.mkdir(parents=True, exist_ok=True)

    results = {}

    # -------------------------------------------------------------------------
    # Calculate every element
    # -------------------------------------------------------------------------
    for elem in elements:
        atom_indices = np.where(species == elem)[0]
        nat = len(atom_indices)
        mass = ELEMENT_MASSES_AMU[elem]

        print("\n" + "-" * 78)
        print(f"Element: {elem}")
        print(f"Atoms  : {nat}")
        print(f"Mass   : {mass:.6f} amu")

        vacf, vacf_raw = compute_species_vacf_fft(
            vel,
            atom_indices,
            max_lag=max_lag,
            chunk_atoms=VACF_CHUNK_ATOMS,
        )

        freq, vdos, window = vacf_to_vdos_cosine(vacf, DT_FS)

        results[elem] = {
            "natoms": nat,
            "mass_amu": mass,
            "vacf": vacf,
            "vacf_raw": vacf_raw,
            "wavenumber_cm1": freq,
            "vdos": vdos,
            "window": window,
        }

        print(f"VACF(0) raw       : {vacf_raw[0]:.8e} (Angstrom/fs)^2")
        print(f"VDOS Nyquist      : {freq[-1]:.3f} cm^-1")
        print(f"VDOS resolution   : {freq[1] - freq[0]:.6f} cm^-1")

    # All elements share the same VACF time/frequency grids.
    time_fs = np.arange(max_lag + 1) * DT_FS
    wavenumber_cm1 = next(iter(results.values()))["wavenumber_cm1"]
    nyquist_cm1 = wavenumber_cm1[-1]

    # -------------------------------------------------------------------------
    # Cutoff convergence
    # -------------------------------------------------------------------------
    cutoff_max = min(FINAL_CUTOFF_CM1, nyquist_cm1)

    if cutoff_max < FINAL_CUTOFF_CM1:
        print(
            f"\nWARNING: FINAL_CUTOFF_CM1={FINAL_CUTOFF_CM1:g} cm^-1 exceeds "
            f"the Nyquist limit ({nyquist_cm1:.3f} cm^-1). "
            f"Using {cutoff_max:.3f} cm^-1."
        )

    if cutoff_max <= CUTOFF_MIN_CM1:
        cutoffs = np.linspace(
            max(0.1, 0.1 * cutoff_max),
            cutoff_max,
            20,
        )
    else:
        cutoffs = np.arange(
            CUTOFF_MIN_CM1,
            cutoff_max + 0.5 * CUTOFF_STEP_CM1,
            CUTOFF_STEP_CM1,
        )
        if cutoffs[-1] < cutoff_max:
            cutoffs = np.append(cutoffs, cutoff_max)

    summary_rows = []

    for elem in elements:
        r = results[elem]

        F_conv = np.array([
            force_constant_from_vdos(
                r["wavenumber_cm1"],
                r["vdos"],
                c,
                r["mass_amu"],
            )
            for c in cutoffs
        ])

        F_final = force_constant_from_vdos(
            r["wavenumber_cm1"],
            r["vdos"],
            cutoff_max,
            r["mass_amu"],
        )

        r["F_conv"] = F_conv
        r["F_final"] = F_final

        summary_rows.append([
            elem,
            r["natoms"],
            r["mass_amu"],
            TEMPERATURE_K,
            DT_FS,
            actual_max_corr_fs,
            cutoff_max,
            F_final,
        ])

    # -------------------------------------------------------------------------
    # Print summary
    # -------------------------------------------------------------------------
    print("\n" + "=" * 92)
    print("FINAL MEAN FORCE CONSTANTS")
    print("VACF -> half-cosine taper -> DCT-I cosine transform")
    print("No empirical baseline subtraction; no negative-value clipping")
    print(f"Temperature         : {TEMPERATURE_K:g} K")
    print(f"dt                  : {DT_FS:g} fs")
    print(f"VACF max lag        : {actual_max_corr_fs:g} fs")
    print(f"Final cutoff        : {cutoff_max:g} cm^-1")
    print("=" * 92)
    print(f"{'Element':<12}{'N atoms':>12}{'Mass (amu)':>16}{'F (N/m)':>18}")
    print("-" * 58)

    for elem in elements:
        r = results[elem]
        print(
            f"{elem:<12}"
            f"{r['natoms']:>12d}"
            f"{r['mass_amu']:>16.6f}"
            f"{r['F_final']:>18.4f}"
        )

    print("=" * 92)
    print(
        "IMPORTANT: the reported final F uses the cutoff above. "
        "Use the convergence plot to verify that F has reached a stable plateau."
    )

    # -------------------------------------------------------------------------
    # CSV output
    # -------------------------------------------------------------------------
    write_column_csv(
        outdir / "VACF_all_elements.csv",
        "time_fs",
        time_fs,
        {e: results[e]["vacf"] for e in elements},
    )

    write_column_csv(
        outdir / "VDOS_all_elements.csv",
        "wavenumber_cm^-1",
        wavenumber_cm1,
        {e: results[e]["vdos"] for e in elements},
    )

    write_column_csv(
        outdir / "ForceConstant_convergence_all_elements.csv",
        "cutoff_cm^-1",
        cutoffs,
        {e: results[e]["F_conv"] for e in elements},
    )

    with open(outdir / "ForceConstant_summary.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "element",
            "n_atoms",
            "mass_amu",
            "temperature_K",
            "dt_fs",
            "max_corr_fs",
            "final_cutoff_cm-1",
            "mean_force_constant_N_per_m",
        ])
        writer.writerows(summary_rows)

    # -------------------------------------------------------------------------
    # Figure 1: VACF
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(9, 6))

    for elem in elements:
        ax.plot(
            time_fs,
            results[elem]["vacf"],
            label=elem,
            linewidth=1.5,
        )

    ax.axhline(0.0, linewidth=0.8)
    ax.set_xlabel("Correlation time (fs)")
    ax.set_ylabel("Normalized VACF")
    ax.set_title(f"Velocity autocorrelation functions ({TEMPERATURE_K:g} K)")
    ax.set_xlim(0, actual_max_corr_fs)
    ax.legend()
    fig.tight_layout()

    save_figure(fig, outdir / "VACF_all_elements")

    # -------------------------------------------------------------------------
    # Figure 2: VDOS
    # Normalize each curve by its own maximum positive amplitude for plotting only.
    # Raw VDOS is retained in the CSV and used for force constants.
    # -------------------------------------------------------------------------
    if VDOS_PLOT_MAX_CM1 is None:
        plot_max = min(
            nyquist_cm1,
            max(5000.0, 1.2 * cutoff_max),
        )
    else:
        plot_max = min(float(VDOS_PLOT_MAX_CM1), nyquist_cm1)

    fig, ax = plt.subplots(figsize=(9, 6))

    for elem in elements:
        g = results[elem]["vdos"]
        mask_plot = (wavenumber_cm1 >= 0) & (wavenumber_cm1 <= plot_max)

        positive_max = np.max(g[mask_plot])
        if positive_max > 0:
            g_plot = g / positive_max
        else:
            g_plot = g

        ax.plot(
            wavenumber_cm1,
            g_plot,
            label=elem,
            linewidth=1.5,
        )

    ax.axhline(0.0, linewidth=0.8)
    ax.set_xlabel(r"Wavenumber (cm$^{-1}$)")
    ax.set_ylabel("Normalized VDOS (plot only)")
    ax.set_title(
        f"Element-resolved VDOS from VACF cosine transform ({TEMPERATURE_K:g} K)"
    )
    ax.set_xlim(0, plot_max)
    ax.legend()
    fig.tight_layout()

    save_figure(fig, outdir / "VDOS_all_elements")

    # -------------------------------------------------------------------------
    # Figure 3: F convergence
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(9, 6))

    for elem in elements:
        ax.plot(
            cutoffs,
            results[elem]["F_conv"],
            marker="o",
            markersize=3.5,
            label=elem,
            linewidth=1.3,
        )

    ax.set_xlabel(r"Upper integration wavenumber (cm$^{-1}$)")
    ax.set_ylabel("Mean force constant (N/m)")
    ax.set_title(
        f"Force-constant convergence from VACF-derived VDOS ({TEMPERATURE_K:g} K)"
    )
    ax.legend()
    fig.tight_layout()

    save_figure(fig, outdir / "ForceConstant_convergence_all_elements")

    # -------------------------------------------------------------------------
    # Figure 4: final F summary
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5.5))

    F_values = [results[e]["F_final"] for e in elements]
    x = np.arange(len(elements))

    ax.bar(x, F_values)
    ax.set_xticks(x)
    ax.set_xticklabels(elements)
    ax.set_ylabel("Mean force constant (N/m)")
    ax.set_title(
        f"Mean force constants at cutoff = {cutoff_max:g} cm^-1 "
        f"({TEMPERATURE_K:g} K)"
    )

    fig.tight_layout()
    save_figure(fig, outdir / "ForceConstant_final_all_elements")

    print(f"\nAll outputs written to:\n  {outdir.resolve()}")


if __name__ == "__main__":
    main()
