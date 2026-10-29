#!/usr/bin/env python3
"""Make figures of orbital elements and absolute magnitude."""

import argparse
from pathlib import Path
import warnings

import numpy as np
import matplotlib.pyplot as plt
from erfa import ErfaWarning

from minor_planet_painter import MPCORB, mycolor

warnings.simplefilter("ignore", ErfaWarning)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot minor planets from MPCORB.DAT"
    )
    parser.add_argument(
        "--MPCORB", default=None,
        help="Path to MPCORB.DAT"
    )
    parser.add_argument(
        "--onlyNEA", action="store_true", default=False,
        help="Plot only NEA"
    )
    parser.add_argument(
        "--out", type=str, default="orbelem.jpg",
        help="Output filename for orbital-element figure"
    )
    parser.add_argument(
        "--col_MBA", type=str, default=mycolor[1],
        help="Color for MBA (default: mycolor[1])"
    )
    parser.add_argument(
        "--col_NEA", type=str, default=mycolor[0],
        help="Color for NEA (default: mycolor[0])"
    )

    args = parser.parse_args()

    fi = MPCORB if args.MPCORB is None else args.MPCORB

    # Extract orbital elements ================================================
    M, omega, Omega, i, e, a = [], [], [], [], [], []
    H, n, epoch = [], [], []

    with open(fi, "r") as f:
        for line in f:
            if line.startswith("---"):
                break

        for line in f:
            if not line.strip():
                continue

            epoch.append(line[20:25].strip())
            M.append(float(line[26:35].strip()))
            omega.append(float(line[37:46].strip()))
            Omega.append(float(line[48:57].strip()))
            i.append(float(line[59:68].strip()))
            e.append(float(line[70:79].strip()))
            n.append(float(line[80:91].strip()))
            a.append(float(line[92:103].strip()))

            H_value = line[8:13].strip()
            H.append(float(H_value) if H_value else -99.0)

    print(f"  N_sssbs = {len(M)}")

    # Convert to arrays =======================================================
    e = np.asarray(e)
    i = np.asarray(i)
    a = np.asarray(a)
    H = np.asarray(H)

    q = a * (1.0 - e)

    # Masks for classification
    mask_nea = q < 1.3
    mask_mba = (q >= 1.3) & (a >= 1.8) & (a <= 3.3)
    mask_hilda = (
        (q >= 1.3) & (a >= 3.7) & (a <= 4.0)
        & (e >= 0.07) & (e <= 0.3)
    )
    mask_trojan = (q >= 1.3) & (a >= 5.0) & (a <= 5.4)
    mask_tno = (q >= 1.3) & (a >= 30.0)

    # Count numbers
    N_all = len(a)
    N_nea = np.sum(mask_nea)
    N_mba = np.sum(mask_mba)
    N_hilda = np.sum(mask_hilda)
    N_trojan = np.sum(mask_trojan)
    N_tno = np.sum(mask_tno)
    N_others = N_all - (
        N_nea + N_mba + N_hilda + N_trojan + N_tno
    )

    # Colors
    col_MBA = args.col_MBA
    col_NEA = args.col_NEA

    # Output filenames
    out_path = Path(args.out)
    out_H = out_path.with_name(
        f"{out_path.stem}_H{out_path.suffix}"
    )
    # Absolute magnitude validity mask ========================================
    # Exclude dummy and unreasonably large H values
    valid_H = np.isfinite(H) & (H >= 0) & (H < 35)

    print(f"  N with valid H (0 <= H < 35) = {np.sum(valid_H)}")
    print(f"  NEAs with valid H = {np.sum(mask_nea & valid_H)}")
    print(f"  MBAs with valid H = {np.sum(mask_mba & valid_H)}")

    # Diameter conversion (assuming pV = 0.15)
    pv = 0.15

    def H_to_D(H):
        return 1329.0 / np.sqrt(pv) * 10.0 ** (-np.asarray(H) / 5.0)

    def D_to_H(D):
        D = np.asarray(D, dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            return 5.0 * np.log10(
                1329.0 / (D * np.sqrt(pv))
            )

    # Figure 2: Absolute magnitude H ==========================================
    # Histogram helper
    def plot_histograms(ax, x, bins, mask, label, color, ls):
        ax.hist(
            x[mask],
            histtype="step",
            bins=bins,
            color=color,
            ls=ls,
            linewidth=1.5,
            density=not args.onlyNEA,
            label=label
        )

    fig2, ax_H = plt.subplots(figsize=(8, 6))

    mask_nea_H = mask_nea & valid_H
    mask_mba_H = mask_mba & valid_H

    # H bins: 0.5 mag
    H_valid = H[valid_H]
    bin_width_H = 0.5

    bin_H = np.arange(
        np.floor(np.min(H_valid) / bin_width_H) * bin_width_H,
        np.ceil(np.max(H_valid) / bin_width_H) * bin_width_H
        + bin_width_H * 0.5,
        bin_width_H
    )

    if args.onlyNEA:
        plot_histograms(
            ax_H, H, bin_H, mask_nea_H,
            f"NEA N={np.sum(mask_nea_H)}", "black", "solid"
        )
        ax_H.set_ylabel("N")

    else:
        plot_histograms(
            ax_H, H, bin_H, mask_mba_H,
            f"MBA N={np.sum(mask_mba_H)}", col_MBA, "solid"
        )
        plot_histograms(
            ax_H, H, bin_H, mask_nea_H,
            f"NEA N={np.sum(mask_nea_H)}", col_NEA, "dashed"
        )
        ax_H.set_ylabel("Normalized fraction")

    ax_H.set_xlabel("Absolute magnitude H [mag]")

    # Put large H (small diameter) on the left
    ax_H.set_xlim(np.max(bin_H), np.min(bin_H))
    ax_H.legend()
    

    from matplotlib.ticker import FixedLocator, FixedFormatter, AutoMinorLocator
    
    # Secondary x-axis: diameter in km
    ax_top = ax_H.secondary_xaxis(
        "top",
        functions=(H_to_D, D_to_H)
    )
    ax_top.set_xlabel("D [km]")
    
    # Major ticks: powers of 10 in diameter [km]
    D_major = np.array([0.001, 0.01, 0.1, 1, 10, 100, 1000])
    H_major = D_to_H(D_major)
    
    ax_top.xaxis.set_major_locator(FixedLocator(D_major))
    ax_top.xaxis.set_major_formatter(
        FixedFormatter([f"{d:g}" for d in D_major])
    )
    
    # Minor ticks between powers of 10
    D_minor = np.concatenate([
        np.arange(2, 10) * decade
        for decade in [1e-3, 1e-2, 1e-1, 1, 10, 100]
    ])
    ax_top.xaxis.set_minor_locator(FixedLocator(D_minor))
    ax_top.tick_params(axis="x", which="minor", length=3)

    fig2.savefig(out_H, dpi=300, bbox_inches="tight")
    print(f"  Figure is saved as {out_H}")
