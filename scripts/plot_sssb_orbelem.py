#!/usr/bin/env python3
"""Make a figure of orbital elements.
"""
import argparse
import re
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import warnings
from erfa import ErfaWarning
# To suppress ErfaWarning: ERFA function "dtf2d" yielded 1 of "dubious year (Note 6)"
warnings.simplefilter('ignore', ErfaWarning)

from minor_planet_painter import MPCORB, mycolor


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot minor planets from MPCORB.DAT")
    parser.add_argument(
        "--MPCORB", default=None, 
        help="Path to MPCORB.DAT")
    parser.add_argument(
        "--onlyNEA", action="store_true", default=False, 
        help="Plot only NEA")
    parser.add_argument(
        "--out", type=str, default="orbelem.jpg", 
        help="Figure name")
    args = parser.parse_args()
    
    if args.MPCORB is None:
        fi = MPCORB
    else:
        fi = args.MPCORB


    # Extract orbital elements ================================================
    M, omega, Omega, i, e, a = [], [], [], [], [], []
    n = []
    epoch = []
    with open(fi, 'r') as f:
        for line in f:
            # End of header
            if line.startswith('---'):
                break 
        for line in f:
            if len(line.strip()) == 0:
                continue
            #if float(line[92:103].strip()) < 5.5:
            #    continue
            epoch.append(line[20:25].strip())
            M.append(float(line[26:35].strip()))
            omega.append(float(line[37:46].strip()))
            Omega.append(float(line[48:57].strip()))
            i.append(float(line[59:68].strip()))
            e.append(float(line[70:79].strip()))
            a.append(float(line[92:103].strip()))
            n.append(float(line[80:91].strip()))
    print(f"  N_sssbs = {len(M)}")
    # Extract orbital elements ================================================

    
    # Plot ====================================================================
    # Masks for classification
    e = np.array(e)
    i = np.array(i)
    a = np.array(a)
    q = a * (1 - e) 
    mask_nea = (q < 1.3)
    # q >= 1.3 is necessary to choose exclusively
    mask_mba = (q >= 1.3) & (a >= 1.8) & (a <= 3.3)
    mask_hilda = (q >= 1.3) & (a >= 3.7) & (a <= 4.0) & (e >= 0.07) & (e <= 0.3)
    mask_trojan = (q >= 1.3) & (a >= 5.0) & (a <= 5.4)
    mask_tno = (q >= 1.3) & (a >= 30.0)

    # Colors
    cols = ["red", "green", "orange", "blue", "skyblue", "gray"]

    # Default: others
    colors = np.full(a.shape, cols[-1], dtype=object) 
    colors[mask_nea] = cols[0]
    colors[mask_mba] = cols[1]
    colors[mask_hilda] = cols[2]
    colors[mask_trojan] = cols[3]
    colors[mask_tno] = cols[4]
    
    # Count numbers 
    N_all = len(a)
    N_nea = np.sum(mask_nea)
    N_mba = np.sum(mask_mba)
    N_hilda = np.sum(mask_hilda)
    N_trojan = np.sum(mask_trojan)
    N_tno = np.sum(mask_tno)
    N_others = N_all - (N_nea + N_mba + N_hilda + N_trojan + N_tno)
    
    # Labels
    legend_labels = [
        f"NEA (q<1.3): N={N_nea}",
        f"MBA (q>=1.3, a=1.8–3.3): N={N_mba}",
        f"Hilda (q$\geq$1.3, a=3.7-4.0, e=0.07-0.30): N={N_hilda}",
        f"Trojan (q$\geq$1.3, a=5.0-5.4): N={N_trojan}",
        f"TNO (q$\geq$1.3, a$\geq$30.0): N={N_tno}",
        f"Others: N={N_others}"
    ]

    # 1. eccentricity
    bin_e = np.arange(0, 1.01, 0.01)
    # Do not plot retrograde bodies
    bin_i = np.arange(0, 91, 1)

    fig = plt.figure(figsize=(12, 8))
    ax_e = fig.add_axes([0.15, 0.60, 0.82, 0.35])
    ax_i = fig.add_axes([0.15, 0.12, 0.82, 0.35])


    if args.onlyNEA:
        ax_e.hist(
            e[mask_nea], histtype="step", bins=bin_e, color="black", 
            label=f"NEA N={N_nea}")
        ax_i.hist(
            i[mask_nea], histtype="step", bins=bin_i, color="black", 
            label=f"NEA N={N_nea}")
        ax_e.set_ylabel('N')
        ax_i.set_ylabel('N')
    else:
        col_mba = mycolor[0]
        col_nea = mycolor[1]
        ax_e.hist(
            e[mask_mba], histtype="step", ls="solid", bins=bin_e, color=col_mba, 
            label=f"MBA N={N_mba}", density=True)
        ax_i.hist(
            i[mask_mba], histtype="step", ls="solid", bins=bin_i, color=col_mba, 
            label=f"MBA N={N_mba}", density=True)
        ax_e.hist(
            e[mask_nea], histtype="step", ls="dashed", bins=bin_e, color=col_nea, 
            label=f"NEA N={N_nea}", density=True)
        ax_i.hist(
            i[mask_nea], histtype="step", ls="dashed", bins=bin_i, color=col_nea, 
            label=f"NEA N={N_nea}", density=True)
        ax_e.set_ylabel('Normalized fraction')
        ax_i.set_ylabel('Normalized fraction')

    ax_e.set_xlabel('Eccentricity')
    ax_e.legend(fontsize=12)

    ax_i.set_xlabel('Inclination [deg]')
    ax_i.legend(fontsize=12)

    # Align 
    x, y = -0.08, 0.5
    ax_e.yaxis.set_label_coords(x, y)
    ax_i.yaxis.set_label_coords(x, y)

    plt.show(block=False)
    
    ans = input("  Save figure? (y/n): ").strip().lower()
    if ans != 'y':
        plt.close()
    else:
        try:
            plt.savefig(args.out)
            print(f"  Figure is saved as {args.out}")
        except ValueError:
            print("  Not saved. Exiting.")
        plt.close()

    # 2. Inclination
    # 3. and so on......
    # Plot ====================================================================
