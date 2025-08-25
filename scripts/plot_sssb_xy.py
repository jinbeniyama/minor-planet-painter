#!/usr/bin/env python3
"""Make a figure of spatial distribution of minor bodies.

See about Hildas.
http://www.solexorb.it/Animated/Rotframe.html
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

from minor_planet_painter.common import (
    MPCORB, jd2utc, utc2jd, get_planet_positions, get_planet_orbits,
    solve_kepler_eq, mpcepoch2jd
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot minor planets from MPCORB.DAT")
    parser.add_argument(
        "epoch", 
        help="UTC epoch, format YYYY-MM-DDTHH:MM:SS")
    parser.add_argument(
        "--Nobj", type=int, 
        help="Number of asteroids to be plotted")
    parser.add_argument(
        "--range", type=float, default=6.5, 
        help="Plot range in AU (square from -range to +range)")
    parser.add_argument(
        "--MPCORB", default=None, 
        help="Path to MPCORB.DAT")
    parser.add_argument(
        "--out", type=str, default="MPCORB.jpg", 
        help="Figure name")
    args = parser.parse_args()
    
    if args.MPCORB is None:
        fi = MPCORB
    else:
        fi = args.MPCORB

    # parse UTC
    t_utc_dt = datetime.fromisoformat(args.epoch)
    t_utc_iso = t_utc_dt.isoformat()


    # Calculate locations =====================================================
    M, omega, Omega, i, e, a = [], [], [], [], [], []
    n = []
    epoch = []
    pattern = re.compile(r"[-+]?\d*\.\d+|\d+")

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

            if args.Nobj is not None and len(a) >= args.Nobj:
                break  

    print(f"  N_sssbs = {len(M)}")

    # convert to numpy arrays and radians
    M = np.deg2rad(np.array(M))
    omega = np.deg2rad(np.array(omega))
    Omega = np.deg2rad(np.array(Omega))
    i = np.deg2rad(np.array(i))
    e = np.array(e)
    a = np.array(a)
    n = np.deg2rad(n)
    epoch = np.array(epoch)

    jd_now = utc2jd(t_utc_iso)

    epoch_jd = np.array([mpcepoch2jd(e) for e in epoch])
    dt_days = jd_now - epoch_jd
    # Linear assumption
    M_t = M + n*dt_days

    ejd_min = float(np.min(epoch_jd))
    ejd_max = float(np.max(epoch_jd))
    #print("jd_now:", jd_now)
    #print(f"epoch_jd min: {ejd_min}")
    #print(f"epoch_jd max: {ejd_max}")
    print("  Check epochs:")
    print(f"    Minimum (earliest): {jd2utc(ejd_min)}")
    print(f"    Maximum (latest)  : {jd2utc(ejd_max)}")
    
    # Calculate eccentric anomaly
    E = solve_kepler_eq(M_t, e)
    # Calculate nu and r 
    nu = 2*np.arctan2(np.sqrt(1+e)*np.sin(E/2), np.sqrt(1-e)*np.cos(E/2))
    nu_total = np.mod(omega + nu, 2*np.pi)
    r = a * (1 - e*np.cos(E))
    
    cos_Omega = np.cos(Omega)
    sin_Omega = np.sin(Omega)
    cos_i = np.cos(i)
    sin_i = np.sin(i)
    cos_omega = np.cos(omega)
    sin_omega = np.sin(omega)

    # Finally calculate x and y
    x = r * (cos_Omega * np.cos(nu_total) - sin_Omega * np.sin(nu_total) * cos_i)
    y = r * (sin_Omega * np.cos(nu_total) + cos_Omega * np.sin(nu_total) * cos_i)
    # Calculate locations =====================================================
     


    # Plot ====================================================================
    # Masks for classification
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

    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_axes([0.1, 0.1, 0.85, 0.85])

    ax.scatter(x, y, s=0.5, color=colors, alpha=0.9, label='_nolegend_')
    for label, color in zip(legend_labels, cols):
        ax.scatter([], [], color=color, label=label, s=20)

    ax.set_xlabel('X [AU]')
    ax.set_ylabel('Y [AU]')

    planet_positions = get_planet_positions(t_utc_iso)
    planet_orbits = get_planet_orbits(t_utc_iso)
    
    # Planet orbit 
    for name, (ox, oy) in planet_orbits.items():
        ax.plot(ox, oy, linestyle='-', linewidth=0.8, alpha=0.6, label=None)
    
    # Planet
    for name, (px, py, pz) in planet_positions.items():
        ax.scatter(px, py, s=40, ec="black", label=None)

    # The Sun
    ax.scatter(0, 0, color='yellow', s=50, ec="black", label=None, zorder=10)

    ax.set_title(f'{t_utc_iso} UTC')
    r = args.range
    ax.set_xlim(-r, r)
    ax.set_ylim(-r, r)
    ax.set_aspect('equal')
    ax.legend(fontsize=12)
    plt.show(block=False)
    
    print()
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
    # Plot ====================================================================
