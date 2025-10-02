#!/usr/bin/env python3
"""Make a figure of sky motion of minor bodies.
"""
import argparse
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
import warnings
from erfa import ErfaWarning
from astroquery.jplhorizons import Horizons
# To suppress ErfaWarning: ERFA function "dtf2d" yielded 1 of "dubious year (Note 6)"
warnings.simplefilter('ignore', ErfaWarning)

from minor_planet_painter.common import MPCORB


def fetch_ephemeris(obj_name, obscode, epoch_jd):
    obj = Horizons(id=obj_name, location=obscode, id_type="smallbody", epochs=epoch_jd)
    eph = obj.ephemerides()
    eph["RA_rate_arcsec_s"] = eph["RA_rate"]/3600.
    eph["DEC_rate_arcsec_s"] = eph["DEC_rate"]/3600.
    eph["vel"] = (eph["RA_rate_arcsec_s"]**2 + eph["DEC_rate_arcsec_s"]**2)**0.5
    eph["targetname"] = obj_name
    return eph.to_pandas()


def D_from_Hp(H, p, Herr=0, perr=0):
    """Calculate diameter from absolute magnitude H and albedo.

    Parameters
    ----------
    H : float or np.ndarray
        Absolute magnitude(s)
    p : float or np.ndarray
        Visual geometric albedo(s)
    Herr : float or np.ndarray, optional
        Uncertainty in H (default=0)
    perr : float or np.ndarray, optional
        Uncertainty in p (default=0)

    Returns
    -------
    D : float or np.ndarray
        Diameter(s) in km
    Derr : float or np.ndarray
        Uncertainty in diameter
    """
    H = np.asanyarray(H, dtype=float)
    p = np.asanyarray(p, dtype=float)
    Herr = np.asanyarray(Herr, dtype=float)
    perr = np.asanyarray(perr, dtype=float)

    # Diameter
    D = 1329.0 * p**(-0.5) * 10.0**(-0.2 * H)

    # Error propagation
    errH_2 = Herr**2 * (1329.0 * p**(-0.5) / 5.0 * 10.0**(-0.2 * H) * np.log(10))**2
    errp_2 = perr**2 * (1329.0 * p**(-1.5) / 2.0 * 10.0**(-0.2 * H))**2

    Derr = np.sqrt(errp_2 + errH_2)

    return D, Derr


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot sky motion of minor planets from MPCORB.DAT")
    parser.add_argument(
        "--Nobj", type=int, 
        help="Number of asteroids to be plotted")
    parser.add_argument(
        "--MPCORB", default=None, 
        help="Path to MPCORB.DAT")
    parser.add_argument(
        "--out", type=str, default="skymotion.jpg", 
        help="Figure name")
    args = parser.parse_args()
    
    if args.MPCORB is None:
        fi = MPCORB
    else:
        fi = args.MPCORB
    
    # Apophis close approach
    epoch_utc = "2029-04-13T21:52"
    epoch_jd = 2462240.4111111113

    # Extract object name and H ================================================
    obj_list = []
    H_list = []
    i, e, a = [], [], []

    with open(fi, 'r') as f:
        for line in f:
            # End of header
            if line.startswith('---'):
                break 
        for line in f:
            if len(line.strip()) == 0:
                continue
            H = line[7:14].strip()
            # N=18
            obj = line[175:193].strip()

            # Some asteroids have no H (e.g., WISE discovery)
            if H == "":
                continue
            
            # Do not plot H=99.99 for this purpose
            if H == "99.99":
                continue
            H_list.append(H)
            obj_list.append(obj)
            i.append(float(line[59:68].strip()))
            e.append(float(line[70:79].strip()))
            a.append(float(line[92:103].strip()))

            if args.Nobj is not None and len(a) >= args.Nobj:
                break  

    print(f"  N_sssbs = {len(obj_list)}")

    H_list = np.array(H_list)
    obj_list = np.array(obj_list)
    i = np.deg2rad(np.array(i))
    e = np.array(e)
    a = np.array(a)
    H = np.array(H)

    # Extract object name and H ================================================


    # Query sky motion ========================================================
    # t_elapse ~ 1 m when N = 100
    # From Kiso, Japan
    obscode = "381"
    results = []
    for obj in obj_list:
        res = fetch_ephemeris(obj, obscode, epoch_jd)
        results.append(res)
    df = pd.concat(results)
    vel = df["vel"].to_numpy()
    # Query sky motion ========================================================


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
    label_list = []
    if N_nea > 0:
        label_list.append(f"NEA (q<1.3): N={N_nea}")
    else:
        del cols[-1]
    if N_mba > 0:
        label_list.append(f"MBA (q>=1.3, a=1.8–3.3): N={N_mba}")
    else:
        del cols[-1]
    if N_hilda > 0:
        label_list.append(f"Hilda (q$\geq$1.3, a=3.7-4.0, e=0.07-0.30): N={N_hilda}")
    else:
        del cols[-1]
    if N_trojan > 0:
        label_list.append(f"Trojan (q$\geq$1.3, a=5.0-5.4): N={N_trojan}")
    else:
        del cols[-1]
    if N_tno > 0:
        label_list.append(f"TNO (q$\geq$1.3, a$\geq$30.0): N={N_tno}")
    else:
        del cols[-1]
    if N_others > 0:
        label_list.append(f"Others: N={N_others}")
    else:
        del cols[-1]

    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_axes([0.1, 0.15, 0.85, 0.75])

    ax.scatter(
        a, vel, s=1.5, color=colors, alpha=0.5, label='_nolegend_')
    # Dummy to show legends
    for label, color in zip(label_list, cols):
        ax.scatter([], [], color=color, label=label, s=20)

    ax.set_xlabel('Semimajor axis [au]', fontsize=20)
    ax.set_ylabel('Sky motion [arcsec/s]', fontsize=20)

    ax.set_title(f'{epoch_utc} UTC')
    ax.legend(fontsize=12)
    ax.set_yscale("log")
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
