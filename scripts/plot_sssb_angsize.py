#!/usr/bin/env python3
"""Make a figure of angular size of minor bodies.

TODO:
Plot angular size as seen from the Earth.
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
from asteropy.constants import au_km

from minor_planet_painter.common import (
    MPCORB, jd2utc, utc2jd, 
    solve_kepler_eq, mpcepoch2jd
    )


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



def ang_size(D, d):
    """Calculate angular size.

    Parameters
    --------------
    D : float
        diameter in km
    d : float
        distance between target and observer in au

    Return
    --------
    theta : float
        angular distance in arcsec
    """
    theta = np.tan(D/(d*au_km))*(180/np.pi)*3600.
    return theta


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
        "--MPCORB", default=None, 
        help="Path to MPCORB.DAT")
    parser.add_argument(
        "--out", type=str, default="angsize.jpg", 
        help="Figure name")
    args = parser.parse_args()
    
    if args.MPCORB is None:
        fi = MPCORB
    else:
        fi = args.MPCORB

    # parse UTC
    t_utc_dt = datetime.fromisoformat(args.epoch)
    t_utc_iso = t_utc_dt.isoformat()


    # Calculate angular distance ==============================================
    M, omega, Omega, i, e, a = [], [], [], [], [], []
    n = []
    H = []
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
            H_obj = line[7:14].strip()
            # Some asteroids have no H (e.g., WISE discovery)
            if H_obj == "":
                continue
            
            # Do not plot H=99.99 for this purpose
            if H_obj == "99.99":
                continue
            H.append(H_obj)
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
    H = np.array(H)

    jd_now = utc2jd(t_utc_iso)

    epoch_jd = np.array([mpcepoch2jd(e) for e in epoch])
    dt_days = jd_now - epoch_jd
    # Linear assumption
    M_t = M + n*dt_days

    ejd_min = float(np.min(epoch_jd))
    ejd_max = float(np.max(epoch_jd))
    print("  Check epochs:")
    print(f"    Minimum (earliest): {jd2utc(ejd_min)}")
    print(f"    Maximum (latest)  : {jd2utc(ejd_max)}")
    
    # Calculate eccentric anomaly
    E = solve_kepler_eq(M_t, e)
    # Calculate nu and r 
    nu = 2*np.arctan2(np.sqrt(1+e)*np.sin(E/2), np.sqrt(1-e)*np.cos(E/2))
    nu_total = np.mod(omega + nu, 2*np.pi)
    r = a * (1 - e*np.cos(E))

    # Calculate angular size with H (absolute mag) and r (distance)
    # Assume pv = 0.100
    pv = 0.100
    D, _ = D_from_Hp(H, pv, Herr=0, perr=0)

    angsize = ang_size(D, r)
    # Calculate angular distance ==============================================
     

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

    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_axes([0.1, 0.15, 0.85, 0.75])

    ax.scatter(r, angsize, s=1.5, color=colors, alpha=0.5, label='_nolegend_')
    for label, color in zip(legend_labels, cols):
        ax.scatter([], [], color=color, label=label, s=20)

    ax.set_xlabel('Distance from the Sun [au]', fontsize=20)
    # To be updated to from the Earth
    ax.set_ylabel('Angular size as seen from the Sun [arcsec]', fontsize=20)

    ax.set_title(f'{t_utc_iso} UTC')
    ax.legend(fontsize=12)
    ax.set_yscale("log")
    plt.show()
    
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
