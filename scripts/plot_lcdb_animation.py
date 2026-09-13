#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Access LCDB text file and make animation.
"""
from argparse import ArgumentParser as ap
import pandas as pd
import re
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np


def handle_LCDB(lcdb, tumbler=False):
    """
    Obtain all period from lc from lc_summary_pub.txt.
    Note:
        Can be used for LCDB2023Apr, LCDB2023Oct, etc.
  
    Parameters
    ----------
    lcdb : str
        path to lcdb file
    tumbler : bool, optional
        tumbler only option
  
    Return
    ------
    df : pandas.DataFrame
        including object name, absolute magnitude H, Diameter[km], Period[h]
    """
    # Number of headers
    n_header = 5
    # Number of skipped objects
    n_skip = 0
  
    name_list, stype_list, s_stype_list, fam_list, D_list = [], [], [], [], []
    num_list = []
    s_D_list, H_list, s_H_list, G_list, s_G_list = [], [], [], [], []
    pv_list, s_pv_list, f_pv_list, rotP_list, f_rotP_list = [], [], [], [], []
    amin_list, amax_list, f_a_list, U_list, bin_list = [], [], [], [], []
    pol_list, sur_list, pri_list = [], [], []
  
    with open(lcdb, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()[n_header:]
        for line in lines:
            # Read p24, 4.1.3 in readme.pdf
            # Use NAME, not DESIG, since numbered objects have no designation. 
            obj = line[10:10+30].strip()
            obj = obj.replace(" ", "")
            num = line[0:7].strip()

            fam = line[62:62+8].strip()
            s_stype = line[71].strip()
            stype = line[73:73+10].strip()
            s_D = line[86].strip()

            D = line[88:88+8].strip()
            s_H = line[97].strip()
            H = line[99:99+6].strip()
            s_G = line[109].strip()
            G = line[111:111+6].strip()
            G1 = line[118:118+6].strip()
            G2 = line[125:125+6].strip()
            
            # Dummies (for e.g., 67P)
            if D == "":
                D = -99
            if H == "":
                H = -99

            s_pv = line[132].strip()
            f_pv = line[134].strip()
            pv = line[136:136+10].strip()

            f_rotP = line[143].strip()
            rotP = line[145:145+13].strip()

            f_a = line[175].strip()
            amin = line[177:177+4].strip()
            amax = line[182:182+4].strip()

            # Substitute 0 when no amin/amax
            if amin == "":
                amin = 0
            if amax == "":
                amax = 0
            
            # U (quality) flag
            U = line[187:187+2].strip()

            sign_binary = line[196:196+3].strip()
            if sign_binary == "":
                sign_binary = 0
            sign_pole = line[200:200+3].strip()
            if sign_pole == "":
                sign_pole = 0
            sign_survey = line[204:204+3].strip()
            if sign_survey == "":
                sign_survey = 0
            sign_private = line[214:214+3].strip()
            if sign_private == "":
                sign_private = 0
                
            if not (obj and fam and rotP and U):
                n_skip += 1
                continue
  
            name_list.append(obj)
            num_list.append(num)
            stype_list.append(stype)
            s_stype_list.append(s_stype)
            fam_list.append(fam)
            D_list.append(D)

            s_D_list.append(s_D)
            H_list.append(H)
            s_H_list.append(s_H)
            G_list.append(G)
            s_G_list.append(s_G)

            pv_list.append(pv)
            s_pv_list.append(s_pv)
            f_pv_list.append(f_pv)
            rotP_list.append(rotP)
            f_rotP_list.append(f_rotP)

            amin_list.append(amin)
            amax_list.append(amax)
            f_a_list.append(f_a)
            U_list.append(U)
            bin_list.append(sign_binary)

            pol_list.append(sign_pole)
            sur_list.append(sign_survey)
            pri_list.append(sign_private)
    
    info = dict(
        obj=name_list, num=num_list, stype=stype_list, stypesource=s_stype_list, 
        family=fam_list, H=H_list, 
        Hsource=s_H_list, G=G_list, Gsource=s_G_list, pv=pv_list, pvsource=s_pv_list, 
        pvflag=f_pv_list, D=D_list, Dsource=s_D_list, P=rotP_list, Pflag=f_rotP_list, 
        amin=amin_list, amax=amax_list, aflag=f_a_list, U=U_list, Binary=bin_list, 
        Pole=pol_list, Survey=sur_list, Private=pri_list)
    
    df = pd.DataFrame(info.values(), index=info.keys()).T
    print(f"n_df, nskip = {len(df)}, {n_skip}")
    return df, n_skip


if __name__ == "__main__":
    parser = ap(description="Handling LCDB text file")
    parser.add_argument(
        "file", type=str,
        help="path to text file"
    )
    args = parser.parse_args()

    df, n_skip = handle_LCDB(args.file)

    df['D'] = pd.to_numeric(df['D'], errors='coerce')
    df['P'] = pd.to_numeric(df['P'], errors='coerce')
    df['U'] = pd.to_numeric(df['U'], errors='coerce')

    # 1. Unit conversion
    # Diameter D [km] -> [m]
    df['D_m'] = df['D'] * 1000.0
    # Rotation period P [hour] -> [s]
    df['P_s'] = df['P'] * 3600.0


    # 2. Function to extract discovery year from object name
    def get_discovery_year(obj_name):
        if not isinstance(obj_name, str):
            return 1800  # Default for initial display

        # Search for a 4-digit year near the beginning regardless of spacing
        match = re.search(r'^(?:[A-Za-z0-9]+\s*)?(1[89]\d{2}|20\d{2})', obj_name)
        if not match:
            # Fallback: search for the first 4-digit year anywhere in the string
            match = re.search(r'\b(1[89]\d{2}|20\d{2})\b', obj_name)

        if match:
            return int(match.group(1))

        # Names without years (e.g., Ceres) are treated as pre-1990
        return 1800


    # Add discovery year column
    df['DiscoveryYear'] = df['obj'].apply(get_discovery_year)

    # 3. Animation configuration
    fig, ax = plt.subplots(figsize=(15, 7.5))

    # Adjust margins to eliminate excess top space
    plt.subplots_adjust(left=0.07, right=0.97, top=0.92, bottom=0.12)

    start_year = 1990
    end_year = 2023
    years = range(start_year, end_year + 1)

    # Threshold for quality flag U
    U_th = 2

    def update(frame_year):
        ax.clear()

        # Ensure grid lines are drawn below data points
        ax.set_axisbelow(True)

        # Filter by discovery year and quality flag U >= U_th
        u_numeric = pd.to_numeric(df['U'], errors='coerce')
        subset = df[
            ((df['DiscoveryYear'] <= frame_year) | (df['DiscoveryYear'] < 1990)) & 
            (u_numeric >= U_th)
        ]

        # Drop rows with missing values in D_m or P_s (also drop dummy -99000m)
        subset = subset.dropna(subset=['D_m', 'P_s'])
        subset = subset[(subset['D_m'] > 0) & (subset['P_s'] > 0)]

        # Split into NEA (family 9101) and others
        is_nea = (subset['family'] == '9101') | (subset['family'] == 9101)
        nea_df = subset[is_nea]
        other_df = subset[~is_nea]

        # 1. Plot other asteroids (gray crosses)
        ax.scatter(
            other_df['D_m'],
            other_df['P_s'],
            marker='x',
            color='gray',
            alpha=0.5,
            s=30,
            label='Other asteroids'
        )

        # 2. Plot NEA (dark red open circles) - family = 9101
        ax.scatter(
            nea_df['D_m'],
            nea_df['P_s'],
            marker='o',
            facecolors='none',
            edgecolors='darkred',
            linewidths=1.5,
            alpha=0.8,
            s=45,
            label='NEA (family=9101)'
        )

        # Set log scale for both axes
        ax.set_xscale('log')
        ax.set_yscale('log')

        # Fix axis limits (X: 1 m ~ 10^6 m, Y: 1 s ~ 10^6 s, inverted Y)
        ax.set_xlim(1, 1e6)
        ax.set_ylim(1e6, 1)

        # Horizontal dashed line at P = 2 hours (7200 seconds)
        p_2h_sec = 2.0 * 3600.0
        ax.axhline(p_2h_sec, color='green', linestyle='--', linewidth=2, alpha=1.0, label='P = 2 h')

        # Axis labels and title
        ax.set_xlabel('Diameter [m]', fontsize=18)
        ax.set_ylabel('Rotation period [s]', fontsize=18)
        ax.set_title(
            f'Year: {frame_year} (N={len(subset)})',
            fontsize=20,
            fontweight='bold'
        )

        # Tick parameters
        ax.tick_params(axis='both', which='major', labelsize=15)
        ax.tick_params(axis='both', which='minor', labelsize=12)

        # Legend configuration
        ax.legend(loc='upper right', frameon=True, fontsize=15)

        # Grid lines
        ax.grid(True, which='both', linestyle='--', alpha=0.5)


    # Create animation with 1 frame per second (interval = 1000 ms)
    ani = animation.FuncAnimation(
        fig, update, frames=years, interval=1000, repeat=False
    )

    # 4. Save as MP4 (DPI = 200)
    output_file = 'asteroid_rotation_animation.mp4'
    print(f'Saving animation to {output_file}...')
    ani.save(output_file, writer='ffmpeg', fps=1, dpi=200)
    print('Done!')

    plt.close()
