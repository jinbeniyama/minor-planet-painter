#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot location of a SSSB.
"""

import os
from argparse import ArgumentParser
import numpy as np
import pandas as pd
from astropy.time import Time
from astroquery.jplhorizons import Horizons
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors


MYMARK = ["o", "^", "s", "p", "h", "<", ">", "D", "v", "P"]


def plot_ast_coord_single(
    asteroid,
    df,
    date0,
    date1,
    wid=3.5,
    quad=None,
    out=None,
):

    # Get heliocentric asteroid positions (full period)
    ast_full = Horizons(
        location="500@10",
        id=asteroid,
        epochs={
            "start": date0,
            "stop": date1,
            "step": "1d",
        },
    ).vectors()

    # Get heliocentric Earth positions (full period)
    earth_full = Horizons(
        location="500@10",
        id="399",
        epochs={
            "start": date0,
            "stop": date1,
            "step": "1d",
        },
    ).vectors()

    # Interpolation functions
    f_ast_x = interpolate.interp1d(
        ast_full["datetime_jd"],
        ast_full["x"],
        fill_value="extrapolate",
    )
    f_ast_y = interpolate.interp1d(
        ast_full["datetime_jd"],
        ast_full["y"],
        fill_value="extrapolate",
    )
    f_earth_x = interpolate.interp1d(
        earth_full["datetime_jd"],
        earth_full["x"],
        fill_value="extrapolate",
    )
    f_earth_y = interpolate.interp1d(
        earth_full["datetime_jd"],
        earth_full["y"],
        fill_value="extrapolate",
    )

    # Figure
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_axes([0.10, 0.15, 0.75, 0.75])
    cax = fig.add_axes([0.88, 0.15, 0.02, 0.75])

    # Sun
    ax.scatter(
        0, 0,
        marker="x",
        color="black",
        s=100,
        label=None,
    )

    # Axis limits
    if quad is None:
        ax.set_xlim(-wid, wid)
        ax.set_ylim(-wid, wid)

    elif quad == "Q1":
        ax.set_xlim(0, wid)
        ax.set_ylim(0, wid)

    elif quad == "Q2":
        ax.set_xlim(-wid, 0)
        ax.set_ylim(0, wid)

    elif quad == "Q3":
        ax.set_xlim(-wid, 0)
        ax.set_ylim(-wid, 0)

    elif quad == "Q4":
        ax.set_xlim(0, wid)
        ax.set_ylim(-wid, 0)

    ax.set_aspect("equal")
    ax.set_xlabel("x [au]")
    ax.set_ylabel("y [au]")

    # Calculate phase angles
    phase_list = []

    for _, row in df.iterrows():

        eph = Horizons(
            location=row["loc"],
            id=asteroid,
            epochs=row["jd"],
        ).ephemerides()

        phase = eph["alpha"][0]
        phase_list.append(phase)

        print(
            f"{row['lab']}: "
            f"phase angle = {phase:.1f} deg"
        )

    # Color map
    amin, amax = np.min(phase_list), np.max(phase_list)
    norm = colors.Normalize(vmin=amin, vmax=amax)

    mappable = cm.ScalarMappable(
        norm=norm,
        cmap=cm.viridis,
    )

    # Asteroid orbit (full period)
    ax.plot(
        ast_full["x"],
        ast_full["y"],
        lw=1.0,
        color="gray",
        alpha=0.5,
        label="Asteroid orbit",
    )

    # Earth orbit (plot only the first year)
    earth_start_jd = Time(
        date0,
        format="iso",
        scale="utc",
    ).jd

    earth_end_jd = earth_start_jd + 365.25

    earth_mask = (
        (earth_full["datetime_jd"] >= earth_start_jd)
        & (earth_full["datetime_jd"] <= earth_end_jd)
    )

    ax.plot(
        earth_full["x"][earth_mask],
        earth_full["y"][earth_mask],
        lw=1.0,
        color="gray",
        ls="dashed",
        alpha=0.5,
        label="Earth orbit",
    )

    # Observation positions
    for i, (_, row) in enumerate(df.iterrows()):

        jd = row["jd"]
        lab = row["lab"]

        phase = phase_list[i]
        col = mappable.to_rgba(phase)

        marker = MYMARK[i % len(MYMARK)]

        # Earth-asteroid connecting line
        ax.plot(
            [f_earth_x(jd), f_ast_x(jd)],
            [f_earth_y(jd), f_ast_y(jd)],
            color=col,
            lw=1.0,
            linestyle="-",
            alpha=0.6,
            zorder=8,
        )

        # Asteroid position
        ax.scatter(
            f_ast_x(jd),
            f_ast_y(jd),
            color=col,
            marker=marker,
            s=150,
            ec="k",
            zorder=10,
            label=lab,
        )

        # Earth position
        ax.scatter(
            f_earth_x(jd),
            f_earth_y(jd),
            facecolor="none",
            marker=marker,
            s=60,
            ec=col,
            lw=1.5,
            alpha=1.0,
            zorder=9,
        )

    # Legend
    ax.legend(
        ncol=3,
        fontsize=10,
        columnspacing=1.0,
    )

    # Colorbar
    fig.colorbar(
        mappable,
        cax=cax,
        label="Phase angle [deg]",
    )

    if out:
        plt.savefig(out, bbox_inches="tight")

    plt.close(fig)


def main():

    parser = ArgumentParser(
        description="Plot heliocentric asteroid orbit and observation positions."
    )

    parser.add_argument(
        "--asteroid",
        type=str,
        required=True,
        help="Target asteroid name or Horizons ID.",
    )

    parser.add_argument(
        "--date0",
        type=str,
        required=True,
        help="Start date for orbital calculation (YYYY-MM-DD).",
    )

    parser.add_argument(
        "--date1",
        type=str,
        required=True,
        help="End date for orbital calculation (YYYY-MM-DD).",
    )

    parser.add_argument(
        "--width",
        type=float,
        default=3.5,
        help="Plot half-width in AU (default: 3.5).",
    )

    parser.add_argument(
        "--obsdate",
        type=str,
        required=True,
        help="Comma-separated observation dates (YYYY-MM-DD).",
    )

    parser.add_argument(
        "--obsloc",
        type=str,
        default=None,
        help="Comma-separated Horizons observer locations (default: 500).",
    )

    parser.add_argument(
        "--quad",
        type=str,
        choices=["Q1", "Q2", "Q3", "Q4"],
        default=None,
        help="Display only the specified quadrant (Q1-Q4).",
    )

    parser.add_argument(
        "--outdir",
        type=str,
        default="fig",
        help="Output directory.",
    )

    parser.add_argument(
        "--outtype",
        type=str,
        default="pdf",
        help="Output figure format.",
    )

    args = parser.parse_args()

    # Observation dates
    obsdates = [
        x.strip() for x in args.obsdate.split(",")
    ]

    # Observer locations
    if args.obsloc is None:
        obslocs = ["500"] * len(obsdates)
    else:
        obslocs = [
            x.strip() for x in args.obsloc.split(",")
        ]

    if len(obsdates) != len(obslocs):
        parser.error(
            "The number of --obsdate and --obsloc entries must be the same."
        )

    # Convert UTC dates to Julian dates
    jd_list = Time(
        obsdates,
        format="iso",
        scale="utc",
    ).jd

    # Create DataFrame
    df = pd.DataFrame({
        "jd": jd_list,
        "lab": obsdates,
        "loc": obslocs,
    })

    # Output directory
    os.makedirs(args.outdir, exist_ok=True)

    # Output filename
    asteroid_name = args.asteroid.replace(" ", "")
    suffix = f"_{args.quad}" if args.quad else ""

    out = os.path.join(
        args.outdir,
        f"{asteroid_name}_loc{suffix}.{args.outtype}",
    )

    # Plot
    plot_ast_coord_single(
        asteroid=args.asteroid,
        df=df,
        date0=args.date0,
        date1=args.date1,
        wid=args.width,
        quad=args.quad,
        out=out,
    )

    print(f"Saved: {out}")


if __name__ == "__main__":
    main()

