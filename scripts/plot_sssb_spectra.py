#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import pandas as pd
from argparse import ArgumentParser as ap
import matplotlib.pyplot as plt

import classy


def plot_spec(irtf_data_list, stypes=None, plotcolors=None, xnorm=1.0, xrange=(0.7, 2.5), yrange=(0.8, 1.4), out=None):
    """
    Plot reflectance spectra (IRTF) with taxonomic templates
    
    Parameters
    ----------
    irtf_data_list : list of tuple
        List containing (df_irtf, label) pairs.
    stypes : list of str, optional
        Taxonomic types to plot (e.g., ["S", "Q", "B"]).
    plotcolors : list of str, optional
        Colors for each observation data.
    xnorm : float
        Normalization wavelength.
    xrange : tuple of float
        Wavelength range (xmin, xmax).
    yrange : tuple of float
        Reflectance range (ymin, ymax).
    out : str, optional
        Output filename (pdf).
    """
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_axes([0.15, 0.19, 0.8, 0.76])

    # -------------------------------------------------------------------------
    # Taxonomic Templates
    # -------------------------------------------------------------------------
    if stypes:
        temp = classy.taxonomies.mahlke.load_templates()
        w_norm_vis = xnorm  
        
        col_stype = {
            "A": "green", 
            "S": "red", 
            "Q": "orange", 
            "D": "red", 
            "L": "orange", 
            "C": "blue", 
            "X": "gray", 
            "K": "orange", 
            "V": "purple",
            "B": "dodgerblue"
        }

        for st in stypes:
            st_upper = st.upper()
            if st_upper in temp:
                sp = temp[st_upper]
                col = col_stype.get(st_upper, "gray")
                label = f"{st_upper}-type"
                sp.normalize(at=w_norm_vis)
                ax.plot(sp.wave, sp.refl, lw=1, alpha=0.1, color="gray")
                ax.fill_between(
                    sp.wave, sp.refl - sp.refl_err, sp.refl + sp.refl_err,
                    color=col, alpha=0.2, label=label
                )

    # -------------------------------------------------------------------------
    # IRTF Observation Data (Multiple)
    # -------------------------------------------------------------------------
    default_colors = ["black", "brown", "gray", "magenta", "cyan"]

    for i, (df_irtf, objlabel) in enumerate(irtf_data_list):
        if plotcolors and i < len(plotcolors):
            col_irtf = plotcolors[i]
        else:
            col_irtf = default_colors[i % len(default_colors)]

        ax.errorbar(
            df_irtf["w"],
            df_irtf["ref"],
            yerr=df_irtf["referr"],
            fmt="o",
            markersize=3,
            linewidth=0.5,
            capsize=0,
            label=objlabel,
            color=col_irtf,
            zorder=3,
        )

    # Labels and limits
    ax.set_xlabel(r"Wavelength [$\mu$m]")
    ax.set_ylabel("Reflectance")
    ax.set_ylim(yrange)

    ax.set_xlim(xrange)

    ax.legend(loc="upper left", fontsize=8)
    ax.grid(alpha=0.3)

    if out is not None:
        plt.savefig(out)
        plt.close(fig)
    else:
        plt.show()


if __name__ == "__main__":
    parser = ap(description="Plot reflectance spectra from IRTF observations with templates")
    parser.add_argument(
        "f_irtf", type=str, nargs="+",
        help="IRTF observation file(s)"
    )
    parser.add_argument(
        "--label", type=str, nargs="*", default=None,
        help="Label(s) for the observation file(s)"
    )
    parser.add_argument(
        "--plotcol", type=str, nargs="*", default=None,
        help="Color(s) for the observation file(s) (e.g., black black)"
    )
    parser.add_argument(
        "--stype", type=str, nargs="*", default=None,
        help="Taxonomic templates to plot (e.g., S Q B)"
    )
    parser.add_argument(
        "--xnorm", type=float, default=1.0,
        help="Normalization wavelength (default: 1.0 um)"
    )
    parser.add_argument(
        "--xrange", type=float, nargs=2, default=[0.7, 2.5],
        help="Wavelength range [xmin, xmax] (default: 0.7 2.5)"
    )
    parser.add_argument(
        "--yrange", type=float, nargs=2, default=[0.8, 1.4],
        help="Reflectance range [ymin, ymax] (default: 0.8 1.4)"
    )
    parser.add_argument(
        "--out", type=str, default=None,
        help="Output filename"
    )
    parser.add_argument(
        "--outdir", type=str, default="fig",
        help="Output directory name"
    )
    args = parser.parse_args()

    outdir = args.outdir
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    cols = ["w", "ref", "referr", "flag1", "flag2", "flag3", "flag4"]

    # Handle multiple labels
    labels = args.label if args.label is not None else []
    if len(labels) < len(args.f_irtf):
        labels = list(labels) + [None] * (len(args.f_irtf) - len(labels))

    irtf_data_list = []

    for f_path, lbl in zip(args.f_irtf, labels):
        # Load IRTF reflectance
        df_irtf = pd.read_csv(
            f_path,
            sep=r'\s+',
            names=cols,
            header=None,
            engine="python"
        )
        
        # Normalize at specified wavelength (xnorm ± 0.01 um)
        xnorm = args.xnorm
        norm_value = df_irtf.loc[
            (df_irtf["w"] > xnorm - 0.01) & (df_irtf["w"] < xnorm + 0.01),
            "ref"
        ].mean()
        df_irtf["ref"] /= norm_value
        df_irtf["referr"] /= norm_value

        if lbl is None:
            lbl = os.path.basename(f_path)

        irtf_data_list.append((df_irtf, lbl))

    # Plot output filename
    if args.out:
        out = args.out
    else:
        out = "sssb_spec.pdf"
    out = os.path.join(outdir, out)

    plot_spec(
        irtf_data_list, 
        stypes=args.stype, 
        plotcolors=args.plotcol,
        xnorm=args.xnorm, 
        xrange=tuple(args.xrange), 
        yrange=tuple(args.yrange), 
        out=out
    )
