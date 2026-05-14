#!/usr/bin/env python3
"""Make a figure of Psid vs. Psyn.

Note on Period Analysis:
The period derived from long-term lightcurve observations (e.g., via Lomb-Scargle)
is the synodic period (Psyn), not the intrinsic sidereal period (Psid).

When the spin pole is nearly perpendicular to the ecliptic plane (high ecliptic
latitude), the viewing geometry remains relatively stable. In such cases,
it is possible to obtain a 'clean' folded lightcurve using Psyn, as it
compensates for the phase shift caused by the Earth's orbital motion.
"""
import argparse
import numpy as np
import matplotlib.pyplot as plt

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=str, default="Psid_vs_Psyn_final_fixed.png")
    return parser.parse_args()

def main(args=None):
    if args is None:
        args = get_args()

    P_orb_earth = 1.0
    P_orb_ast = 5.0
    P_sid = 0.01
    aE, aA = 1.0, 2.5

    earth_phases_deg = [0, 90, 180, 270, 360, 500, 600, 700, 800]
    earth_phases_deg = [0, 90, 180, 270]
    num_epochs = len(earth_phases_deg)
    t_list = [(deg / 360.0) * P_orb_earth for deg in earth_phases_deg]

    cmap = plt.get_cmap('tab10')
    colors = [cmap(i % 10) for i in range(num_epochs)]

    obs_t = []
    obs_flux = []

    for i, t_start in enumerate(t_list):
        t_obs = np.linspace(t_start, t_start + 0.04, 200)
        thetaE = 2 * np.pi * t_obs / P_orb_earth
        thetaA = 2 * np.pi * t_obs / P_orb_ast
        alpha = np.arctan2(aE*np.sin(thetaE) - aA*np.sin(thetaA),
                           aE*np.cos(thetaE) - aA*np.cos(thetaA))

        flux = 1.0 + 0.25 * np.cos(2 * (2 * np.pi * t_obs / P_sid - alpha))
        flux += np.random.normal(0, 0.005, size=len(flux))
        obs_t.append(t_obs)
        obs_flux.append(flux)

    # Periodic analysis
    all_t = np.concatenate(obs_t)
    all_f = np.concatenate(obs_flux)

    def find_best_p(p_min, p_max, steps):
        ps = np.linspace(p_min, p_max, steps)
        min_d = float('inf')
        best = ps[0]
        for p in ps:
            phase = (all_t / p) % 1.0
            idx = np.argsort(phase)
            diff = np.sum(np.diff(all_f[idx])**2)
            if diff < min_d:
                min_d = diff
                best = p
        return best
    
    # Should find local minimum
    p_coarse = find_best_p(P_sid*0.96, P_sid * 1.04, 10000)
    p_coarse2 = find_best_p(p_coarse*0.98, p_coarse * 1.02, 10000)
    best_p = find_best_p(p_coarse2 * 0.98, p_coarse2 * 1.02, 1000000)

    plt.rcParams.update({'font.size': 12})
    fig = plt.figure(figsize=(18, 8))

    ax_lc = fig.add_axes([0.07, 0.15, 0.26, 0.70])
    ax_orb = fig.add_axes([0.38, 0.15, 0.26, 0.70])
    ax_fold = fig.add_axes([0.70, 0.15, 0.26, 0.70])

    for i in range(num_epochs):
        ax_lc.scatter((obs_t[i] / P_sid) % 1.0, obs_flux[i], s=3, color=colors[i], alpha=0.5)
    ax_lc.set_title(f"Folded by Sidereal ($P_{{sid}}$)")
    ax_lc.set_xlabel("Phase")
    ax_lc.set_ylabel("Relative Brightness")
    ax_lc.grid(alpha=0.2)
    ax_lc.text(0.05, 0.92, f"$P_{{sid}}={P_sid:.6f}$", transform=ax_lc.transAxes,
               bbox=dict(facecolor='white', alpha=0.8), fontsize=12)

    t_bg = np.linspace(0, max(t_list), 1000)
    ax_orb.plot(aE*np.cos(2*np.pi*t_bg), aE*np.sin(2*np.pi*t_bg), "--b", alpha=0.1, label="Earth Orbit")
    ax_orb.plot(aA*np.cos(2*np.pi*t_bg/P_orb_ast), aA*np.sin(2*np.pi*t_bg/P_orb_ast), "--r", alpha=0.1, label="Asteroid Orbit")
    ax_orb.scatter(0, 0, c="gold", s=200, edgecolors="orange", zorder=10, label="Sun")

    for i, t in enumerate(t_list):
        xE, yE = aE*np.cos(2*np.pi*t), aE*np.sin(2*np.pi*t)
        xA, yA = aA*np.cos(2*np.pi*t/P_orb_ast), aA*np.sin(2*np.pi*t/P_orb_ast)
        l_e = "Earth" if i == 0 else ""
        l_a = "Asteroid" if i == 0 else ""
        l_los = "Line of Sight" if i == 0 else ""

        ax_orb.scatter(xE, yE, color=colors[i], s=40, edgecolors="k", zorder=6, label=l_e)
        ax_orb.scatter(xA, yA, color=colors[i], s=80, marker="^", edgecolors="k", zorder=6, label=l_a)
        ax_orb.plot([xE, xA], [yE, yA], color=colors[i], alpha=0.3, lw=1, label=l_los)

    ax_orb.set_title("Orbital Configuration")
    ax_orb.set_aspect("equal")
    ax_orb.legend(loc='upper right', fontsize=9, framealpha=0.9)
    ax_orb.set_xlabel("X (AU)")
    ax_orb.set_ylabel("Y (AU)")

    for i in range(num_epochs):
        ax_fold.scatter((obs_t[i] / best_p) % 1.0, obs_flux[i], s=3, color=colors[i], alpha=0.5)
    ax_fold.set_title(f"Folded by Synodic ($P_{{syn}}$)")
    ax_fold.set_xlabel("Phase")
    ax_fold.grid(alpha=0.2)
    ax_fold.text(0.05, 0.92, f"$P_{{syn}}={best_p:.8f}$", transform=ax_fold.transAxes,
                 bbox=dict(facecolor='orange', alpha=0.2), fontsize=12, fontweight='bold')

    plt.savefig(args.out, dpi=200, bbox_inches='tight')
    print(f"Saved: {args.out}")
    print(f"Final Psyn detected: {best_p:.8f}")

if __name__ == "__main__":
    main()
