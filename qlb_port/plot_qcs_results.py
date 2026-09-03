"""Plot QPU vs emulator vs exact for the 2-site Dirac-QLB hardware run.

Reads qcs_results.json (written by qcs_quil/run_qcs.py after a Cepheus-1 run) and
draws the velocity and density trajectories in the same colour code as the other
hardware figures: exact grey, emulator blue (C0), hardware red (C3).
"""

import argparse
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results", help="qcs_results.json from run_qcs.py")
    ap.add_argument("--out", default="qlb_port/hw_2site_qcs.png")
    args = ap.parse_args()

    d = json.load(open(args.results))
    jobs = d["jobs"]
    alpha = sorted((j for j in jobs if j["kind"] == "alpha"), key=lambda j: j["t"])
    dens = sorted((j for j in jobs if j["kind"] == "density"), key=lambda j: j["t"])
    ta = [j["t"] for j in alpha]
    td = [j["t"] for j in dens]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))

    ax1.plot(ta, [j["exact"] for j in alpha], "-o", color="0.5", lw=1.5, ms=6,
             label="exact", zorder=2)
    ax1.plot(ta, [j["emulator"] for j in alpha], "s", color="C0", mfc="none", mew=1.7,
             ms=9, label="emulator (QVM)", zorder=3)
    ax1.plot(ta, [j["qpu"] for j in alpha], "D", color="C3", ms=8,
             label="Cepheus-1 QPU", zorder=4)
    ax1.axhline(0, color="k", lw=0.6, ls=":")
    ax1.set_xlabel("step  $t$")
    ax1.set_ylabel(r"$\langle\alpha_x\rangle(t)$")
    ax1.set_title("Velocity (Zitterbewegung trembling)")
    ax1.set_xticks(ta)
    ax1.set_ylim(-1.2, 1.2)
    ax1.legend(frameon=False, loc="lower left")

    ax2.plot(td, [j["exact"] for j in dens], "-o", color="0.5", lw=1.5, ms=6,
             label="exact", zorder=2)
    ax2.plot(td, [j["emulator"] for j in dens], "s", color="C0", mfc="none", mew=1.7,
             ms=9, label="emulator (QVM)", zorder=3)
    ax2.plot(td, [j["qpu"] for j in dens], "D", color="C3", ms=8,
             label="Cepheus-1 QPU", zorder=4)
    ax2.set_xlabel("step  $t$")
    ax2.set_ylabel(r"$\langle x\rangle(t)$")
    ax2.set_title("Position (density sloshing)")
    ax2.set_xticks(td)
    ax2.set_ylim(0, 1)
    ax2.legend(frameon=False, loc="upper right")

    dev = d.get("device", "QPU")
    us = d.get("total_execution_us", 0.0)
    fig.suptitle(f"Two-site Dirac QLB on {dev}  "
                 f"($N=2$, $\\tilde m={d['mass']}$, {d.get('shots', '?')} shots, "
                 f"{us / 1e6:.2f}s QPU)", y=1.02)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print("wrote", args.out)


if __name__ == "__main__":
    main()
