"""
Assemble the four-site (n_pos=2) Fourier-streaming hardware runs into the paper
figures and print the numbers used in the tables of the adder-streaming paper.

Reads the saved data (no cloud, no credits):
    qlb_port/hw_4site_fourier_density.npz      packet density rho(x,t), pinned qubits
    qlb_port/hw_4site_fourier_trembling.npz    velocity <alpha_x>(t), pinned qubits
    qlb_port/hw_4site_fourier_q19_density.npz  first density run (automatic placement)
    qlb_port/hw_4site_t0.npz                   t=0 preparation/readout check

Writes, in the style of plot_2site_results.py (exact grey, emulator blue, hardware red):
    <outdir>/hw_4site_results.png   (a) mean position, (b) velocity trembling
    <outdir>/hw_4site_density.png   per-site density rho(x) at t=0..4

Usage:
    PYTHONPATH=. python3 -m qlb_port.plot_4site_results
"""

import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def tvd(p, q):
    return 0.5 * float(np.abs(np.asarray(p) - np.asarray(q)).sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default="proposal/quairk/ancillaStreaming/figures")
    args = ap.parse_args()

    DD = np.load("qlb_port/hw_4site_fourier_density.npz")
    DT = np.load("qlb_port/hw_4site_fourier_trembling.npz")
    DQ = np.load("qlb_port/hw_4site_fourier_q19_density.npz")
    D0 = np.load("qlb_port/hw_4site_t0.npz")

    ts = DD["t"]
    N = DD["rho_exact"].shape[1]
    xs = np.arange(N)

    # ---- tables ----
    print("density rho(x,t)  (pinned run):")
    for i, t in enumerate(ts):
        print(f"  t={t}  exact={np.round(DD['rho_exact'][i], 3).tolist()}  "
              f"emu={np.round(DD['rho_aer'][i], 3).tolist()}  hw={np.round(DD['rho_hw'][i], 3).tolist()}"
              f"  TVD(emu)={tvd(DD['rho_exact'][i], DD['rho_aer'][i]):.3f}"
              f"  TVD(hw)={tvd(DD['rho_exact'][i], DD['rho_hw'][i]):.3f}"
              f"  <x> ex/emu/hw={DD['x_exact'][i]:.3f}/{DD['x_aer'][i]:.3f}/{DD['x_hw'][i]:.3f}")
    print("density, automatic placement (first run):")
    for i, t in enumerate(DQ["t"]):
        r = DQ["rho_hw"][i]
        print(f"  t={t}  hw={np.round(r, 3).tolist()}  TVD={tvd(DQ['rho_exact'][i], r):.3f}  "
              f"P(x0=1) exact={DQ['rho_exact'][i][1::2].sum():.3f} hw={r[1::2].sum():.3f}")
    print(f"t=0 check: rho hw={np.round(D0['rho_hw'], 3).tolist()}  "
          f"TVD={tvd(D0['rho_exact'], D0['rho_hw']):.3f}  <alpha_x> hw={float(D0['av_hw']):+.3f}")

    E = float(DT["E"])
    ex, hw = DT["av_exact"], DT["av_hw"]
    A = np.column_stack([np.ones_like(ex), ex])
    (off, scale), *_ = np.linalg.lstsq(A, hw, rcond=None)
    print(f"\nvelocity <alpha_x>(t)  2E={2 * E:.3f}  period={np.pi / E:.2f} steps:")
    for i, t in enumerate(DT["t"]):
        print(f"  t={t}  exact={ex[i]:+.3f}  emu={DT['av_aer'][i]:+.3f}  hw={hw[i]:+.3f}")
    print(f"  fit hw = offset + scale*exact:  offset={off:+.3f}  scale={scale:.3f}")

    # ---- figure 1: (a) mean position, (b) velocity trembling ----
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.3))
    axA.plot(ts, DD["x_exact"], "-o", color="0.5", lw=1.3, ms=5, zorder=2,
             label=r"exact  $\langle x\rangle$")
    axA.plot(ts, DD["x_aer"], "s", color="C0", ms=8, mfc="none", mew=1.7, zorder=3,
             label="emulator")
    axA.plot(ts, DD["x_hw"], "D", color="C3", ms=8, zorder=4, label="hardware")
    axA.axhline((N - 1) / 2, color="k", lw=0.6, ls=":")
    axA.set_xticks(ts)
    axA.set_xlabel("step  $t$")
    axA.set_ylabel(r"mean position  $\langle x\rangle$")
    axA.set_ylim(-0.05, N - 0.95)
    axA.set_title("(a) packet position")
    axA.legend(frameon=False, loc="lower left", ncol=1)

    tv = DT["t"]
    tsm = np.linspace(0, float(tv.max()), 400)
    M = np.column_stack([np.ones_like(tv, float), np.cos(2 * E * tv), np.sin(2 * E * tv)])
    coef, *_ = np.linalg.lstsq(M, ex, rcond=None)
    avsm = coef[0] + coef[1] * np.cos(2 * E * tsm) + coef[2] * np.sin(2 * E * tsm)
    axB.plot(tsm, avsm, color="0.6", lw=1.6, label=r"exact  $\sin 2Et$")
    axB.plot(tv, ex, "o", color="0.4", ms=5)
    axB.plot(tv, DT["av_aer"], "s", color="C0", ms=8, mfc="none", mew=1.7, label="emulator")
    axB.plot(tv, hw, "D", color="C3", ms=8, label="hardware")
    axB.axhline(0, color="k", lw=0.6, ls=":")
    axB.set_xticks(tv)
    axB.set_xlabel("step  $t$")
    axB.set_ylabel(r"$\langle\alpha_x\rangle(t)$")
    axB.set_ylim(-1.18, 1.18)
    axB.set_title("(b) velocity trembling")
    axB.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    out1 = f"{args.outdir}/hw_4site_results.png"
    fig.savefig(out1, dpi=150, bbox_inches="tight")

    # ---- figure 2: per-site density at every recorded step ----
    fig, axs = plt.subplots(1, len(ts), figsize=(13, 3.1), sharey=True)
    w = 0.27
    for i, (ax, t) in enumerate(zip(axs, ts)):
        ax.bar(xs - w, DD["rho_exact"][i], w, color="0.65", label="exact")
        ax.bar(xs, DD["rho_aer"][i], w, color="C0", label="emulator")
        ax.bar(xs + w, DD["rho_hw"][i], w, color="C3", label="hardware")
        ax.set_xticks(xs)
        ax.set_xlabel("site  $x$")
        ax.set_title(rf"$t={t}$" + f"   (TVD {tvd(DD['rho_exact'][i], DD['rho_hw'][i]):.2f})",
                     fontsize=10)
        ax.set_ylim(0, 0.72)
    axs[0].set_ylabel(r"$\rho(x)$")
    axs[0].legend(frameon=False, loc="upper left", fontsize=8)
    fig.tight_layout()
    out2 = f"{args.outdir}/hw_4site_density.png"
    fig.savefig(out2, dpi=150, bbox_inches="tight")
    print(f"\nwrote {out1} and {out2}")


if __name__ == "__main__":
    main()
