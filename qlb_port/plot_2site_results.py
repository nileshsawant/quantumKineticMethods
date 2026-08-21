"""
Assemble the two two-site hardware runs into the paper figure
(figures/hw_2site_results.png) and print the position and velocity tables.

Reads the saved data (no cloud, no credits):
    qlb_port/hw_2site_density.npz    position sloshing: t, x_exact / x_aer / x_hw, rho_*
    qlb_port/hw_2site_trembling.npz  velocity trembling: t, av_exact / av_aer / av_hw, E

Usage:
    PYTHONPATH=. python3 -m qlb_port.plot_2site_results
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DD = np.load("qlb_port/hw_2site_density.npz")
DT = np.load("qlb_port/hw_2site_trembling.npz")
OUT = "proposal/quairk/qlbDiracCircuits/figures/hw_2site_results.png"


def tvd(p, q):
    return 0.5 * float(np.abs(np.asarray(p) - np.asarray(q)).sum())


# ---- tables (the numbers that go into the paper) ----
tsd = DD["t"]
print("position  <x>(t) = rho(1):")
for i, t in enumerate(tsd):
    print(f"  t={t}  exact={DD['x_exact'][i]:+.3f}  emu={DD['x_aer'][i]:+.3f}  "
          f"hw={DD['x_hw'][i]:+.3f}")
print(f"  t=0 density TVD exact-hw = {tvd(DD['rho_exact'][0], DD['rho_hw'][0]):.3f}")

E = float(DT["E"])
tsv = DT["t"]
print(f"\nvelocity  <alpha_x>(t)   2E={2 * E:.3f}  period={2 * np.pi / (2 * E):.2f} steps:")
for i, t in enumerate(tsv):
    print(f"  t={t}  exact={DT['av_exact'][i]:+.3f}  emu={DT['av_aer'][i]:+.3f}  "
          f"hw={DT['av_hw'][i]:+.3f}")

# ---- figure: (a) position sloshing, (b) velocity trembling ----
fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.3))

# (a) mean position <x>(t): period-two sloshing of a localized packet
axA.plot(tsd, DD["x_exact"], "-o", color="0.5", lw=1.3, ms=5, zorder=2,
         label=r"exact  $\langle x\rangle$")
axA.plot(tsd, DD["x_aer"], "s", color="C0", ms=8, mfc="none", mew=1.7, zorder=3,
         label="emulator")
axA.plot(tsd, DD["x_hw"], "D", color="C3", ms=8, zorder=4, label="hardware")
axA.axhline(0.5, color="k", lw=0.6, ls=":")
axA.set_xticks(tsd)
axA.set_xlabel("step  $t$")
axA.set_ylabel(r"mean position  $\langle x\rangle=\rho(1)$")
axA.set_ylim(-0.03, 1.03)
axA.set_title("(a) position sloshing")
axA.legend(frameon=False, loc="lower center", ncol=3)

# (b) mean velocity <alpha_x>(t): the trembling, one full period
tsm = np.linspace(0, float(tsv.max()), 400)
Amat = np.column_stack([np.ones_like(tsv, float), np.cos(2 * E * tsv), np.sin(2 * E * tsv)])
coef, *_ = np.linalg.lstsq(Amat, DT["av_exact"], rcond=None)
avsm = coef[0] + coef[1] * np.cos(2 * E * tsm) + coef[2] * np.sin(2 * E * tsm)
axB.plot(tsm, avsm, color="0.6", lw=1.6, label=r"exact  $\sin 2Et$")
axB.plot(tsv, DT["av_exact"], "o", color="0.4", ms=5)
axB.plot(tsv, DT["av_aer"], "s", color="C0", ms=8, mfc="none", mew=1.7, label="emulator")
axB.plot(tsv, DT["av_hw"], "D", color="C3", ms=8, label="hardware")
axB.axhline(0, color="k", lw=0.6, ls=":")
axB.set_xticks(tsv)
axB.set_xlabel("step  $t$")
axB.set_ylabel(r"$\langle\alpha_x\rangle(t)$")
axB.set_ylim(-1.18, 1.18)
axB.set_title("(b) velocity trembling")
axB.legend(frameon=False, loc="upper right")

fig.tight_layout()
fig.savefig(OUT, dpi=150, bbox_inches="tight")
print(f"\nwrote {OUT}")
