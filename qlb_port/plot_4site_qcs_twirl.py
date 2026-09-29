"""
Pauli twirling (randomized compiling) of the four-site t=4 velocity circuits on Rigetti QCS.

Reads qlb_port/qcs_quil_n2_twirl/qcs_results_twirl_{A,B}.json (two identical rounds, each with the
untwirled native program of the even-parity and the step-1 velocity circuit, 2000 shots, and 16
Pauli-twirled variants of each, 125 shots) and the twirled programs themselves.

Writes <outdir>/hw_4site_qcs_twirl.png and prints the numbers used in the paper.

Usage:
    PYTHONPATH=. python3 -m qlb_port.plot_4site_qcs_twirl
"""

import argparse
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = "qlb_port/qcs_quil_n2_twirl/"
GROUPS = (("even_all_t4", "even-parity state"), ("step1_q1_t4", "original state, $q_1$ read"))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default="proposal/quairk/ancillaStreaming/figures")
    args = ap.parse_args()

    rounds = {r: json.load(open(D + f"qcs_results_twirl_{r}.json")) for r in "AB"}
    for g, _ in GROUPS:
        base = [l for l in open(D + f"{g}_native.quil") if l.split() and l.split()[0][:2] in ("RX", "RZ", "CZ")]
        extra = []
        for k in range(16):
            tw = [l for l in open(D + f"{g}_tw{k:02d}.quil") if l.split() and l.split()[0][:2] in ("RX", "RZ", "CZ")]
            extra.append((sum(l.startswith("RX") for l in tw) - sum(l.startswith("RX") for l in base),
                          sum(l.startswith("RZ") for l in tw) - sum(l.startswith("RZ") for l in base)))
        e = np.array(extra)
        print(f"{g}: {sum(l.startswith('CZ') for l in base)} CZ, {len(base)} native gates; twirling adds "
              f"{e[:, 0].mean():.1f} RX(pi) and {e[:, 1].mean():.1f} RZ(pi) on average")
    for r, R in rounds.items():
        print(f"round {r}: {R['total_execution_us'] / 1e6:.2f} s QPU execution")

    stats = {}
    for r, R in rounds.items():
        J = R["jobs"]
        for g, _ in GROUPS:
            tw = [j for j in J if j["group"] == g and j["variant"] >= 0]
            ctl = [j for j in J if j["group"] == g and j["variant"] < 0][0]
            v = np.array([j["qpu"] for j in tw])
            n = np.array([sum(j["counts"].values()) for j in tw])
            pooled = (v * n).sum() / n.sum()
            stats[(g, r)] = dict(v=v, ctl=ctl["qpu"], pooled=pooled, exact=ctl["exact"],
                                 err=np.sqrt((1 - pooled ** 2) / n.sum()), sd=v.std(ddof=1))
            print(f"  {r} {g}: control {ctl['qpu']:+.3f}  pooled {pooled:+.3f} +- {stats[(g, r)]['err']:.3f} "
                  f"({n.sum()} shots)  variant std {v.std(ddof=1):.3f} (shot noise {np.sqrt(1 / n[0]):.3f})  "
                  f"min {v.min():+.3f} max {v.max():+.3f}")
    allp = [s["pooled"] for s in stats.values()]
    allc = [s["ctl"] for s in stats.values()]
    print(f"pooled twirled: range {min(allp):+.3f}..{max(allp):+.3f}, mean {np.mean(allp):+.3f}; "
          f"controls: range {min(allc):+.3f}..{max(allc):+.3f}")

    fig, ax = plt.subplots(figsize=(8.5, 4.4))
    exact = stats[(GROUPS[0][0], "A")]["exact"]
    xt, xl = [], []
    for k, ((g, lbl), r) in enumerate([(gr, r) for gr in GROUPS for r in "AB"]):
        s = stats[(g, r)]
        jit = np.linspace(-0.18, 0.18, len(s["v"]))
        ax.plot(k + jit, s["v"], "o", color="C4", ms=5, alpha=0.45,
                label="twirled variants (125 shots each)" if k == 0 else None)
        ax.errorbar(k, s["pooled"], yerr=s["err"], fmt="o", color="C4", ms=10, capsize=4,
                    label="twirled, pooled (2000 shots)" if k == 0 else None)
        ax.plot(k, s["ctl"], "D", color="C3", ms=9, label="untwirled (2000 shots)" if k == 0 else None)
        xt.append(k); xl.append(f"{lbl}\nround {r}")
    ax.axhline(exact, color="0.4", lw=1.2, label=f"exact  {exact:+.3f}")
    ax.axhline(0, color="k", lw=0.6, ls=":")
    ax.set_xticks(xt)
    ax.set_xticklabels(xl, fontsize=8.5)
    ax.set_ylabel(r"$\langle\alpha_x\rangle(4)$, raw")
    ax.set_ylim(-0.5, 0.45)
    ax.legend(frameon=False, loc="upper left", fontsize=8.5, ncol=2)
    fig.tight_layout()
    out = f"{args.outdir}/hw_4site_qcs_twirl.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
