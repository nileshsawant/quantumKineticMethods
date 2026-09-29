"""
Dynamical-decoupling test of the four-site t=4 velocity circuits on Rigetti QCS.

Reads qlb_port/qcs_quil_n2_dd/qcs_results_dd.json: two circuits (original state with q1 read, and the
even-parity state with all qubits read), each in three versions -- as compiled (orig), with every CZ
in its own FENCE-bounded layer (fence), and fenced with paired RX(pi) echo pulses on idle qubits (dd)
-- interleaved over four rounds of 2000 shots.

Writes <outdir>/hw_4site_qcs_dd.png and prints the numbers used in the paper.

Usage:
    PYTHONPATH=. python3 -m qlb_port.plot_4site_qcs_dd
"""

import argparse
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = "qlb_port/qcs_quil_n2_dd/qcs_results_dd.json"
CIRCUITS = (("step1_q1_t4", "original state, $q_1$ read"), ("even_all_t4", "even-parity state"))
VERSIONS = (("orig", "as compiled"), ("fence", "fenced layers"), ("dd", "fenced + echoes"))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default="proposal/quairk/ancillaStreaming/figures")
    args = ap.parse_args()

    R = json.load(open(RES))
    J = R["jobs"]
    exact = J[0]["exact"]
    print(f"DD batch: {len(J)} jobs, {R['total_execution_us'] / 1e6:.1f} s QPU execution; exact {exact:+.3f}")
    stats = {}
    for c, _ in CIRCUITS:
        for v, _ in VERSIONS:
            vals = np.array([j["qpu"] for j in sorted(J, key=lambda j: j["variant"]) if j["group"] == f"{c}_{v}"])
            stats[(c, v)] = vals
            print(f"  {c:12s} {v:6s} {np.round(vals, 3).tolist()}  mean {vals.mean():+.3f}  sd {vals.std(ddof=1):.3f}"
                  f"  |mean - exact| {abs(vals.mean() - exact):.3f}")
    for v, _ in VERSIONS:
        err = np.mean([abs(stats[(c, v)].mean() - exact) for c, _ in CIRCUITS])
        print(f"  mean |error| over both circuits, {v}: {err:.3f}")

    fig, ax = plt.subplots(figsize=(8.5, 4.4))
    cols = {"orig": "C3", "fence": "C0", "dd": "C2"}
    xt, xl = [], []
    for i, (c, clbl) in enumerate(CIRCUITS):
        for k, (v, vlbl) in enumerate(VERSIONS):
            x = i * 4 + k
            vals = stats[(c, v)]
            ax.plot(x + np.linspace(-0.15, 0.15, len(vals)), vals, "o", color=cols[v], ms=7,
                    label=vlbl if i == 0 else None)
            ax.plot([x - 0.3, x + 0.3], [vals.mean()] * 2, "-", color=cols[v], lw=2)
            xt.append(x); xl.append(f"{vlbl}")
        ax.text(i * 4 + 1, 0.3, clbl, ha="center", fontsize=9)
    ax.axhline(exact, color="0.4", lw=1.2, label=f"exact  {exact:+.3f}")
    ax.axhline(0, color="k", lw=0.6, ls=":")
    ax.set_xticks(xt)
    ax.set_xticklabels(xl, fontsize=8, rotation=20)
    ax.set_ylabel(r"$\langle\alpha_x\rangle(4)$, raw")
    ax.set_ylim(-0.45, 0.38)
    ax.legend(frameon=False, loc="lower right", fontsize=8.5, ncol=2)
    fig.tight_layout()
    out = f"{args.outdir}/hw_4site_qcs_dd.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
