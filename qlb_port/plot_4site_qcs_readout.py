"""
Readout-error mitigation of the four-site (n_pos=2) Fourier-streaming run on Rigetti QCS.

Reads (no cloud, no credits):
    qlb_port/qcs_quil_n2_fourier/qcs_results_n2f_pinned.json   QCS run: 10 QLB + 16 calibration jobs
    qlb_port/hw_4site_fourier_{trembling,density}.npz           Open Quantum run, for comparison

The 16 calibration jobs prepare each computational basis state |b> and measure all four
qubits, giving the confusion matrix A[i, b] = P(read i | prepared b).  A measured
distribution p is mitigated by the non-negative least-squares solution of A q = p,
renormalized.  The velocity jobs read only q1, so they are corrected with the 2x2
marginal of A for that qubit.

Writes <outdir>/hw_4site_qcs_readout.png and prints the numbers used in the paper.

Usage:
    PYTHONPATH=. python3 -m qlb_port.plot_4site_qcs_readout
"""

import argparse
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import nnls

RESULTS = "qlb_port/qcs_quil_n2_fourier/qcs_results_n2f_pinned.json"
QUBITS = ["$q_0$", "$q_1$", "$p_0$", "$p_1$"]


def tvd(p, q):
    return 0.5 * float(np.abs(np.asarray(p) - np.asarray(q)).sum())


def confusion(jobs, shots, n=4):
    A = np.zeros((2 ** n, 2 ** n))
    for b in range(2 ** n):
        for i, v in jobs[f"cal_{b:0{n}b}"]["counts"].items():
            A[int(i), b] = v / shots
    return A


def marginal(A, k):
    """2x2 confusion of qubit k, averaged over the prepared states of the others."""
    n = A.shape[0]
    M = np.zeros((2, 2))
    for a in (0, 1):
        cols = [b for b in range(n) if (b >> k) & 1 == a]
        for r in (0, 1):
            M[r, a] = np.mean([A[[i for i in range(n) if (i >> k) & 1 == r], b].sum() for b in cols])
    return M


def mitigate(A, p):
    q, _ = nnls(A, p)
    return q / q.sum()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default="proposal/quairk/ancillaStreaming/figures")
    args = ap.parse_args()

    R = json.load(open(RESULTS))
    jobs = {j["name"]: j for j in R["jobs"]}
    S = R["shots"]
    DT = np.load("qlb_port/hw_4site_fourier_trembling.npz")
    DD = np.load("qlb_port/hw_4site_fourier_density.npz")
    ts = np.arange(5)
    xs = np.arange(4)

    A = confusion(jobs, S)
    print(f"QCS run: qubits {R['qubits']}, {S} shots, "
          f"{R['total_execution_us'] / 1e6:.1f} s QPU execution")
    print("P(read b | prepared b):", np.round(np.diag(A), 3).tolist())
    for k, nm in enumerate(["q0", "q1", "p0", "p1"]):
        M = marginal(A, k)
        print(f"  {nm}: P(1|0)={M[1, 0]:.3f}  P(0|1)={M[0, 1]:.3f}")
    for p1 in (0, 1):
        cols = [b for b in range(16) if not (b >> 2) & 1 and (b >> 3) & 1 == p1]
        e = np.mean([A[[i for i in range(16) if (i >> 2) & 1], b].sum() for b in cols])
        print(f"  P(p0 reads 1 | p0=0, p1={p1}) = {e:.3f}")

    M1 = marginal(A, 1)
    av_ex, av_raw, av_mit = [], [], []
    for t in ts:
        j = jobs[f"alpha_t{t}"]
        p = np.array([j["counts"].get("0", 0), j["counts"].get("1", 0)]) / S
        q = np.clip(np.linalg.solve(M1, p), 0, None)
        q /= q.sum()
        av_ex.append(j["exact"])
        av_raw.append(p[1] - p[0])
        av_mit.append(q[1] - q[0])
    av_ex, av_raw, av_mit = map(np.array, (av_ex, av_raw, av_mit))

    rho_ex, rho_raw, rho_mit = [], [], []
    for t in ts:
        j = jobs[f"dens_t{t}"]
        p = np.zeros(16)
        for i, v in j["counts"].items():
            p[int(i)] = v / S
        rho_ex.append(j["exact_rho"])
        rho_raw.append(p.reshape(4, 4).sum(1))
        rho_mit.append(mitigate(A, p).reshape(4, 4).sum(1))
    rho_ex, rho_raw, rho_mit = map(np.array, (rho_ex, rho_raw, rho_mit))

    fit = lambda y: np.linalg.lstsq(np.column_stack([np.ones(5), av_ex]), y, rcond=None)[0]
    print("\nvelocity:  t  exact  OpenQuantum  QCS-raw  QCS-mitigated")
    for t in ts:
        print(f"  {t}  {av_ex[t]:+.3f}  {DT['av_hw'][t]:+.3f}  {av_raw[t]:+.3f}  {av_mit[t]:+.3f}")
    for nm, y in (("OpenQuantum", DT["av_hw"]), ("QCS raw", av_raw), ("QCS mitigated", av_mit)):
        a, b = fit(y)
        print(f"  fit {nm:14s} a={a:+.3f} b={b:.3f}")
    print("\ndensity:  t  TVD OpenQuantum / QCS-raw / QCS-mitigated   <x> raw/mitigated   "
          "rho raw | rho mitigated")
    for t in ts:
        print(f"  {t}  {tvd(rho_ex[t], DD['rho_hw'][t]):.3f} / {tvd(rho_ex[t], rho_raw[t]):.3f} / "
              f"{tvd(rho_ex[t], rho_mit[t]):.3f}   {rho_raw[t] @ xs:.3f}/{rho_mit[t] @ xs:.3f}   "
              f"{np.round(rho_raw[t], 3).tolist()} | {np.round(rho_mit[t], 3).tolist()}")

    # ---- figure: (a) velocity, (b) density TVD, (c) confusion matrix ----
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(16, 4.4),
                                        gridspec_kw={"width_ratios": [1.1, 1.0, 0.95]})
    E = float(DT["E"])
    tsm = np.linspace(0, 4, 400)
    axA.plot(tsm, np.sin(2 * E * tsm), color="0.6", lw=1.6, label=r"exact  $\sin 2Et$")
    axA.plot(ts, av_ex, "o", color="0.4", ms=5)
    axA.plot(ts, DT["av_hw"], "D", color="C3", ms=7, label="Open Quantum (qubits 0,1,11,10)")
    axA.plot(ts, av_raw, "^", color="C2", ms=9, mfc="none", mew=1.7,
             label="QCS raw (qubits 103,102,100,101)")
    axA.plot(ts, av_mit, "^", color="C2", ms=8, label="QCS readout-mitigated")
    axA.axhline(0, color="k", lw=0.6, ls=":")
    axA.set_xticks(ts)
    axA.set_xlabel("step  $t$")
    axA.set_ylabel(r"$\langle\alpha_x\rangle(t)$")
    axA.set_ylim(-1.18, 1.18)
    axA.set_title("(a) velocity trembling")
    axA.legend(frameon=False, loc="lower left", fontsize=8.5)

    w = 0.27
    axB.bar(ts - w, [tvd(rho_ex[t], DD["rho_hw"][t]) for t in ts], w, color="C3",
            label="Open Quantum")
    axB.bar(ts, [tvd(rho_ex[t], rho_raw[t]) for t in ts], w, color="C2", alpha=0.45,
            label="QCS raw")
    axB.bar(ts + w, [tvd(rho_ex[t], rho_mit[t]) for t in ts], w, color="C2",
            label="QCS readout-mitigated")
    axB.axhline(0.32, color="0.4", lw=0.8, ls="--")
    axB.text(4.45, 0.325, "uniform at $t=0,1,4$", ha="right", va="bottom", fontsize=8, color="0.3")
    axB.set_xticks(ts)
    axB.set_xlabel("step  $t$")
    axB.set_ylabel("TVD from exact density")
    axB.set_ylim(0, 0.36)
    axB.set_title("(b) packet density error")
    axB.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.86), fontsize=8.5)

    im = axC.imshow(np.log10(np.maximum(A, 1e-4)), cmap="viridis", vmin=-3.3, vmax=0,
                    origin="upper")
    lab = [f"{b:04b}"[::-1] for b in range(16)]
    axC.set_xticks(range(16))
    axC.set_xticklabels(lab, rotation=90, fontsize=6.5)
    axC.set_yticks(range(16))
    axC.set_yticklabels(lab, fontsize=6.5)
    axC.set_xlabel(r"prepared  $q_0q_1p_0p_1$")
    axC.set_ylabel(r"read  $q_0q_1p_0p_1$")
    axC.set_title("(c) readout confusion matrix")
    cb = fig.colorbar(im, ax=axC, fraction=0.046, pad=0.03)
    cb.set_label(r"$\log_{10} P(\mathrm{read}\,|\,\mathrm{prepared})$")
    fig.tight_layout()
    out = f"{args.outdir}/hw_4site_qcs_readout.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
