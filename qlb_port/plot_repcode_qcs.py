"""
Analysis of the repetition-code memory experiments on Rigetti QCS (see qcs_quil/repcode_qcs.py).

For every code (bit flip, phase flip), distance d and number of rounds R, the logical error rate is
estimated in three ways from the same shots:
  * unencoded: the final value of data qubit 0 alone, which sees the same gates and waits;
  * majority vote over the final values of the d data qubits (ignores the ancillas);
  * matching: detection events from the ancilla checks (measured at the end, each the parity of its
    check summed over the rounds) and the final data parities, decoded by minimum-weight perfect
    matching (PyMatching) with the detector error model of a
    simple circuit noise model (the weights only need to be roughly right).
With --sim the same analysis is run on shots sampled from that noise model instead of device data.
With --populations it instead compiles the programs with quilc (as for the device) and prints the
fraction of the circuit, sampled at every CZ, that each data qubit spends in |1> in the ideal
evolution (needs pyquil and a quilc server).

Needs stim, pymatching, matplotlib (the project's pyquil venv).

Usage:
    python plot_repcode_qcs.py [--results file] [--sim] [--populations] [--outdir dir]
"""

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "qcs_quil"))
from repcode_qcs import layout, to_stim                   # noqa: E402

RES = "qlb_port/qcs_quil_repcode/qcs_results_repcode.json"
NOISE = {"p1": 0.001, "p2": 0.008, "pm": 0.015, "pz": 0.01, "px": 0.003}


def majority_fail_independent(p, d):
    """Probability that more than d/2 of d independent bits flip, bit i with probability p[i]."""
    dist = np.zeros(d + 1)
    dist[0] = 1.0
    for pi in p:
        dist[1:] = dist[1:] * (1 - pi) + dist[:-1] * pi
        dist[0] *= 1 - pi
    return float(dist[d // 2 + 1:].sum())


def analyse(code, d, R, state, bits):
    import pymatching
    data, anc = layout(d)
    na = len(anc)
    final = bits[:, na:]
    wrong = final != state
    p_i = wrong.mean(axis=0)
    circ = to_stim(code, d, R, state)
    det, obs = circ.compile_m2d_converter().convert(measurements=bits.astype(bool),
                                                     separate_observables=True)
    dem = to_stim(code, d, R, state, NOISE).detector_error_model(decompose_errors=True)
    pred = pymatching.Matching.from_detector_error_model(dem).decode_batch(det)
    return {"unenc": float(wrong[:, 0].mean()), "maj": float((wrong.sum(axis=1) > d / 2).mean()),
            "mwpm": float((pred[:, 0] != obs[:, 0]).mean()), "det": float(det.mean()),
            "pbar": float(p_i.mean()), "indep": majority_fail_independent(p_i, d), "n": len(bits)}


def excited_fraction(quil_text, chain, d):
    """Mean P(|1>) of each data qubit over the CZ steps of the quilc-compiled native program."""
    import math
    from pyquil import Program, get_qc
    from run_qcs import _pin_text
    from cdr_qcs import parse
    qc = get_qc("Cepheus-1-108Q", compiler_timeout=300)
    nat = str(qc.compiler.quil_to_native_quil(
        Program(_pin_text(quil_text, dict(enumerate(chain)))), protoquil=True))
    gates, meas = parse(nat)
    qs = sorted({q for _, _, qq in gates for q in qq} | set(meas))
    idx = {q: k for k, q in enumerate(qs)}
    n = len(qs)
    psi = np.zeros([2] * n, complex)
    psi[(0,) * n] = 1
    data = [chain[i] for i in layout(d)[0]]
    samples = []
    for g, ang, q in gates:
        if g == "CZ":
            sl = [slice(None)] * n
            sl[idx[q[0]]], sl[idx[q[1]]] = 1, 1
            psi[tuple(sl)] *= -1
            p = np.abs(psi) ** 2
            samples.append([p.take(1, axis=idx[dq]).sum() for dq in data])
            continue
        c, s = math.cos(ang / 2), math.sin(ang / 2)
        U = (np.array([[c, -1j * s], [-1j * s, c]]) if g == "RX"
             else np.diag([np.exp(-1j * ang / 2), np.exp(1j * ang / 2)]))
        k = idx[q[0]]
        psi = np.moveaxis(np.tensordot(U, psi, axes=([1], [k])), 0, k)
    return np.mean(samples, axis=0)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default=RES)
    ap.add_argument("--sim", action="store_true", help="use shots sampled from the noise model")
    ap.add_argument("--populations", action="store_true")
    ap.add_argument("--outdir", default="proposal/quairk/ancillaStreaming/figures")
    args = ap.parse_args()

    if args.populations:
        chain = [100, 101, 102, 93, 84, 83, 74, 65, 56, 57, 58, 67, 68, 77, 86, 95, 104]
        base = os.path.dirname(RES)
        print("time-averaged P(1) of the data qubits in the compiled programs (R = 3)")
        for code in ("bit", "phase"):
            for s in (0, 1):
                f = [excited_fraction(open(os.path.join(base, f"rep_{code}_d{d}_r3_s{s}.quil")).read(),
                                      chain, d).mean() for d in (3, 5, 7, 9)]
                print(f"  {code:5s} state {s}: d = 3, 5, 7, 9 -> " + ", ".join(f"{x:.2f}" for x in f))
        return

    R_ = json.load(open(args.results))
    jobs = R_["jobs"]
    print(f"{len(jobs)} jobs, {R_.get('total_execution_us', 0) / 1e6:.1f} s QPU execution"
          + ("  [SIMULATED with the noise model]" if args.sim else ""))

    res = {}
    for j in jobs:
        if j.get("role") != "rep":
            continue
        code, d, R, s = j["code"], j["d"], j["rounds"], j["state"]
        if args.sim:
            bits = to_stim(code, d, R, s, NOISE).compile_sampler(seed=7).sample(j["shots"]).astype(int)
        elif "bitstrings" in j:
            bits = np.array([[int(c) for c in x] for x in j["bitstrings"]])
        else:
            continue
        res[(code, d, R, s)] = analyse(code, d, R, s, bits)

    print("\nper logical state: majority-vote logical error (observed | predicted for independent errors"
          " with the measured per-qubit rates), unencoded error, mean data-qubit error")
    print("code   R  state   d=3 obs|ind       d=5 obs|ind       d=7 obs|ind       d=9 obs|ind"
          "      unencoded (d=3..9)          mean data-qubit error (d=3..9)")
    keys = sorted({(k[0], k[2], k[3]) for k in res})
    for code, R, s in keys:
        ds = sorted(k[1] for k in res if k[0] == code and k[2] == R and k[3] == s)
        v = [res[(code, d, R, s)] for d in ds]
        print(f"{code:5s} {R:2d}   {s}     " + "  ".join(f"{x['maj']:.4f}|{x['indep']:.4f}" for x in v)
              + "    " + " ".join(f"{x['unenc']:.3f}" for x in v)
              + "    " + " ".join(f"{x['pbar']:.3f}" for x in v))

    print("\nmean of the two logical states")
    print("code   d  R   unencoded  majority  matching   detection-event rate")
    table = {}
    for code, d, R in sorted({k[:3] for k in res}):
        v = [res[(code, d, R, s)] for s in (0, 1) if (code, d, R, s) in res]
        m = {k: float(np.mean([x[k] for x in v])) for k in ("unenc", "maj", "mwpm", "det")}
        table[(code, d, R)] = m
        print(f"{code:5s} {d:2d} {R:2d}   {m['unenc']:.4f}     {m['maj']:.4f}    {m['mwpm']:.4f}"
              f"     {m['det']:.3f}")
    print("\nunencoded range over the four codes, mean of the two states")
    for code, R in sorted({(k[0], k[2]) for k in table}):
        u = [table[k]["unenc"] for k in table if k[0] == code and k[2] == R]
        print(f"  {code:5s} R={R:2d}: {min(u):.3f} to {max(u):.3f}   (mean {np.mean(u):.3f})")
    same = all(abs(t["maj"] - t["mwpm"]) < 1e-12 for t in table.values())
    print(f"matching decoder identical to majority vote in every case: {same}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.6), sharex=True, sharey=True)
    labels = {("bit", 0): r"(a) bit-flip code, $|0_L\rangle$", ("bit", 1): r"(b) bit-flip code, $|1_L\rangle$",
              ("phase", 0): r"(c) phase-flip code, $|+_L\rangle$", ("phase", 1): r"(d) phase-flip code, $|-_L\rangle$"}
    floor = 2.5e-4
    for (row, code), (col, s) in [((r, c), (k, s)) for r, c in enumerate(("bit", "phase"))
                                  for k, s in enumerate((0, 1))]:
        ax = axes[row, col]
        Rs = sorted({k[2] for k in res if k[0] == code})
        for R, colr in zip(Rs, ("C0", "C1", "C2", "C3")):
            ds = sorted(k[1] for k in res if k[0] == code and k[2] == R and k[3] == s)
            unenc = np.mean([res[(code, d, R, s)]["unenc"] for d in ds])
            xs = [1] + ds
            ys = [unenc] + [res[(code, d, R, s)]["maj"] for d in ds]
            ax.semilogy(xs, [max(y, floor) for y in ys], "-", color=colr, lw=1.3, label=f"$R={R}$")
            for x, y in zip(xs, ys):
                ax.semilogy(x, max(y, floor), "v" if y == 0 else "o", color=colr, ms=6,
                            mfc="none" if y == 0 else colr)
        ax.axhline(floor, color="0.6", lw=0.6, ls=":")
        ax.set_title(labels[(code, s)], fontsize=10)
        ax.set_xticks([1, 3, 5, 7, 9])
        ax.set_xticklabels(["1\n(unencoded)", "3", "5", "7", "9"])
    for ax in axes[1]:
        ax.set_xlabel("code distance $d$")
    for ax in axes[:, 0]:
        ax.set_ylabel("logical error rate")
    axes[0, 0].legend(frameon=False, fontsize=8.5, loc="upper right")
    fig.tight_layout()
    out = os.path.join(args.outdir, "hw_repcode" + ("_sim" if args.sim else "") + ".png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
