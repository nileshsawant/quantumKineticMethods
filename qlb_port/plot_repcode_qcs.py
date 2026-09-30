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

Needs stim, pymatching, matplotlib (the project's pyquil venv).

Usage:
    python plot_repcode_qcs.py [--results file] [--sim] [--outdir dir]
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


def analyse(code, d, R, state, bits):
    import pymatching
    data, anc = layout(d)
    na = len(anc)
    final = bits[:, na:]
    unenc = float((final[:, 0] != state).mean())
    maj = float(((final != state).sum(axis=1) > d / 2).mean())
    circ = to_stim(code, d, R, state)
    det, obs = circ.compile_m2d_converter().convert(measurements=bits.astype(bool),
                                                     separate_observables=True)
    dem = to_stim(code, d, R, state, NOISE).detector_error_model(decompose_errors=True)
    pred = pymatching.Matching.from_detector_error_model(dem).decode_batch(det)
    mwpm = float((pred[:, 0] != obs[:, 0]).mean())
    return unenc, maj, mwpm, float(det.mean())


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default=RES)
    ap.add_argument("--sim", action="store_true", help="use shots sampled from the noise model")
    ap.add_argument("--outdir", default="proposal/quairk/ancillaStreaming/figures")
    args = ap.parse_args()

    R_ = json.load(open(args.results))
    jobs = R_["jobs"]
    print(f"{len(jobs)} jobs, {R_.get('total_execution_us', 0) / 1e6:.1f} s QPU execution"
          + ("  [SIMULATED with the noise model]" if args.sim else ""))
    for j in jobs:
        if j.get("role") == "mcm" and "bitstrings" in j:
            b = np.array([[int(c) for c in s] for s in j["bitstrings"]])
            print(f"  {j['name']}: P(1) per measurement {np.round(b.mean(axis=0), 3).tolist()}"
                  f" (ideal {j['expected']}); all as ideal {float((b == j['expected']).all(axis=1).mean()):.3f}")

    rows = {}
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
        rows.setdefault((code, d, R), []).append((analyse(code, d, R, s, bits), bits.shape[0]))

    print("\ncode   d  R   unencoded  majority  matching   detection-event rate   (mean of the two states)")
    table = {}
    for key in sorted(rows):
        v = np.mean([r[0] for r in rows[key]], axis=0)
        n = sum(r[1] for r in rows[key])
        table[key] = v
        se = np.sqrt(v[2] * (1 - v[2]) / n)
        print(f"{key[0]:5s} {key[1]:2d} {key[2]:2d}   {v[0]:.4f}     {v[1]:.4f}    {v[2]:.4f}"
              f" ({se:.4f})   {v[3]:.3f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2), sharey=True)
    for ax, code in zip(axes, ("bit", "phase")):
        Rs = sorted({k[2] for k in table if k[0] == code})
        for R, col in zip(Rs, ("C0", "C1", "C2", "C3", "C4")):
            ds = sorted({k[1] for k in table if k[0] == code and k[2] == R})
            ax.semilogy(ds, [max(table[(code, d, R)][2], 1e-4) for d in ds], "o-", color=col,
                        label=f"$R={R}$, matching")
            ax.semilogy(ds, [max(table[(code, d, R)][1], 1e-4) for d in ds], "s--", color=col,
                        mfc="none", alpha=0.7, label=f"$R={R}$, majority")
        ax.set_xlabel("code distance $d$")
        ax.set_title(f"({'a' if code == 'bit' else 'b'}) {code}-flip code")
        ax.set_xticks(sorted({k[1] for k in table}))
    axes[0].set_ylabel("logical error rate")
    axes[1].legend(frameon=False, fontsize=7.5, ncol=2)
    fig.tight_layout()
    out = os.path.join(args.outdir, "hw_repcode" + ("_sim" if args.sim else "") + ".png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
