"""
Noisy simulation (Stim + PyMatching) of the repetition-code memory with mid-circuit measurement of
the ancillas, compared with the end-only scheme that was run on Cepheus-1 (repcode_qcs.py).

Schemes, all with the same gates and the circuit noise model NOISE of plot_repcode_qcs.py:
  end:    ancillas measured once at the end (summed parities), as on the device;
  reset:  ancillas measured every round and reset by a conditional X (feed-forward);
  noreset: ancillas measured every round, no reset (running parity, no control flow in a shot).
km scales the data-qubit idle errors during each mid-circuit readout (relative to one round's idle
errors pz, px); readout on superconducting devices is slow compared with gates, so km = 1-3.

Logical error rate = fraction of shots where the MWPM-corrected value of data qubit 0 is wrong,
averaged over the two logical states.  The unencoded reference is data qubit 0 without decoding.

Usage (on a compute node):
    python sim_repcode_mid.py [--shots 100000] [--km 0,1,3] [--outdir dir]
"""

import argparse
import os
import sys
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "qcs_quil"))
import repcode_qcs as endq                                  # noqa: E402
import repcode_mid_qcs as midq                              # noqa: E402
from plot_repcode_qcs import NOISE                          # noqa: E402

SCHEMES = ("end", "reset", "noreset")


def circuit(scheme, code, d, R, s, noise):
    if scheme == "end":
        return endq.to_stim(code, d, R, s, noise)
    return midq.to_stim(code, d, R, s, scheme == "reset", noise)


def run(task):
    import pymatching
    scheme, code, d, R, s, km, shots = task
    noise = dict(NOISE, km=km)
    ideal = circuit(scheme, code, d, R, s, None)
    det0, obs0 = ideal.compile_detector_sampler(seed=1).sample(64, separate_observables=True)
    assert not det0.any() and not obs0.any(), (scheme, code, d, R, s)
    c = circuit(scheme, code, d, R, s, noise)
    meas = c.compile_sampler(seed=11).sample(shots)
    det, obs = c.compile_m2d_converter().convert(measurements=meas, separate_observables=True)
    dem = c.detector_error_model(decompose_errors=True)
    pred = pymatching.Matching.from_detector_error_model(dem).decode_batch(det)
    i0 = R * (d - 1) if scheme != "end" else d - 1          # record index of data qubit 0
    return task, float((pred[:, 0] != obs[:, 0]).mean()), float((meas[:, i0] != s).mean())


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--shots", type=int, default=100000)
    ap.add_argument("--distances", default="3,5,7,9")
    ap.add_argument("--rounds", default="1,3,6,10")
    ap.add_argument("--km", default="0,1,3")
    ap.add_argument("--procs", type=int, default=os.cpu_count())
    ap.add_argument("--outdir", default="proposal/quairk/ancillaStreaming/figures")
    args = ap.parse_args()
    ds = [int(x) for x in args.distances.split(",")]
    Rs = [int(x) for x in args.rounds.split(",")]
    kms = [float(x) for x in args.km.split(",")]

    tasks = [(sc, code, d, R, s, km, args.shots) for sc in SCHEMES for code in endq.CODES
             for d in ds for R in Rs for s in (0, 1) for km in (kms if sc != "end" else [0.0])]
    with Pool(args.procs) as pool:
        out = pool.map(run, tasks)
    res = {}
    for (sc, code, d, R, s, km, _), pl, pu in out:
        res.setdefault((sc, code, d, R, km), []).append((pl, pu))
    tab = {k: np.mean(v, axis=0) for k, v in res.items()}

    print(f"noise {NOISE}, {args.shots} shots per state; logical error (MWPM), mean of 2 states")
    for code in endq.CODES:
        print(f"\n{code}-flip code          " + "".join(f"   d={d:<5d}" for d in ds) + "  unencoded")
        for R in Rs:
            rows = [("end", 0.0)] + [(sc, km) for sc in ("reset", "noreset") for km in kms]
            for sc, km in rows:
                v = [tab[(sc, code, d, R, km)] for d in ds]
                lab = f"R={R:2d} {sc:7s}" + (f" km={km:g}" if sc != "end" else "     ")
                print(f"  {lab:20s}" + "".join(f"  {x[0]:.5f}" for x in v) + f"   {v[0][1]:.4f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    km_fig = 1.0 if 1.0 in kms else kms[0]
    fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True, sharey=True)
    floor = 0.5 / args.shots
    titles = {"end": "ancillas read at the end", "noreset": "ancillas read every round"}
    labels = iter("abcd")
    for row, code in enumerate(endq.CODES):
        for col, sc in enumerate(("end", "noreset")):        # reset gives the same rates as noreset
            ax = axes[row, col]
            km = 0.0 if sc == "end" else km_fig
            for d, colr in zip(ds, ("C0", "C1", "C2", "C3", "C4")):
                y = [max(tab[(sc, code, d, R, km)][0], floor) for R in Rs]
                ax.semilogy(Rs, y, "o-", color=colr, label=f"$d={d}$")
            u = [tab[(sc, code, ds[0], R, km)][1] for R in Rs]
            ax.semilogy(Rs, u, "k--", lw=1, label="unencoded")
            ax.set_title(f"({next(labels)}) {code}-flip code, {titles[sc]}", fontsize=9.5)
    for ax in axes[1]:
        ax.set_xlabel("rounds $R$")
    for ax in axes[:, 0]:
        ax.set_ylabel("logical error rate")
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    os.makedirs(args.outdir, exist_ok=True)
    path = os.path.join(args.outdir, "sim_repcode_mid.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"\nwrote {path} (mid-circuit panels with km={km_fig:g})")


if __name__ == "__main__":
    main()
