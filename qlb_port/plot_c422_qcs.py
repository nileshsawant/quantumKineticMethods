"""
Analysis of the [[4,2,2]] error-detection experiment on Rigetti QCS (see qcs_quil/c422_qcs.py).

For every number of rounds R the Bell-state fidelity F = (1 + C_X + C_Y + C_Z) / 4 is formed from the
correlations C_B of the measured pair in the X, Y and Z bases (sign-corrected by their ideal values):
  * bare: the physical Bell pair on qubits (0,1);
  * encoded, no detection: the logical Bell state, read from pair (2,3) in every shot;
  * encoded, detection: the same, keeping only shots whose four-qubit parity (XXXX, YYYY or ZZZZ)
    has its ideal value; the kept fraction is reported.
Uncertainties are binomial.  With --sim the shots are sampled from a simple Stim noise model.

Usage:
    python plot_c422_qcs.py [--results file] [--sim] [--outdir dir]
"""

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "qcs_quil"))
from c422_qcs import to_stim                               # noqa: E402

RES = "qlb_port/qcs_quil_c422/qcs_results_c422.json"
NOISE = {"p1": 0.001, "p2": 0.008, "pm": 0.015}


def corr(b, i, j, s):
    """Sign-corrected correlation <(-1)^(b_i xor b_j)> and its binomial variance."""
    if len(b) == 0:
        return float("nan"), float("nan")
    c = s * float(np.mean(1 - 2 * (b[:, i] ^ b[:, j])))
    return c, (1 - c * c) / len(b)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default=RES)
    ap.add_argument("--sim", action="store_true")
    ap.add_argument("--outdir", default="proposal/quairk/ancillaStreaming/figures")
    args = ap.parse_args()
    R_ = json.load(open(args.results))
    print(f"{len(R_['jobs'])} jobs, {R_.get('total_execution_us', 0) / 1e6:.1f} s QPU execution"
          + ("  [SIMULATED with the noise model]" if args.sim else ""))

    data = {}
    for j in R_["jobs"]:
        kind, R, B = j["code"], j["rounds"], j["basis"]
        if args.sim:
            b = to_stim(kind, R, B, NOISE).compile_sampler(seed=5).sample(j["shots"]).astype(int)
        elif "bitstrings" in j:
            b = np.array([[int(c) for c in x] for x in j["bitstrings"]])
        else:
            continue
        ideal = to_stim(kind, R, B).compile_sampler(seed=1).sample(1)[0].astype(int)
        data[(kind, R, B)] = (b, ideal)

    rows = {}
    print("\n  R  basis   bare C     enc C (all)   enc C (kept)   kept   [C sign-corrected]")
    for R in sorted({k[1] for k in data}):
        acc = {"bare": [], "raw": [], "post": []}
        kept = []
        for B in ("X", "Y", "Z"):
            if ("enc", R, B) not in data or ("bare", R, B) not in data:
                continue
            bb, ib = data[("bare", R, B)]
            be, ie = data[("enc", R, B)]
            sb = 1 - 2 * (ib[0] ^ ib[1])
            se = 1 - 2 * (ie[2] ^ ie[3])
            ok = (be.sum(axis=1) % 2) == (ie.sum() % 2)
            acc["bare"].append(corr(bb, 0, 1, sb))
            acc["raw"].append(corr(be, 2, 3, se))
            acc["post"].append(corr(be[ok], 2, 3, se))
            kept.append(ok.mean())
            print(f"{R:3d}   {B}      {acc['bare'][-1][0]:+.3f}      {acc['raw'][-1][0]:+.3f}"
                  f"        {acc['post'][-1][0]:+.3f}       {ok.mean():.3f}")
        if len(kept) < 3:
            continue
        F = {k: ((1 + sum(c for c, _ in v)) / 4, np.sqrt(sum(var for _, var in v)) / 4) for k, v in acc.items()}
        rows[R] = (F, float(np.mean(kept)))

    print("\n  R   Bell fidelity:  bare            encoded, no detection   encoded, detection    kept")
    for R, (F, k) in sorted(rows.items()):
        print(f"{R:3d}                  {F['bare'][0]:.3f} ({F['bare'][1]:.3f})     {F['raw'][0]:.3f} ({F['raw'][1]:.3f})"
              f"           {F['post'][0]:.3f} ({F['post'][1]:.3f})        {k:.3f}")
    if not rows:
        return
    print("\n  R   infidelity ratio (1 - F_detection) / (1 - F_bare)")
    for R, (F, k) in sorted(rows.items()):
        print(f"{R:3d}   {(1 - F['post'][0]) / (1 - F['bare'][0]):.3f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    Rs = sorted(rows)
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.0), gridspec_kw={"width_ratios": [1.4, 1]})
    for key, lab, st in (("bare", "physical Bell pair", "o-"), ("raw", "logical, no detection", "s--"),
                         ("post", "logical, with detection", "D-")):
        axA.errorbar(Rs, [rows[R][0][key][0] for R in Rs], yerr=[rows[R][0][key][1] for R in Rs],
                     fmt=st, capsize=3, ms=6, label=lab)
    axA.set_xlabel("rounds $R$ of CZ$\\cdot$CZ on each pair")
    axA.set_ylabel("Bell-state fidelity")
    axA.set_title("(a) fidelity")
    axA.legend(frameon=False, fontsize=9)
    axB.plot(Rs, [rows[R][1] for R in Rs], "D-", color="C2")
    axB.set_xlabel("rounds $R$")
    axB.set_ylabel("fraction of shots kept")
    axB.set_title("(b) shots passing the check")
    fig.tight_layout()
    out = os.path.join(args.outdir, "hw_c422" + ("_sim" if args.sim else "") + ".png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
