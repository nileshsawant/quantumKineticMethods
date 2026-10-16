"""
SI1000 noise model (Gidney et al., Quantum 5, 605 (2021)) for the repetition codes run on Cepheus-1,
fitted to the device results of plot_repmid_qcs.py (d = 3, R = 1, 2, 3, qubits 104, 95, 86, 77, 68).

SI1000 with strength p: two-qubit gate DEPOLARIZE2(p) (each CNOT counts as one CZ), single-qubit
gate DEPOLARIZE1(p/10), idle qubit in a gate layer DEPOLARIZE1(p/10), initialization X_ERROR(2p),
measurement result flipped with 5p, and every qubit that waits while others are measured
DEPOLARIZE1(2p).  Here the last factor is a parameter f (SI1000: f = 2), and the readout flip can
be fixed to a measured value instead of 5p.  The circuits are the same op lists and detectors as
repcode_qcs.py (ancillas read at the end) and repcode_mid_qcs.py (read every round, no reset).

Steps:
  1. fit p to the end-only device results (unencoded error, majority-vote error, detection rate);
  2. at that p, predict the every-round results and compare with the device;
  3. fit f, the error of a waiting qubit per neighboring readout, to the every-round results;
  4. break-even: the p below which reading every round beats reading at the end (d = 3, 5, 7).
Steps 1-3 are done twice: with the SI1000 readout flip 5p, and with the measured readout errors.

Usage (compute node):
    python sim_si1000_repcode.py [--shots 20000] [--procs 32]
"""

import argparse
import json
import os
import sys
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "qcs_quil"))
import repcode_mid_qcs as midq  # noqa: E402
import repcode_qcs as endq  # noqa: E402

RES = os.path.join(HERE, "qcs_repmid_results")
CODES, ROUNDS = ("bit", "phase"), (1, 2, 3)


def ops_for(scheme, code, d, R, s):
    return endq.ops(code, d, R, s) if scheme == "end" else midq.ops_mid(code, d, R, s, False)


def layers(op_list):
    """Group ops into time steps: gates on disjoint qubits share a layer; consecutive measurements
    form one measurement layer."""
    out, cur, busy, kind = [], [], set(), None
    for o in op_list:
        if o[0] in ("IDLE", "IDLE_M", "ROUND_IDLE"):
            continue
        k = "meas" if o[0] in ("M", "MM") else "gate"
        qs = set(o[1:3]) if o[0] == "CX" else {o[1]}
        if cur and (k != kind or (k == "gate" and qs & busy)):
            out.append((kind, cur))
            cur, busy = [], set()
        cur.append(o)
        busy |= qs
        kind = k
    if cur:
        out.append((kind, cur))
    return out


def si1000(scheme, code, d, R, s, p, f=2.0, pm=None):
    import stim
    nq = 2 * d - 1
    c = stim.Circuit()
    c.append("X_ERROR", range(nq), 2 * p)
    lays = layers(ops_for(scheme, code, d, R, s))
    for n, (kind, lay) in enumerate(lays):
        used = set()
        if kind == "gate":
            for o in lay:
                if o[0] == "CX":
                    c.append("CX", [o[1], o[2]])
                    c.append("DEPOLARIZE2", [o[1], o[2]], p)
                    used |= {o[1], o[2]}
                else:
                    c.append(o[0], [o[1]])
                    c.append("DEPOLARIZE1", [o[1]], p / 10)
                    used.add(o[1])
            idle = [q for q in range(nq) if q not in used]
            if idle:
                c.append("DEPOLARIZE1", idle, p / 10)
        else:
            mq = [o[1] for o in lay]
            idle = [q for q in range(nq) if q not in mq]
            c.append("M", mq, 5 * p if pm is None else pm)
            if idle and n < len(lays) - 1:
                c.append("DEPOLARIZE1", idle, f * p)
    ideal = endq.to_stim(code, d, R, s) if scheme == "end" else midq.to_stim(code, d, R, s, False)
    for inst in ideal:
        if inst.name in ("DETECTOR", "OBSERVABLE_INCLUDE"):
            c.append(inst)
    return c, ideal


def metrics(bits, ideal, noisy, d, R_na):
    import pymatching
    det, obs = ideal.compile_m2d_converter().convert(measurements=bits.astype(bool), separate_observables=True)
    pred = pymatching.Matching.from_detector_error_model(
        noisy.detector_error_model(decompose_errors=True, approximate_disjoint_errors=True)).decode_batch(det)
    ref = ideal.compile_sampler().sample(1)[0].astype(int)
    dq = bits[:, R_na:]
    return np.array([np.mean(dq[:, 0] != ref[R_na]),
                     np.mean((dq.sum(1) > d / 2).astype(int) != int(ref[R_na:].sum() > d / 2)),
                     np.mean(pred[:, 0] != obs[:, 0]), det.mean()])


def rna(scheme, d, R):
    return d - 1 if scheme == "end" else R * (d - 1)


def device():
    out = {}
    for scheme, f in (("end", "repcode_end.json"), ("mid", "repcode_mid.json")):
        rows = {}
        for j in json.load(open(os.path.join(RES, f)))["jobs"]:
            bits = np.array([[int(c) for c in s] for s in j["bitstrings"]])
            code, d, R, s = j["code"], j["d"], j["rounds"], j["state"]
            ideal = (endq.to_stim(code, d, R, s) if scheme == "end" else midq.to_stim(code, d, R, s, False))
            noisy, _ = si1000(scheme, code, d, R, s, 0.005)
            rows[(code, R, s)] = metrics(bits, ideal, noisy, d, rna(scheme, d, R))
        for code in CODES:
            for R in ROUNDS:
                out[(scheme, code, R)] = np.mean([rows[(code, R, s)] for s in (0, 1)], axis=0)
    return out


def sim_task(t):
    scheme, code, d, R, s, p, f, pm, shots = t
    noisy, ideal = si1000(scheme, code, d, R, s, p, f, pm)
    bits = noisy.compile_sampler(seed=3).sample(shots).astype(int)
    return t, metrics(bits, ideal, noisy, d, rna(scheme, d, R))


def predict(pool, scheme, ps, shots, f=2.0, pm=None, d=3, rounds=ROUNDS):
    tasks = [(scheme, code, d, R, s, p, f, pm, shots) for p in ps for code in CODES for R in rounds
             for s in (0, 1)]
    acc = {}
    for t, m in pool.map(sim_task, tasks):
        acc.setdefault((t[5], t[0], t[1], t[3]), []).append(m)
    return {k: np.mean(v, axis=0) for k, v in acc.items()}


def loss(pred, dev, p, scheme, cols=(0, 1, 3)):
    e = 0.0
    for code in CODES:
        for R in ROUNDS:
            a, b = pred[(p, scheme, code, R)], dev[(scheme, code, R)]
            e += sum((np.log(a[i] + 2e-3) - np.log(b[i] + 2e-3)) ** 2 for i in cols)
    return e


def table(title, dev, pred, p, scheme):
    print(f"   {title}")
    print("   code   R   unencoded dev/model   majority dev/model   detection dev/model")
    for code in CODES:
        for R in ROUNDS:
            a, b = dev[(scheme, code, R)], pred[(p, scheme, code, R)]
            print(f"   {code:5s} {R}    {a[0]:.3f} / {b[0]:.3f}        {a[1]:.3f} / {b[1]:.3f}"
                  f"        {a[3]:.3f} / {b[3]:.3f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--shots", type=int, default=20000)
    ap.add_argument("--procs", type=int, default=os.cpu_count())
    args = ap.parse_args()
    dev = device()
    ps = list(np.round(np.geomspace(1e-3, 3e-2, 22), 5))
    fs = [2, 5, 10, 20, 40, 80]
    # measured readout errors: standard readout at the end, readout without 1->2 pulse in the middle
    pm_meas = {"end": 1 - 0.967, "mid": 1 - 0.938}
    with Pool(args.procs) as pool:
        for label, pmode in (("SI1000 (readout flip 5p)", None), ("SI1000 with measured readout errors", pm_meas)):
            print(f"\n=== {label}")
            pm_end = None if pmode is None else pmode["end"]
            pm_mid = None if pmode is None else pmode["mid"]
            pe = predict(pool, "end", ps, args.shots, pm=pm_end)
            p_fit = min(ps, key=lambda p: loss(pe, dev, p, "end"))
            print(f"1. p fitted to the end-only results: p = {p_fit:.4f}")
            table("ancillas read at the end", dev, pe, p_fit, "end")
            pmid = predict(pool, "mid", [p_fit], args.shots, pm=pm_mid)
            print(f"2. ancillas read every round, same p, waiting-qubit error f*p with f = 2:")
            table("ancillas read every round", dev, pmid, p_fit, "mid")
            fits = {}
            for f in [f for f in fs if f * p_fit <= 0.5]:
                pf = predict(pool, "mid", [p_fit], args.shots, f=f, pm=pm_mid)
                fits[f] = (loss(pf, dev, p_fit, "mid"), pf)
            f_fit = min(fits, key=lambda f: fits[f][0])
            print(f"3. waiting-qubit error per neighboring readout fitted to the every-round results: "
                  f"f = {f_fit} (depolarizing {f_fit * p_fit:.3f} per readout; SI1000: {2 * p_fit:.4f})")
            table(f"ancillas read every round, f = {f_fit}", dev, fits[f_fit][1], p_fit, "mid")

        print("\n=== 4. break-even under SI1000 (f = 2, readout flip 5p): logical error (matching), mean of states")
        pb = list(np.round(np.geomspace(5e-4, 2e-2, 12), 5))
        for d, rounds in ((3, (3, 10)), (5, (3, 10)), (7, (10,))):
            for R in rounds:
                e = predict(pool, "end", pb, args.shots * 5, d=d, rounds=(R,))
                m = predict(pool, "mid", pb, args.shots * 5, d=d, rounds=(R,))
                for code in CODES:
                    row = [(p, e[(p, "end", code, R)][2], m[(p, "mid", code, R)][2]) for p in pb]
                    better = [p for p, a, b in row if b < a]
                    print(f"   d={d} R={R:2d} {code:5s}: every round better for p <= "
                          f"{max(better) if better else 'none':}   "
                          + "  ".join(f"p={p:.4f}: {a:.4f}/{b:.4f}" for p, a, b in row[::3]))


if __name__ == "__main__":
    main()
