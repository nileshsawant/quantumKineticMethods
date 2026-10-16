"""
Analysis of the mid-circuit measurement tests on Rigetti Cepheus-1 (results in qcs_repmid_results/).

Prints the numbers quoted in the paper:
  1. a gate right after a mid-circuit measurement, by readout calibration and prepared state
     (qcs_quil/mcm_char2.py on 20 qubits: device MEASURE, its 1->2 pulse repeated right after
     the readout or after the pause, MEASURE redefined without the 1->2 pulse);
  2. the pause after a measurement (qcs_quil/mcm_probe2.py, qubit 84; qcs_quil/mcm_char.py, 20 qubits);
  3. the readout length, with the 0/1 threshold refitted from the raw IQ points (mcm_char2.py,
     10 qubits, MEASURE without the 1->2 pulse);
  4. the repetition codes (d = 3, R = 1, 2, 3, qubits 104, 95, 86, 77, 68) measured at the end
     (qcs_quil/repcode_qcs.py, run_qcs.py) and every round, without and with an echo on the data
     qubits during the readouts (qcs_quil/repcode_mid_qcs.py, run_repmid_qcs.py), same session.

Needs stim and pymatching (the project's pyquil venv).

Usage:
    python plot_repmid_qcs.py [--results qcs_repmid_results]
"""

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "qcs_quil"))


def load(res, name):
    return json.load(open(os.path.join(res, name)))


def gate_after(res):
    print("1. gate after a mid-circuit measurement (20 qubits, 100 shots, x alternates 0, 1)")
    print("   readout                         P(read = x)   gate worked after 0   after 1")
    names = {"none": "device MEASURE", "early": "1->2 pulse repeated at once",
             "late": "1->2 pulse repeated after 1 us", "noshelve": "MEASURE without 1->2 pulse"}
    for v, label in names.items():
        d = load(res, f"mcm_char2_{v}.json")
        K = len(d["streams"][f"q{d['qubits'][0]}_classified"]) // 2
        x = np.arange(K) % 2
        a, g0, g1, nos = [], [], [], []
        for q in d["qubits"]:
            b = np.array(d["streams"][f"q{q}_classified"]).reshape(K, 2)
            a.append(np.mean(b[:, 0] == x))
            g0.append(np.mean(b[x == 0, 1] == 1))
            g1.append(np.mean(b[x == 1, 1] == 0))
            if q not in d["f12"]:
                nos.append((q, g1[-1]))
        print(f"   {label:32s} {np.mean(a):.3f}         {np.mean(g0):.3f}                {np.mean(g1):.3f}"
              f"   (range after 1: {min(g1):.2f} to {max(g1):.2f})")
        if v == "none":
            print("      qubits without a 1->2 pulse in their MEASURE, after 1: "
                  + ", ".join(f"q{q} {g:.2f}" for q, g in nos))


def pause(res):
    print("\n2. pause after a measurement")
    d = load(res, "mcm_probe2.json")
    for k in ("A q84 flips", "B q84 flips+delay"):
        b = np.array(d[k]["bits"])
        print(f"   {k:20s} (50 x [RX(pi); MEASURE]): readout changed after {np.mean(b[1:] != b[:-1]):.2f} of the flips")
    print("   20 qubits, device MEASURE, gate after the second of two readouts worked (mean over qubits):")
    for t in ("5e-7", "1e-6", "2e-6", "4e-6", "8e-6"):
        d = load(res, f"mcm_char_t{t}.json")
        K = len(d["streams"][f"q{d['qubits'][0]}_classified"]) // 3
        ok = [np.mean(np.diff(np.array(d["streams"][f"q{q}_classified"]).reshape(K, 3)[:, 1:], axis=1) != 0)
              for q in d["qubits"]]
        print(f"      pause {float(t) * 1e6:4.1f} us: {np.mean(ok):.3f}")


def readout_length(res):
    print("\n3. readout length (MEASURE without 1->2 pulse, 10 qubits, threshold refitted, 2-fold)")
    print("   length    P(read = x)   gate worked after 0   after 1")
    for t in ("full", "8e-7", "6e-7", "4e-7"):
        d = load(res, f"mcm_char2_noshelve_{t}.json")
        K = len(d["streams"][f"q{d['qubits'][0]}_classified"]) // 2
        x = np.arange(K) % 2
        a, g0, g1 = [], [], []
        for q in d["qubits"]:
            iq = np.array(d["iq"][f"q{q}_unclassified"]).reshape(K, 2, 2)
            p1, p2 = np.zeros(K, int), np.zeros(K, int)
            for tr in (np.arange(K) % 4 < 2, np.arange(K) % 4 >= 2):
                m0, m1 = iq[tr & (x == 0), 0].mean(0), iq[tr & (x == 1), 0].mean(0)
                w, c = m1 - m0, (m1 - m0) @ (m0 + m1) / 2
                p1[~tr], p2[~tr] = iq[~tr, 0] @ w > c, iq[~tr, 1] @ w > c
            a.append(np.mean(p1 == x))
            g0.append(np.mean(p2[x == 0] == 1))
            g1.append(np.mean(p2[x == 1] == 0))
        lab = "full" if t == "full" else f"{float(t) * 1e6:.1f} us"
        print(f"   {lab:8s}  {np.mean(a):.3f}         {np.mean(g0):.3f}                {np.mean(g1):.3f}")


def repcode(res):
    import pymatching
    import repcode_mid_qcs as midq
    import repcode_qcs as endq
    sys.path.insert(0, HERE)
    from plot_repcode_qcs import NOISE

    def logical(bits, ci, cn, d, R_na):
        det, obs = ci.compile_m2d_converter().convert(measurements=bits.astype(bool), separate_observables=True)
        pred = pymatching.Matching.from_detector_error_model(
            cn.detector_error_model(decompose_errors=True)).decode_batch(det)
        ref = ci.compile_sampler().sample(1)[0].astype(int)
        dq = bits[:, R_na:]
        return (np.mean(dq[:, 0] != ref[R_na]),
                np.mean((dq.sum(1) > d / 2).astype(int) != int(ref[R_na:].sum() > d / 2)),
                np.mean(pred[:, 0] != obs[:, 0]), det.mean(), len(bits), (dq != ref[R_na:]).mean(0))

    rows = {}
    for scheme, f in (("end", "repcode_end.json"), ("mid", "repcode_mid.json"), ("echo", "repcode_mid_echo.json")):
        for j in load(res, f)["jobs"]:
            bits = np.array([[int(c) for c in s] for s in j["bitstrings"]])
            code, d, R, s = j["code"], j["d"], j["rounds"], j["state"]
            if scheme == "end":
                ci, cn, rna = endq.to_stim(code, d, R, s), endq.to_stim(code, d, R, s, NOISE), d - 1
            else:
                ci = midq.to_stim(code, d, R, s, False)
                cn, rna = midq.to_stim(code, d, R, s, False, dict(NOISE, km=1)), R * (d - 1)
            rows[(scheme, code, R, s)] = logical(bits, ci, cn, d, rna)
    print("\n4. repetition codes, d = 3, qubits 104, 95, 86, 77, 68, mean of the two logical states")
    print("   end: all qubits measured at the end (1000 shots); mid: ancillas measured every round (300 shots);")
    print("   echo: as mid with two RX(pi) on every data qubit during each readout (300 shots)")
    print("   code   R   scheme  unencoded  majority  matching  detection rate")
    for code in ("bit", "phase"):
        for R in (1, 2, 3):
            for scheme in ("end", "mid", "echo"):
                v = np.mean([rows[(scheme, code, R, s)][:4] for s in (0, 1)], axis=0)
                print(f"   {code:5s} {R}   {scheme:5s}   {v[0]:.3f}      {v[1]:.3f}     {v[2]:.3f}     {v[3]:.3f}")
    print("\n5. by logical state (data qubits 104, 86, 68): majority-vote error and error of each data qubit")
    print("   code  state  R    end: majority  data qubits          every round: majority  data qubits")
    for code in ("bit", "phase"):
        for s in (0, 1):
            for R in (1, 2, 3):
                e, m = rows[("end", code, R, s)], rows[("mid", code, R, s)]
                print(f"   {code:5s}  {s}     {R}       {e[1]:.3f}     {' '.join(f'{x:.3f}' for x in e[5])}"
                      f"          {m[1]:.3f}     {' '.join(f'{x:.3f}' for x in m[5])}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default=os.path.join(HERE, "qcs_repmid_results"))
    args = ap.parse_args()
    gate_after(args.results)
    pause(args.results)
    readout_length(args.results)
    repcode(args.results)


if __name__ == "__main__":
    main()
