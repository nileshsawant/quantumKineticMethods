"""
Relaxation-aware model of the repetition codes run on Cepheus-1 (qubits 104, 95, 86, 77, 68),
compared with the device results in qcs_repmid_results/ (ancillas read at the end, and every round).

The native programs that ran on the device (qcs_quil_repend_B/native, the first shot of
qcs_quil_repmid_gapB30/native) are simulated on a density matrix (Qiskit Aer) with:
  * every qubit decaying with its measured T1 and dephasing with its measured Ramsey T2 while it
    waits, during gates, and while it or another qubit is read out (Kraus thermal relaxation);
  * gates timed as on the device: RZ virtual, RX 40 ns (assumed), CZ from its pulse length; a
    mid-circuit readout from its pulse length plus the 1 us pause, final readouts at the end;
  * depolarizing gate errors from the randomized-benchmarking (1q) and CZ fidelities;
  * a symmetric readout error per qubit: 1 - fRO for the standard readout (end-only programs),
    and 1 - P(read = x) of the readout without the 1->2 pulse (mcm_char2_noshelve.json) for the
    every-round programs, which used that readout.
Unlike Pauli models such as SI1000, T1 decay moves |1> to |0> only.  Leakage, crosstalk, and
coherent gate errors are not included.

Usage (aer-gpu venv):
    python sim_t1t2_repcode.py [--shots 4000]
"""

import argparse
import json
import math
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "qcs_repmid_results")
T_RX = 40e-9
PAUSE = 1e-6


def parse(native, n_slots=None):
    """[(op, qubits, angle, slot)] of a native program; stops at the first DELAY (end of shot 0)."""
    out = []
    for ln in native.splitlines():
        ln = ln.strip()
        if not ln or ln.startswith(("DECLARE", "PRAGMA", "HALT")):
            continue
        if ln.startswith("DELAY"):
            break
        if ln.startswith("MEASURE"):
            t = ln.split()
            k = int(re.search(r"\[(\d+)\]", t[2]).group(1))
            if n_slots is not None and k >= n_slots:
                break
            out.append(("M", [int(t[1])], None, k))
            continue
        m = re.match(r"([A-Z]+)(?:\((.*)\))?\s+([\d\s]+)$", ln)
        g, ang, qs = m.group(1), m.group(2), [int(q) for q in m.group(3).split()]
        a = eval(ang, {"pi": math.pi}) if ang else None
        out.append((g, qs, a, None))
    return out


def build(ops, cal, ro_err, mid, noisy=True):
    from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister
    from qiskit.quantum_info import Kraus
    from qiskit_aer.noise import depolarizing_error, thermal_relaxation_error

    def chan(err):
        return Kraus(err.to_quantumchannel()).to_instruction()
    qs = cal["qubits"]
    idx = {q: i for i, q in enumerate(qs)}
    n_slots = 1 + max(o[3] for o in ops if o[0] == "M")
    qc = QuantumCircuit(QuantumRegister(len(qs)), ClassicalRegister(n_slots))
    t = {q: 0.0 for q in qs}

    def relax(q, dt):
        if noisy and dt > 1e-12:
            T1, T2 = cal["T1"][str(q)], min(cal["T2_ramsey"][str(q)], 2 * cal["T1"][str(q)])
            qc.append(chan(thermal_relaxation_error(T1, T2, dt)), [idx[q]])

    def edge(a, b):
        return f"{a}-{b}" if f"{a}-{b}" in cal["cz_s"] else f"{b}-{a}"

    last_meas = {}
    for g, qq, a, k in ops:
        if g == "RZ":
            qc.rz(a, idx[qq[0]])
            continue
        start = max(t[q] for q in qq)
        for q in qq:
            relax(q, start - t[q])
        if g == "RX":
            qc.rx(a, idx[qq[0]])
            relax(qq[0], T_RX)
            if noisy:
                qc.append(chan(depolarizing_error(2 * (1 - cal["fRB_1q"][str(qq[0])]), 1)),
                          [idx[qq[0]]])
            t[qq[0]] = start + T_RX
        elif g == "CZ":
            e = edge(*qq)
            qc.cz(idx[qq[0]], idx[qq[1]])
            for q in qq:
                relax(q, cal["cz_s"][e])
            if noisy:
                qc.append(chan(depolarizing_error(4 / 3 * (1 - cal["fCZ"][e]), 2)),
                          [idx[qq[0]], idx[qq[1]]])
            for q in qq:
                t[q] = start + cal["cz_s"][e]
        elif g == "M":
            q = qq[0]
            qc.measure(idx[q], k)
            last_meas[q] = k
            t[q] = start + cal["readout_s"][str(q)] + (PAUSE if mid else 0.0)
            relax(q, t[q] - start)
    return qc


def run_aer(qc, shots, seed):
    from qiskit_aer import AerSimulator
    sim = AerSimulator(method="density_matrix", seed_simulator=seed)
    counts = sim.run(qc, shots=shots).result().get_counts()
    rows = []
    for key, c in counts.items():
        bits = [int(b) for b in key.replace(" ", "")[::-1]]
        rows += [bits] * c
    return np.array(rows)


def detectors(bits, scheme, d, R):
    na = d - 1
    if scheme == "end":
        anc, dat = bits[:, :na], bits[:, na:]
        return np.hstack([anc, dat[:, :-1] ^ dat[:, 1:]])
    m = bits[:, :R * na].reshape(len(bits), R, na)
    dat = bits[:, R * na:]
    z = np.zeros_like(m[:, 0])
    dets = [m[:, r] ^ (m[:, r - 2] if r >= 2 else z) for r in range(R)]
    prev = m[:, R - 1] ^ (m[:, R - 2] if R >= 2 else z)
    dets.append(dat[:, :-1] ^ dat[:, 1:] ^ prev)
    return np.hstack(dets)


def metrics(bits, ref, scheme, d, R):
    rna = d - 1 if scheme == "end" else R * (d - 1)
    dq, rd = bits[:, rna:], ref[rna:]
    return {"unenc": np.mean(dq[:, 0] != rd[0]),
            "maj": np.mean((dq.sum(1) > d / 2).astype(int) != int(rd.sum() > d / 2)),
            "data": (dq != rd).mean(0), "det": detectors(bits, scheme, d, R).mean()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--shots", type=int, default=4000)
    args = ap.parse_args()
    cal = json.load(open(os.path.join(RES, "calibration_104_95_86_77_68.json")))
    ch = json.load(open(os.path.join(RES, "mcm_char2_noshelve.json")))
    K = len(ch["streams"][f"q{ch['qubits'][0]}_classified"]) // 2
    x = np.arange(K) % 2
    ro_mid = {q: 1 - np.mean(np.array(ch["streams"][f"q{q}_classified"]).reshape(K, 2)[:, 0] == x)
              for q in cal["qubits"]}
    ro_end = {q: 1 - cal["fRO"][str(q)] for q in cal["qubits"]}
    rng = np.random.default_rng(5)
    jobs = {"end": (json.load(open(os.path.join(RES, "repcode_end.json")))["jobs"],
                    os.path.join(HERE, "qcs_quil_repend_B", "native")),
            "mid": (json.load(open(os.path.join(RES, "repcode_mid.json")))["jobs"],
                    os.path.join(HERE, "qcs_quil_repmid_gapB30", "native"))}
    print("model: T1/T2 relaxation with device timing, depolarizing gate errors, symmetric readout errors")
    print("dev = Cepheus-1, model = this simulation;  majority-vote error and error of each data qubit (104, 86, 68)")
    print("scheme code  state R   majority dev/model   data qubits dev            data qubits model"
          "          detection dev/model")
    for scheme, (jl, ndir) in jobs.items():
        for j in sorted(jl, key=lambda j: (j["code"], j["state"], j["rounds"])):
            code, d, R, s = j["code"], j["d"], j["rounds"], j["state"]
            n_slots = (d - 1) + d if scheme == "end" else R * (d - 1) + d
            ops = parse(open(os.path.join(ndir, j["quil"])).read(), n_slots)
            ro = ro_end if scheme == "end" else ro_mid
            ref = run_aer(build(ops, cal, ro, scheme == "mid", noisy=False), 1, 1)[0]
            sim = run_aer(build(ops, cal, ro, scheme == "mid"), args.shots, 7)
            slot_q = {o[3]: o[1][0] for o in ops if o[0] == "M"}
            flip = np.array([rng.random(len(sim)) < ro[slot_q[k]] for k in range(n_slots)]).T
            sim = sim ^ flip
            dev = np.array([[int(c) for c in b] for b in j["bitstrings"]])
            a, b = metrics(dev, ref, scheme, d, R), metrics(sim, ref, scheme, d, R)
            print(f"{scheme:5s} {code:5s}  {s}    {R}    {a['maj']:.3f} / {b['maj']:.3f}       "
                  f"{' '.join(f'{v:.3f}' for v in a['data'])}      {' '.join(f'{v:.3f}' for v in b['data'])}"
                  f"      {a['det']:.3f} / {b['det']:.3f}")


if __name__ == "__main__":
    main()
