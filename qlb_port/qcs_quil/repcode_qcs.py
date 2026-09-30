"""
Repetition-code experiments on Rigetti QCS (Cepheus-1) with all measurements at the end.

A distance-d repetition code stores one logical bit in d data qubits on a chain, with d-1 ancilla
qubits between them.  Each of R rounds applies the parity-check gates: every ancilla collects the
parity of its two neighbors (Z Z for the bit-flip code, X X for the phase-flip code).  The device
cannot measure in the middle of a circuit (only ProtoQuil runs with many shots), so the ancillas are
not measured between rounds; they are measured once at the end together with the data qubits.  An
ancilla then holds the parity of its check summed over all rounds, which is ideally 0.  The rounds
thus act as a realistic source of gate noise, and the final measurement gives d-1 ancilla checks
plus the d data values for decoding.  R = 0 is preparation and readout only.

The unencoded reference is data qubit 0 of each code, which sees the same gates.

The same abstract circuit is emitted as Quil (for the device) and as a Stim circuit (for noisy
simulation, the detector error model, and matching decoders), with identical measurement records.

Writes <out>/<name>.quil and <out>/manifest.json.  Free (no QPU time).

Usage:
    python repcode_qcs.py --distances 3,5,7,9 --rounds 0,1,3,6 --out ../qcs_quil_repcode
"""

import argparse
import json
import os

CODES = ("bit", "phase")


def layout(d):
    """Logical chain indices: data on even, ancillas on odd positions."""
    return list(range(0, 2 * d - 1, 2)), list(range(1, 2 * d - 2, 2))


def ops(code, d, R, state):
    """Abstract op list: ('X',q) ('H',q) ('CX',c,t) ('M',q) ('ROUND_IDLE',[data])."""
    data, anc = layout(d)
    out = []
    for q in data:
        if state:
            out.append(("X", q))
        if code == "phase":
            out.append(("H", q))
    for _ in range(R):
        if code == "bit":                               # ancilla collects Z parity: CX data -> anc
            for k, a in enumerate(anc):
                out.append(("CX", data[k], a))
            for k, a in enumerate(anc):
                out.append(("CX", data[k + 1], a))
        else:                                           # X parity: H(anc), CX anc -> data, H(anc)
            for a in anc:
                out.append(("H", a))
            for k, a in enumerate(anc):
                out.append(("CX", a, data[k]))
            for k, a in enumerate(anc):
                out.append(("CX", a, data[k + 1]))
            for a in anc:
                out.append(("H", a))
        out.append(("ROUND_IDLE", list(data)))
    if code == "phase":                                 # all gates before any MEASURE (ProtoQuil)
        for q in data:
            out.append(("H", q))
    for a in anc:
        out.append(("M", a))
    for q in data:
        out.append(("M", q))
    return out


def to_quil(op_list):
    nm = sum(1 for o in op_list if o[0] == "M")
    lines = [f"DECLARE ro BIT[{nm}]"]
    k = 0
    for o in op_list:
        if o[0] == "X":
            lines.append(f"X {o[1]}")
        elif o[0] == "H":
            lines.append(f"H {o[1]}")
        elif o[0] == "CX":
            lines.append(f"CNOT {o[1]} {o[2]}")
        elif o[0] == "M":
            lines.append(f"MEASURE {o[1]} ro[{k}]")
            k += 1
    return "\n".join(lines) + "\n"


def to_stim(code, d, R, state, noise=None):
    """Stim circuit with detectors and one logical observable (final value of data qubit 0).
    Record: d-1 ancilla results, then d data results.  Detectors: each ancilla (summed parity,
    ideally 0) and each final data parity f_j xor f_{j+1}.
    noise: dict p1, p2 (depolarizing after 1q/2q gates), pm (measurement flip), pz, px (data errors
    per round)."""
    import stim
    nz = noise or {}
    c = stim.Circuit()
    data, anc = layout(d)
    for o in ops(code, d, R, state):
        if o[0] in ("X", "H"):
            c.append(o[0], [o[1]])
            if nz.get("p1"):
                c.append("DEPOLARIZE1", [o[1]], nz["p1"])
        elif o[0] == "CX":
            c.append("CX", [o[1], o[2]])
            if nz.get("p2"):
                c.append("DEPOLARIZE2", [o[1], o[2]], nz["p2"])
        elif o[0] == "ROUND_IDLE":
            if nz.get("pz"):
                c.append("Z_ERROR", o[1], nz["pz"])
            if nz.get("px"):
                c.append("X_ERROR", o[1], nz["px"])
        elif o[0] == "M":
            if nz.get("pm"):
                c.append("X_ERROR", [o[1]], nz["pm"])
            c.append("M", [o[1]])
    na, nd = len(anc), len(data)
    total = na + nd

    def rec(idx):
        return stim.target_rec(idx - total)

    for j in range(na):
        c.append("DETECTOR", [rec(j)], [2 * j + 1, 0])
    for j in range(na):
        c.append("DETECTOR", [rec(na + j), rec(na + j + 1)], [2 * j + 1, 1])
    c.append("OBSERVABLE_INCLUDE", [rec(na)], 0)
    return c


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--distances", default="3,5,7,9")
    ap.add_argument("--rounds", default="0,1,3,6")
    ap.add_argument("--shots", type=int, default=2000)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    jobs = []
    for code in CODES:
        for d in map(int, args.distances.split(",")):
            for R in map(int, args.rounds.split(",")):
                for state in (0, 1):
                    name = f"rep_{code}_d{d}_r{R}_s{state}"
                    with open(os.path.join(args.out, name + ".quil"), "w") as fh:
                        fh.write(to_quil(ops(code, d, R, state)))
                    data, anc = layout(d)
                    jobs.append({"name": name, "quil": name + ".quil", "kind": "raw", "role": "rep",
                                 "code": code, "d": d, "rounds": R, "state": state,
                                 "n_qubits": len(data) + len(anc), "shots": args.shots})
    manifest = {"description": "repetition-code experiments, all measurements at the end",
                "jobs": jobs}
    with open(os.path.join(args.out, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"wrote {len(jobs)} programs to {args.out}")


if __name__ == "__main__":
    main()
