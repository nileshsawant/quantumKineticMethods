"""
[[4,2,2]] error-detection experiment on Rigetti QCS (Cepheus-1): a logical Bell state with error
detection against a physical Bell state, all measurements at the end.

The [[4,2,2]] code stores two logical qubits in four physical qubits with stabilizers XXXX and ZZZZ;
any single-qubit error flips at least one of them and is detected.  With logical operators
Z1 = Z0 Z2, Z2 = Z0 Z3, X1 = X0 X3, X2 = X0 X2 (physical qubits 0..3 on a chain), the logical Bell
state (|00>+|11>)_L is the product of two physical Bell pairs on (0,1) and (2,3), which needs only
nearest-neighbor gates, and its preparation is fault-tolerant (every single fault is either detected or
harmless).  The reference is one physical Bell pair on qubits (0,1).

Each program prepares the state, applies R rounds of CZ.CZ (= identity) on every pair inside a
PRESERVE_BLOCK, so quilc keeps them as a controlled source of two-qubit gate noise, rotates to the
X, Y or Z basis and measures all its qubits.

The same abstract circuit is emitted as Quil and as Stim (noise model, ideal parities).
Writes <out>/<name>.quil and <out>/manifest.json.  Free.

Usage:
    python c422_qcs.py --rounds 0,2,5,10 --out ../qcs_quil_c422
"""

import argparse
import json
import os

BASES = ("X", "Y", "Z")


def ops(kind, R, basis):
    qubits = [0, 1, 2, 3] if kind == "enc" else [0, 1]
    pairs = [(0, 1), (2, 3)] if kind == "enc" else [(0, 1)]
    out = []
    for a, b in pairs:
        out += [("H", a), ("CX", a, b)]
    if R:
        out.append(("PRESERVE_START",))
        for _ in range(R):
            for a, b in pairs:
                out += [("CZ", a, b), ("CZ", a, b)]
        out.append(("PRESERVE_END",))
    for q in qubits:
        if basis == "X":
            out.append(("H", q))
        elif basis == "Y":
            out += [("SDG", q), ("H", q)]
    for q in qubits:
        out.append(("M", q))
    return out


def to_quil(op_list):
    nm = sum(1 for o in op_list if o[0] == "M")
    lines, k = [f"DECLARE ro BIT[{nm}]"], 0
    for o in op_list:
        if o[0] == "H":
            lines.append(f"H {o[1]}")
        elif o[0] == "SDG":
            lines.append(f"RZ(-pi/2) {o[1]}")
        elif o[0] == "CX":
            lines.append(f"CNOT {o[1]} {o[2]}")
        elif o[0] == "CZ":
            lines.append(f"CZ {o[1]} {o[2]}")
        elif o[0] == "PRESERVE_START":
            lines.append("PRAGMA PRESERVE_BLOCK")
        elif o[0] == "PRESERVE_END":
            lines.append("PRAGMA END_PRESERVE_BLOCK")
        elif o[0] == "M":
            lines.append(f"MEASURE {o[1]} ro[{k}]")
            k += 1
    return "\n".join(lines) + "\n"


def to_stim(kind, R, basis, noise=None):
    """noise: p1, p2 (depolarizing after 1q/2q gates), pm (readout flip)."""
    import stim
    nz = noise or {}
    c = stim.Circuit()
    for o in ops(kind, R, basis):
        if o[0] in ("H", "SDG"):
            c.append("H" if o[0] == "H" else "S_DAG", [o[1]])
            if nz.get("p1"):
                c.append("DEPOLARIZE1", [o[1]], nz["p1"])
        elif o[0] in ("CX", "CZ"):
            c.append(o[0], [o[1], o[2]])
            if nz.get("p2"):
                c.append("DEPOLARIZE2", [o[1], o[2]], nz["p2"])
        elif o[0] == "M":
            if nz.get("pm"):
                c.append("X_ERROR", [o[1]], nz["pm"])
            c.append("M", [o[1]])
    return c


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rounds", default="0,2,5,10")
    ap.add_argument("--shots", type=int, default=2000)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    jobs = []
    for kind in ("enc", "bare"):
        for R in map(int, args.rounds.split(",")):
            for basis in BASES:
                name = f"c422_{kind}_r{R}_{basis}"
                with open(os.path.join(args.out, name + ".quil"), "w") as fh:
                    fh.write(to_quil(ops(kind, R, basis)))
                jobs.append({"name": name, "quil": name + ".quil", "kind": "raw", "role": "c422",
                             "code": kind, "rounds": R, "basis": basis, "shots": args.shots})
    with open(os.path.join(args.out, "manifest.json"), "w") as fh:
        json.dump({"description": "[[4,2,2]] logical Bell state vs physical Bell state", "jobs": jobs},
                  fh, indent=2)
    print(f"wrote {len(jobs)} programs to {args.out}")


if __name__ == "__main__":
    main()
