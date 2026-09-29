"""
Clifford data regression (CDR) for the twirled four-site velocity circuits on Rigetti QCS.

For each target program the pinned native program (RX in {+-pi/2, pi}, RZ(theta), CZ) is taken from
quilc once.  Training circuits keep every RX and CZ and replace a random fraction of the
non-Clifford RZ angles by random multiples of pi/2.  RZ is a virtual frame change on the device, so a
training circuit runs the same physical pulse sequence as the target and differs only in phases.  The
ideal value of every training circuit is computed here exactly (4-qubit state vector).  Every
training circuit is Pauli-twirled once (a fresh twirl each), the target is run as 16 twirled
variants; a linear fit ideal = a + b * measured over the training set is then applied to the target.

Writes <out>/<name>_tw<k>.quil (target variants), <out>/<name>_cdr<k>.quil (training circuits) and
<out>/manifest.json with the exact training values.  pyQuil environment with quilc; free.

Usage:
    python cdr_qcs.py --src ../qcs_quil_n2_fourier --jobs alpha_t1,alpha_t2,alpha_t3,alpha_t4 \
        --qubits 103,102,100,101 --train 32 --train-shots 250 --variants 16 --shots 125 --out ../qcs_quil_n2_cdr
"""

import argparse
import json
import math
import os
import random

import numpy as np

from run_qcs import _pin_text
from twirl_qcs import preserve, twirl


def _angle(text):
    return eval(text, {"__builtins__": {}}, {"pi": math.pi})    # quilc angle expressions, e.g. (3*pi)/4


def parse(native_text):
    """[(name, angle or None, [qubits])] for the gate lines, and the measured qubits in ro order."""
    gates, meas = [], []
    for line in native_text.splitlines():
        s = line.strip()
        if s.startswith(("RX(", "RZ(")):
            name, rest = s.split("(", 1)
            ang, q = rest.rsplit(")", 1)
            gates.append((name, _angle(ang), [int(q)]))
        elif s.startswith("CZ "):
            gates.append(("CZ", None, [int(x) for x in s.split()[1:3]]))
        elif s.startswith("MEASURE"):
            meas.append(int(s.split()[1]))
    return gates, meas


def emit(gates):
    out = []
    for name, ang, q in gates:
        out.append(f"CZ {q[0]} {q[1]}" if name == "CZ" else f"{name}({ang:.17g}) {q[0]}")
    return out


def ideal_value(gates, meas_qubit):
    """P(1) - P(0) of the measured qubit for the ideal circuit from |0000>."""
    qs = sorted({q for _, _, qq in gates for q in qq} | {meas_qubit})
    idx = {q: k for k, q in enumerate(qs)}
    n = len(qs)
    psi = np.zeros([2] * n, complex)
    psi[(0,) * n] = 1
    for name, ang, q in gates:
        if name == "CZ":
            a, b = idx[q[0]], idx[q[1]]
            sl = [slice(None)] * n
            sl[a], sl[b] = 1, 1
            psi[tuple(sl)] *= -1
            continue
        c, s = math.cos(ang / 2), math.sin(ang / 2)
        U = (np.array([[c, -1j * s], [-1j * s, c]]) if name == "RX"
             else np.array([[np.exp(-1j * ang / 2), 0], [0, np.exp(1j * ang / 2)]]))
        k = idx[q[0]]
        psi = np.moveaxis(np.tensordot(U, psi, axes=([1], [k])), 0, k)
    p = np.abs(psi) ** 2
    p1 = p.take(1, axis=idx[meas_qubit]).sum()
    return float(2 * p1 - 1)


def is_clifford(ang):
    return abs((ang / (math.pi / 2)) - round(ang / (math.pi / 2))) < 1e-9


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", required=True)
    ap.add_argument("--jobs", required=True)
    ap.add_argument("--qubits", required=True)
    ap.add_argument("--train", type=int, default=32)
    ap.add_argument("--train-shots", type=int, default=250)
    ap.add_argument("--replace", type=float, default=0.7, help="fraction of non-Clifford RZ replaced")
    ap.add_argument("--variants", type=int, default=16)
    ap.add_argument("--shots", type=int, default=125)
    ap.add_argument("--seed", type=int, default=2027)
    ap.add_argument("--device", default="Cepheus-1-108Q")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    from pyquil import Program, get_qc
    qc = get_qc(args.device, compiler_timeout=120)
    qmap = {i: int(p) for i, p in enumerate(args.qubits.split(","))}
    src = json.load(open(os.path.join(args.src, "manifest.json")))
    ref = {j["quil"][:-5]: j for j in src["jobs"]}
    rng = random.Random(args.seed)
    os.makedirs(args.out, exist_ok=True)

    jobs = []
    for name in args.jobs.split(","):
        base = dict(ref[name])
        raw = _pin_text(open(os.path.join(args.src, name + ".quil")).read(), qmap)
        native = str(qc.compiler.quil_to_native_quil(Program(raw)))
        head = [l.strip() for l in native.splitlines() if l.strip().startswith(("DECLARE", "PRAGMA"))]
        mlines = [l.strip() for l in native.splitlines() if l.strip().startswith("MEASURE")]
        gates, meas = parse(native)
        assert all(is_clifford(a) for n, a, _ in gates if n == "RX"), "non-Clifford RX in native program"
        nonclif = [i for i, (n, a, _) in enumerate(gates) if n == "RZ" and not is_clifford(a)]
        target_ideal = ideal_value(gates, meas[0])
        print(f"{name}: {sum(g[0] == 'CZ' for g in gates)} CZ, {len(nonclif)} non-Clifford RZ; "
              f"ideal {target_ideal:+.4f} (exact {base['exact']:+.4f})")

        def program(gl):
            return "\n".join(head + emit(gl) + mlines) + "\n"

        for k in range(args.variants):
            fname = f"{name}_tw{k:02d}.quil"
            with open(os.path.join(args.out, fname), "w") as fh:
                fh.write(preserve(twirl(program(gates), rng)))
            jobs.append({**base, "name": f"{name}_tw{k:02d}", "quil": fname, "group": name, "variant": k,
                         "native": True, "shots": args.shots, "role": "target", "ideal": target_ideal})
        ideals = []
        for k in range(args.train):
            g2 = list(gates)
            for i in nonclif:
                if rng.random() < args.replace:
                    g2[i] = ("RZ", rng.choice((0.0, math.pi / 2, math.pi, -math.pi / 2)), g2[i][2])
            val = ideal_value(g2, meas[0])
            ideals.append(val)
            fname = f"{name}_cdr{k:02d}.quil"
            with open(os.path.join(args.out, fname), "w") as fh:
                fh.write(preserve(twirl(program(g2), rng)))
            jobs.append({**base, "name": f"{name}_cdr{k:02d}", "quil": fname, "group": f"{name}_cdr",
                         "variant": k, "native": True, "shots": args.train_shots, "role": "train",
                         "ideal": val, "exact": val})
        print(f"   training ideals: min {min(ideals):+.3f} max {max(ideals):+.3f} "
              f"mean {np.mean(ideals):+.3f} sd {np.std(ideals):.3f}")

    manifest = {k: v for k, v in src.items() if k != "jobs"}
    manifest.update({"description": "Clifford data regression: twirled targets + twirled near-Clifford training",
                     "qubits_native": args.qubits, "seed": args.seed, "replace": args.replace, "jobs": jobs})
    with open(os.path.join(args.out, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"wrote {len(jobs)} programs to {args.out}")


if __name__ == "__main__":
    main()
