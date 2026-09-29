"""
Pauli twirling (randomized compiling) of Quil programs for the Rigetti QCS runs.

Each selected program is pinned to physical qubits and compiled once by quilc to native Quil
{RX, RZ, CZ}.  Every CZ is then dressed with random Paulis: P_a (x) P_b before it and
CZ (P_a (x) P_b) CZ^dagger after it, so each variant implements the same unitary up to a global
phase while coherent errors of the CZ layers are turned into stochastic Pauli errors when the
variants are averaged.  The native gates are wrapped in PRAGMA PRESERVE_BLOCK, so quilc passes them
through verbatim and the inserted Paulis are not optimized away (run_qcs.py skips pinning for jobs
flagged "native").

Writes <out>/<name>_native.quil (untwirled control), <out>/<name>_tw<k>.quil (variants), and
<out>/manifest.json.  Run in the pyQuil environment with a quilc server (compilation is free).

Options for zero-noise extrapolation and rescaling:
  --scales 1,3,5   also fold every CZ into 3 or 5 CZs (CZ^2 = I) before twirling: same unitary, more
                   two-qubit noise; groups <name>_x3, <name>_x5
  --mirror         twirled <name>_mirror = native program followed by its inverse (ideal: all |0>)
  --cal            basis-state readout calibration of the measured qubits (native, untwirled)
  --control-shots 0 skips the untwirled control

Usage:
    python twirl_qcs.py --src ../qcs_quil_n2_t4check --jobs even_all_t4,step1_q1_t4 \
        --qubits 103,102,100,101 --variants 16 --shots 125 --control-shots 2000 --out ../qcs_quil_n2_twirl
"""

import argparse
import json
import math
import os
import random

from run_qcs import _pin_text

# CZ (P_a x P_b) CZ^dagger = (P_a Z^[P_b has X]) x (P_b Z^[P_a has X]), up to sign
_HAS_X = {"I": False, "X": True, "Y": True, "Z": False}
_TIMES_Z = {"I": "Z", "X": "Y", "Y": "X", "Z": "I"}
_NATIVE = {"I": [], "X": ["RX(pi) {q}"], "Z": ["RZ(pi) {q}"], "Y": ["RZ(pi) {q}", "RX(pi) {q}"]}


def _pauli_lines(p, q):
    return [s.format(q=q) for s in _NATIVE[p]]


def twirl(native_text, rng):
    out = []
    for line in native_text.splitlines():
        toks = line.split()
        if toks and toks[0] == "CZ":
            a, b = toks[1], toks[2]
            pa, pb = rng.choice("IXYZ"), rng.choice("IXYZ")
            qa = _TIMES_Z[pa] if _HAS_X[pb] else pa
            qb = _TIMES_Z[pb] if _HAS_X[pa] else pb
            out += _pauli_lines(pa, a) + _pauli_lines(pb, b) + [line]
            out += _pauli_lines(qa, a) + _pauli_lines(qb, b)
        else:
            out.append(line)
    return "\n".join(out) + "\n"


def preserve(native_text):
    """Wrap the gates in PRESERVE_BLOCK so quilc submits them verbatim; declarations and
    measurements stay outside the block."""
    head, body, meas = [], [], []
    for line in native_text.splitlines():
        s = line.strip()
        if not s or s == "HALT":
            continue
        if s.startswith(("DECLARE", "PRAGMA")):
            head.append(s)
        elif s.startswith("MEASURE"):
            meas.append(s)
        else:
            body.append(s)
    return "\n".join(head + ["PRAGMA PRESERVE_BLOCK"] + body + ["PRAGMA END_PRESERVE_BLOCK"] + meas) + "\n"


def fold(native_text, scale):
    """Replace every CZ by `scale` (odd) consecutive CZs."""
    out = []
    for line in native_text.splitlines():
        out += [line] * scale if line.split()[:1] == ["CZ"] else [line]
    return "\n".join(out) + "\n"


def _inverse_native(line):
    """Native inverse: RX(+-pi/2) -> RX(-+pi/2), RX(pi) unchanged (up to phase), RZ(a) -> RZ(-a), CZ."""
    name, rest = line.split("(", 1) if "(" in line else (line.split()[0], None)
    if rest is None:
        return line
    angle, qubit = rest.rsplit(")", 1)
    val = eval(angle, {"__builtins__": {}}, {"pi": math.pi})    # quilc angle expressions, e.g. (3*pi)/4
    if name == "RX":
        if abs(abs(val) - math.pi) < 1e-9:
            return line
        return f"RX({'-' if val > 0 else ''}pi/2){qubit}"
    return f"RZ({-val:.17g}){qubit}"


def mirror(native_text):
    """Gates followed by their inverse, then the original measurements: ideally all qubits in |0>."""
    lines = [l.strip() for l in native_text.splitlines() if l.strip() and l.strip() != "HALT"]
    head = [l for l in lines if l.startswith(("DECLARE", "PRAGMA"))]
    meas = [l for l in lines if l.startswith("MEASURE")]
    gates = [l for l in lines if not l.startswith(("DECLARE", "PRAGMA", "MEASURE"))]
    inv = [_inverse_native(l) for l in reversed(gates)]
    return "\n".join(head + gates + inv + meas) + "\n"


def calibration(native_text):
    """Native programs preparing every basis state of the measured qubits (ro order)."""
    lines = [l.strip() for l in native_text.splitlines()]
    meas = [l for l in lines if l.startswith("MEASURE")]
    decl = [l for l in lines if l.startswith("DECLARE")]
    qubits = [l.split()[1] for l in meas]
    progs = []
    for b in range(2 ** len(qubits)):
        xs = [f"RX(pi) {q}" for k, q in enumerate(qubits) if (b >> k) & 1]
        body = ["PRAGMA PRESERVE_BLOCK"] + xs + ["PRAGMA END_PRESERVE_BLOCK"] if xs else []
        progs.append((b, "\n".join(decl + body + meas) + "\n"))
    return progs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", required=True, help="folder with manifest.json and the .quil files")
    ap.add_argument("--jobs", required=True, help="comma-separated quil basenames (without .quil)")
    ap.add_argument("--qubits", required=True)
    ap.add_argument("--variants", type=int, default=16)
    ap.add_argument("--shots", type=int, default=125, help="shots per twirled variant")
    ap.add_argument("--control-shots", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--scales", default="1", help="comma-separated odd CZ fold factors")
    ap.add_argument("--mirror", action="store_true")
    ap.add_argument("--cal", action="store_true")
    ap.add_argument("--cal-shots", type=int, default=2000)
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
    scales = [int(s) for s in args.scales.split(",")]
    for name in args.jobs.split(","):
        base = dict(ref[name])
        raw = _pin_text(open(os.path.join(args.src, name + ".quil")).read(), qmap)
        native = str(qc.compiler.quil_to_native_quil(Program(raw)))
        ncz = sum(1 for l in native.splitlines() if l.startswith("CZ "))
        if args.control_shots > 0:
            with open(os.path.join(args.out, f"{name}_native.quil"), "w") as fh:
                fh.write(preserve(native))
            jobs.append({**base, "name": f"{name}_native", "quil": f"{name}_native.quil", "group": name,
                         "variant": -1, "native": True, "shots": args.control_shots, "scale": 1})
        variants = [(s, name if s == 1 else f"{name}_x{s}", fold(native, s), {}) for s in scales]
        if args.mirror:
            variants.append((1, f"{name}_mirror", mirror(native),
                             {"exact": -1.0, "exact_alpha": -1.0, "mirror": True}))
        for s, group, prog, extra in variants:
            for k in range(args.variants):
                fname = f"{group}_tw{k:02d}.quil"
                with open(os.path.join(args.out, fname), "w") as fh:
                    fh.write(preserve(twirl(prog, rng)))
                jobs.append({**base, **extra, "name": f"{group}_tw{k:02d}", "quil": fname, "group": group,
                             "variant": k, "native": True, "shots": args.shots, "scale": s})
        print(f"{name}: {ncz} CZ; groups {[v[1] for v in variants]}, {args.variants} twirled variants each")
    if args.cal:
        for b, prog in calibration(native):
            fname = f"cal_ro{b}.quil"
            with open(os.path.join(args.out, fname), "w") as fh:
                fh.write(prog)
            jobs.append({"name": f"cal_ro{b}", "kind": "cal", "t": 0, "quil": fname, "npos": base["npos"],
                         "ro_bits": base["ro_bits"], "prepared": b, "exact": 1.0, "group": "cal",
                         "variant": -1, "native": True, "shots": args.cal_shots})

    manifest = {k: v for k, v in src.items() if k != "jobs"}
    manifest.update({"description": "Pauli-twirled native programs (verbatim via PRESERVE_BLOCK)",
                     "qubits_native": args.qubits, "twirl_seed": args.seed, "jobs": jobs})
    with open(os.path.join(args.out, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"wrote {len(jobs)} programs to {args.out}")


if __name__ == "__main__":
    main()
