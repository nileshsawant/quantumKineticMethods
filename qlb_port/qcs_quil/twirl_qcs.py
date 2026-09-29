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

Usage:
    python twirl_qcs.py --src ../qcs_quil_n2_t4check --jobs even_all_t4,step1_q1_t4 \
        --qubits 103,102,100,101 --variants 16 --shots 125 --control-shots 2000 --out ../qcs_quil_n2_twirl
"""

import argparse
import json
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
        ncz = sum(1 for l in native.splitlines() if l.startswith("CZ "))
        with open(os.path.join(args.out, f"{name}_native.quil"), "w") as fh:
            fh.write(preserve(native))
        jobs.append({**base, "name": f"{name}_native", "quil": f"{name}_native.quil", "group": name,
                     "variant": -1, "native": True, "shots": args.control_shots})
        for k in range(args.variants):
            fname = f"{name}_tw{k:02d}.quil"
            with open(os.path.join(args.out, fname), "w") as fh:
                fh.write(preserve(twirl(native, rng)))
            jobs.append({**base, "name": f"{name}_tw{k:02d}", "quil": fname, "group": name,
                         "variant": k, "native": True, "shots": args.shots})
        print(f"{name}: {ncz} CZ, {args.variants} twirled variants")

    manifest = {k: v for k, v in src.items() if k != "jobs"}
    manifest.update({"description": "Pauli-twirled native programs (quilc bypassed at submission)",
                     "qubits_native": args.qubits, "twirl_seed": args.seed, "jobs": jobs})
    with open(os.path.join(args.out, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"wrote {len(jobs)} programs to {args.out}")


if __name__ == "__main__":
    main()
