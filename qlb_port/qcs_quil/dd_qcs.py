"""
Dynamical decoupling of native Quil programs for the Rigetti QCS runs, using FENCE for alignment.

Each selected program is pinned and compiled once by quilc to native Quil.  The gate list is then
split into layers, one per CZ, separated by FENCE on all used qubits, so that every CZ runs while the
other qubits wait.  A qubit that has no gate at all during a run of consecutive CZ layers receives an
RX(pi) at the start of each of those layers, in pairs: X U1 X U2 = U2 U1^dagger, which cancels a
quasi-static phase (e.g. a detuning) accumulated while idle if the two layers last equally long.  An
unpaired last idle layer of an odd run is left alone.  Quil offers no finer timing control
(DELAY is not accepted by the compiler path used here), so this is the simplest echo that can be
expressed.  Three versions of each program are written, all wrapped in PRESERVE_BLOCK:
    <name>_orig   the native program as compiled (no fences)
    <name>_fence  fenced layers, no echo pulses (control for the cost of the fences)
    <name>_dd     fenced layers with the paired echo pulses
and repeated --rounds times in the manifest, interleaved, to sample the drift.

Usage (pyQuil environment with a quilc server; free):
    python dd_qcs.py --src ../qcs_quil_n2_t4check --jobs step1_q1_t4,even_all_t4 \
        --qubits 103,102,100,101 --rounds 4 --shots 2000 --out ../qcs_quil_n2_dd
"""

import argparse
import json
import os

from run_qcs import _pin_text
from twirl_qcs import preserve


def _split(native_text):
    lines = [l.strip() for l in native_text.splitlines() if l.strip() and l.strip() != "HALT"]
    head = [l for l in lines if l.startswith(("DECLARE", "PRAGMA"))]
    meas = [l for l in lines if l.startswith("MEASURE")]
    gates = [l for l in lines if not l.startswith(("DECLARE", "PRAGMA", "MEASURE"))]
    return head, gates, meas


def _qubits(gate):
    toks = gate.split()
    return [int(t) for t in (toks[1:3] if toks[0] == "CZ" else toks[-1:])]


def layers(gates):
    """Segments ending with a CZ (the tail after the last CZ is its own segment)."""
    segs, cur = [], []
    for g in gates:
        cur.append(g)
        if g.startswith("CZ "):
            segs.append(cur)
            cur = []
    if cur:
        segs.append(cur)
    return segs


def build(native_text, echo):
    head, gates, meas = _split(native_text)
    segs = layers(gates)
    used = sorted({q for g in gates for q in _qubits(g)})
    fence = "FENCE " + " ".join(map(str, used))
    busy = [{q for g in s for q in _qubits(g)} for s in segs]
    is_cz = [s[-1].startswith("CZ ") for s in segs]
    pulses = [[] for _ in segs]
    n_echo = 0
    if echo:
        for q in used:
            run = []
            for k in range(len(segs) + 1):
                if k < len(segs) and is_cz[k] and q not in busy[k]:
                    run.append(k)
                    continue
                for k2 in run[: len(run) - len(run) % 2]:
                    pulses[k2].append(f"RX(pi) {q}")
                    n_echo += 1
                run = []
    body = []
    for k, s in enumerate(segs):
        body.append(fence)
        body += s[:-1] + pulses[k] + [s[-1]] if is_cz[k] else s
    body.append(fence)
    return "\n".join(head + body + meas) + "\n", n_echo


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", required=True)
    ap.add_argument("--jobs", required=True)
    ap.add_argument("--qubits", required=True)
    ap.add_argument("--rounds", type=int, default=4)
    ap.add_argument("--shots", type=int, default=2000)
    ap.add_argument("--device", default="Cepheus-1-108Q")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    from pyquil import Program, get_qc
    qc = get_qc(args.device, compiler_timeout=120)
    qmap = {i: int(p) for i, p in enumerate(args.qubits.split(","))}
    src = json.load(open(os.path.join(args.src, "manifest.json")))
    ref = {j["quil"][:-5]: j for j in src["jobs"]}
    os.makedirs(args.out, exist_ok=True)

    progs = []
    for name in args.jobs.split(","):
        raw = _pin_text(open(os.path.join(args.src, name + ".quil")).read(), qmap)
        native = str(qc.compiler.quil_to_native_quil(Program(raw)))
        fenced, _ = build(native, echo=False)
        dd, n_echo = build(native, echo=True)
        for tag, text in (("orig", native), ("fence", fenced), ("dd", dd)):
            fname = f"{name}_{tag}.quil"
            with open(os.path.join(args.out, fname), "w") as fh:
                fh.write(preserve(text) if tag == "orig" else _preserve_fenced(text))
            progs.append((name, tag, fname))
        print(f"{name}: {sum(l.startswith('CZ ') for l in native.splitlines())} CZ layers, "
              f"{n_echo} echo pulses inserted")

    jobs = []
    for r in range(args.rounds):
        for name, tag, fname in progs:
            jobs.append({**ref[name], "name": f"{name}_{tag}_r{r}", "quil": fname, "group": f"{name}_{tag}",
                         "variant": r, "native": True, "shots": args.shots})
    manifest = {k: v for k, v in src.items() if k != "jobs"}
    manifest.update({"description": "dynamical decoupling test: orig / fenced / fenced+echo, interleaved",
                     "qubits_native": args.qubits, "jobs": jobs})
    with open(os.path.join(args.out, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"wrote {len(progs)} programs, {len(jobs)} jobs to {args.out}")


def _preserve_fenced(text):
    """PRESERVE_BLOCK wrapper that keeps FENCE lines inside the block."""
    head, body, meas = [], [], []
    for line in text.splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith(("DECLARE", "PRAGMA")):
            head.append(s)
        elif s.startswith("MEASURE"):
            meas.append(s)
        else:
            body.append(s)
    return "\n".join(head + ["PRAGMA PRESERVE_BLOCK"] + body + ["PRAGMA END_PRESERVE_BLOCK"] + meas) + "\n"


if __name__ == "__main__":
    main()
