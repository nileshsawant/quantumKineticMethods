"""Where does |1> go on data qubits 86 and 68 under CZ gates (Cepheus-1)?

Each program repeats one case K times with a passive reset (DELAY gap) between shots and reads
all five chain qubits 104, 95, 86, 77, 68 once per shot.  Cases:
    ref s       all five qubits in |s>, s = 0, 1, 2 (|2> = RX(pi) then the 1->2 shelving pulse)
    cz g p n    pair group g (pairs target<-partner), targets and spectators in |1>,
                partners in |p> (p = 0, 1, +), then n CZs per pair; "+" partners are rotated
                back with RX(-pi/2), so ideally they read n mod 2
    idle g      as cz g 0 0 with a DELAY of 8 CZ lengths instead
Groups (pairs act in parallel): A = 86<-77, 104<-95 (spectator 68);  C = 86<-95, 68<-77
(spectator 104).  The CZ orientations are those quilc emits for the repetition-code programs.
--readout noshelve redefines MEASURE without the 1->2 pulse, so |2> can be told from |1>.
Usage on QCS JupyterLab:  python leak_probe.py translate|run --readout default|noshelve
"""
import argparse
import json
import os
import sys

import numpy as np
from pyquil import Program, get_qc

sys.path.insert(0, os.path.expanduser("~/repmid/qcs_quil"))
from run_repmid_qcs import measure_defcals_without_shelving, shelving_pulses  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("mode", choices=("translate", "run"))
ap.add_argument("--readout", choices=("default", "noshelve"), default="noshelve")
ap.add_argument("--K", type=int, default=50, help="shots per program")
ap.add_argument("--reps", type=int, default=4, help="runs per program")
ap.add_argument("--gap", type=float, default=2e-4)
ap.add_argument("--idle", type=float, default=8e-7)
ap.add_argument("--ns", default="0,1,2,4,8")
ap.add_argument("--tag", default="")
args = ap.parse_args()

Q = [104, 95, 86, 77, 68]
GROUPS = {"A": ([(86, 77), (104, 95)], [68]), "C": ([(86, 95), (68, 77)], [104])}
CZ = {frozenset(p): f"CZ {a} {b}" for p, (a, b) in
      {(86, 77): (77, 86), (68, 77): (77, 68), (104, 95): (95, 104), (86, 95): (95, 86)}.items()}

qc = get_qc("Cepheus-1-108Q", compiler_timeout=300, execution_timeout=300)
f12 = shelving_pulses(qc)
missing = [q for q in Q if q not in f12]
if missing:
    print(f"no 1->2 pulse for {missing}: their 'ref 2' is |1>", flush=True)


def shot(case):
    """Native Quil for one shot (without MEASURE) and the ideal outcome per qubit."""
    kind = case[0]
    if kind == "ref":
        s = case[1]
        body = [f"RX(pi) {q}" for q in Q] if s else []
        body += [f12[q] for q in Q if q in f12] if s == 2 else []
        return body, {q: min(s, 1) for q in Q}
    g, p, n = case[1], case[2], case[3]
    pairs, spect = GROUPS[g]
    ideal = {q: 0 for q in Q}
    body = [f"RX(pi) {t}" for t, _ in pairs] + [f"RX(pi) {q}" for q in spect]
    ideal.update({t: 1 for t, _ in pairs})
    ideal.update({q: 1 for q in spect})
    for _, b in pairs:
        if p == "1":
            body.append(f"RX(pi) {b}")
        elif p == "+":
            body.append(f"RX(pi/2) {b}")
        ideal[b] = {"0": 0, "1": 1, "+": n % 2}[p]
    if kind == "idle":
        body += [f"DELAY {q} {args.idle:g}" for q in Q]
    else:
        for _ in range(n):
            body += [CZ[frozenset(pr)] for pr in pairs]
    if p == "+":
        body += [f"RX(-pi/2) {b}" for _, b in pairs]
    return body, ideal


cases = [("ref", s) for s in (0, 1, 2)]
for g in GROUPS:
    cases.append(("idle", g, "0", 0))
    cases += [("cz", g, p, n) for p in ("0", "1", "+") for n in map(int, args.ns.split(","))]

head = measure_defcals_without_shelving(qc, Q) if args.readout == "noshelve" else []
progs = []
for case in cases:
    body, ideal = shot(case)
    lines = head + [f"DECLARE ro BIT[{len(Q) * args.K}]"]
    for k in range(args.K):
        lines += body + [f"MEASURE {q} ro[{len(Q) * k + i}]" for i, q in enumerate(Q)]
        lines += [f"DELAY {q} {args.gap:g}" for q in Q]
    p = Program("\n".join(lines) + "\n")
    p.wrap_in_numshots_loop(1)
    progs.append((case, ideal, qc.compiler.native_quil_to_executable(p)))
print(f"translated {len(progs)} programs, readout {args.readout}, K {args.K}", flush=True)
if args.mode != "run":
    sys.exit(0)

out, total = [], 0.0
for case, ideal, exe in progs:
    streams, iq, us = {}, {}, 0.0
    for _ in range(args.reps):
        r = qc.run(exe)
        us += float(getattr(r, "execution_duration_microseconds", 0) or 0)
        qd = r.data.result_data.to_qpu()
        for k, s in qd.readout_values.items():
            v = np.asarray(s.inner()).ravel()
            if k.endswith("_classified"):
                streams.setdefault(k, []).extend(int(x) for x in v)
            elif k.endswith("_unclassified"):
                iq.setdefault(k, []).extend([round(float(z.real), 4), round(float(z.imag), 4)] for z in v)
    total += us
    err = {q: float(np.mean(np.array(streams[f"q{q}_classified"]) != ideal[q])) for q in Q}
    out.append({"case": list(case), "ideal": ideal, "execution_us": us, "streams": streams, "iq": iq})
    print(f"  {' '.join(map(str, case)):14s} " + "  ".join(f"q{q} {err[q]:.2f}" for q in Q), flush=True)
    json.dump({"qubits": Q, "readout": args.readout, "K": args.K, "reps": args.reps, "gap": args.gap,
               "groups": {g: [list(map(list, v[0])), v[1]] for g, v in GROUPS.items()},
               "cz": {" ".join(map(str, sorted(k))): v for k, v in CZ.items()},
               "execution_us": total, "jobs": out},
              open(os.path.expanduser(f"~/repmid/leak_probe_{args.readout}{args.tag}.json"), "w"))
print(f"total execution {total / 1e3:.1f} ms")
