"""Mid-circuit readout characterization on many qubits at once (Cepheus-1).

Each shot, on every listed qubit in parallel:  prepare x (x alternates 0,1 by shot); MEASURE (m1);
DELAY t; MEASURE (m2); DELAY t; RX(pi); MEASURE (m3); then DELAY gap (passive reset).
Ideal m1 = m2 = x, m3 = 1 - x.  Per qubit this gives the readout fidelity P(m1 = x), whether a
measurement leaves the qubit in the state it reported P(m2 = m1), and whether a gate right after a
measurement works P(m3 != m2).  m1 of the next shot also shows whether the gap reset the qubit.

Usage on QCS JupyterLab:  python mcm_char.py translate|run --qubits 100,101 [--t 1e-6] [--gap 2e-4]
                                                       [--shots 100] [--tag name]"""
import argparse
import json
import os
import sys

import numpy as np
from pyquil import Program, get_qc

sys.path.insert(0, os.path.expanduser("~/repmid/qcs_quil"))
from repcode_mid_qcs import readout_stream_to_bits  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("mode", choices=("translate", "run"))
ap.add_argument("--qubits", required=True)
ap.add_argument("--t", type=float, default=1e-6)
ap.add_argument("--gap", type=float, default=2e-4)
ap.add_argument("--shots", type=int, default=100)
ap.add_argument("--tag", default="char")
args = ap.parse_args()

Q = [int(x) for x in args.qubits.split(",")]
nq, K = len(Q), args.shots
lines = [f"DECLARE ro BIT[{3 * nq * K}]"]
for k in range(K):
    base = 3 * nq * k
    if k % 2:
        lines += [f"RX(pi) {q}" for q in Q]
    lines += [f"MEASURE {q} ro[{base + i}]" for i, q in enumerate(Q)]
    lines += [f"DELAY {q} {args.t:g}" for q in Q]
    lines += [f"MEASURE {q} ro[{base + nq + i}]" for i, q in enumerate(Q)]
    lines += [f"DELAY {q} {args.t:g}" for q in Q]
    lines += [f"RX(pi) {q}" for q in Q]
    lines += [f"MEASURE {q} ro[{base + 2 * nq + i}]" for i, q in enumerate(Q)]
    lines += [f"DELAY {q} {args.gap:g}" for q in Q]
p = Program("\n".join(lines) + "\n")
p.wrap_in_numshots_loop(1)
qc = get_qc("Cepheus-1-108Q", compiler_timeout=300, execution_timeout=300)
exe = qc.compiler.native_quil_to_executable(p)
print(f"translated: {nq} qubits, {K} shots, t {args.t:g} s, gap {args.gap:g} s", flush=True)
if args.mode != "run":
    sys.exit(0)

r = qc.run(exe)
us = float(getattr(r, "execution_duration_microseconds", 0) or 0)
qd = r.data.result_data.to_qpu()
streams = {k: [int(v) for v in np.asarray(s.inner()).ravel()]
           for k, s in qd.readout_values.items() if k.endswith("_classified")}
path = os.path.expanduser(f"~/repmid/mcm_char_{args.tag}.json")
json.dump({"qubits": Q, "t": args.t, "gap": args.gap, "execution_us": us,
           "mappings": dict(qd.mappings), "streams": streams}, open(path, "w"))
lens = {k: len(v) for k, v in streams.items()}
if set(lens.values()) != {3 * K}:
    sys.exit(f"{us / 1e3:.1f} ms; expected {3 * K} values per qubit, got {lens}")
bits = readout_stream_to_bits(qd, 3 * nq * K).reshape(K, 3, nq)
x = np.arange(K) % 2
m1, m2, m3 = bits[:, 0], bits[:, 1], bits[:, 2]
print(f"{us / 1e3:.1f} ms   qubit  P(m1=x)  P(m2=m1)  P(m3!=m2)  P(m1=1|x=0)  P(m1=0|x=1)")
for i, q in enumerate(Q):
    a = np.mean(m1[:, i] == x)
    b = np.mean(m2[:, i] == m1[:, i])
    c = np.mean(m3[:, i] != m2[:, i])
    e01 = np.mean(m1[x == 0, i] == 1)
    e10 = np.mean(m1[x == 1, i] == 0)
    print(f"          {q:4d}   {a:.2f}     {b:.2f}      {c:.2f}       {e01:.2f}         {e10:.2f}")
