"""Does a gate work after a mid-circuit measurement, depending on the prepared state and on how the
qubit is returned from |2> (Cepheus-1)?  Per shot, on all listed qubits in parallel:
    prepare x (alternates 0,1);  MEASURE (m1);  [variant];  RX(pi);  MEASURE (m2);  DELAY gap
variant  none:  DELAY t
         early: 1->2 pulse right after MEASURE, then DELAY t
         late:  DELAY t, then the 1->2 pulse
         noshelve: DELAY t, with MEASURE redefined (DEFCAL) without the 1->2 pulse
Ideal m1 = x, m2 = 1 - x.  The 1->2 pulse is the one each MEASURE calibration plays before readout.
The raw IQ points are saved too, so a threshold can be refitted for the noshelve readout.
Usage on QCS JupyterLab:  python mcm_char2.py translate|run --variant none|early|late --qubits ...
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
ap.add_argument("--variant", choices=("none", "early", "late", "noshelve"), default="none")
ap.add_argument("--qubits", required=True)
ap.add_argument("--t", type=float, default=1e-6)
ap.add_argument("--gap", type=float, default=2e-4)
ap.add_argument("--shots", type=int, default=100)
ap.add_argument("--ro-duration", type=float, default=None, help="noshelve: readout length (s)")
ap.add_argument("--tag", default="")
args = ap.parse_args()

qc = get_qc("Cepheus-1-108Q", compiler_timeout=300, execution_timeout=300)
f12 = shelving_pulses(qc)
Q = [int(x) for x in args.qubits.split(",")]
nq, K = len(Q), args.shots
lines = [f"DECLARE ro BIT[{2 * nq * K}]"]
if args.variant == "noshelve":
    lines = measure_defcals_without_shelving(qc, Q, args.ro_duration) + lines
for k in range(K):
    base = 2 * nq * k
    if k % 2:
        lines += [f"RX(pi) {q}" for q in Q]
    lines += [f"MEASURE {q} ro[{base + i}]" for i, q in enumerate(Q)]
    if args.variant == "early":
        lines += [f12[q] for q in Q if q in f12]
    lines += [f"DELAY {q} {args.t:g}" for q in Q]
    if args.variant == "late":
        lines += [f12[q] for q in Q if q in f12]
    lines += [f"RX(pi) {q}" for q in Q]
    lines += [f"MEASURE {q} ro[{base + nq + i}]" for i, q in enumerate(Q)]
    lines += [f"DELAY {q} {args.gap:g}" for q in Q]
p = Program("\n".join(lines) + "\n")
p.wrap_in_numshots_loop(1)
exe = qc.compiler.native_quil_to_executable(p)
print(f"translated: variant {args.variant}, {nq} qubits, {K} shots", flush=True)
if args.mode != "run":
    sys.exit(0)

r = qc.run(exe)
us = float(getattr(r, "execution_duration_microseconds", 0) or 0)
qd = r.data.result_data.to_qpu()
streams = {k: [int(v) for v in np.asarray(s.inner()).ravel()]
           for k, s in qd.readout_values.items() if k.endswith("_classified")}
iq = {k: [[float(z.real), float(z.imag)] for z in np.asarray(s.inner()).ravel()]
      for k, s in qd.readout_values.items() if k.endswith("_unclassified")}
json.dump({"qubits": Q, "variant": args.variant, "t": args.t, "gap": args.gap, "execution_us": us,
           "f12": [q for q in Q if q in f12], "streams": streams, "iq": iq},
          open(os.path.expanduser(f"~/repmid/mcm_char2_{args.variant}{args.tag}.json"), "w"))
x = np.arange(K) % 2
print(f"{us / 1e3:.1f} ms  qubit  f12  P(m1=x)  gate worked: x=0  x=1")
w0, w1 = [], []
for q in Q:
    v = np.array(streams[f"q{q}_classified"])
    if len(v) != 2 * K:
        print(f"  q{q}: {len(v)} values, expected {2 * K}")
        continue
    v = v.reshape(K, 2)
    a, b0, b1 = np.mean(v[:, 0] == x), np.mean(v[x == 0, 1] == 1), np.mean(v[x == 1, 1] == 0)
    w0.append(b0)
    w1.append(b1)
    print(f"        {q:4d}   {'y' if q in f12 else 'n'}    {a:.2f}        {b0:.2f}  {b1:.2f}")
print(f"  mean gate worked: x=0 {np.mean(w0):.2f}   x=1 {np.mean(w1):.2f}")
