"""Mid-circuit readout probes, one program per case (50 measurements each), on Cepheus-1.
  A  q84  (readout without 1->2 shelving)    50 x [RX(pi); MEASURE]                 ideal 1010...
  B  q84                                     50 x [RX(pi); MEASURE; DELAY 2us]      ideal 1010...
  C  q100 (shelving readout), device default 50 x [RX(pi); MEASURE; DELAY 2us]      ideal 1010...
  D  q100 with 1->2 pulse repeated after MEASURE, + DELAY 2us                       ideal 1010...
  E  q100 device default 25 x [RX(pi); MEASURE; DELAY 2us; MEASURE; DELAY 2us]      ideal 11001100...
Raw IQ readout is kept.  Usage on QCS JupyterLab:  python mcm_probe2.py translate | run"""
import json
import os
import sys

import numpy as np
from pyquil import Program, get_qc

sys.path.insert(0, os.path.expanduser("~/repmid/qcs_quil"))
from run_repmid_qcs import add_unshelve, shelving_pulses  # noqa: E402
from repcode_mid_qcs import readout_stream_to_bits  # noqa: E402

K, DLY = 50, "2e-06"


def flips(q, delay):
    body = [f"RX(pi) {q}\nMEASURE {q} ro[{k}]" + (f"\nDELAY {q} {DLY}" if delay else "")
            for k in range(K)]
    return "\n".join([f"DECLARE ro BIT[{K}]"] + body) + "\n", [(k + 1) % 2 for k in range(K)]


def qnd(q):
    body = [f"RX(pi) {q}\nMEASURE {q} ro[{2 * k}]\nDELAY {q} {DLY}\nMEASURE {q} ro[{2 * k + 1}]\nDELAY {q} {DLY}"
            for k in range(K // 2)]
    return "\n".join([f"DECLARE ro BIT[{K}]"] + body) + "\n", [(k // 2 + 1) % 2 for k in range(K)]


qc = get_qc("Cepheus-1-108Q", compiler_timeout=300, execution_timeout=300)
f12 = shelving_pulses(qc)
print("shelving pulse on q84:", 84 in f12, " on q100:", 100 in f12)
cases = {"A q84 flips": (84, *flips(84, False)), "B q84 flips+delay": (84, *flips(84, True)),
         "C q100 flips+delay": (100, *flips(100, True))}
txt, ideal = flips(100, True)
cases["D q100 flips+delay+unshelve"] = (100, add_unshelve(txt, f12), ideal)
cases["E q100 measure twice"] = (100, *qnd(100))

exes = {}
for name, (q, text, ideal) in cases.items():
    p = Program(text)
    p.wrap_in_numshots_loop(1)
    try:
        exes[name] = qc.compiler.native_quil_to_executable(p)
        print("translated", name, flush=True)
    except Exception as e:
        print("FAILED", name, str(e)[:300], flush=True)
if sys.argv[1] != "run" or len(exes) != len(cases):
    sys.exit(0)

out = {}
for name, (q, text, ideal) in cases.items():
    r = qc.run(exes[name])
    qd = r.data.result_data.to_qpu()
    bits = readout_stream_to_bits(qd, K)[0]
    iq = np.asarray(qd.readout_values[f"q{q}_unclassified"].inner()).ravel()
    us = float(getattr(r, "execution_duration_microseconds", 0) or 0)
    agree = sum(int(b) == i for b, i in zip(bits, ideal))
    print(f"{name:30s} {''.join(map(str, bits))}  {agree}/{K}  {us / 1e3:.1f} ms", flush=True)
    out[name] = {"bits": [int(b) for b in bits], "ideal": ideal,
                 "iq": [[float(z.real), float(z.imag)] for z in iq]}
json.dump(out, open(os.path.expanduser("~/repmid/mcm_probe2.json"), "w"))
