"""Which part of the one-round bit-flip program (logical |1>) makes data qubit 68 lose |1>?
Runs the device-native program from qcs_quil_repend_B/native/rep_bit_d3_r1_s1.quil and variants
with parts removed, without recompiling, twice in interleaved order.  Data qubits 104, 86, 68
ideally read 1 in every variant; each ancilla ideally reads the parity of the data qubits it
still has a CZ with (0 if it has none).
Usage on QCS JupyterLab:  python ablate_r1.py translate|run [--shots 1000] [--passes 2]
"""
import argparse
import json
import os

import numpy as np
from pyquil import Program, get_qc

ap = argparse.ArgumentParser()
ap.add_argument("mode", choices=("translate", "run"))
ap.add_argument("--native", default=os.path.expanduser("~/repmid/qcs_quil_repend_B/native/rep_bit_d3_r1_s1.quil"))
ap.add_argument("--shots", type=int, default=1000)
ap.add_argument("--passes", type=int, default=2)
ap.add_argument("--qvm", action="store_true", help="run on a noiseless QVM to check the ideal outcomes")
args = ap.parse_args()

src = [ln for ln in open(args.native).read().splitlines() if ln.strip() and not ln.startswith("PRAGMA")]
head = [ln for ln in src if ln.startswith("DECLARE")]
meas = [ln for ln in src if ln.startswith("MEASURE")]
body = [ln for ln in src if ln not in head and ln not in meas]
DATA, ANC = (104, 86, 68), (95, 77)
is_cz = lambda ln: ln.startswith("CZ")
on = lambda ln, q: int(ln.split()[-1]) == q and not is_cz(ln)    # single-qubit gate on q


def drop(pred):
    return [ln for ln in body if not pred(ln)]


V = {
    "full": body,
    "no CZ 77-86": drop(lambda l: l == "CZ 77 86"),
    "no CZ 77-68": drop(lambda l: l == "CZ 77 68"),
    "no CZ 95-*": drop(lambda l: l in ("CZ 95 104", "CZ 95 86")),
    "only CNOT 68->77": drop(lambda l: l in ("CZ 77 86", "CZ 95 104", "CZ 95 86") or on(l, 95)),
    "no ancilla gates": drop(lambda l: is_cz(l) or on(l, 95) or on(l, 77)),
    "68 alone": [l for l in body if on(l, 68)],
    "full, 68 prep RX(pi)": ["RX(pi) 68"] + drop(lambda l: on(l, 68)),
}


def ideal(lines):
    cz = [tuple(map(int, l.split()[1:])) for l in lines if is_cz(l)]
    out = {q: int(any(on(l, q) for l in lines)) for q in DATA}
    for a in ANC:
        out[a] = sum(out[d] for p in cz if a in p for d in p if d != a) % 2
    return out


qc = get_qc("Cepheus-1-108Q", as_qvm=args.qvm, compiler_timeout=300, execution_timeout=300)
exes = {}
for name, lines in V.items():
    p = Program("\n".join(head + lines + meas) + "\n")
    p.wrap_in_numshots_loop(args.shots)
    exes[name] = p if args.qvm else qc.compiler.native_quil_to_executable(p)
    print(f"translated {name:22s} {len(lines)} gates  ideal {ideal(lines)}", flush=True)
if args.mode != "run":
    raise SystemExit(0)

order = [m.split()[1] for m in meas]
res, total = [], 0.0
print("variant                 pass  " + "  ".join(f"q{q}" for q in order))
for k in range(args.passes):
    for name, exe in (exes.items() if k % 2 == 0 else reversed(list(exes.items()))):
        r = qc.run(exe)
        us = float(getattr(r, "execution_duration_microseconds", 0) or 0)
        total += us
        bits = np.asarray(r.get_register_map()["ro"])
        idl = ideal(V[name])
        err = {int(q): float(np.mean(bits[:, i] != idl[int(q)])) for i, q in enumerate(order)}
        res.append({"variant": name, "pass": k, "lines": V[name], "ideal": idl, "error": err,
                    "execution_us": us, "bitstrings": ["".join(map(str, b)) for b in bits]})
        print(f"{name:22s}  {k}     " + "  ".join(f"{err[int(q)]:.2f}" for q in order), flush=True)
        json.dump({"native": args.native, "order": order, "shots": args.shots, "execution_us": total,
                   "runs": res}, open(os.path.expanduser(f"~/repmid/ablate_r1{'_qvm' if args.qvm else ''}.json"), "w"))
print(f"total execution {total / 1e3:.1f} ms")
