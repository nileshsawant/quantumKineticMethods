"""
Diagnostic programs for the failing in-circuit parity checks (Section 5.5): where do the 17-25%
ancilla errors come from?

On the same layout as export_parity_checks.py (logical 0 q0, 1 q1, 2 p0, 3 p1, 4-6 check ancillas):
  rest      nothing applied, all seven qubits measured (readout of |0>)
  anc_cal   X on the three ancillas (readout of |1>)
  idle_t4   the t = 4 QLB circuit without checks, ancillas idle but measured (spectator effects)
  chk_00    three parity checks on q0 q1 = |00>, no QLB gates
  chk_11    three parity checks on q0 q1 = |11>
  chk_bell  three parity checks on (|00> + |11>)/sqrt2 of q0 q1
  chk_t4    the t = 4 QLB circuit with checks (as in Section 5.5), for a same-batch reference
All ancillas ideally read 0 except in anc_cal.  Free; needs qiskit.

Usage:
    PYTHONPATH=. python -m qlb_port.export_check_diagnostics --outdir qlb_port/qcs_quil_n2_pdiag
"""

import argparse
import json
import os

from qiskit import QuantumCircuit

from . import run_hardware_zitterbewegung as H
from .export_parity_checks import check_circuit
from .export_quil import circuit_to_quil

N = 7


def checks(qc):
    for a in (4, 5, 6):
        qc.cx(1, 0)
        qc.cx(0, a)
        qc.cx(1, 0)


def cz_repetition_programs(outdir, shots, nmax=4):
    """n repetitions (n = 0..nmax) of CNOT(q0 -> a) = H_a CZ H_a on q0 = |0>, for the ancillas 4 and 5,
    written in native gates inside a PRESERVE_BLOCK so that nothing cancels.  Ideally the ancillas stay
    in |0>.  A coherent phase error theta per gate gives P(1) = sin^2(n theta / 2); an incoherent flip
    probability p gives (1 - (1 - 2p)^n) / 2."""
    h = lambda q: [f"RZ(pi/2) {q}", f"RX(pi/2) {q}", f"RZ(pi/2) {q}"]
    jobs = []
    for n in range(nmax + 1):
        body = []
        for a in (4, 5):
            for _ in range(n):
                body += h(a) + [f"CZ 0 {a}"] + h(a)
        lines = ["DECLARE ro BIT[3]"]
        if body:
            lines += ["PRAGMA PRESERVE_BLOCK"] + body + ["PRAGMA END_PRESERVE_BLOCK"]
        lines += ["MEASURE 0 ro[0]", "MEASURE 4 ro[1]", "MEASURE 5 ro[2]"]
        name = f"czrep_n{n}"
        with open(os.path.join(outdir, name + ".quil"), "w") as fh:
            fh.write("\n".join(lines) + "\n")
        jobs.append({"name": name, "quil": name + ".quil", "kind": "raw", "role": "czrep", "n": n,
                     "shots": shots})
    controls = {"ctl_h_only": lambda a: h(a) + h(a),              # the two Hadamards without the CZ
                "ctl_cz_only": lambda a: [f"CZ 0 {a}"]}           # the CZ alone, ancilla left in |0>
    for name, gates in controls.items():
        body = [g for a in (4, 5) for g in gates(a)]
        lines = (["DECLARE ro BIT[3]", "PRAGMA PRESERVE_BLOCK"] + body + ["PRAGMA END_PRESERVE_BLOCK",
                 "MEASURE 0 ro[0]", "MEASURE 4 ro[1]", "MEASURE 5 ro[2]"])
        with open(os.path.join(outdir, name + ".quil"), "w") as fh:
            fh.write("\n".join(lines) + "\n")
        jobs.append({"name": name, "quil": name + ".quil", "kind": "raw", "role": "czrep", "n": -1,
                     "shots": shots})
    return jobs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--shots", type=int, default=2000)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--cz-repetition", action="store_true",
                    help="write the CNOT-repetition programs instead (coherent vs incoherent test)")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    if args.cz_repetition:
        jobs = cz_repetition_programs(args.outdir, args.shots)
        with open(os.path.join(args.outdir, "manifest.json"), "w") as fh:
            json.dump({"description": "repeated CNOT(q0 -> ancilla) on q0 = |0>", "jobs": jobs}, fh, indent=2)
        print(f"wrote {len(jobs)} programs to {args.outdir}")
        return
    psi = H.even_sector(H.mode_state(2, 0.0, 0.8)[0], 2)

    progs = {}
    qc = QuantumCircuit(N, N)
    progs["rest"] = qc
    qc = QuantumCircuit(N, N)
    qc.x([4, 5, 6])
    progs["anc_cal"] = qc
    base = check_circuit(2, psi, 0.8, 4, False).remove_final_measurements(inplace=False)
    qc = QuantumCircuit(N, N)
    qc.compose(base, qubits=range(4), inplace=True)
    progs["idle_t4"] = qc
    for name, prep in (("chk_00", []), ("chk_11", ["x0", "x1"]), ("chk_bell", ["h0", "cx01"])):
        qc = QuantumCircuit(N, N)
        for g in prep:
            if g == "x0":
                qc.x(0)
            elif g == "x1":
                qc.x(1)
            elif g == "h0":
                qc.h(0)
            elif g == "cx01":
                qc.cx(0, 1)
        checks(qc)
        progs[name] = qc
    progs["chk_t4"] = check_circuit(2, psi, 0.8, 4, True).remove_final_measurements(inplace=False)

    jobs = []
    for name, qc in progs.items():
        if qc.num_clbits == 0:
            full = QuantumCircuit(N, N)
            full.compose(qc, qubits=range(N), inplace=True)
            qc = full
        qc.measure(range(N), range(N))
        quil, measured = circuit_to_quil(qc)
        assert measured == list(range(N))
        with open(os.path.join(args.outdir, name + ".quil"), "w") as fh:
            fh.write(quil)
        n_cx = sum(1 for l in quil.splitlines() if l.startswith("CNOT"))
        jobs.append({"name": name, "quil": name + ".quil", "kind": "raw", "role": "pdiag",
                     "shots": args.shots, "anc_ideal": 1 if name == "anc_cal" else 0})
        print(f"{name}: {n_cx} CNOT")
    with open(os.path.join(args.outdir, "manifest.json"), "w") as fh:
        json.dump({"description": "diagnostics of the in-circuit parity checks", "jobs": jobs}, fh, indent=2)
    print(f"wrote {len(jobs)} programs to {args.outdir}")


if __name__ == "__main__":
    main()
