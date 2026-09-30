"""
Parity checks inside the QLB circuits: the even-parity velocity circuits of the symmetry-verification
test, with the conserved spinor parity q0 xor q1 copied onto a fresh ancilla after every QLB step except
the last, and all qubits measured at the end (no mid-circuit measurement).

In the characteristic frame one step is collision + streaming, both of which conserve q0 xor q1; the
rotations R^-1 R between steps cancel, so the circuit is prep, R^-1, t x (Q_char, stream), measure.
A check computes the parity onto q0 (CX q1->q0), copies it to the ancilla (CX q0->a), and restores q0
(CX q1->q0), so the ancilla needs to be coupled to q0 only.  Ideally every ancilla reads 0.

Writes <outdir>/<name>.quil + manifest.json (kind raw).  Needs qiskit (not pyquil); free.

Usage:
    PYTHONPATH=. python -m qlb_port.export_parity_checks --outdir qlb_port/qcs_quil_n2_pcheck
"""

import argparse
import json
import os

import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit.library import StatePreparation, UnitaryGate
from qiskit.quantum_info import Statevector

from . import operators as ops
from . import run_hardware_zitterbewegung as H
from . import streaming as st
from .export_quil import circuit_to_quil


def check_circuit(npos, psi0, m, t, checks):
    n_anc = max(t - 1, 0) if checks else 0
    n = 2 + npos + n_anc
    qc = QuantumCircuit(n, n)
    qc.append(StatePreparation(psi0), range(2 + npos))
    qc.append(UnitaryGate(H._RINV, label="Rinv"), [0, 1])
    stream = st.streaming_circuit("x", npos, method="fourier")
    Q = ops.collision_operator_char("x", m, 0.0)
    for k in range(1, t + 1):
        qc.append(UnitaryGate(Q, label="Qchar"), [0, 1])
        qc.compose(stream, qubits=list(range(2 + npos)), inplace=True)
        if checks and k < t:
            a = 2 + npos + k - 1
            qc.cx(1, 0)
            qc.cx(0, a)
            qc.cx(1, 0)
    qc.measure(range(n), range(n))
    return qc


def ideal(qc, npos):
    """Exact P(q1=1)-P(q1=0), P(even parity), P(all ancillas 0) of a measured circuit."""
    sv = Statevector(qc.remove_final_measurements(inplace=False))
    p = sv.probabilities()
    idx = np.arange(len(p))
    b = lambda k: (idx >> k) & 1
    anc_zero = np.ones(len(p), bool)
    for k in range(2 + npos, qc.num_qubits):
        anc_zero &= b(k) == 0
    return (float(p @ (2 * b(1) - 1)), float(p[(b(0) ^ b(1)) == 0].sum()), float(p[anc_zero].sum()))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--npos", type=int, default=2)
    ap.add_argument("--mass", type=float, default=0.8)
    ap.add_argument("--k0", type=float, default=0.0)
    ap.add_argument("--tmax", type=int, default=4)
    ap.add_argument("--rounds", type=int, default=2)
    ap.add_argument("--shots", type=int, default=2000)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    npos, m = args.npos, args.mass
    k0 = H.nearest_lattice_k(npos, args.k0)
    psi, E = H.mode_state(npos, k0, m)
    psi = H.even_sector(psi, npos)
    traj = H.classical_evolution(psi, npos, m, args.tmax)

    progs = []
    for t in range(1, args.tmax + 1):
        for checks in ((False, True) if t > 1 else (False,)):
            qc = check_circuit(npos, psi, m, t, checks)
            quil, measured = circuit_to_quil(qc)
            assert measured == list(range(qc.num_qubits))
            v, peven, pacc = ideal(qc, npos)
            exact = H.alpha_x_expectation(traj[t], npos)
            assert abs(v - exact) < 1e-9 and abs(peven - 1) < 1e-9 and abs(pacc - 1) < 1e-9
            name = f"pc_t{t}_{'chk' if checks else 'none'}"
            with open(os.path.join(args.outdir, name + ".quil"), "w") as fh:
                fh.write(quil)
            n_cx = sum(1 for l in quil.splitlines() if l.startswith("CNOT"))
            progs.append({"base": name, "quil": name + ".quil", "t": t, "checks": int(checks),
                          "n_anc": qc.num_qubits - 2 - npos, "exact": exact, "cnot": n_cx})
            print(f"{name}: {qc.num_qubits} qubits, {n_cx} CNOT, exact <alpha_x> {exact:+.4f}, "
                  f"ideal ancillas all 0: {pacc:.6f}")
    jobs = []
    for r in range(args.rounds):                      # interleaved rounds for run-to-run spread
        for p in progs:
            jobs.append({**p, "name": f"{p['base']}_r{r}", "round": r, "kind": "raw",
                         "role": "pcheck", "shots": args.shots})
    manifest = {"description": "QLB even-parity velocity circuits with in-circuit parity checks",
                "npos": npos, "mass": m, "k0": k0, "E": float(E), "jobs": jobs,
                "qubit_order": "0 q0, 1 q1, 2 p0, 3 p1, 4.. check ancillas (step 1, 2, ...)"}
    with open(os.path.join(args.outdir, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"wrote {len(jobs)} jobs ({len(progs)} programs x {args.rounds} rounds) to {args.outdir}")


if __name__ == "__main__":
    main()
