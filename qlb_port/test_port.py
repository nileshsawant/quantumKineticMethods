"""
Validation for the QLB -> circuit porting layer.
================================================

Checks that every ported operator reproduces its classical target on the
qiskit-aer (GPU) backend, and that the porting-layer operators stay consistent
with the reference solver (``dirac_qlb_solver``).

Run from the repository root:

    module load qiskit/aer-gpu
    python -m qlb_port.test_port
"""

import os
import sys
import numpy as np
from qiskit import QuantumCircuit

# Allow running as a bare script (python qlb_port/test_port.py) as well as -m.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from qlb_port import operators as ops
from qlb_port import backend as bk
from qlb_port import port
from qlb_port import streaming

_failures = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        _failures.append(name)
    return ok


def test_consistency_with_solver():
    """Porting-layer operators must equal the validated solver operators."""
    print("Test A: consistency with dirac_qlb_solver (single source of truth)")
    try:
        import dirac_qlb_solver as s
    except Exception as e:
        check("import dirac_qlb_solver", False, f"{type(e).__name__}: {e}")
        return
    check("X_ROTATION == solver X_MATRIX", np.allclose(ops.X_ROTATION, s.X_MATRIX))
    check("Z_ROTATION == solver Z_MATRIX", np.allclose(ops.Z_ROTATION, s.Z_MATRIX))
    check("ALPHA_Y matches solver", np.allclose(ops.ALPHA_Y, s.ALPHA_Y))
    for ax in ("x", "y", "z"):
        check(f"streaming signs {ax} match SWEEP_CONFIG",
              list(ops.streaming_signs(ax)) == list(s.SWEEP_CONFIG[ax]["stream"]))
    # collision operator equals a_hat I - i b_hat ALPHA_Y from the solver convention
    m, g = 0.3, 0.2
    a, b = ops.collision_coefficients(m, g)
    check("collision unitary |a|^2+|b|^2 = 1", abs(abs(a) ** 2 + abs(b) ** 2 - 1) < 1e-12)


def test_rotations():
    print("Test B: rotation circuits reproduce X and Z (fidelity = 1)")
    for ax in ("x", "z"):
        r = port.verify(ops.ROTATIONS[ax])
        check(f"rotation {ax}", r["ok"], f"fid {r['fidelity']:.12f}, cx {r['cx']}, depth {r['depth']}")


def test_collision():
    print("Test C: collision circuits reproduce Q_hat over a (m,g) grid (fidelity = 1)")
    worst = 1.0
    for m in (0.0, 0.1, 0.3, 0.6):
        for g in (0.0, 0.2, 0.5):
            Q = ops.collision_operator(m, g)
            # unitarity of the target itself
            uok = np.allclose(Q.conj().T @ Q, np.eye(4))
            r = port.verify(Q)
            worst = min(worst, r["fidelity"])
            if not (uok and r["ok"]):
                check(f"collision m={m} g={g}", False, f"unitary={uok}, fid {r['fidelity']:.2e}")
    check("all (m,g) collisions ported", worst > 1 - 1e-9, f"worst fidelity {worst:.12f}")


def test_statevector_application():
    print("Test D: GPU statevector application matches classical operator")
    rng = np.random.default_rng(0)
    psi = rng.standard_normal(4) + 1j * rng.standard_normal(4)
    psi /= np.linalg.norm(psi)
    for label, U in (("X", ops.X_ROTATION),
                     ("Z", ops.Z_ROTATION),
                     ("Qhat(0.3,0.2)", ops.collision_operator(0.3, 0.2))):
        qc = port.port_unitary(U, label=label)
        out = bk.apply_to_statevector(qc, psi)
        fid = bk.state_fidelity(out, U @ psi)
        check(f"{label} |circuit psi> == U|psi>", abs(fid - 1.0) < 1e-9, f"fid {fid:.12f}")


def test_composed_sweep_spinor():
    print("Test E: composed rotate->collide->rotate-back matches classical (spinor part)")
    # One sweep's SPINOR-space action (no streaming): R (a I - i b R^-1 ALPHA_Y R) R^-1
    #   = a I - i b ALPHA_Y = Q_hat  (standard frame).  Verify the composed circuit.
    m, g = 0.4, 0.1
    for ax in ("x", "z"):
        R = ops.ROTATIONS[ax]
        Qc = ops.collision_operator_char(ax, m, g)          # collision in char frame
        composed = R @ Qc @ R.conj().T                       # back to standard frame
        target = ops.collision_operator(m, g)                # should equal Q_hat
        check(f"{ax}: R Qc R^-1 == Q_hat (math)", np.allclose(composed, target))
        r = port.verify(composed)
        check(f"{ax}: composed circuit fidelity", r["ok"], f"fid {r['fidelity']:.12f}")


def test_increment():
    print("Test F: increment circuit adds 1 (mod 2^n)")
    for n in (2, 3, 4):
        qc = QuantumCircuit(n)
        qc.append(streaming.increment_gate(n), list(range(n)))
        U = bk.circuit_unitary(qc)
        N = 2 ** n
        P = np.zeros((N, N), dtype=complex)
        for x in range(N):
            P[(x + 1) % N, x] = 1.0
        check(f"increment n={n}", np.allclose(U, P))


def test_streaming():
    print("Test G: streaming circuits reproduce the +/-1 shift permutation (fidelity = 1)")
    for axis in ("x", "y", "z"):
        for n_pos in (2, 3):
            qc = streaming.streaming_circuit(axis, n_pos)
            U = bk.circuit_unitary(qc)
            P = ops.streaming_reference(axis, n_pos)
            fid = bk.gate_fidelity(P, U)
            check(f"streaming {axis} n_pos={n_pos}", np.allclose(U, P),
                  f"fid {fid:.12f}")


def test_streaming_statevector():
    print("Test H: streaming applied to a localized packet moves it by one site (GPU)")
    axis, n_pos = "x", 4
    N = 2 ** n_pos
    qc = streaming.streaming_circuit(axis, n_pos)
    # +x mover (spinor component 2, qubit1=1) localized at site x0 -> should move to x0+1
    x0, c = 5, 2
    psi = np.zeros(4 * N, dtype=complex)
    psi[x0 * 4 + c] = 1.0
    out = bk.apply_to_statevector(qc, psi)
    expect = np.zeros_like(psi); expect[((x0 + 1) % N) * 4 + c] = 1.0
    check("+x component moves +1 site", bk.state_fidelity(out, expect) > 1 - 1e-9)
    # -x mover (component 0, qubit1=0) should move to x0-1
    c = 0
    psi = np.zeros(4 * N, dtype=complex); psi[x0 * 4 + c] = 1.0
    out = bk.apply_to_statevector(qc, psi)
    expect = np.zeros_like(psi); expect[((x0 - 1) % N) * 4 + c] = 1.0
    check("-x component moves -1 site", bk.state_fidelity(out, expect) > 1 - 1e-9)


if __name__ == "__main__":
    print("=" * 70)
    print("QLB -> circuit porting validation   (device: %s)" % bk.device_report())
    print("=" * 70)
    test_consistency_with_solver()
    test_rotations()
    test_collision()
    test_statevector_application()
    test_composed_sweep_spinor()
    test_increment()
    test_streaming()
    test_streaming_statevector()
    print("=" * 70)
    if _failures:
        print(f"FAILED ({len(_failures)}): " + ", ".join(_failures))
        raise SystemExit(1)
    print("All porting validation tests passed.")
