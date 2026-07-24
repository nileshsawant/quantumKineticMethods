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
from qlb_port import sweep
from qlb_port import potential

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


def test_sweep_assembly():
    print("Test I: full single-axis sweep circuit == classical sub-step operator")
    for axis in ("x", "y", "z"):
        for m in (0.0, 0.3):
            qc = sweep.sweep_circuit(axis, 3, m_tilde=m)
            U = bk.circuit_unitary(qc)
            M = sweep.sweep_operator(axis, 3, m_tilde=m)
            check(f"sweep {axis} m={m}", np.allclose(U, M), f"fid {bk.gate_fidelity(M, U):.12f}")


def test_free_particle_evolution():
    print("Test J: multi-step free-particle circuit matches classical solver operator (GPU)")
    axis, n_pos, T = "x", 4, 6
    N = 2 ** n_pos
    Op = sweep.sweep_operator(axis, n_pos, m_tilde=0.0)
    # physical +x mover spinor, Gaussian packet on the ring
    sp = ops.X_ROTATION @ (np.array([0, 0, 1, 1], dtype=complex) / np.sqrt(2))
    x = np.arange(N)
    env = np.exp(-((x - N // 2) ** 2) / (2 * 2.0 ** 2)) * np.exp(1j * 0.6 * x)
    psi0 = np.zeros(4 * N, dtype=complex)
    for xi in range(N):
        for c in range(4):
            psi0[xi * 4 + c] = env[xi] * sp[c]
    psi0 /= np.linalg.norm(psi0)
    pc = psi0.copy()
    for _ in range(T):
        pc = Op @ pc
    out = bk.apply_to_statevector(sweep.evolution_circuit(axis, n_pos, T), psi0)
    check(f"{T}-step circuit vs classical", bk.state_fidelity(out, pc) > 1 - 1e-9,
          f"fid {bk.state_fidelity(out, pc):.12f}")


def test_free_particle_propagation():
    print("Test K: one sweep moves a localized +x / -x packet by +1 / -1 site")
    axis, n_pos = "x", 4
    N = 2 ** n_pos
    x0 = 6
    step = sweep.sweep_circuit(axis, n_pos, m_tilde=0.0)
    for label, char, direction in (("+x", np.array([0, 0, 1, 1]), +1),
                                   ("-x", np.array([1, 1, 0, 0]), -1)):
        sp = ops.X_ROTATION @ (char.astype(complex) / np.sqrt(2))
        psi = np.zeros(4 * N, dtype=complex)
        for c in range(4):
            psi[x0 * 4 + c] = sp[c]
        out = bk.apply_to_statevector(step, psi)
        expect = np.zeros(4 * N, dtype=complex)
        for c in range(4):
            expect[((x0 + direction) % N) * 4 + c] = sp[c]
        check(f"{label} packet moves {direction:+d} site",
              bk.state_fidelity(out, expect) > 1 - 1e-9,
              f"fid {bk.state_fidelity(out, expect):.12f}")


def test_potential_oracle():
    print("Test L: massless sweep-with-potential circuit == classical operator")
    rng = np.random.default_rng(1)
    n_pos = 3
    for axis in ("x", "y", "z"):
        V = rng.uniform(0, 0.5, 2 ** n_pos)
        qc = potential.sweep_circuit_potential(axis, n_pos, V)
        U = bk.circuit_unitary(qc)
        M = potential.sweep_operator_potential(axis, n_pos, V)
        check(f"sweepV {axis}", np.allclose(U, M), f"fid {bk.gate_fidelity(M, U):.12f}")


def test_barrier_scattering():
    print("Test M: massless packet scattering off a potential barrier (circuit vs classical, GPU)")
    axis, n_pos, T = "x", 5, 8
    N = 2 ** n_pos
    # barrier of a few sites in the middle of the ring
    V = potential.impurity_field(n_pos, range(N // 2, N // 2 + 3), g_value=0.8)
    Op = potential.sweep_operator_potential(axis, n_pos, V)
    # +x packet approaching the barrier
    sp = ops.X_ROTATION @ (np.array([0, 0, 1, 1], dtype=complex) / np.sqrt(2))
    x = np.arange(N)
    env = np.exp(-((x - N // 4) ** 2) / (2 * 3.0 ** 2)) * np.exp(1j * 0.5 * x)
    psi0 = (env[:, None] * sp[None, :]).reshape(-1)
    psi0 /= np.linalg.norm(psi0)
    pc = psi0.copy()
    for _ in range(T):
        pc = Op @ pc
    out = bk.apply_to_statevector(potential.evolution_circuit_potential(axis, n_pos, T, V), psi0)
    check(f"{T}-step barrier scattering matches classical",
          bk.state_fidelity(out, pc) > 1 - 1e-9, f"fid {bk.state_fidelity(out, pc):.12f}")
    check("probability conserved", abs(np.linalg.norm(out) - 1) < 1e-9)


def test_massive_barrier():
    print("Test N: massive (gapped) packet reflects off a barrier (circuit vs classical, GPU)")
    axis, n_pos, T, m = "x", 5, 10, 0.5
    N = 2 ** n_pos
    V = potential.impurity_field(n_pos, range(N // 2, N // 2 + 3), g_value=0.8)
    # circuit == classical operator for the massive multiplexed sweep
    U = bk.circuit_unitary(potential.sweep_circuit_potential(axis, 3,
                           potential.impurity_field(3, [3, 4], 0.6), m_tilde=m))
    M = potential.sweep_operator_potential(axis, 3,
                           potential.impurity_field(3, [3, 4], 0.6), m_tilde=m)
    check("massive sweepV circuit == classical", np.allclose(U, M),
          f"fid {bk.gate_fidelity(M, U):.12f}")
    # scattering: circuit matches classical, and reflection is non-zero (barrier acts)
    sp = ops.X_ROTATION @ (np.array([0, 0, 1, 1], dtype=complex) / np.sqrt(2))
    x = np.arange(N)
    env = np.exp(-((x - N // 4) ** 2) / (2 * 3.0 ** 2)) * np.exp(1j * 0.5 * x)
    psi0 = (env[:, None] * sp[None, :]).reshape(-1)
    psi0 /= np.linalg.norm(psi0)
    Op = potential.sweep_operator_potential(axis, n_pos, V, m_tilde=m)
    pc = psi0.copy()
    for _ in range(T):
        pc = Op @ pc
    out = bk.apply_to_statevector(
        potential.evolution_circuit_potential(axis, n_pos, T, V, m_tilde=m), psi0)
    check("massive barrier: circuit matches classical",
          bk.state_fidelity(out, pc) > 1 - 1e-9, f"fid {bk.state_fidelity(out, pc):.12f}")
    reflected = (np.abs(out.reshape(N, 4)[:N // 2 - 2]) ** 2).sum()
    check("massive barrier reflects (non-zero back-flux)", reflected > 0.02,
          f"reflected fraction {reflected:.3f}")


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
    test_sweep_assembly()
    test_free_particle_evolution()
    test_free_particle_propagation()
    test_potential_oracle()
    test_barrier_scattering()
    test_massive_barrier()
    print("=" * 70)
    if _failures:
        print(f"FAILED ({len(_failures)}): " + ", ".join(_failures))
        raise SystemExit(1)
    print("All porting validation tests passed.")
