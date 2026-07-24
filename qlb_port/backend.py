"""
qiskit-aer backend helpers (GPU-first) for the QLB porting layer.
=================================================================

Two Aer simulators are used:

  * a ``unitary`` simulator to extract the full unitary a circuit implements
    (used to verify a ported gate against its target matrix), and
  * a ``statevector`` simulator to apply a circuit to an input spinor/state
    (used for state-level checks and, later, full circuit emulation).

Both prefer the GPU (H100) and fall back to CPU if a GPU device is unavailable,
so the module also runs on a login node without a GPU.
"""

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator

_UNITARY_SIM = None
_STATEVECTOR_SIM = None


def _make(method):
    """Build an AerSimulator for `method`, preferring GPU, falling back to CPU."""
    try:
        sim = AerSimulator(method=method, device="GPU")
        probe = QuantumCircuit(1)
        probe.h(0)
        if method == "unitary":
            probe.save_unitary()
        else:
            probe.save_statevector()
        sim.run(probe).result()          # raises if GPU method unsupported
        return sim, "GPU"
    except Exception:
        return AerSimulator(method=method, device="CPU"), "CPU"


def unitary_simulator():
    """Cached unitary-method AerSimulator (GPU if available)."""
    global _UNITARY_SIM
    if _UNITARY_SIM is None:
        _UNITARY_SIM, _ = _make("unitary")
    return _UNITARY_SIM


def statevector_simulator():
    """Cached statevector-method AerSimulator (GPU if available)."""
    global _STATEVECTOR_SIM
    if _STATEVECTOR_SIM is None:
        _STATEVECTOR_SIM, _ = _make("statevector")
    return _STATEVECTOR_SIM


def device_report():
    """Return which device each simulator resolved to (for diagnostics)."""
    _, ud = _make("unitary")
    _, sd = _make("statevector")
    return {"unitary": ud, "statevector": sd}


def circuit_unitary(qc):
    """Return the 2^n x 2^n unitary implemented by circuit `qc` (n = qc.num_qubits)."""
    sim = unitary_simulator()
    work = qc.copy()
    work.save_unitary()
    res = sim.run(transpile(work, sim)).result()
    return np.asarray(res.get_unitary())


def apply_to_statevector(qc, psi_in):
    """Apply circuit `qc` to input statevector `psi_in`; return the output statevector."""
    sim = statevector_simulator()
    n = qc.num_qubits
    psi_in = np.asarray(psi_in, dtype=complex).reshape(-1)
    work = QuantumCircuit(n)
    work.initialize(psi_in / np.linalg.norm(psi_in), list(range(n)))
    work.compose(qc, inplace=True)
    work.save_statevector()
    res = sim.run(transpile(work, sim)).result()
    return np.asarray(res.get_statevector())


def gate_fidelity(A, B):
    """Phase-invariant gate fidelity |tr(A^dag B)| / dim for two same-size unitaries."""
    A = np.asarray(A); B = np.asarray(B)
    return float(abs(np.trace(A.conj().T @ B)) / A.shape[0])


def state_fidelity(a, b):
    """|<a|b>| for two (unnormalized) state vectors."""
    a = np.asarray(a).reshape(-1); b = np.asarray(b).reshape(-1)
    return float(abs(np.vdot(a, b)) / (np.linalg.norm(a) * np.linalg.norm(b)))


def evolve_snapshots(step_circuit, psi0, snapshots):
    """
    Apply a single-step circuit repeatedly, saving the statevector at each snapshot time.

    Efficient for time evolution: the full circuit (initialize + repeated steps with
    labelled statevector saves) is transpiled once and simulated in a single run, rather
    than rebuilding and re-transpiling a fresh circuit per snapshot.

    Parameters
    ----------
    step_circuit : QuantumCircuit for one time step (n qubits).
    psi0         : initial statevector (length 2**n).
    snapshots    : iterable of step indices at which to record the statevector.

    Returns
    -------
    dict {t: statevector} for each t in snapshots.
    """
    sim = statevector_simulator()
    n = step_circuit.num_qubits
    psi0 = np.asarray(psi0, dtype=complex).reshape(-1)
    snaps = sorted(set(int(t) for t in snapshots))

    qc = QuantumCircuit(n)
    qc.initialize(psi0 / np.linalg.norm(psi0), list(range(n)))
    if 0 in snaps:
        qc.save_statevector(label="t0")
    for t in range(1, max(snaps) + 1):
        qc.compose(step_circuit, inplace=True)
        if t in snaps:
            qc.save_statevector(label=f"t{t}")
    data = sim.run(transpile(qc, sim)).result().data()
    return {t: np.asarray(data[f"t{t}"]) for t in snaps}
