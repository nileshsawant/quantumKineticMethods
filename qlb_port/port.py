"""
General "unitary -> verified circuit" porting harness, plus QLB-specific builders.
==================================================================================

The core primitive is :func:`port_unitary`, which compiles any 2^n x 2^n unitary into
a gate circuit over a chosen basis (KAK / Shannon decomposition via the transpiler),
and :func:`verify`, which checks the compiled circuit reproduces the target unitary on
the qiskit-aer backend.  :func:`port_and_verify` combines them and reports gate metrics.

These are reused to port each QLB unit operator (rotations, collision).  The same
harness will serve the AI-driven synthesis work later: any target unitary in, a
verified circuit + fidelity out.
"""

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import UnitaryGate

from . import operators as ops
from . import backend as bk

# CNOT + single-qubit rotations (the KAK target of Task 1).
DEFAULT_BASIS = ("rz", "ry", "rx", "cx")
# IBM Heron-native set (kept for later hardware-aware ports).
HERON_BASIS = ("cz", "rz", "sx", "x")

# Two-qubit entangling gates whose counts we report.
_TWO_QUBIT = ("cx", "cz", "ecr")


def num_qubits_for(U):
    """Number of qubits for a 2^n x 2^n unitary (raises if not a power-of-two square)."""
    d = np.asarray(U).shape[0]
    n = int(round(np.log2(d)))
    if 2 ** n != d or np.asarray(U).shape != (d, d):
        raise ValueError(f"operator shape {np.asarray(U).shape} is not 2^n x 2^n")
    return n


def port_unitary(U, basis=DEFAULT_BASIS, optimization_level=3, label=None):
    """
    Compile a unitary matrix `U` into a gate circuit over `basis`.

    Parameters
    ----------
    U : (2^n, 2^n) array — the target unitary.
    basis : tuple[str] — basis gate set to compile to (default: CNOT + rotations).
    optimization_level : int — transpiler optimization level (0-3).
    label : str | None — optional gate label.

    Returns
    -------
    QuantumCircuit implementing `U` (up to global phase) over `basis`.
    """
    n = num_qubits_for(U)
    qc = QuantumCircuit(n)
    qc.append(UnitaryGate(np.asarray(U, dtype=complex), label=label), list(range(n)))
    return transpile(qc, basis_gates=list(basis), optimization_level=optimization_level)


def two_qubit_count(qc):
    """Total number of two-qubit entangling gates in `qc`."""
    ops_count = qc.count_ops()
    return sum(ops_count.get(g, 0) for g in _TWO_QUBIT)


def verify(U, qc=None, basis=DEFAULT_BASIS, tol=1e-9):
    """
    Port `U` (if `qc` not given) and verify the circuit reproduces it.

    Returns a dict with the circuit, phase-invariant gate fidelity, two-qubit gate
    count, depth, basis op counts, and a boolean ``ok`` (fidelity within `tol` of 1).
    """
    if qc is None:
        qc = port_unitary(U, basis=basis)
    U_circ = bk.circuit_unitary(qc)
    fid = bk.gate_fidelity(U, U_circ)
    return {
        "circuit": qc,
        "fidelity": fid,
        "ok": abs(fid - 1.0) < tol,
        "cx": two_qubit_count(qc),
        "depth": qc.depth(),
        "ops": dict(qc.count_ops()),
    }


def port_and_verify(U, basis=DEFAULT_BASIS, tol=1e-9, label=None):
    """Convenience: port `U` and verify in one call (see :func:`verify`)."""
    qc = port_unitary(U, basis=basis, label=label)
    return verify(U, qc=qc, basis=basis, tol=tol)


# --- QLB-specific convenience builders -----------------------------------------
def rotation_circuit(axis, basis=DEFAULT_BASIS):
    """Ported circuit for the QLB rotation R along `axis` ('x' or 'z'; 'y' is identity)."""
    R = ops.ROTATIONS[axis]
    return port_unitary(R, basis=basis, label=f"R_{axis}")


def collision_circuit(m_tilde, g_tilde=0.0, basis=DEFAULT_BASIS):
    """Ported circuit for the standard-frame collision Q_hat(m_tilde, g_tilde)."""
    Q = ops.collision_operator(m_tilde, g_tilde)
    return port_unitary(Q, basis=basis, label="Qhat")
