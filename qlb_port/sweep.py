"""
Assemble the ported unit operators into a full single-axis QLB sweep circuit.
=============================================================================

One QLB sub-step along an axis is

    rotate (R^-1)  ->  collide (Q_char)  ->  stream (+/-1 shift)  ->  rotate back (R)

acting on a 4-component spinor (qubits 0,1) at every lattice site (position qubits
2..(1+n_pos)).  For a **free particle** the mass/potential couplings are uniform, so
the collision Q_char is a single fixed 2-qubit gate; the potential-dependent
(position-varying) collision is handled separately as a later phase-oracle stage.

The classical reference :func:`sweep_operator` is the exact linear operator of one
solver sub-step on the (spinor x position) space (periodic streaming), and
:func:`sweep_circuit` is its gate-level implementation.  They agree to machine
precision (see test_port.py::test_sweep_assembly).

Qubit / index layout (matches operators.streaming_reference): i = x*4 + c, with
spinor = qubits 0,1 (LSB) and position = qubits 2..(1+n_pos).
"""

import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit.library import UnitaryGate

from . import operators as ops
from . import streaming as st


def _embed_spinor(A, n_pos):
    """Embed a 4x4 spinor operator into the (spinor x position) space (i = x*4 + c)."""
    return np.kron(np.eye(2 ** n_pos, dtype=complex), np.asarray(A, dtype=complex))


def sweep_operator(axis, n_pos, m_tilde=0.0, g_tilde=0.0):
    """
    Exact classical operator of one QLB sub-step along `axis` (periodic streaming).

        Sub = (R) . Stream . (Q_char) . (R^-1)

    Returns a (4*N, 4*N) unitary, N = 2**n_pos.
    """
    R = ops.ROTATIONS[axis]
    R_inv = R.conj().T
    Q_char = ops.collision_operator_char(axis, m_tilde, g_tilde)
    S = ops.streaming_reference(axis, n_pos)
    return _embed_spinor(R, n_pos) @ S @ _embed_spinor(Q_char, n_pos) @ _embed_spinor(R_inv, n_pos)


def sweep_circuit(axis, n_pos, m_tilde=0.0, g_tilde=0.0):
    """
    Gate-level circuit for one QLB sub-step along `axis`.

    Parameters
    ----------
    axis    : 'x', 'y', or 'z'
    n_pos   : number of position qubits (lattice N = 2**n_pos)
    m_tilde : (per-sweep) dimensionless mass coupling (0 => massless)
    g_tilde : (per-sweep) dimensionless potential coupling (0 => free particle)

    Returns
    -------
    QuantumCircuit on (2 + n_pos) qubits.
    """
    R = ops.ROTATIONS[axis]
    R_inv = R.conj().T
    Q_char = ops.collision_operator_char(axis, m_tilde, g_tilde)

    qc = QuantumCircuit(2 + n_pos, name=f"sweep_{axis}")
    qc.append(UnitaryGate(R_inv, label="Rinv"), [0, 1])       # rotate into char frame
    qc.append(UnitaryGate(Q_char, label="Qchar"), [0, 1])     # collide
    qc.compose(st.streaming_circuit(axis, n_pos), inplace=True)  # stream +/-1
    qc.append(UnitaryGate(R, label="R"), [0, 1])              # rotate back
    return qc


def evolution_circuit(axis, n_pos, n_steps, m_tilde=0.0, g_tilde=0.0):
    """Circuit for `n_steps` repeated single-axis sub-steps (free-particle evolution)."""
    qc = QuantumCircuit(2 + n_pos, name=f"evolve_{axis}_{n_steps}")
    step = sweep_circuit(axis, n_pos, m_tilde, g_tilde)
    for _ in range(n_steps):
        qc.compose(step, inplace=True)
    return qc
