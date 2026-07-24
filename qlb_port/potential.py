"""
Position-dependent potential (impurity scattering) as a circuit.
================================================================

In the QLB collision Q_hat(m, g) = a_hat I - i b_hat ALPHA_Y the potential enters
through g (the local coupling g_tilde = V_energy * DT / (N_sweeps * hbar)).  For the
**massless** Dirac / graphene case (m = 0) the mass coupling vanishes (b_hat = 0) and
the collision reduces to a position-dependent global phase

    Q_hat(0, g(x)) = a_hat(g(x)) * I ,   |a_hat| = 1 ,

i.e. a diagonal phase oracle diag(a_hat(x)) acting on the position register only
(the same phase for all four spinor components).  This is exactly how a potential
landscape - such as the random impurity barriers of the Klein-tunneling benchmark -
scatters the wave packet.

This module builds that phase oracle and the massless sweep circuit that includes it,
plus a general classical reference (valid for massive potentials too) for validation.

Layout matches operators.streaming_reference: position = qubits 2..(1+n_pos), and the
oracle's diagonal entry x multiplies position basis state |x>.
"""

import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit.library import UnitaryGate, DiagonalGate

from . import operators as ops
from . import streaming as st


def collision_phases(V_tilde):
    """
    Massless collision phases a_hat(x) for a per-site potential coupling array.

    Parameters
    ----------
    V_tilde : array (N,) of per-site g_tilde (dimensionless potential coupling).

    Returns
    -------
    a_hat : complex array (N,), unit modulus (a_hat = collision_coefficients(0, g)[0]).
    """
    V_tilde = np.asarray(V_tilde, dtype=float)
    return np.array([ops.collision_coefficients(0.0, float(g))[0] for g in V_tilde])


def potential_oracle(V_tilde):
    """
    Diagonal phase-oracle gate diag(a_hat(x)) for the massless potential collision.

    Acts on n_pos = log2(len(V_tilde)) position qubits.
    """
    a_hat = collision_phases(V_tilde)
    n = int(round(np.log2(len(a_hat))))
    if 2 ** n != len(a_hat):
        raise ValueError("len(V_tilde) must be a power of two (2**n_pos)")
    return DiagonalGate(list(a_hat))


def impurity_field(n_pos, positions, g_value):
    """
    Build a per-site potential coupling array with `g_value` at `positions`, else 0.

    Parameters
    ----------
    n_pos     : number of position qubits (lattice N = 2**n_pos).
    positions : iterable of lattice sites carrying a barrier.
    g_value   : dimensionless potential coupling g_tilde at those sites.
    """
    V = np.zeros(2 ** n_pos, dtype=float)
    V[np.asarray(list(positions), dtype=int)] = g_value
    return V


def sweep_circuit_potential(axis, n_pos, V_tilde, m_tilde=0.0):
    """
    Massless QLB sub-step with a position-dependent potential:

        rotate R^-1  ->  potential phase oracle diag(a_hat(x))  ->  stream  ->  rotate R

    Only the massless case (m_tilde = 0) is supported here, since then the collision is
    a diagonal phase on the position register.  (A massive potential requires a
    position-multiplexed 2-qubit collision; see sweep_operator_potential for the exact
    classical operator.)
    """
    if m_tilde != 0.0:
        raise NotImplementedError(
            "massive + potential collision needs a position-multiplexed gate; "
            "only m_tilde=0 (graphene/Klein) is ported here.")
    R = ops.ROTATIONS[axis]
    R_inv = R.conj().T
    qc = QuantumCircuit(2 + n_pos, name=f"sweepV_{axis}")
    qc.append(UnitaryGate(R_inv, label="Rinv"), [0, 1])
    qc.append(potential_oracle(V_tilde), list(range(2, 2 + n_pos)))
    qc.compose(st.streaming_circuit(axis, n_pos), inplace=True)
    qc.append(UnitaryGate(R, label="R"), [0, 1])
    return qc


def sweep_operator_potential(axis, n_pos, V_tilde, m_tilde=0.0):
    """
    Exact classical operator of one QLB sub-step with a per-site potential (general;
    valid for massive potentials too).  Collision is block-diagonal over position,
    each block Q_char(m_tilde, V_tilde[x]).
    """
    N = 2 ** n_pos
    R = ops.ROTATIONS[axis]
    R_inv = R.conj().T
    Coll = np.zeros((4 * N, 4 * N), dtype=complex)
    for x in range(N):
        Coll[x * 4:x * 4 + 4, x * 4:x * 4 + 4] = ops.collision_operator_char(
            axis, m_tilde, float(V_tilde[x]))
    S = ops.streaming_reference(axis, n_pos)
    IN = np.eye(N, dtype=complex)
    return np.kron(IN, R) @ S @ Coll @ np.kron(IN, R_inv)


def evolution_circuit_potential(axis, n_pos, n_steps, V_tilde):
    """Circuit for `n_steps` massless sub-steps through a fixed potential landscape."""
    qc = QuantumCircuit(2 + n_pos, name=f"evolveV_{axis}_{n_steps}")
    step = sweep_circuit_potential(axis, n_pos, V_tilde)
    for _ in range(n_steps):
        qc.compose(step, inplace=True)
    return qc
