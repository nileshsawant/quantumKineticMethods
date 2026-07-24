"""
2D QLB: compose the x- and y-sweeps on two position registers.
==============================================================

A 2D Dirac time step is operator-split into an x-sweep followed by a y-sweep
(Lie-Trotter).  On qubits:

    spinor      = qubits 0, 1
    x-position  = qubits 2 .. (1 + nx)
    y-position  = qubits (2 + nx) .. (1 + nx + ny)

The x-sweep acts on (spinor, x-register) and the y-sweep on (spinor, y-register),
sharing the spinor qubits.  Each is the validated 1D :func:`sweep.sweep_circuit`
mapped onto the appropriate register, so the 2D circuit inherits the 1D correctness.

Full-state index layout (matches classical_step_2d): i = (y*Nx + x)*4 + c, with
c the spinor component (qubits 0,1), x the x-register, y the y-register.
"""

import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit.library import UnitaryGate

from . import operators as ops
from . import sweep
from . import streaming as st
from . import potential as pot


def _xy_qubits(nx, ny):
    xq = [0, 1] + list(range(2, 2 + nx))
    yq = [0, 1] + list(range(2 + nx, 2 + nx + ny))
    return xq, yq


def sweep2d_circuit(nx, ny, m_tilde=0.0, g_tilde=0.0):
    """One 2D QLB step (x-sweep then y-sweep) for uniform mass/potential."""
    qc = QuantumCircuit(2 + nx + ny, name="sweep2d")
    xq, yq = _xy_qubits(nx, ny)
    qc.compose(sweep.sweep_circuit("x", nx, m_tilde=m_tilde, g_tilde=g_tilde),
               qubits=xq, inplace=True)
    qc.compose(sweep.sweep_circuit("y", ny, m_tilde=m_tilde, g_tilde=g_tilde),
               qubits=yq, inplace=True)
    return qc


def evolution2d_circuit(nx, ny, n_steps, m_tilde=0.0, g_tilde=0.0):
    """`n_steps` repeated 2D steps."""
    qc = QuantumCircuit(2 + nx + ny, name=f"evolve2d_{n_steps}")
    step = sweep2d_circuit(nx, ny, m_tilde=m_tilde, g_tilde=g_tilde)
    for _ in range(n_steps):
        qc.compose(step, inplace=True)
    return qc


# -----------------------------------------------------------------------------
# Classical reference (array-based, mirrors the circuit exactly)
# -----------------------------------------------------------------------------
def _line_step(field, axis, m_tilde, g_tilde):
    """One 1D QLB sub-step on a (N,4) line: rotate -> collide -> stream(periodic) -> rotate.

    g_tilde may be a scalar (uniform) or a per-site array of length N.
    """
    R = ops.ROTATIONS[axis]
    R_inv = R.conj().T
    signs = ops.streaming_signs(axis)
    pr = field @ R_inv.T
    if np.isscalar(g_tilde):
        pr = pr @ ops.collision_operator_char(axis, m_tilde, g_tilde).T
    elif m_tilde == 0.0:                       # massless: per-site scalar phase
        pr = pr * pot.collision_phases(g_tilde)[:, None]
    else:                                      # massive: per-site 2-qubit collision
        for i in range(pr.shape[0]):
            pr[i] = ops.collision_operator_char(axis, m_tilde, float(g_tilde[i])) @ pr[i]
    out = np.empty_like(pr)
    for c in range(4):
        out[:, c] = np.roll(pr[:, c], int(signs[c]))
    return out @ R.T


def classical_step_2d(psi_flat, nx, ny, m_tilde=0.0, g_tilde=0.0):
    """One classical 2D QLB step on a flat statevector (index i = (y*Nx + x)*4 + c)."""
    Nx, Ny = 2 ** nx, 2 ** ny
    psi = np.asarray(psi_flat, dtype=complex).reshape(Ny, Nx, 4).copy()
    for y in range(Ny):                       # x-sweep for each y-row
        psi[y] = _line_step(psi[y], "x", m_tilde, g_tilde)
    for xi in range(Nx):                      # y-sweep for each x-column
        psi[:, xi, :] = _line_step(psi[:, xi, :], "y", m_tilde, g_tilde)
    return psi.reshape(-1)


def sweep2d_operator(nx, ny, m_tilde=0.0, g_tilde=0.0):
    """Exact 2D step operator (built column-by-column; use only for small nx, ny)."""
    dim = 4 * (2 ** nx) * (2 ** ny)
    M = np.zeros((dim, dim), dtype=complex)
    for j in range(dim):
        e = np.zeros(dim, dtype=complex)
        e[j] = 1.0
        M[:, j] = classical_step_2d(e, nx, ny, m_tilde, g_tilde)
    return M


# -----------------------------------------------------------------------------
# 2D with a position-dependent potential V(x, y)  (massless / graphene)
# -----------------------------------------------------------------------------
def _oracle_diag(V2d):
    """Phase-oracle diagonal a_hat(V(x,y)) ordered by register index (y*Nx + x)."""
    V2d = np.asarray(V2d, dtype=float)           # shape (Nx, Ny), V2d[x, y]
    return pot.collision_phases(V2d.T.reshape(-1))   # index = y*Nx + x


def sweep2d_circuit_potential(nx, ny, V2d, m_tilde=0.0):
    """
    Massless 2D QLB step through a potential landscape V(x, y).

    The collision is the diagonal phase oracle diag(a_hat(V(x,y))) on the full
    (x, y) position register, applied in each sweep:

        x-sweep: rotate_x -> oracle(x,y) -> stream_x -> rotate_x^-1
        y-sweep:            oracle(x,y) -> stream_y

    (Only m_tilde = 0 is ported as a circuit; the massive case would need a
    position-multiplexed collision.)
    """
    if m_tilde != 0.0:
        raise NotImplementedError("2D potential circuit is massless-only (graphene/Klein).")
    from qiskit.circuit.library import DiagonalGate
    oracle = DiagonalGate(list(_oracle_diag(V2d)))
    posq = list(range(2, 2 + nx + ny))
    xq, yq = _xy_qubits(nx, ny)
    Rx = ops.ROTATIONS["x"]

    qc = QuantumCircuit(2 + nx + ny, name="sweep2dV")
    # x-sweep
    qc.append(UnitaryGate(Rx.conj().T, label="Rx_inv"), [0, 1])
    qc.append(oracle, posq)
    qc.compose(st.streaming_circuit("x", nx), qubits=xq, inplace=True)
    qc.append(UnitaryGate(Rx, label="Rx"), [0, 1])
    # y-sweep (Ry = I)
    qc.append(oracle, posq)
    qc.compose(st.streaming_circuit("y", ny), qubits=yq, inplace=True)
    return qc


def evolution2d_circuit_potential(nx, ny, n_steps, V2d, m_tilde=0.0):
    """`n_steps` repeated massless 2D steps through a fixed potential landscape."""
    qc = QuantumCircuit(2 + nx + ny, name=f"evolve2dV_{n_steps}")
    step = sweep2d_circuit_potential(nx, ny, V2d, m_tilde=m_tilde)
    for _ in range(n_steps):
        qc.compose(step, inplace=True)
    return qc


def classical_step_2d_potential(psi_flat, nx, ny, V2d, m_tilde=0.0):
    """Classical 2D QLB step with a per-cell potential V2d (shape (Nx, Ny), V2d[x,y])."""
    Nx, Ny = 2 ** nx, 2 ** ny
    V2d = np.asarray(V2d, dtype=float)
    psi = np.asarray(psi_flat, dtype=complex).reshape(Ny, Nx, 4).copy()
    for y in range(Ny):                       # x-sweep: g along x is V2d[:, y]
        psi[y] = _line_step(psi[y], "x", m_tilde, V2d[:, y])
    for xi in range(Nx):                      # y-sweep: g along y is V2d[xi, :]
        psi[:, xi, :] = _line_step(psi[:, xi, :], "y", m_tilde, V2d[xi, :])
    return psi.reshape(-1)


def barrier_field_2d(nx, ny, x_sites, g_value):
    """2D potential (shape (Nx, Ny)) with a vertical barrier at the given x-sites."""
    V = np.zeros((2 ** nx, 2 ** ny), dtype=float)
    V[np.asarray(list(x_sites), dtype=int), :] = g_value
    return V
