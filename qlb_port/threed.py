"""
3D QLB: compose the x-, y-, and z-sweeps on three position registers.
=====================================================================

A 3D Dirac time step is operator-split into x-, y-, then z-sweeps (Lie-Trotter).
On qubits:

    spinor      = qubits 0, 1
    x-position  = qubits 2 .. (1 + nx)
    y-position  = next ny qubits
    z-position  = next nz qubits

Each sweep is the validated 1D :func:`sweep.sweep_circuit` mapped onto the matching
register (spinor shared), so the 3D circuit inherits the 1D correctness.  Dellar's
representation makes this clean: the three spatial streaming matrices are
{ALPHA_X, BETA, ALPHA_Z}, so the z-sweep uses the Z rotation exactly as the 1D case.

Full-state index layout (matches classical_step_3d):
    i = ((z*Ny + y)*Nx + x)*4 + c
with c the spinor component and x/y/z the three registers (x lowest).
"""

import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit.library import UnitaryGate, DiagonalGate

from . import operators as ops
from . import sweep
from . import streaming as st
from . import potential as pot
from .twod import _line_step   # shared 1D line step (rotate->collide->stream->rotate)

AXES = ("x", "y", "z")


def _registers(nx, ny, nz):
    """Qubit lists [spinor(0,1) + position register] for each axis."""
    q = {}
    start = 2
    for axis, n in zip(AXES, (nx, ny, nz)):
        q[axis] = [0, 1] + list(range(start, start + n))
        start += n
    return q


def sweep3d_circuit(nx, ny, nz, m_tilde=0.0, g_tilde=0.0):
    """One 3D QLB step (x-, y-, z-sweeps) for uniform mass/potential."""
    qc = QuantumCircuit(2 + nx + ny + nz, name="sweep3d")
    regs = _registers(nx, ny, nz)
    for axis, n in zip(AXES, (nx, ny, nz)):
        qc.compose(sweep.sweep_circuit(axis, n, m_tilde=m_tilde, g_tilde=g_tilde),
                   qubits=regs[axis], inplace=True)
    return qc


def evolution3d_circuit(nx, ny, nz, n_steps, m_tilde=0.0, g_tilde=0.0):
    """`n_steps` repeated 3D steps."""
    qc = QuantumCircuit(2 + nx + ny + nz, name=f"evolve3d_{n_steps}")
    step = sweep3d_circuit(nx, ny, nz, m_tilde=m_tilde, g_tilde=g_tilde)
    for _ in range(n_steps):
        qc.compose(step, inplace=True)
    return qc


def classical_step_3d(psi_flat, nx, ny, nz, m_tilde=0.0, g_tilde=0.0):
    """One classical 3D QLB step on a flat statevector (i = ((z*Ny+y)*Nx+x)*4+c)."""
    Nx, Ny, Nz = 2 ** nx, 2 ** ny, 2 ** nz
    psi = np.asarray(psi_flat, dtype=complex).reshape(Nz, Ny, Nx, 4).copy()
    for z in range(Nz):                              # x-sweep
        for y in range(Ny):
            psi[z, y] = _line_step(psi[z, y], "x", m_tilde, g_tilde)
    for z in range(Nz):                              # y-sweep
        for x in range(Nx):
            psi[z, :, x, :] = _line_step(psi[z, :, x, :], "y", m_tilde, g_tilde)
    for y in range(Ny):                              # z-sweep
        for x in range(Nx):
            psi[:, y, x, :] = _line_step(psi[:, y, x, :], "z", m_tilde, g_tilde)
    return psi.reshape(-1)


def sweep3d_operator(nx, ny, nz, m_tilde=0.0, g_tilde=0.0):
    """Exact 3D step operator (built column-by-column; use only for small sizes)."""
    dim = 4 * (2 ** nx) * (2 ** ny) * (2 ** nz)
    M = np.zeros((dim, dim), dtype=complex)
    for j in range(dim):
        e = np.zeros(dim, dtype=complex)
        e[j] = 1.0
        M[:, j] = classical_step_3d(e, nx, ny, nz, m_tilde, g_tilde)
    return M


# -----------------------------------------------------------------------------
# 3D with a position-dependent potential V(x, y, z)  (massless / graphene)
# -----------------------------------------------------------------------------
def _oracle_diag_3d(V3d):
    """Phase-oracle diagonal a_hat(V(x,y,z)) ordered by register index x + Nx*y + Nx*Ny*z."""
    V3d = np.asarray(V3d, dtype=float)               # shape (Nx, Ny, Nz), V3d[x,y,z]
    return pot.collision_phases(V3d.flatten("F"))    # F-order: x fastest -> matches register


def sweep3d_circuit_potential(nx, ny, nz, V3d, m_tilde=0.0):
    """Massless 3D QLB step through a potential landscape V(x,y,z) (phase oracle per sweep)."""
    if m_tilde != 0.0:
        raise NotImplementedError("3D potential circuit is massless-only (graphene/Klein).")
    oracle = DiagonalGate(list(_oracle_diag_3d(V3d)))
    posq = list(range(2, 2 + nx + ny + nz))
    regs = _registers(nx, ny, nz)
    qc = QuantumCircuit(2 + nx + ny + nz, name="sweep3dV")
    for axis, n in zip(AXES, (nx, ny, nz)):
        R = ops.ROTATIONS[axis]
        if axis != "y":
            qc.append(UnitaryGate(R.conj().T, label=f"R{axis}_inv"), [0, 1])
        qc.append(oracle, posq)
        qc.compose(st.streaming_circuit(axis, n), qubits=regs[axis], inplace=True)
        if axis != "y":
            qc.append(UnitaryGate(R, label=f"R{axis}"), [0, 1])
    return qc


def evolution3d_circuit_potential(nx, ny, nz, n_steps, V3d, m_tilde=0.0):
    """`n_steps` repeated massless 3D steps through a fixed potential landscape."""
    qc = QuantumCircuit(2 + nx + ny + nz, name=f"evolve3dV_{n_steps}")
    step = sweep3d_circuit_potential(nx, ny, nz, V3d, m_tilde=m_tilde)
    for _ in range(n_steps):
        qc.compose(step, inplace=True)
    return qc


def classical_step_3d_potential(psi_flat, nx, ny, nz, V3d, m_tilde=0.0):
    """Classical 3D QLB step with a per-cell potential V3d (shape (Nx,Ny,Nz), V3d[x,y,z])."""
    Nx, Ny, Nz = 2 ** nx, 2 ** ny, 2 ** nz
    V3d = np.asarray(V3d, dtype=float)
    psi = np.asarray(psi_flat, dtype=complex).reshape(Nz, Ny, Nx, 4).copy()
    for z in range(Nz):                              # x-sweep: g along x = V3d[:, y, z]
        for y in range(Ny):
            psi[z, y] = _line_step(psi[z, y], "x", m_tilde, V3d[:, y, z])
    for z in range(Nz):                              # y-sweep: g along y = V3d[x, :, z]
        for x in range(Nx):
            psi[z, :, x, :] = _line_step(psi[z, :, x, :], "y", m_tilde, V3d[x, :, z])
    for y in range(Ny):                              # z-sweep: g along z = V3d[x, y, :]
        for x in range(Nx):
            psi[:, y, x, :] = _line_step(psi[:, y, x, :], "z", m_tilde, V3d[x, y, :])
    return psi.reshape(-1)


def planar_barrier_field_3d(nx, ny, nz, x_sites, g_value):
    """3D potential (shape (Nx,Ny,Nz)) with a planar barrier at the given x-sites (all y,z)."""
    V = np.zeros((2 ** nx, 2 ** ny, 2 ** nz), dtype=float)
    V[np.asarray(list(x_sites), dtype=int), :, :] = g_value
    return V
