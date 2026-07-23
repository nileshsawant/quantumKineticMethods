"""
QLB unit operators as NumPy matrices (single source of truth for the porting layer).
====================================================================================

The QLB Dirac time step is a sequence of per-axis sub-steps, each of the form

    rotate (R^-1)  ->  collide (Q_hat)  ->  stream (+/-1 shift)  ->  rotate back (R)

acting on a 4-component spinor (= 2 qubits).  The portable, gate-level unit operators
are therefore small and exact:

    * ROTATIONS  R in {X, Z}   — fixed 2-qubit unitaries (Y = I).
    * COLLISION  Q_hat(m,g)    — a parameterized 2-qubit unitary, a_hat I - i b_hat ALPHA_Y.
    * STREAMING                — a per-component +/-1 lattice shift (position-register
                                 operator; handled in a later module).

These definitions match ``dirac_qlb_solver`` exactly (see test_port.py::test_consistency).
Representation (Dellar 2011): the spatial streaming matrices are {ALPHA_X, BETA, ALPHA_Z}
and the mass/collision term acts along ALPHA_Y; these four form a mutually anticommuting
Dirac set.  The x-rotation X = S @ Z, with S = exp(-i (pi/4) Sigma_y) the spinor rotation
mapping the z streaming axis onto x, diagonalizes ALPHA_X and preserves the collision
structure (X^-1 Q_hat X = Q).
"""

import numpy as np

# --- Pauli matrices -------------------------------------------------------------
SIGMA_X = np.array([[0, 1], [1, 0]], dtype=np.complex128)
SIGMA_Y = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
SIGMA_Z = np.array([[1, 0], [0, -1]], dtype=np.complex128)
I2 = np.eye(2, dtype=np.complex128)
I4 = np.eye(4, dtype=np.complex128)

# --- Dirac matrices (standard representation) -----------------------------------
_Z2 = np.zeros((2, 2), dtype=np.complex128)
ALPHA_X = np.block([[_Z2, SIGMA_X], [SIGMA_X, _Z2]])
ALPHA_Y = np.block([[_Z2, SIGMA_Y], [SIGMA_Y, _Z2]])   # = SIGMA_X (x) SIGMA_Y ; the mass/collision generator
ALPHA_Z = np.block([[_Z2, SIGMA_Z], [SIGMA_Z, _Z2]])
BETA = np.block([[I2, _Z2], [_Z2, -I2]])

# --- Rotations ------------------------------------------------------------------
# Z diagonalizes ALPHA_Z (Dellar 2011, Eq. 10): Z^-1 ALPHA_Z Z = diag(-1,-1,1,1).
Z_ROTATION = (1.0 / np.sqrt(2)) * np.array([[0, -1, 0, 1],
                                            [1, 0, -1, 0],
                                            [0, 1, 0, 1],
                                            [1, 0, 1, 0]], dtype=np.complex128)

# Spinor rotation S = exp(-i (pi/4) Sigma_y), Sigma_y = blockdiag(sigma_y, sigma_y),
# mapping the z streaming axis onto x (a rotation about y).
_SIGMA_Y_BLOCK = np.block([[SIGMA_Y, _Z2], [_Z2, SIGMA_Y]])
_S_ZTOX = np.cos(np.pi / 4) * I4 - 1j * np.sin(np.pi / 4) * _SIGMA_Y_BLOCK

# X diagonalizes ALPHA_X: X^-1 ALPHA_X X = diag(-1,-1,1,1).
X_ROTATION = _S_ZTOX @ Z_ROTATION
Y_ROTATION = I4.copy()   # BETA is already diagonal along y

ROTATIONS = {"x": X_ROTATION, "y": Y_ROTATION, "z": Z_ROTATION}

# Streaming matrix per axis (Dellar representation) and the per-component +/-1 shift.
STREAM_MATRIX = {"x": ALPHA_X, "y": BETA, "z": ALPHA_Z}


def streaming_signs(axis):
    """Per-spinor-component lattice shift (+/-1) for a sweep along `axis`."""
    R = ROTATIONS[axis]
    diag = np.diag(R.conj().T @ STREAM_MATRIX[axis] @ R).real
    return np.rint(diag).astype(int)


# --- Collision ------------------------------------------------------------------
def collision_coefficients(m_tilde, g_tilde=0.0):
    """
    Return the QLB collision coefficients (a_hat, b_hat).

        Omega = m_tilde^2 - g_tilde^2
        D     = 1 + Omega/4 - i g_tilde
        a_hat = (1 - Omega/4) / D          (diagonal / phase factor)
        b_hat = m_tilde / D                (mass coupling; 0 for massless)

    They satisfy |a_hat|^2 + |b_hat|^2 = 1, so the collision is unitary.  Here
    m_tilde and g_tilde are the (already per-sweep-rescaled) dimensionless mass and
    potential couplings.
    """
    Omega = m_tilde ** 2 - g_tilde ** 2
    D = 1.0 + Omega / 4.0 - 1j * g_tilde
    a_hat = (1.0 - Omega / 4.0) / D
    b_hat = m_tilde / D
    return a_hat, b_hat


def collision_operator(m_tilde, g_tilde=0.0):
    """
    Standard-frame collision unitary  Q_hat = a_hat I - i b_hat ALPHA_Y.

    This is the operator applied on the spinor (2 qubits) for a site with the given
    local mass and potential couplings.  In the characteristic frame of a sweep it
    becomes  a_hat I - i b_hat (R^-1 ALPHA_Y R)  (see collision_operator_char).
    """
    a_hat, b_hat = collision_coefficients(m_tilde, g_tilde)
    return a_hat * I4 - 1j * b_hat * ALPHA_Y


def collision_operator_char(axis, m_tilde, g_tilde=0.0):
    """Collision in the characteristic frame of `axis`: R^-1 Q_hat R."""
    R = ROTATIONS[axis]
    return R.conj().T @ collision_operator(m_tilde, g_tilde) @ R


# --- Streaming ------------------------------------------------------------------
# Qubit/index convention for the combined spinor(x)position state:
#   * spinor  = qubits 0,1  (qubit 0 = least significant); component c = index & 3
#   * position = qubits 2..(1+n_pos); position x = index >> 2
# so the full basis index is  i = x*4 + c  on a lattice of N = 2**n_pos sites.
# Streaming shifts each component c by streaming_signs(axis)[c] sites (periodic /
# modular wrap-around), which matches the periodic y- and z-sweeps and the bulk of
# the x-sweep; the open / bounce-back x boundary is a separate boundary operator.
def streaming_reference(axis, n_pos):
    """
    Classical streaming permutation on the (spinor x position) space for `axis`.

    Returns a (4*N, 4*N) permutation matrix, N = 2**n_pos, with i = x*4 + c layout.
    """
    signs = streaming_signs(axis)
    N = 2 ** n_pos
    dim = 4 * N
    P = np.zeros((dim, dim), dtype=np.complex128)
    for c in range(4):
        s = int(signs[c])
        for x in range(N):
            P[((x + s) % N) * 4 + c, x * 4 + c] = 1.0
    return P


# Convenient aliases used by the porting harness / tests
X_ROT = X_ROTATION
Z_ROT = Z_ROTATION

