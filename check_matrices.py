
import numpy as np
import sys

# Pauli Matrices
SIGMA_X = np.array([[0, 1], [1, 0]], dtype=complex)
SIGMA_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
SIGMA_Z = np.array([[1, 0], [0, -1]], dtype=complex)
IDENTITY = np.eye(2, dtype=complex)

# Alpha Matrices
# ZERO = np.zeros((2, 2), dtype=complex)
# ALPHA_X = np.block([[ZERO, SIGMA_X], [SIGMA_X, ZERO]])
# ALPHA_Y = np.block([[ZERO, SIGMA_Y], [SIGMA_Y, ZERO]])
# ALPHA_Z = np.block([[ZERO, SIGMA_Z], [SIGMA_Z, ZERO]])
# BETA = np.block([[IDENTITY, ZERO], [ZERO, -IDENTITY]])

ALPHA_X = np.array([[0, 0, 0, 1],
                    [0, 0, 1, 0],
                    [0, 1, 0, 0],
                    [1, 0, 0, 0]], dtype=complex)

# X Matrix from file
X_MATRIX = (1.0 / np.sqrt(2)) * np.array([[-1, 0, 1, 0],
                                           [0, 1, 0, -1],
                                           [1, 0, 1, 0],
                                           [0, 1, 0, 1]], dtype=complex)
X_INV = X_MATRIX.conj().T

# Check products
A_rot = X_MATRIX @ ALPHA_X @ X_INV
print("X @ ALPHA_X @ X_INV:")
print(np.round(A_rot.real, 2) + 1j*np.round(A_rot.imag, 2))

# Check eigenvalues of A_rot
eig = np.linalg.eigvals(A_rot)
print("Eigenvalues:", eig)

# Check if X_MATRIX acts as expected
# Expected: A_rot should be diagonal (1, 1, -1, -1) or similar.
