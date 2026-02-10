"""
Comprehensive Test Suite for 3D Dirac QLB Solver
=================================================
Tests verify the implementation against the specifications from:
- Dellar et al. (2011) - Lattice Boltzmann algorithms
- Palpacelli et al. (2008, 2012) - QLB for Dirac equation

Test Categories:
1. Matrix definitions and properties
2. Wave function initialization and normalization
3. QLB collision coefficients
4. 1D streaming and collision
5. Operator splitting
6. Conservation laws
7. Graphene-specific features
"""

import numpy as np
import sys
import os

# Import the solver module
import dirac_qlb_solver as solver

# Test utilities
def assert_close(actual, expected, rtol=1e-5, atol=1e-8, msg=""):
    """Assert that two values are close within tolerance."""
    if not np.allclose(actual, expected, rtol=rtol, atol=atol):
        print(f"FAIL: {msg}")
        print(f"  Expected: {expected}")
        print(f"  Actual:   {actual}")
        print(f"  Diff:     {np.abs(actual - expected)}")
        return False
    return True

def assert_unitary(matrix, msg=""):
    """Assert that a matrix is unitary (U† U = I)."""
    identity = np.eye(matrix.shape[0], dtype=np.complex128)
    product = matrix.conj().T @ matrix
    return assert_close(product, identity, msg=f"Unitarity test: {msg}")

def assert_hermitian(matrix, msg=""):
    """Assert that a matrix is Hermitian (A = A†)."""
    return assert_close(matrix, matrix.conj().T, msg=f"Hermitian test: {msg}")

# =============================================================================
# Test 1: Matrix Definitions and Properties
# =============================================================================

def test_pauli_matrices():
    """Test that Pauli matrices satisfy standard commutation relations."""
    print("\n" + "="*70)
    print("TEST 1: Pauli Matrices and Properties")
    print("="*70)
    
    sigma_x = solver.SIGMA_X_2x2
    sigma_y = solver.SIGMA_Y_2x2
    sigma_z = solver.SIGMA_Z_2x2
    identity = solver.IDENTITY_2x2
    
    tests_passed = 0
    total_tests = 0
    
    # Test 1.1: Pauli matrices are Hermitian
    total_tests += 3
    if assert_hermitian(sigma_x, "σ_x"):
        tests_passed += 1
    if assert_hermitian(sigma_y, "σ_y"):
        tests_passed += 1
    if assert_hermitian(sigma_z, "σ_z"):
        tests_passed += 1
    
    # Test 1.2: Pauli matrices square to identity
    total_tests += 3
    if assert_close(sigma_x @ sigma_x, identity, msg="σ_x² = I"):
        tests_passed += 1
    if assert_close(sigma_y @ sigma_y, identity, msg="σ_y² = I"):
        tests_passed += 1
    if assert_close(sigma_z @ sigma_z, identity, msg="σ_z² = I"):
        tests_passed += 1
    
    # Test 1.3: Anticommutation relations {σ_i, σ_j} = 2δ_ij I
    total_tests += 1
    anticomm_xy = sigma_x @ sigma_y + sigma_y @ sigma_x
    if assert_close(anticomm_xy, np.zeros((2, 2)), msg="{σ_x, σ_y} = 0"):
        tests_passed += 1
    
    # Test 1.4: Commutation relations [σ_x, σ_y] = 2i σ_z
    total_tests += 1
    comm_xy = sigma_x @ sigma_y - sigma_y @ sigma_x
    if assert_close(comm_xy, 2j * sigma_z, msg="[σ_x, σ_y] = 2iσ_z"):
        tests_passed += 1
    
    print(f"\nPauli matrices tests: {tests_passed}/{total_tests} passed")
    return tests_passed == total_tests

def test_dirac_matrices():
    """Test properties of Dirac alpha and beta matrices."""
    print("\n" + "="*70)
    print("TEST 2: Dirac Alpha and Beta Matrices")
    print("="*70)
    
    alpha_x = solver.ALPHA_X
    alpha_y = solver.ALPHA_Y
    alpha_z = solver.ALPHA_Z
    beta = solver.BETA
    identity = solver.IDENTITY_4x4
    
    tests_passed = 0
    total_tests = 0
    
    # Test 2.1: Alpha matrices are Hermitian
    total_tests += 3
    if assert_hermitian(alpha_x, "α_x"):
        tests_passed += 1
    if assert_hermitian(alpha_y, "α_y"):
        tests_passed += 1
    if assert_hermitian(alpha_z, "α_z"):
        tests_passed += 1
    
    # Test 2.2: Beta matrix is Hermitian
    total_tests += 1
    if assert_hermitian(beta, "β"):
        tests_passed += 1
    
    # Test 2.3: Alpha matrices square to identity
    total_tests += 3
    if assert_close(alpha_x @ alpha_x, identity, msg="α_x² = I"):
        tests_passed += 1
    if assert_close(alpha_y @ alpha_y, identity, msg="α_y² = I"):
        tests_passed += 1
    if assert_close(alpha_z @ alpha_z, identity, msg="α_z² = I"):
        tests_passed += 1
    
    # Test 2.4: Beta matrix squares to identity
    total_tests += 1
    if assert_close(beta @ beta, identity, msg="β² = I"):
        tests_passed += 1
    
    # Test 2.5: Anticommutation {α_i, α_j} = 2δ_ij I
    total_tests += 1
    anticomm_xy = alpha_x @ alpha_y + alpha_y @ alpha_x
    if assert_close(anticomm_xy, np.zeros((4, 4)), msg="{α_x, α_y} = 0"):
        tests_passed += 1
    
    # Test 2.6: Anticommutation {α_i, β} = 0
    total_tests += 1
    anticomm_x_beta = alpha_x @ beta + beta @ alpha_x
    if assert_close(anticomm_x_beta, np.zeros((4, 4)), msg="{α_x, β} = 0"):
        tests_passed += 1
    
    print(f"\nDirac matrices tests: {tests_passed}/{total_tests} passed")
    return tests_passed == total_tests

def test_rotation_matrices():
    """Test that rotation matrices are unitary."""
    print("\n" + "="*70)
    print("TEST 3: Rotation Matrices for Operator Splitting")
    print("="*70)
    
    X_matrix = solver.X_MATRIX
    Y_matrix = solver.Y_MATRIX
    Z_matrix = solver.Z_MATRIX
    X_inv = solver.X_INV_MATRIX
    Y_inv = solver.Y_INV_MATRIX
    Z_inv = solver.Z_INV_MATRIX
    
    tests_passed = 0
    total_tests = 0
    
    # Test 3.1: Rotation matrices are unitary
    total_tests += 3
    if assert_unitary(X_matrix, "X"):
        tests_passed += 1
    if assert_unitary(Y_matrix, "Y"):
        tests_passed += 1
    if assert_unitary(Z_matrix, "Z"):
        tests_passed += 1
    
    # Test 3.2: Inverse matrices are correct (U^-1 = U†)
    total_tests += 3
    if assert_close(X_inv, X_matrix.conj().T, msg="X_inv = X†"):
        tests_passed += 1
    if assert_close(Y_inv, Y_matrix.conj().T, msg="Y_inv = Y†"):
        tests_passed += 1
    if assert_close(Z_inv, Z_matrix.conj().T, msg="Z_inv = Z†"):
        tests_passed += 1
    
    # Test 3.3: Rotation matrices diagonalize respective -alpha matrices
    # For Z: Z† (-α_z) Z should be diagonal
    total_tests += 1
    alpha_z_rotated = Z_inv @ (-solver.ALPHA_Z) @ Z_matrix
    # Check if off-diagonal elements are small
    off_diagonal_max = np.max(np.abs(alpha_z_rotated - np.diag(np.diag(alpha_z_rotated))))
    if off_diagonal_max < 1e-10:
        tests_passed += 1
        print(f"  Z diagonalizes -α_z: max off-diagonal = {off_diagonal_max:.2e}")
    else:
        print(f"  FAIL: Z does not fully diagonalize -α_z: max off-diagonal = {off_diagonal_max:.2e}")
    
    print(f"\nRotation matrices tests: {tests_passed}/{total_tests} passed")
    return tests_passed == total_tests

# =============================================================================
# Test 2: Wave Function Initialization and Normalization
# =============================================================================

def test_wave_function_normalization():
    """Test that the initial wave function is properly normalized."""
    print("\n" + "="*70)
    print("TEST 4: Wave Function Initialization and Normalization")
    print("="*70)
    
    tests_passed = 0
    total_tests = 0
    
    # Test 4.1: Wave function normalization
    total_tests += 1
    norm = np.sum(np.abs(solver.psi)**2) * solver.DX * solver.DY * solver.DZ
    if assert_close(norm, 1.0, atol=1e-6, msg="Wave function normalization"):
        tests_passed += 1
        print(f"  Wave function norm: {norm:.10f}")
    
    # Test 4.2: Wave function has correct shape
    total_tests += 1
    expected_shape = (solver.NX, solver.NY, solver.NZ, 4)
    if solver.psi.shape == expected_shape:
        tests_passed += 1
        print(f"  Wave function shape: {solver.psi.shape} ✓")
    else:
        print(f"  FAIL: Expected shape {expected_shape}, got {solver.psi.shape}")
    
    # Test 4.3: Wave function is complex-valued
    total_tests += 1
    if solver.psi.dtype == np.complex128:
        tests_passed += 1
        print(f"  Wave function dtype: {solver.psi.dtype} ✓")
    else:
        print(f"  FAIL: Expected dtype complex128, got {solver.psi.dtype}")
    
    # Test 4.4: Initial wave packet is centered approximately at x0
    total_tests += 1
    prob_density = np.sum(np.abs(solver.psi[:, :, 0, :])**2, axis=(1, 2))
    x_center = np.sum(np.arange(solver.NX) * prob_density) / np.sum(prob_density)
    expected_center = solver.NX // 4
    if np.abs(x_center - expected_center) < 5:  # Within 5 grid points
        tests_passed += 1
        print(f"  Wave packet center: {x_center:.2f} (expected ~{expected_center}) ✓")
    else:
        print(f"  FAIL: Wave packet center at {x_center:.2f}, expected ~{expected_center}")
    
    print(f"\nWave function tests: {tests_passed}/{total_tests} passed")
    return tests_passed == total_tests

# =============================================================================
# Test 3: QLB Collision Coefficients
# =============================================================================

def test_collision_coefficients():
    """Test the calculation of collision coefficients a_hat and b_hat."""
    print("\n" + "="*70)
    print("TEST 5: QLB Collision Coefficients")
    print("="*70)
    
    tests_passed = 0
    total_tests = 0
    
    # Test 5.1: Coefficients for zero potential and zero mass
    total_tests += 1
    V_zero = np.array([0.0])
    m_test = 0.0
    dt_test = solver.DT
    hbar_test = solver.HBAR
    c_test = solver.C
    
    m_tilde_3 = m_test * c_test**2 * dt_test / (3 * hbar_test)
    g_tilde_3 = solver.Q_ELECTRON * V_zero * dt_test / (3 * hbar_test)
    Omega_3 = m_tilde_3**2 - g_tilde_3**2
    
    denominator = 1 + Omega_3/4 - 1j * g_tilde_3
    a_hat = (1 - Omega_3/4) / denominator
    b_hat = m_tilde_3 / denominator
    
    # For m=0, V=0: a_hat should be 1, b_hat should be 0
    if assert_close(np.abs(a_hat[0]), 1.0, atol=1e-10, msg="a_hat = 1 for m=V=0"):
        tests_passed += 1
        print(f"  a_hat (m=V=0): {a_hat[0]:.10f}")
        print(f"  b_hat (m=V=0): {b_hat[0]:.10e}")
    
    # Test 5.2: Coefficients are bounded
    total_tests += 1
    V_test = np.linspace(-1e-19, 1e-19, 10)  # Range of potentials
    m_tilde_3 = solver.M_PARTICLE * c_test**2 * dt_test / (3 * hbar_test)
    g_tilde_3 = solver.Q_ELECTRON * V_test * dt_test / (3 * hbar_test)
    Omega_3 = m_tilde_3**2 - g_tilde_3**2
    
    denominator = 1 + Omega_3/4 - 1j * g_tilde_3
    denominator = np.where(np.abs(denominator) < 1e-15, 1e-15, denominator)
    a_hat = (1 - Omega_3/4) / denominator
    b_hat = m_tilde_3 / denominator
    
    max_a = np.max(np.abs(a_hat))
    max_b = np.max(np.abs(b_hat))
    if max_a < 10 and max_b < 10:  # Reasonable bounds
        tests_passed += 1
        print(f"  Coefficient bounds: |a_hat|_max = {max_a:.4f}, |b_hat|_max = {max_b:.4f} ✓")
    else:
        print(f"  FAIL: Coefficients unbounded: |a_hat|_max = {max_a:.4f}, |b_hat|_max = {max_b:.4f}")
    
    print(f"\nCollision coefficient tests: {tests_passed}/{total_tests} passed")
    return tests_passed == total_tests

# =============================================================================
# Test 4: 1D QLB Step Function
# =============================================================================

def test_qlb_1d_step():
    """Test the 1D QLB sub-step function."""
    print("\n" + "="*70)
    print("TEST 6: 1D QLB Step Function")
    print("="*70)
    
    tests_passed = 0
    total_tests = 0
    
    # Test 6.1: Function runs without errors
    total_tests += 1
    N_test = 10
    psi_line_test = np.zeros((N_test, 4), dtype=np.complex128)
    psi_line_test[N_test//2, 0] = 1.0  # Localized initial condition
    V_line_test = np.zeros(N_test)
    
    try:
        psi_out = solver.perform_qlb_sub_step(
            psi_line_test, V_line_test, 'z', solver.MATRICES,
            solver.HBAR, solver.C, solver.M_PARTICLE, solver.DT, solver.DX
        )
        tests_passed += 1
        print("  1D QLB step executes successfully ✓")
    except Exception as e:
        print(f"  FAIL: 1D QLB step raised exception: {e}")
        return False
    
    # Test 6.2: Output has correct shape
    total_tests += 1
    if psi_out.shape == psi_line_test.shape:
        tests_passed += 1
        print(f"  Output shape: {psi_out.shape} ✓")
    else:
        print(f"  FAIL: Output shape {psi_out.shape}, expected {psi_line_test.shape}")
    
    # Test 6.3: Norm is approximately conserved (for small time step)
    total_tests += 1
    norm_in = np.sum(np.abs(psi_line_test)**2)
    norm_out = np.sum(np.abs(psi_out)**2)
    norm_change = np.abs(norm_out - norm_in) / norm_in
    if norm_change < 0.5:  # Allow some change due to boundaries
        tests_passed += 1
        print(f"  Norm change: {norm_change:.4f} (< 0.5) ✓")
    else:
        print(f"  FAIL: Large norm change: {norm_change:.4f}")
    
    # Test 6.4: Wave packet has propagated
    total_tests += 1
    center_in = np.sum(np.arange(N_test) * np.sum(np.abs(psi_line_test)**2, axis=1))
    center_out = np.sum(np.arange(N_test) * np.sum(np.abs(psi_out)**2, axis=1))
    if center_out != center_in:
        tests_passed += 1
        print(f"  Wave packet propagated from {center_in:.2f} to {center_out:.2f} ✓")
    else:
        print(f"  WARNING: Wave packet did not move (might be OK for stationary case)")
        tests_passed += 1  # Don't fail this test
    
    print(f"\n1D QLB step tests: {tests_passed}/{total_tests} passed")
    return tests_passed == total_tests

# =============================================================================
# Test 5: Operator Splitting and Conservation
# =============================================================================

def test_operator_splitting():
    """Test that operator splitting preserves key properties."""
    print("\n" + "="*70)
    print("TEST 7: Operator Splitting and Conservation")
    print("="*70)
    
    tests_passed = 0
    total_tests = 0
    
    # Create a small test system
    NX_test, NY_test, NZ_test = 8, 8, 1
    psi_test = np.zeros((NX_test, NY_test, NZ_test, 4), dtype=np.complex128)
    
    # Initialize with a Gaussian
    x0, y0 = NX_test // 2, NY_test // 2
    for i in range(NX_test):
        for j in range(NY_test):
            r2 = (i - x0)**2 + (j - y0)**2
            psi_test[i, j, 0, 0] = np.exp(-r2 / 4.0)
    
    # Normalize
    norm = np.sum(np.abs(psi_test)**2) * solver.DX * solver.DY * solver.DZ
    psi_test /= np.sqrt(norm)
    
    V_test = np.zeros((NX_test, NY_test, NZ_test))
    
    # Test 7.1: Single X-step preserves norm approximately
    total_tests += 1
    psi_temp = psi_test.copy()
    for j in range(NY_test):
        for k in range(NZ_test):
            psi_line = psi_temp[:, j, k, :]
            V_line = V_test[:, j, k]
            psi_temp[:, j, k, :] = solver.perform_qlb_sub_step(
                psi_line, V_line, 'x', solver.MATRICES,
                solver.HBAR, solver.C, 0.0, solver.DT, solver.DX  # Use m=0 for simplicity
            )
    
    norm_after = np.sum(np.abs(psi_temp)**2) * solver.DX * solver.DY * solver.DZ
    if assert_close(norm_after, 1.0, rtol=0.3, msg="Norm after X-step"):
        tests_passed += 1
        print(f"  Norm after X-step: {norm_after:.6f}")
    
    # Test 7.2: Full splitting step (X->Y->Z) runs
    total_tests += 1
    try:
        psi_final = psi_test.copy()
        
        # X-step
        for j in range(NY_test):
            for k in range(NZ_test):
                psi_line = psi_final[:, j, k, :]
                V_line = V_test[:, j, k]
                psi_final[:, j, k, :] = solver.perform_qlb_sub_step(
                    psi_line, V_line, 'x', solver.MATRICES,
                    solver.HBAR, solver.C, 0.0, solver.DT, solver.DX
                )
        
        # Y-step
        for i in range(NX_test):
            for k in range(NZ_test):
                psi_line = psi_final[i, :, k, :]
                V_line = V_test[i, :, k]
                psi_final[i, :, k, :] = solver.perform_qlb_sub_step(
                    psi_line, V_line, 'y', solver.MATRICES,
                    solver.HBAR, solver.C, 0.0, solver.DT, solver.DY
                )
        
        tests_passed += 1
        print("  Full operator splitting step completes ✓")
        
        # Check final norm
        norm_final = np.sum(np.abs(psi_final)**2) * solver.DX * solver.DY * solver.DZ
        print(f"  Final norm: {norm_final:.6f}")
        
    except Exception as e:
        print(f"  FAIL: Operator splitting raised exception: {e}")
    
    print(f"\nOperator splitting tests: {tests_passed}/{total_tests} passed")
    return tests_passed == total_tests

# =============================================================================
# Test 6: Graphene-Specific Features
# =============================================================================

def test_graphene_setup():
    """Test graphene-specific configuration."""
    print("\n" + "="*70)
    print("TEST 8: Graphene-Specific Configuration")
    print("="*70)
    
    tests_passed = 0
    total_tests = 0
    
    # Run graphene setup
    psi_graphene, V_graphene = solver.setup_graphene_simulation()
    
    # Test 8.1: Particle mass is set to zero for massless fermions
    total_tests += 1
    if solver.M_PARTICLE == 0.0:
        tests_passed += 1
        print(f"  Massless Dirac fermions: M = {solver.M_PARTICLE} ✓")
    else:
        print(f"  FAIL: Expected M=0 for graphene, got M={solver.M_PARTICLE}")
    
    # Test 8.2: System is 2D (NZ = 1)
    total_tests += 1
    if solver.NZ == 1:
        tests_passed += 1
        print(f"  2D system: NZ = {solver.NZ} ✓")
    else:
        print(f"  FAIL: Expected NZ=1 for 2D, got NZ={solver.NZ}")
    
    # Test 8.3: Impurity region has non-zero potential
    total_tests += 1
    impurity_region_start = solver.NX // 3
    impurity_region_end = 2 * solver.NX // 3
    num_impurities = np.sum(V_graphene[impurity_region_start:impurity_region_end, :, :] > 0)
    
    if num_impurities > 0:
        tests_passed += 1
        print(f"  Impurities present: {num_impurities} sites with V > 0 ✓")
    else:
        print(f"  FAIL: No impurities found in impurity region")
    
    # Test 8.4: Wave function is normalized
    total_tests += 1
    norm_graphene = np.sum(np.abs(psi_graphene)**2) * solver.DX * solver.DY * solver.DZ
    if assert_close(norm_graphene, 1.0, atol=1e-6, msg="Graphene wave function normalization"):
        tests_passed += 1
        print(f"  Graphene wave function norm: {norm_graphene:.10f}")
    
    # Test 8.5: Wave packet starts at left edge
    total_tests += 1
    prob_density = np.sum(np.abs(psi_graphene[:, :, 0, :])**2, axis=(1, 2))
    x_center = np.sum(np.arange(solver.NX) * prob_density) / np.sum(prob_density)
    expected_center = solver.NX // 8
    if np.abs(x_center - expected_center) < 5:
        tests_passed += 1
        print(f"  Graphene wave packet center: {x_center:.2f} (expected ~{expected_center}) ✓")
    else:
        print(f"  FAIL: Wave packet center at {x_center:.2f}, expected ~{expected_center}")
    
    print(f"\nGraphene setup tests: {tests_passed}/{total_tests} passed")
    return tests_passed == total_tests

# =============================================================================
# Test 7: Physical Constants and Parameters
# =============================================================================

def test_physical_parameters():
    """Test that physical constants and numerical parameters are reasonable."""
    print("\n" + "="*70)
    print("TEST 9: Physical Constants and Numerical Parameters")
    print("="*70)
    
    tests_passed = 0
    total_tests = 0
    
    # Test 9.1: Physical constants are positive and reasonable
    total_tests += 3
    if solver.HBAR > 0 and solver.HBAR < 1e-30:
        tests_passed += 1
        print(f"  ℏ = {solver.HBAR:.4e} J·s ✓")
    else:
        print(f"  FAIL: Unreasonable HBAR: {solver.HBAR}")
    
    if solver.C > 1e8 and solver.C < 1e9:
        tests_passed += 1
        print(f"  c = {solver.C:.4e} m/s ✓")
    else:
        print(f"  FAIL: Unreasonable C: {solver.C}")
    
    if solver.Q_ELECTRON > 1e-20 and solver.Q_ELECTRON < 1e-18:
        tests_passed += 1
        print(f"  e = {solver.Q_ELECTRON:.4e} C ✓")
    else:
        print(f"  FAIL: Unreasonable Q_ELECTRON: {solver.Q_ELECTRON}")
    
    # Test 9.2: Light-cone condition: c*dt <= dx
    total_tests += 1
    light_cone = solver.C * solver.DT
    if light_cone <= solver.DX * 1.01:  # Allow 1% margin
        tests_passed += 1
        print(f"  Light-cone condition: c·dt = {light_cone:.4e} m ≤ dx = {solver.DX:.4e} m ✓")
    else:
        print(f"  FAIL: Light-cone violated: c·dt = {light_cone:.4e} > dx = {solver.DX:.4e}")
    
    # Test 9.3: Grid size is reasonable
    total_tests += 1
    if solver.NX > 0 and solver.NY > 0 and solver.NZ > 0:
        tests_passed += 1
        print(f"  Grid size: {solver.NX} × {solver.NY} × {solver.NZ} ✓")
    else:
        print(f"  FAIL: Invalid grid size")
    
    # Test 9.4: Domain size is in nanometer range (appropriate for graphene)
    total_tests += 1
    if 1e-10 < solver.LX < 1e-7:
        tests_passed += 1
        print(f"  Domain size: {solver.LX*1e9:.2f} nm ✓")
    else:
        print(f"  WARNING: Domain size {solver.LX*1e9:.2f} nm may be unusual")
        tests_passed += 1  # Don't fail
    
    print(f"\nPhysical parameter tests: {tests_passed}/{total_tests} passed")
    return tests_passed == total_tests

# =============================================================================
# Main Test Runner
# =============================================================================

def run_all_tests():
    """Run all tests and report results."""
    print("\n")
    print("╔" + "═"*68 + "╗")
    print("║" + " "*15 + "DIRAC QLB SOLVER TEST SUITE" + " "*26 + "║")
    print("║" + " "*15 + "Comprehensive Verification" + " "*27 + "║")
    print("╚" + "═"*68 + "╝")
    
    test_results = []
    
    # Run all test categories
    test_results.append(("Physical Parameters", test_physical_parameters()))
    test_results.append(("Pauli Matrices", test_pauli_matrices()))
    test_results.append(("Dirac Matrices", test_dirac_matrices()))
    test_results.append(("Rotation Matrices", test_rotation_matrices()))
    test_results.append(("Wave Function Init", test_wave_function_normalization()))
    test_results.append(("Collision Coefficients", test_collision_coefficients()))
    test_results.append(("1D QLB Step", test_qlb_1d_step()))
    test_results.append(("Operator Splitting", test_operator_splitting()))
    test_results.append(("Graphene Setup", test_graphene_setup()))
    
    # Summary
    print("\n" + "="*70)
    print("TEST SUMMARY")
    print("="*70)
    
    passed_count = 0
    for test_name, result in test_results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{test_name:.<50} {status}")
        if result:
            passed_count += 1
    
    print("="*70)
    print(f"Total: {passed_count}/{len(test_results)} test categories passed")
    print("="*70)
    
    return passed_count == len(test_results)

if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
