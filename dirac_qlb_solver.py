"""
3D Dirac Equation Solver using Quantum Lattice Boltzmann (QLB) Method
=====================================================================
This solver implements the QLB scheme for the Dirac equation with applications
to electron transport in graphene and other relativistic quantum systems.

References:
[1] Dellar, P. J. (2011). "Lattice Boltzmann algorithms without cubic defects..."
[3] Palpacelli, S., et al. (2012). "Quantum lattice Boltzmann simulation..."
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# =============================================================================
# Part 1: Project Setup and Global Constants
# =============================================================================

# 1.1 Physical Constants (SI units for clarity, often scaled in practice)
HBAR = 1.0545718e-34      # Reduced Planck's constant (J·s)
C = 2.99792458e8          # Speed of light (m/s)
Q_ELECTRON = 1.60217663e-19  # Elementary charge (C)

# 1.2 Numerical Parameters - Spatial dimensions
NX = 64   # Number of grid points in x-direction
NY = 64   # Number of grid points in y-direction
NZ = 1    # For 2D graphene example, set NZ=1. For full 3D, set NZ > 1.

# Physical domain size (m)
LX = 1e-9  # 1 nm
LY = 1e-9  # 1 nm
LZ = 1e-9  # 1 nm

# Discretization
DX = LX / NX
DY = LY / NY
DZ = LZ / NZ if NZ > 1 else DX

# Time step (s) - derived from light-cone rule: c * dt = dx
DT = DX / C

# Particle mass (kg)
# For graphene massless Dirac fermions, set M_PARTICLE = 0
# For standard electron: 9.1093837e-31 kg
M_PARTICLE = 9.1093837e-31  # Will be modified for graphene

# Compton frequency: omega_c = m*c^2 / hbar
OMEGA_C = M_PARTICLE * C**2 / HBAR if M_PARTICLE > 0 else 0.0

# Total simulation steps
T_STEPS = 500

# Output frequency for saving data/plots
OUTPUT_FREQ = 50

print("=" * 70)
print("3D Dirac QLB Solver - Simulation Parameters")
print("=" * 70)
print(f"Grid:                {NX} x {NY} x {NZ}")
print(f"Domain size:         {LX*1e9:.2f} x {LY*1e9:.2f} x {LZ*1e9:.2f} nm")
print(f"Spatial step (Dx):   {DX:.2e} m")
print(f"Time step (Dt):      {DT:.2e} s")
print(f"Particle Mass:       {M_PARTICLE:.2e} kg")
print(f"Compton Frequency:   {OMEGA_C:.2e} rad/s")
print("=" * 70)

# =============================================================================
# Part 2: Dirac Equation and Rotation Matrices
# =============================================================================

# 2.1 Pauli Matrices (2x2)
SIGMA_X_2x2 = np.array([[0, 1], 
                        [1, 0]], dtype=np.complex128)
SIGMA_Y_2x2 = np.array([[0, -1j], 
                        [1j, 0]], dtype=np.complex128)
SIGMA_Z_2x2 = np.array([[1, 0], 
                        [0, -1]], dtype=np.complex128)
IDENTITY_2x2 = np.array([[1, 0], 
                         [0, 1]], dtype=np.complex128)

# 2.2 Dirac Matrices (Standard Representation, 4x4)
# Alpha matrices are block-diagonal with Pauli matrices
ALPHA_X = np.block([[np.zeros((2, 2), dtype=np.complex128), SIGMA_X_2x2],
                    [SIGMA_X_2x2, np.zeros((2, 2), dtype=np.complex128)]])

ALPHA_Y = np.block([[np.zeros((2, 2), dtype=np.complex128), SIGMA_Y_2x2],
                    [SIGMA_Y_2x2, np.zeros((2, 2), dtype=np.complex128)]])

ALPHA_Z = np.block([[np.zeros((2, 2), dtype=np.complex128), SIGMA_Z_2x2],
                    [SIGMA_Z_2x2, np.zeros((2, 2), dtype=np.complex128)]])

# Beta matrix (block diagonal)
BETA = np.block([[IDENTITY_2x2, np.zeros((2, 2), dtype=np.complex128)],
                 [np.zeros((2, 2), dtype=np.complex128), -IDENTITY_2x2]])

# 4x4 Identity matrix
IDENTITY_4x4 = np.eye(4, dtype=np.complex128)

# 2.3 Rotation Matrices for Operator Splitting
# These matrices transform the wave function into a basis where the streaming
# operator for a given direction is diagonal

# Z_MATRIX: For diagonalizing -ALPHA_Z (from Dellar et al., Eq. 10)
Z_MATRIX = (1.0 / np.sqrt(2)) * np.array([[0, -1, 0, 1],
                                           [1, 0, -1, 0],
                                           [0, 1, 0, 1],
                                           [1, 0, 1, 0]], dtype=np.complex128)
Z_INV_MATRIX = Z_MATRIX.conj().T  # For unitary matrices, inverse is Hermitian conjugate

# X_MATRIX: For diagonalizing -ALPHA_X (from Dellar et al., Eq. 15)
X_MATRIX = (1.0 / np.sqrt(2)) * np.array([[-1, 0, 1, 0],
                                           [0, 1, 0, -1],
                                           [1, 0, 1, 0],
                                           [0, 1, 0, 1]], dtype=np.complex128)
X_INV_MATRIX = X_MATRIX.conj().T

# Y_MATRIX: For the Majorana form, the y-streaming matrix is already diagonal
# We use Identity for Y
Y_MATRIX = IDENTITY_4x4.copy()
Y_INV_MATRIX = IDENTITY_4x4.copy()

# Store matrices in a dictionary for easier access
MATRICES = {
    'alpha_x': ALPHA_X, 
    'alpha_y': ALPHA_Y, 
    'alpha_z': ALPHA_Z,
    'beta': BETA, 
    'identity': IDENTITY_4x4,
    'X': X_MATRIX, 
    'X_inv': X_INV_MATRIX,
    'Y': Y_MATRIX, 
    'Y_inv': Y_INV_MATRIX,
    'Z': Z_MATRIX, 
    'Z_inv': Z_INV_MATRIX
}

print("\nDirac matrices and rotation matrices initialized.")

# =============================================================================
# Part 3: Wave Function and Potential Initialization
# =============================================================================

# Initialize the 4-component wave function (psi)
# Shape: (NX, NY, NZ, 4) - each grid point holds a 4-component complex vector
psi = np.zeros((NX, NY, NZ, 4), dtype=np.complex128)

# Initialize the scalar potential field V(x,y,z)
V_field = np.zeros((NX, NY, NZ), dtype=np.float64)

# 3.1 Initial Condition: 3D Gaussian Wave Packet
# Parameters for Gaussian wave packet
D0_initial = 8 * DX  # Initial spread of the wave packet
x0, y0, z0 = NX // 4, NY // 2, NZ // 2  # Center of the wave packet

# Initial momentum (px, py, pz). For propagating packet in x-direction:
initial_momentum_x = 0.5 * HBAR / D0_initial
initial_momentum_y = 0.0
initial_momentum_z = 0.0

for i in range(NX):
    for j in range(NY):
        for k in range(NZ):
            # Coordinates relative to center
            rx = (i - x0) * DX
            ry = (j - y0) * DY
            rz = (k - z0) * DZ

            # Gaussian envelope
            gaussian_envelope = np.exp(-(rx**2 + ry**2 + rz**2) / (4 * D0_initial**2))

            # Phase factor for initial momentum
            phase_factor = np.exp(1j * (initial_momentum_x * rx +
                                        initial_momentum_y * ry +
                                        initial_momentum_z * rz) / HBAR)

            # Initialize first component of psi (phi_plus_1)
            psi[i, j, k, 0] = (1 / (np.sqrt(2 * np.pi) * D0_initial)**(3/2)) * gaussian_envelope * phase_factor

# Normalize the wave function: Integral(|psi|^2) = 1
norm_factor = np.sum(np.abs(psi)**2) * DX * DY * DZ
if norm_factor > 0:
    psi /= np.sqrt(norm_factor)

print(f"\nWave function initialized:")
print(f"  Initial spread: {D0_initial:.2e} m")
print(f"  Center position: ({x0}, {y0}, {z0})")
print(f"  Initial momentum (x): {initial_momentum_x:.2e} kg·m/s")

# =============================================================================
# Part 4: Core QLB 1D Step Function
# =============================================================================

def perform_qlb_sub_step(psi_line_in, V_line_in, current_axis, matrices_dict, 
                         hbar, c_light, m_particle, dt, dx):
    """
    Performs one sub-step (Rotate -> Collide -> Stream -> Rotate Back) for a 1D line of psi.
    
    Parameters:
    -----------
    psi_line_in : ndarray, shape (N_points, 4)
        The 4-component wave function along one line
    V_line_in : ndarray, shape (N_points,)
        The scalar potential along the line
    current_axis : str
        'x', 'y', or 'z'
    matrices_dict : dict
        Dictionary containing rotation matrices
    hbar, c_light, m_particle, dt, dx : float
        Physical and numerical parameters
        
    Returns:
    --------
    psi_line_out : ndarray, shape (N_points, 4)
        Updated wave function after QLB step
    """
    N_points = psi_line_in.shape[0]
    psi_line_out = np.zeros_like(psi_line_in, dtype=np.complex128)

    # Determine rotation matrices for the current axis
    if current_axis == 'x':
        R_inv = matrices_dict['X_inv']
        R_op = matrices_dict['X']
    elif current_axis == 'y':
        R_inv = matrices_dict['Y_inv']  # Y is Identity
        R_op = matrices_dict['Y']
    elif current_axis == 'z':
        R_inv = matrices_dict['Z_inv']
        R_op = matrices_dict['Z']
    else:
        raise ValueError("Invalid axis. Must be 'x', 'y', or 'z'.")

    # 4.1 Calculate Collision Coefficients (a_hat, b_hat) for each point
    # These depend on local potential and effective Compton frequency
    # Factor of 1/3 accounts for operator splitting in 3D
    m_tilde_3 = m_particle * c_light**2 * dt / (3 * hbar) if hbar > 0 else 0.0
    g_tilde_3 = Q_ELECTRON * V_line_in * dt / (3 * hbar) if hbar > 0 else 0.0

    Omega_3 = m_tilde_3**2 - g_tilde_3**2
    
    # Collision coefficients (position-dependent)
    denominator = 1 + Omega_3/4 - 1j * g_tilde_3
    # Avoid division by zero
    denominator = np.where(np.abs(denominator) < 1e-15, 1e-15, denominator)
    
    a_hat = (1 - Omega_3/4) / denominator
    b_hat = m_tilde_3 / denominator

    # 4.2 Perform Rotate -> Collide -> Stream -> Rotate Back
    # Create a buffer for the rotated wave function
    psi_rotated = np.zeros_like(psi_line_in, dtype=np.complex128)
    
    # Rotate all points
    for k in range(N_points):
        psi_rotated[k, :] = R_inv @ psi_line_in[k, :]

    # Collide and Stream in the rotated basis
    # In the rotated basis: (u1, u2) propagate forward, (d1, d2) propagate backward
    psi_collided_rotated = np.zeros_like(psi_rotated, dtype=np.complex128)
    
    for k in range(N_points):
        u1, u2, d1, d2 = psi_rotated[k, 0], psi_rotated[k, 1], psi_rotated[k, 2], psi_rotated[k, 3]
        
        # Apply collision (mixing of components)
        collided_u1 = a_hat[k] * u1 + b_hat[k] * d2
        collided_u2 = a_hat[k] * u2 + b_hat[k] * d1
        collided_d1 = a_hat[k] * d1 - b_hat[k] * u2
        collided_d2 = a_hat[k] * d2 - b_hat[k] * u1
        
        # Stream: u components go forward (to k+1), d components go backward (to k-1)
        # Forward propagation (u1, u2)
        k_next = (k + 1) % N_points  # Periodic boundary conditions
        psi_collided_rotated[k_next, 0] += collided_u1
        psi_collided_rotated[k_next, 1] += collided_u2
        
        # Backward propagation (d1, d2)
        k_prev = (k - 1) % N_points  # Periodic boundary conditions
        psi_collided_rotated[k_prev, 2] += collided_d1
        psi_collided_rotated[k_prev, 3] += collided_d2

    # Rotate back to standard basis
    for k in range(N_points):
        psi_line_out[k, :] = R_op @ psi_collided_rotated[k, :]

    return psi_line_out


print("\nCore QLB 1D step function defined.")

# =============================================================================
# Part 5: Main Simulation Loop
# =============================================================================

def run_simulation(psi_initial, V_field_input, output_dir_name="qlb_simulation_output"):
    """
    Run the main QLB simulation loop.
    
    Parameters:
    -----------
    psi_initial : ndarray
        Initial wave function
    V_field_input : ndarray
        Potential field
    output_dir_name : str
        Directory name for output files
    """
    # Create output directory
    os.makedirs(output_dir_name, exist_ok=True)

    # Initialize current wave function
    psi_current = psi_initial.copy()

    # Storage for observables
    densities_over_time = []
    transmission_coefficients = []

    print("\n" + "=" * 70)
    print("Starting QLB Simulation")
    print("=" * 70)

    for t_step in range(T_STEPS):
        psi_next = psi_current.copy()

        # Operator splitting: Apply QLB step for x, then y, then z directions
        
        # 5.1 X-Direction Step
        for j in range(NY):
            for k in range(NZ):
                psi_line_in = psi_current[:, j, k, :]
                V_line_in = V_field_input[:, j, k]
                psi_next[:, j, k, :] = perform_qlb_sub_step(
                    psi_line_in, V_line_in, 'x', MATRICES, HBAR, C, M_PARTICLE, DT, DX
                )
        psi_current = psi_next.copy()

        # 5.2 Y-Direction Step
        for i in range(NX):
            for k in range(NZ):
                psi_line_in = psi_current[i, :, k, :]
                V_line_in = V_field_input[i, :, k]
                psi_next[i, :, k, :] = perform_qlb_sub_step(
                    psi_line_in, V_line_in, 'y', MATRICES, HBAR, C, M_PARTICLE, DT, DY
                )
        psi_current = psi_next.copy()

        # 5.3 Z-Direction Step (only if NZ > 1)
        if NZ > 1:
            for i in range(NX):
                for j in range(NY):
                    psi_line_in = psi_current[i, j, :, :]
                    V_line_in = V_field_input[i, j, :]
                    psi_next[i, j, :, :] = perform_qlb_sub_step(
                        psi_line_in, V_line_in, 'z', MATRICES, HBAR, C, M_PARTICLE, DT, DZ
                    )
            psi_current = psi_next.copy()

        # 5.4 Compute Observables
        total_density = np.sum(np.abs(psi_current)**2) * DX * DY * DZ
        densities_over_time.append(total_density)

        # Calculate current density in x-direction for transmission
        Jx_field = np.zeros((NX, NY, NZ), dtype=np.float64)
        for i in range(NX):
            for j in range(NY):
                for k in range(NZ):
                    psi_vec = psi_current[i, j, k, :]
                    # J_x = Re(psi† · c · α_x · psi)
                    Jx_field[i, j, k] = np.real(np.vdot(psi_vec, C * MATRICES['alpha_x'] @ psi_vec))
        
        # Transmission coefficient: integrate current at outlet plane
        outlet_plane_x = NX - 2
        transmitted_current = np.sum(Jx_field[outlet_plane_x, :, :]) * DY * DZ
        transmission_coefficients.append(transmitted_current)

        # Output and visualization
        if (t_step + 1) % OUTPUT_FREQ == 0:
            print(f"\nTime step {t_step + 1}/{T_STEPS}:")
            print(f"  Total Density: {total_density:.6f}")
            print(f"  Transmitted Current: {transmitted_current:.4e}")
            
            # Visualize probability density for 2D case (NZ=1)
            if NZ == 1:
                prob_density_slice = np.sum(np.abs(psi_current[:, :, 0, :])**2, axis=-1)
                
                plt.figure(figsize=(10, 8))
                plt.imshow(prob_density_slice.T, origin='lower', cmap='hot',
                           extent=[0, LX*1e9, 0, LY*1e9], aspect='auto')
                plt.colorbar(label='Probability Density')
                plt.title(f'Probability Density at t = {(t_step + 1) * DT:.2e} s')
                plt.xlabel('x (nm)')
                plt.ylabel('y (nm)')
                plt.tight_layout()
                plt.savefig(os.path.join(output_dir_name, f'prob_density_xy_t{t_step+1:04d}.png'), dpi=150)
                plt.close()

    print("\n" + "=" * 70)
    print("Simulation Complete")
    print("=" * 70)

    # Plot final observables
    time_array = np.arange(T_STEPS) * DT

    # Total density over time
    plt.figure(figsize=(12, 5))
    
    plt.subplot(1, 2, 1)
    plt.plot(time_array * 1e15, densities_over_time, 'b-', linewidth=2)
    plt.xlabel('Time (fs)')
    plt.ylabel('Total Density')
    plt.title('Total Density vs Time')
    plt.grid(True, alpha=0.3)
    
    plt.subplot(1, 2, 2)
    plt.plot(time_array * 1e15, transmission_coefficients, 'r-', linewidth=2)
    plt.xlabel('Time (fs)')
    plt.ylabel('Transmitted Current (arb. units)')
    plt.title('Transmission Coefficient vs Time')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir_name, 'observables_summary.png'), dpi=150)
    plt.close()

    return psi_current, densities_over_time, transmission_coefficients


# =============================================================================
# Part 6: Adapting for Electrons in Graphene (2D Tunneling Example)
# =============================================================================

def setup_graphene_simulation():
    """
    Configure parameters and initial conditions for graphene electron transport.
    """
    global M_PARTICLE, OMEGA_C, psi, V_field, NZ, DZ
    
    print("\n" + "=" * 70)
    print("Configuring for Graphene Simulation")
    print("=" * 70)
    
    # 6.1 Set to 2D
    NZ = 1
    DZ = DX
    
    # 6.2 Massless Dirac fermions in graphene
    M_PARTICLE = 0.0  # Massless for ideal graphene
    OMEGA_C = 0.0
    
    print(f"Particle Mass set to: {M_PARTICLE} kg (massless Dirac fermions)")
    
    # 6.3 Potential Field: Random Impurities
    impurity_concentration = 0.05  # 5% impurities
    impurity_potential_strength = 0.1 * Q_ELECTRON  # 0.1 eV barrier
    
    # Define impurity region
    impurity_region_start_x = NX // 3
    impurity_region_end_x = 2 * NX // 3
    
    # Reset V_field
    V_field = np.zeros((NX, NY, NZ), dtype=np.float64)
    
    np.random.seed(42)  # For reproducibility
    for i in range(impurity_region_start_x, impurity_region_end_x):
        for j in range(NY):
            for k in range(NZ):
                if np.random.rand() < impurity_concentration:
                    V_field[i, j, k] = impurity_potential_strength
    
    print(f"\nImpurity region: x ∈ [{impurity_region_start_x}, {impurity_region_end_x}]")
    print(f"Impurity concentration: {impurity_concentration * 100}%")
    print(f"Impurity potential: {impurity_potential_strength / Q_ELECTRON:.2f} eV")
    
    # 6.4 Initial Condition: Wave Packet for Graphene
    D0_graphene = 8 * DX
    initial_momentum_x_graphene = 1.0 * HBAR / D0_graphene
    initial_momentum_y_graphene = 0.0
    initial_momentum_z_graphene = 0.0
    
    # Reset psi
    psi = np.zeros((NX, NY, NZ, 4), dtype=np.complex128)
    
    # Position wave packet at left edge
    x0_graphene, y0_graphene, z0_graphene = NX // 8, NY // 2, NZ // 2
    
    for i in range(NX):
        for j in range(NY):
            for k in range(NZ):
                rx = (i - x0_graphene) * DX
                ry = (j - y0_graphene) * DY
                rz = (k - z0_graphene) * DZ
                
                gaussian_envelope = np.exp(-(rx**2 + ry**2 + rz**2) / (4 * D0_graphene**2))
                phase_factor = np.exp(1j * (initial_momentum_x_graphene * rx +
                                            initial_momentum_y_graphene * ry +
                                            initial_momentum_z_graphene * rz) / HBAR)
                
                # Initialize wave packet (primarily in first component)
                psi[i, j, k, 0] = (1 / (np.sqrt(2 * np.pi) * D0_graphene)**(3/2)) * gaussian_envelope * phase_factor
    
    # Normalize
    norm_factor = np.sum(np.abs(psi)**2) * DX * DY * DZ
    if norm_factor > 0:
        psi /= np.sqrt(norm_factor)
    
    print(f"\nGraphene wave packet initialized:")
    print(f"  Initial position: ({x0_graphene}, {y0_graphene})")
    print(f"  Initial spread: {D0_graphene:.2e} m")
    print(f"  Initial momentum: {initial_momentum_x_graphene:.2e} kg·m/s")
    print("=" * 70)
    
    return psi, V_field


# =============================================================================
# Main Execution
# =============================================================================

if __name__ == "__main__":
    print("\n\n")
    print("╔" + "═" * 68 + "╗")
    print("║" + " " * 15 + "3D DIRAC QLB SOLVER" + " " * 34 + "║")
    print("║" + " " * 10 + "Quantum Lattice Boltzmann Method" + " " * 26 + "║")
    print("╚" + "═" * 68 + "╝")
    
    # Choose simulation type
    print("\nSimulation Options:")
    print("  1. Standard 3D Dirac equation (massive particle)")
    print("  2. Graphene electron transport (massless, 2D with impurities)")
    
    # For automated execution, default to graphene simulation
    simulation_choice = 2
    
    if simulation_choice == 2:
        # Setup and run graphene simulation
        psi_graphene, V_graphene = setup_graphene_simulation()
        psi_final, densities, transmission = run_simulation(
            psi_graphene, V_graphene, 
            output_dir_name="qlb_graphene_output"
        )
        print(f"\nGraphene simulation results saved to: qlb_graphene_output/")
    else:
        # Run standard simulation with initialized parameters
        psi_final, densities, transmission = run_simulation(
            psi, V_field,
            output_dir_name="qlb_simulation_output"
        )
        print(f"\nSimulation results saved to: qlb_simulation_output/")
    
    print("\n" + "=" * 70)
    print("All simulations complete!")
    print("=" * 70)
