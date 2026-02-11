"""
Palpacelli 2012 Configuration - Klein Tunneling Study
======================================================
This script configures the 3D Dirac QLB solver to match the parameters from:
"Klein Tunneling in the Presence of Random Impurities"
Palpacelli et al., Int. J. Mod. Phys. C, Vol. 23, No. 12 (2012)

Key parameters from Section 4 of the paper for random media simulations.
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from dirac_qlb_solver import perform_qlb_sub_step, HBAR, C, Q_ELECTRON, IDENTITY_4x4

# =============================================================================
# PALPACELLI 2012 PARAMETERS - Section 4: Klein Paradox in Random Media
# =============================================================================

# Physical constants (use from main solver)
# HBAR = 1.0545718e-34 J·s
# C = 2.99792458e8 m/s
# Q_ELECTRON = 1.60217663e-19 C

# Grid Configuration (from paper)
NX = 2048  # Total domain length
NY = 512   # Domain height
NZ = 1     # 2D simulation

# Lattice spacing (from paper: Δx = 0.96 nm)
DX = 0.96e-9  # m
DY = DX
DZ = DX

# Physical domain size
LX = NX * DX  # ≈ 1.97 μm
LY = NY * DY  # ≈ 0.49 μm
LZ = DZ

# Time step (light-cone condition: c·Δt = Δx)
DT = DX / C  # ≈ 3.2e-18 s

# Energy parameters (from paper)
# For lattice units: ΔE = ℏ/Δt
DELTA_E = HBAR / DT  # Energy unit ≈ 0.33 eV
# Paper uses E = 0.117 in lattice units ≈ 80 meV
ENERGY_LATTICE = 0.117
ENERGY_PHYSICAL = 80e-3 * Q_ELECTRON  # 80 meV in Joules

# Wave packet parameters (from paper: σ = 48 lattice spacings)
SIGMA_LATTICE = 48  # cells
SIGMA_PHYSICAL = SIGMA_LATTICE * DX  # ≈ 46.08 nm
D0_INITIAL = SIGMA_PHYSICAL

# Initial momentum (from paper: k₀ ≈ 1 in lattice units)
K0_LATTICE = 1.0  # lattice units (1/Δx)
K0_PHYSICAL = K0_LATTICE / DX  # 1/m
INITIAL_MOMENTUM_X = HBAR * K0_PHYSICAL  # kg·m/s

# Particle mass (massless case for graphene)
M_PARTICLE = 0.0  # kg
OMEGA_C = 0.0

# Domain regions (from paper)
INLET_START = 0
INLET_END = 512
IMPURITY_START = 512
IMPURITY_END = 1536
OUTLET_START = 1536
OUTLET_END = 2048

# Impurity configuration (from paper)
IMPURITY_SIZE_CELLS = 8  # d = 8 lattice spacings
IMPURITY_SIZE_PHYSICAL = IMPURITY_SIZE_CELLS * DX  # ≈ 7.68 nm

# Impurity concentrations tested in paper: 0.1%, 0.5%, 1%, 5%
IMPURITY_CONCENTRATIONS = [0.001, 0.005, 0.01, 0.05]
DEFAULT_CONCENTRATION = 0.05  # 5%

# Barrier heights tested in paper: 25, 50, 100, 200, 285 meV
BARRIER_HEIGHTS_MEV = [25, 50, 100, 200, 285]
DEFAULT_BARRIER_HEIGHT = 100e-3 * Q_ELECTRON  # 100 meV in Joules

# Total simulation steps (needs to be long enough for wave packet to traverse domain)
# At c·Δt = Δx, need ~2048 steps to cross domain, plus settling time
T_STEPS = 6000

# Output frequency
OUTPUT_FREQ = 500

print("=" * 80)
print("PALPACELLI 2012 CONFIGURATION")
print("Klein Tunneling in the Presence of Random Impurities")
print("=" * 80)
print(f"\nGrid Configuration:")
print(f"  Grid size:           {NX} × {NY} × {NZ}")
print(f"  Domain size:         {LX*1e6:.3f} × {LY*1e6:.3f} μm")
print(f"  Lattice spacing:     {DX*1e9:.3f} nm")
print(f"  Time step:           {DT:.3e} s")
print(f"\nPhysical Parameters:")
print(f"  Particle mass:       {M_PARTICLE} kg (massless)")
print(f"  Energy:              {ENERGY_PHYSICAL/Q_ELECTRON*1e3:.1f} meV")
print(f"  Initial momentum:    {INITIAL_MOMENTUM_X:.3e} kg·m/s")
print(f"\nWave Packet:")
print(f"  Spread σ:            {SIGMA_PHYSICAL*1e9:.2f} nm ({SIGMA_LATTICE} cells)")
print(f"  Initial position:    Inlet region")
print(f"\nDomain Regions:")
print(f"  Inlet:               [{INLET_START}:{INLET_END}] = {INLET_END-INLET_START} cells")
print(f"  Impurities:          [{IMPURITY_START}:{IMPURITY_END}] = {IMPURITY_END-IMPURITY_START} cells")
print(f"  Outlet:              [{OUTLET_START}:{OUTLET_END}] = {OUTLET_END-OUTLET_START} cells")
print(f"\nImpurity Configuration:")
print(f"  Impurity size:       {IMPURITY_SIZE_PHYSICAL*1e9:.2f} nm ({IMPURITY_SIZE_CELLS} cells)")
print(f"  Default concentration: {DEFAULT_CONCENTRATION*100:.1f}%")
print(f"  Default barrier:     {DEFAULT_BARRIER_HEIGHT/Q_ELECTRON*1e3:.0f} meV")
print(f"\nSimulation:")
print(f"  Total steps:         {T_STEPS}")
print(f"  Output frequency:    {OUTPUT_FREQ}")
print("=" * 80)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def calculate_impurity_count(concentration, domain_area, impurity_area):
    """
    Calculate number of impurities for given concentration.
    
    Concentration C = N·d²/A where:
    - N = number of impurities
    - d² = impurity area
    - A = total area
    """
    N = int(concentration * domain_area / impurity_area)
    return N

def generate_random_impurities(concentration, barrier_height, seed=42):
    """
    Generate random impurity distribution matching Palpacelli 2012.
    
    Parameters:
    -----------
    concentration : float
        Impurity concentration (0.001 to 0.05)
    barrier_height : float
        Potential barrier height in Joules
    seed : int
        Random seed for reproducibility
        
    Returns:
    --------
    V_field : ndarray
        Potential field with random square impurities
    impurity_positions : list
        List of (x, y) positions of impurities
    """
    V_field = np.zeros((NX, NY, NZ), dtype=np.float64)
    
    # Calculate impurity region area
    impurity_region_length = IMPURITY_END - IMPURITY_START
    impurity_region_area = impurity_region_length * NY
    impurity_area = IMPURITY_SIZE_CELLS * IMPURITY_SIZE_CELLS
    
    # Calculate number of impurities
    N_impurities = calculate_impurity_count(concentration, impurity_region_area, impurity_area)
    
    print(f"\nGenerating random impurity distribution:")
    print(f"  Concentration:       {concentration*100:.2f}%")
    print(f"  Number of impurities: {N_impurities}")
    print(f"  Barrier height:      {barrier_height/Q_ELECTRON*1e3:.1f} meV")
    
    np.random.seed(seed)
    impurity_positions = []
    
    # Randomly place square impurities in the impurity region
    for _ in range(N_impurities):
        # Random position for top-left corner of impurity
        ix = np.random.randint(IMPURITY_START, IMPURITY_END - IMPURITY_SIZE_CELLS)
        iy = np.random.randint(0, NY - IMPURITY_SIZE_CELLS)
        
        # Fill square impurity region
        for di in range(IMPURITY_SIZE_CELLS):
            for dj in range(IMPURITY_SIZE_CELLS):
                V_field[ix + di, iy + dj, 0] = barrier_height
        
        impurity_positions.append((ix, iy))
    
    return V_field, impurity_positions

def initialize_gaussian_wave_packet():
    """
    Initialize Gaussian wave packet in inlet region matching Palpacelli 2012.
    
    Returns:
    --------
    psi : ndarray
        Initial wave function
    """
    psi = np.zeros((NX, NY, NZ, 4), dtype=np.complex128)
    
    # Position wave packet at center of inlet region
    x0 = INLET_END // 2  # Middle of inlet
    y0 = NY // 2  # Centered vertically
    z0 = 0
    
    # Work in lattice units to avoid numerical overflow
    # Normalize Gaussian such that peak amplitude is O(1)
    sigma_cells = SIGMA_LATTICE  # In lattice units
    
    for i in range(NX):
        for j in range(NY):
            for k in range(NZ):
                # Coordinates in lattice units (cells)
                dx_cells = i - x0
                dy_cells = j - y0
                r_sq_cells = dx_cells**2 + dy_cells**2
                
                # Gaussian envelope in lattice units
                gaussian_envelope = np.exp(-r_sq_cells / (4 * sigma_cells**2))
                
                # Phase factor for initial momentum (k0 = 1 in lattice units)
                phase_factor = np.exp(1j * K0_LATTICE * dx_cells)
                
                # Initialize first component (positive energy, spin up)
                psi[i, j, k, 0] = gaussian_envelope * phase_factor
    
    # Normalize: ∫|ψ|² d³r = 1 (in lattice units: sum over cells)
    # For lattice units, the "volume element" is just 1 (counting cells)
    prob_before = np.sum(np.abs(psi)**2)
    if prob_before > 0:
        psi /= np.sqrt(prob_before)
    
    # Verify normalization
    prob_after = np.sum(np.abs(psi)**2)
    
    print(f"\nWave packet initialized:")
    print(f"  Position:            ({x0}, {y0}, {z0})")
    print(f"  Spread:              {SIGMA_PHYSICAL*1e9:.2f} nm ({sigma_cells} cells)")
    print(f"  Momentum (k):        {K0_LATTICE} (lattice units)")
    print(f"  Norm (before):       {prob_before:.6f}")
    print(f"  Norm (after):        {prob_after:.6f}")
    
    return psi

def visualize_setup(V_field, psi, output_dir="palpacelli2012_setup"):
    """
    Visualize the initial setup.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Plot potential field
    plt.figure(figsize=(16, 6))
    
    plt.subplot(1, 2, 1)
    V_slice = V_field[:, :, 0].T / Q_ELECTRON * 1e3  # Convert to meV
    plt.imshow(V_slice, origin='lower', cmap='viridis', aspect='auto',
               extent=[0, LX*1e9, 0, LY*1e9])
    plt.colorbar(label='Potential (meV)')
    plt.axvline(x=IMPURITY_START*DX*1e9, color='red', linestyle='--', linewidth=2, label='Impurity region')
    plt.axvline(x=IMPURITY_END*DX*1e9, color='red', linestyle='--', linewidth=2)
    plt.xlabel('x (nm)')
    plt.ylabel('y (nm)')
    plt.title('Potential Field with Random Impurities')
    plt.legend()
    
    # Plot initial wave packet
    plt.subplot(1, 2, 2)
    prob_density = np.sum(np.abs(psi[:, :, 0, :])**2, axis=-1).T
    plt.imshow(prob_density, origin='lower', cmap='hot', aspect='auto',
               extent=[0, LX*1e9, 0, LY*1e9])
    plt.colorbar(label='Probability Density')
    plt.axvline(x=IMPURITY_START*DX*1e9, color='cyan', linestyle='--', linewidth=2, label='Impurity region')
    plt.axvline(x=IMPURITY_END*DX*1e9, color='cyan', linestyle='--', linewidth=2)
    plt.xlabel('x (nm)')
    plt.ylabel('y (nm)')
    plt.title('Initial Wave Packet')
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'initial_setup.png'), dpi=150)
    plt.close()
    
    print(f"\nSetup visualization saved to: {output_dir}/initial_setup.png")

# =============================================================================
# EXAMPLE USAGE
# =============================================================================

if __name__ == "__main__":
    print("\n" + "=" * 80)
    print("EXAMPLE: Setting up Palpacelli 2012 simulation")
    print("=" * 80)
    
    # Generate random impurities with 5% concentration and 100 meV barrier
    V_field, impurity_positions = generate_random_impurities(
        concentration=DEFAULT_CONCENTRATION,
        barrier_height=DEFAULT_BARRIER_HEIGHT,
        seed=42
    )
    
    # Initialize wave packet
    psi = initialize_gaussian_wave_packet()
    
    # Visualize setup
    visualize_setup(V_field, psi)
    
    print("\n" + "=" * 80)
    print("Setup complete!")
    print("=" * 80)
    print(f"\nTo run the simulation, use the returned psi and V_field with")
    print(f"the QLB solver adapted for the larger grid size.")
    print(f"\nNote: This configuration requires ~2048×512×4×16 bytes ≈ 67 MB memory")
    print(f"for the wave function alone.")
