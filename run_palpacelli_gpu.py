"""
GPU-Accelerated Palpacelli 2012 Simulation with Reduced Resolution
===================================================================
Uses 1/4 resolution (512×128 instead of 2048×512) and CuPy for GPU acceleration.
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import time

# Try to import cupy, fall back to numpy if not available
try:
    import cupy as cp
    GPU_AVAILABLE = True
    print("GPU (CuPy) detected and will be used for acceleration!")
except ImportError:
    print("CuPy not available, falling back to NumPy (CPU)")
    cp = np
    GPU_AVAILABLE = False

from dirac_qlb_solver import (
    ALPHA_X, ALPHA_Y, ALPHA_Z, BETA, IDENTITY_4x4,
    X_MATRIX, Y_MATRIX, Z_MATRIX,
    X_INV_MATRIX, Y_INV_MATRIX, Z_INV_MATRIX,
    HBAR, C, Q_ELECTRON
)

# =============================================================================
# REDUCED RESOLUTION CONFIGURATION (1/4 of paper)
# =============================================================================

# Grid: 1/4 resolution
NX = 512   # Paper: 2048
NY = 128   # Paper: 512
NZ = 1

# Same physical lattice spacing as paper
DX = 0.96e-9  # m
DY = DX
DZ = DX

# Physical domain (smaller due to fewer cells)
LX = NX * DX  # ~0.49 μm (vs paper's 1.97 μm)
LY = NY * DY  # ~0.12 μm (vs paper's 0.49 μm)
LZ = DZ

# Time step (light-cone condition)
DT = DX / C

# Particle mass (massless for graphene)
M_PARTICLE = 0.0
OMEGA_C = 0.0

# Domain regions (scaled to 1/4)
INLET_START = 0
INLET_END = 128      # Paper: 512
IMPURITY_START = 128  # Paper: 512
IMPURITY_END = 384    # Paper: 1536
OUTLET_START = 384    # Paper: 1536
OUTLET_END = 512      # Paper: 2048

# Wave packet (scaled)
SIGMA_LATTICE = 12    # Paper: 48 (1/4)
K0_LATTICE = 1.0

# Impurity configuration (same as paper)
IMPURITY_SIZE_CELLS = 2  # Paper: 8 (1/4)
DEFAULT_CONCENTRATION = 0.05  # 5%
DEFAULT_BARRIER_HEIGHT = 100e-3 * Q_ELECTRON  # 100 meV

print("\n" + "="*80)
print("GPU-ACCELERATED PALPACELLI 2012 (1/4 RESOLUTION)")
print("="*80)
print(f"\nConfiguration:")
print(f"  Grid: {NX} × {NY} × {NZ} (paper: 2048 × 512)")
print(f"  Domain: {LX*1e9:.1f} × {LY*1e9:.1f} nm (paper: 1966 × 492 nm)")
print(f"  Lattice spacing: {DX*1e9:.3f} nm")
print(f"  Wave packet σ: {SIGMA_LATTICE} cells (paper: 48 cells)")
print(f"  Impurity size: {IMPURITY_SIZE_CELLS}×{IMPURITY_SIZE_CELLS} (paper: 8×8)")
print(f"  GPU acceleration: {'ENABLED' if GPU_AVAILABLE else 'DISABLED'}")
print("="*80)

# Convert matrices to GPU if available
if GPU_AVAILABLE:
    X_MATRIX_GPU = cp.asarray(X_MATRIX)
    Y_MATRIX_GPU = cp.asarray(Y_MATRIX)
    Z_MATRIX_GPU = cp.asarray(Z_MATRIX)
    X_INV_MATRIX_GPU = cp.asarray(X_INV_MATRIX)
    Y_INV_MATRIX_GPU = cp.asarray(Y_INV_MATRIX)
    Z_INV_MATRIX_GPU = cp.asarray(Z_INV_MATRIX)
    BETA_GPU = cp.asarray(BETA)
else:
    X_MATRIX_GPU = X_MATRIX
    Y_MATRIX_GPU = Y_MATRIX
    Z_MATRIX_GPU = Z_MATRIX
    X_INV_MATRIX_GPU = X_INV_MATRIX
    Y_INV_MATRIX_GPU = Y_INV_MATRIX
    Z_INV_MATRIX_GPU = Z_INV_MATRIX
    BETA_GPU = BETA

def generate_random_impurities_gpu(concentration, barrier_height, seed=42):
    """Generate random impurities (uses CPU numpy for random generation)."""
    V_field = np.zeros((NX, NY, NZ), dtype=np.float64)
    
    impurity_region_length = IMPURITY_END - IMPURITY_START
    impurity_region_area = impurity_region_length * NY
    impurity_area = IMPURITY_SIZE_CELLS * IMPURITY_SIZE_CELLS
    
    N_impurities = int(concentration * impurity_region_area / impurity_area)
    
    print(f"\nGenerating {N_impurities} random impurities...")
    
    np.random.seed(seed)
    impurity_positions = []
    
    for _ in range(N_impurities):
        ix = np.random.randint(IMPURITY_START, IMPURITY_END - IMPURITY_SIZE_CELLS)
        iy = np.random.randint(0, NY - IMPURITY_SIZE_CELLS)
        
        for di in range(IMPURITY_SIZE_CELLS):
            for dj in range(IMPURITY_SIZE_CELLS):
                V_field[ix + di, iy + dj, 0] = barrier_height
        
        impurity_positions.append((ix, iy))
    
    # Transfer to GPU if available
    if GPU_AVAILABLE:
        V_field = cp.asarray(V_field)
    
    return V_field, impurity_positions

def initialize_wave_packet_gpu():
    """Initialize Gaussian wave packet on GPU."""
    if GPU_AVAILABLE:
        psi = cp.zeros((NX, NY, NZ, 4), dtype=cp.complex128)
    else:
        psi = np.zeros((NX, NY, NZ, 4), dtype=np.complex128)
    
    x0 = INLET_END // 2
    y0 = NY // 2
    z0 = 0
    
    D0_initial = SIGMA_LATTICE * DX
    initial_momentum_x = HBAR * K0_LATTICE / DX
    
    # Create coordinate grids
    i_grid, j_grid, k_grid = np.meshgrid(
        np.arange(NX), np.arange(NY), np.arange(NZ), indexing='ij'
    )
    
    rx = (i_grid - x0) * DX
    ry = (j_grid - y0) * DY
    rz = (k_grid - z0) * DZ
    
    gaussian_envelope = np.exp(-(rx**2 + ry**2 + rz**2) / (4 * D0_initial**2))
    phase_factor = np.exp(1j * initial_momentum_x * rx / HBAR)
    
    # Initialize first component
    psi_cpu = np.zeros((NX, NY, NZ, 4), dtype=np.complex128)
    psi_cpu[:, :, :, 0] = (1 / (np.sqrt(2 * np.pi) * D0_initial)**(3/2)) * gaussian_envelope * phase_factor
    
    # Normalize
    norm_factor = np.sum(np.abs(psi_cpu)**2) * DX * DY * DZ
    if norm_factor > 0:
        psi_cpu /= np.sqrt(norm_factor)
    
    prob_check = np.sum(np.abs(psi_cpu)**2) * DX * DY * DZ
    
    print(f"\nWave packet initialized:")
    print(f"  Position: ({x0}, {y0}, {z0})")
    print(f"  Spread: {D0_initial*1e9:.2f} nm")
    print(f"  Probability: {prob_check:.6f}")
    
    # Transfer to GPU
    if GPU_AVAILABLE:
        psi = cp.asarray(psi_cpu)
    else:
        psi = psi_cpu
    
    return psi

def perform_qlb_line_gpu(psi_line, V_line, axis, R_inv, R_op):
    """
    Perform QLB step on a line using GPU acceleration.
    Uses the same algorithm as existing solver but vectorized for GPU.
    Boundary conditions are handled externally before calling this function.
    """
    N_points = psi_line.shape[0]
    xp = cp if GPU_AVAILABLE else np
    
    # Collision coefficients (1/3 for operator splitting)
    m_tilde_3 = M_PARTICLE * C**2 * DT / (3 * HBAR) if HBAR > 0 else 0.0
    g_tilde_3 = Q_ELECTRON * V_line * DT / (3 * HBAR) if HBAR > 0 else xp.zeros_like(V_line)
    
    Omega_3 = m_tilde_3**2 - g_tilde_3**2
    denominator = 1 + Omega_3/4 - 1j * g_tilde_3
    denominator = xp.where(xp.abs(denominator) < 1e-15, 1e-15, denominator)
    
    a_hat = (1 - Omega_3/4) / denominator
    b_hat = m_tilde_3 / denominator
    
    # Rotate to characteristic basis
    psi_rotated = xp.einsum('ij,kj->ki', R_inv, psi_line)
    
    # Collide
    u1, u2, d1, d2 = psi_rotated[:, 0], psi_rotated[:, 1], psi_rotated[:, 2], psi_rotated[:, 3]
    
    collided_u1 = a_hat * u1 + b_hat * d2
    collided_u2 = a_hat * u2 + b_hat * d1
    collided_d1 = a_hat * d1 - b_hat * u2
    collided_d2 = a_hat * d2 - b_hat * u1
    
    # Stream in characteristic basis (periodic for interior points)
    psi_collided_rotated = xp.zeros_like(psi_rotated)
    
    # Interior streaming
    psi_collided_rotated[1:, 0] = collided_u1[:-1]      # u1 shifts right
    psi_collided_rotated[1:, 1] = collided_u2[:-1]      # u2 shifts right
    psi_collided_rotated[:-1, 2] = collided_d1[1:]      # d1 shifts left
    psi_collided_rotated[:-1, 3] = collided_d2[1:]      # d2 shifts left
    
    # Note: Boundary conditions are applied externally after X-sweep
    # to avoid interfering with Y-direction propagation
    
    # Rotate back
    psi_out = xp.einsum('ij,kj->ki', R_op, psi_collided_rotated)
    
    return psi_out

def run_simulation_gpu(n_steps=1000, output_freq=100, output_dir="palpacelli_gpu", save_snapshots=True):
    """Run GPU-accelerated simulation."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Create snapshots subdirectory if saving snapshots
    if save_snapshots:
        snapshots_dir = os.path.join(output_dir, 'snapshots')
        os.makedirs(snapshots_dir, exist_ok=True)
    
    print("\n" + "="*80)
    print("STARTING GPU SIMULATION")
    print("="*80)
    
    # Generate impurities
    print("\n[1/4] Generating impurities...")
    V_field, impurity_positions = generate_random_impurities_gpu(
        DEFAULT_CONCENTRATION, DEFAULT_BARRIER_HEIGHT, seed=42
    )
    
    # Initialize wave packet
    print("\n[2/4] Initializing wave packet...")
    psi = initialize_wave_packet_gpu()
    
    # Get initial probability density for consistent colormap scaling
    if GPU_AVAILABLE:
        psi_cpu_temp = cp.asnumpy(psi)
    else:
        psi_cpu_temp = psi
    prob_density_initial = np.sum(np.abs(psi_cpu_temp[:, :, 0, :])**2, axis=-1)
    PROB_DENSITY_MAX = np.max(prob_density_initial) * 1.2  # 20% headroom
    
    # Storage (keep on CPU for efficiency)
    time_points = []
    total_probability = []
    inlet_probability = []
    impurity_probability = []
    outlet_probability = []
    
    print("\n[3/4] Running QLB simulation...")
    print(f"  Steps: {n_steps}, Output freq: {output_freq}")
    
    start_time = time.time()
    
    for step in range(n_steps):
        # X-direction sweep
        psi_new = psi.copy()
        for j in range(NY):
            psi_line = psi[:, j, 0, :]
            V_line = V_field[:, j, 0]
            psi_new[:, j, 0, :] = perform_qlb_line_gpu(
                psi_line, V_line, 'x', X_INV_MATRIX_GPU, X_MATRIX_GPU
            )
        psi = psi_new
        
        # Apply open/absorbing boundary conditions at x=0 and x=NX-1
        # in characteristic basis to avoid reflections
        xp = cp if GPU_AVAILABLE else np
        for j in range(NY):
            # Transform to characteristic basis
            psi_rot_inlet = xp.einsum('ij,j->i', X_INV_MATRIX_GPU, psi[0, j, 0, :])
            psi_rot_outlet = xp.einsum('ij,j->i', X_INV_MATRIX_GPU, psi[-1, j, 0, :])
            
            # Inlet (x=0): zero incoming forward waves (u1, u2)
            psi_rot_inlet[0] = 0.0  # u1
            psi_rot_inlet[1] = 0.0  # u2
            
            # Outlet (x=NX-1): zero incoming backward waves (d1, d2)
            psi_rot_outlet[2] = 0.0  # d1
            psi_rot_outlet[3] = 0.0  # d2
            
            # Transform back to physical basis
            psi[0, j, 0, :] = xp.einsum('ij,j->i', X_MATRIX_GPU, psi_rot_inlet)
            psi[-1, j, 0, :] = xp.einsum('ij,j->i', X_MATRIX_GPU, psi_rot_outlet)
        
        # Y-direction sweep
        psi_new = psi.copy()
        for i in range(NX):
            psi_line = psi[i, :, 0, :]
            V_line = V_field[i, :, 0]
            psi_new[i, :, 0, :] = perform_qlb_line_gpu(
                psi_line, V_line, 'y', Y_INV_MATRIX_GPU, Y_MATRIX_GPU
            )
        psi = psi_new
        
        # Re-apply boundary conditions after Y-sweep with exponential damping
        # Apply stronger absorption near boundaries to prevent reflections
        sponge_width = 10  # cells
        for j in range(NY):
            # Inlet (x=0): zero incoming forward waves
            psi_rot_inlet = xp.einsum('ij,j->i', X_INV_MATRIX_GPU, psi[0, j, 0, :])
            psi_rot_inlet[0] = 0.0  # u1
            psi_rot_inlet[1] = 0.0  # u2
            psi[0, j, 0, :] = xp.einsum('ij,j->i', X_MATRIX_GPU, psi_rot_inlet)
            
            # Outlet (x=NX-1): zero incoming backward waves
            psi_rot_outlet = xp.einsum('ij,j->i', X_INV_MATRIX_GPU, psi[-1, j, 0, :])
            psi_rot_outlet[2] = 0.0  # d1
            psi_rot_outlet[3] = 0.0  # d2
            psi[-1, j, 0, :] = xp.einsum('ij,j->i', X_MATRIX_GPU, psi_rot_outlet)
        
        # Apply exponential damping in sponge layer near outlet
        for i in range(max(0, NX - sponge_width), NX):
            distance_from_outlet = NX - 1 - i
            damping_factor = xp.exp(-0.5 * (sponge_width - distance_from_outlet) / sponge_width)
            for j in range(NY):
                psi[i, j, 0, :] *= damping_factor
        
        # Diagnostics
        if step % output_freq == 0:
            # Transfer to CPU for analysis
            if GPU_AVAILABLE:
                psi_cpu = cp.asnumpy(psi)
            else:
                psi_cpu = psi
            
            prob_density = np.sum(np.abs(psi_cpu)**2, axis=-1) * DX * DY * DZ
            
            total_prob = np.sum(prob_density)
            inlet_prob = np.sum(prob_density[INLET_START:INLET_END, :, :])
            impurity_prob = np.sum(prob_density[IMPURITY_START:IMPURITY_END, :, :])
            outlet_prob = np.sum(prob_density[OUTLET_START:OUTLET_END, :, :])
            
            time_points.append(step * DT)
            total_probability.append(total_prob)
            inlet_probability.append(inlet_prob)
            impurity_probability.append(impurity_prob)
            outlet_probability.append(outlet_prob)
            
            elapsed = time.time() - start_time
            time_per_step = elapsed / (step + 1)
            eta = time_per_step * (n_steps - step - 1)
            
            print(f"  Step {step:5d}/{n_steps} | "
                  f"P={total_prob:.4f} | "
                  f"Inlet={inlet_prob:.3f} | "
                  f"Impurity={impurity_prob:.3f} | "
                  f"Outlet={outlet_prob:.3f} | "
                  f"ETA: {eta:.0f}s")
            
            # Save snapshot
            if save_snapshots:
                # Get V_field on CPU
                if GPU_AVAILABLE:
                    V_field_cpu_local = cp.asnumpy(V_field)
                else:
                    V_field_cpu_local = V_field
                
                save_snapshot(psi_cpu, V_field_cpu_local, step, step * DT,
                            total_prob, inlet_prob, impurity_prob, outlet_prob,
                            snapshots_dir, PROB_DENSITY_MAX)
    
    total_time = time.time() - start_time
    print(f"\n  Completed in {total_time:.1f} seconds")
    print(f"  Performance: {total_time/n_steps:.4f} s/step, {n_steps/total_time:.1f} steps/s")
    
    # Transfer final state to CPU
    if GPU_AVAILABLE:
        psi_final = cp.asnumpy(psi)
        V_field_cpu = cp.asnumpy(V_field)
    else:
        psi_final = psi
        V_field_cpu = V_field
    
    # Plot
    print("\n[4/4] Creating plots...")
    plot_results(time_points, total_probability, inlet_probability,
                 impurity_probability, outlet_probability, psi_final, V_field_cpu, output_dir)
    
    # Create animation if snapshots were saved
    if save_snapshots:
        print("\n[5/5] Creating animation from snapshots...")
        create_animation(snapshots_dir, output_dir, len(time_points))
    
    print("\n" + "="*80)
    print("SIMULATION COMPLETE")
    print("="*80)
    print(f"\nOutput directory: {output_dir}/")
    if save_snapshots:
        print(f"  - snapshots/: {len(time_points)} individual frames")
        print(f"  - animation.gif: Animated evolution")
    print(f"  - results.png: Final analysis plots")
    
    return {
        'time': np.array(time_points),
        'total_prob': np.array(total_probability),
        'inlet_prob': np.array(inlet_probability),
        'impurity_prob': np.array(impurity_probability),
        'outlet_prob': np.array(outlet_probability),
        'psi_final': psi_final,
        'V_field': V_field_cpu
    }

def save_snapshot(psi_cpu, V_field_cpu, step, time_val, total_prob, inlet_prob, 
                  impurity_prob, outlet_prob, output_dir, prob_density_max):
    """Save a snapshot of the current simulation state."""
    
    fig = plt.figure(figsize=(14, 5))
    
    # Probability density
    ax1 = plt.subplot(1, 2, 1)
    # Probability per pixel (dimensionless): integrate |ψ|² over voxel volume
    prob_per_pixel = np.sum(np.abs(psi_cpu[:, :, 0, :])**2, axis=-1).T * DX * DY * DZ
    # Use fixed colormap scale to avoid rescaling when wave exits
    vmax_fixed = prob_density_max * DX * DY * DZ
    im1 = ax1.imshow(prob_per_pixel, origin='lower', cmap='hot', aspect='equal',
                     extent=[0, LX*1e9, 0, LY*1e9], vmin=0, vmax=vmax_fixed)
    ax1.axvline(x=IMPURITY_START*DX*1e9, color='cyan', linestyle='--', alpha=0.5, lw=1)
    ax1.axvline(x=IMPURITY_END*DX*1e9, color='cyan', linestyle='--', alpha=0.5, lw=1)
    ax1.set_xlabel('x (nm)')
    ax1.set_ylabel('y (nm)')
    ax1.set_title(f'Probability Density | t={time_val*1e9:.2f} ns')
    plt.colorbar(im1, ax=ax1, label='Probability per Pixel')
    
    # Potential with overlay
    ax2 = plt.subplot(1, 2, 2)
    V_slice = V_field_cpu[:, :, 0].T / Q_ELECTRON * 1e3
    im2 = ax2.imshow(V_slice, origin='lower', cmap='viridis', aspect='equal',
                     extent=[0, LX*1e9, 0, LY*1e9], alpha=0.7)
    
    # Overlay probability density contours
    X, Y = np.meshgrid(np.linspace(0, LX*1e9, NX), np.linspace(0, LY*1e9, NY))
    ax2.contour(X, Y, prob_per_pixel, levels=5, colors='white', alpha=0.6, linewidths=1)
    
    ax2.set_xlabel('x (nm)')
    ax2.set_ylabel('y (nm)')
    ax2.set_title('Potential + Wave Packet')
    plt.colorbar(im2, ax=ax2, label='Potential (meV)')
    
    # Add text overlay with statistics
    stats_text = f'Step: {step}\nP_total: {total_prob:.4f}\nInlet: {inlet_prob:.3f}\nImpurity: {impurity_prob:.3f}\nOutlet: {outlet_prob:.3f}'
    ax2.text(0.02, 0.98, stats_text, transform=ax2.transAxes,
             fontsize=9, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'snapshot_{step:05d}.png'), dpi=100, bbox_inches='tight')
    plt.close()

def create_animation(snapshots_dir, output_dir, n_frames):
    """Create an animated GIF from saved snapshots."""
    try:
        from PIL import Image
        import glob
        
        # Get all snapshot files
        snapshot_files = sorted(glob.glob(os.path.join(snapshots_dir, 'snapshot_*.png')))
        
        if len(snapshot_files) == 0:
            print("  Warning: No snapshots found to create animation")
            return
        
        # Load images
        images = [Image.open(f) for f in snapshot_files]
        
        # Save as GIF
        gif_path = os.path.join(output_dir, 'animation.gif')
        images[0].save(
            gif_path,
            save_all=True,
            append_images=images[1:],
            duration=200,  # 200ms per frame
            loop=0
        )
        
        print(f"  Animation created: {gif_path}")
        print(f"  Frames: {len(images)}, Duration: {len(images)*0.2:.1f}s")
        
    except ImportError:
        print("  PIL/Pillow not available, skipping animation")
    except Exception as e:
        print(f"  Error creating animation: {e}")

def plot_results(time_pts, total_prob, inlet_prob, impurity_prob, outlet_prob,
                 psi, V_field, output_dir):
    """Create analysis plots."""
    
    time_ns = np.array(time_pts) * 1e9
    
    fig = plt.figure(figsize=(16, 10))
    
    # Probabilities
    ax1 = plt.subplot(2, 3, 1)
    ax1.plot(time_ns, total_prob, 'k-', lw=2, label='Total')
    ax1.plot(time_ns, inlet_prob, 'b-', label='Inlet')
    ax1.plot(time_ns, impurity_prob, 'r-', label='Impurity')
    ax1.plot(time_ns, outlet_prob, 'g-', label='Outlet')
    ax1.set_xlabel('Time (ns)')
    ax1.set_ylabel('Probability')
    ax1.set_title('Regional Probability Evolution')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Transmission
    ax2 = plt.subplot(2, 3, 2)
    transmission = np.array(outlet_prob) / (np.array(outlet_prob) + np.array(inlet_prob) + 1e-10)
    ax2.plot(time_ns, transmission, 'g-', lw=2)
    ax2.set_xlabel('Time (ns)')
    ax2.set_ylabel('Transmission')
    ax2.set_title('Transmission Coefficient')
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim([0, 1.1])
    
    # Conservation
    ax3 = plt.subplot(2, 3, 3)
    error = np.abs(np.array(total_prob) - 1.0)
    ax3.semilogy(time_ns, error, 'k-', lw=2)
    ax3.set_xlabel('Time (ns)')
    ax3.set_ylabel('|P - 1|')
    ax3.set_title('Conservation Error')
    ax3.grid(True, alpha=0.3)
    
    # Final probability
    ax4 = plt.subplot(2, 3, 4)
    prob_final = np.sum(np.abs(psi[:, :, 0, :])**2, axis=-1).T * DX * DY * DZ
    im = ax4.imshow(prob_final, origin='lower', cmap='hot', aspect='equal',
                    extent=[0, LX*1e9, 0, LY*1e9])
    ax4.axvline(x=IMPURITY_START*DX*1e9, color='cyan', linestyle='--', alpha=0.7)
    ax4.axvline(x=IMPURITY_END*DX*1e9, color='cyan', linestyle='--', alpha=0.7)
    ax4.set_xlabel('x (nm)')
    ax4.set_ylabel('y (nm)')
    ax4.set_title('Final Probability Density')
    plt.colorbar(im, ax=ax4)
    
    # Potential
    ax5 = plt.subplot(2, 3, 5)
    V_slice = V_field[:, :, 0].T / Q_ELECTRON * 1e3
    im = ax5.imshow(V_slice, origin='lower', cmap='viridis', aspect='equal',
                    extent=[0, LX*1e9, 0, LY*1e9])
    ax5.set_xlabel('x (nm)')
    ax5.set_ylabel('y (nm)')
    ax5.set_title('Potential Field (meV)')
    plt.colorbar(im, ax=ax5)
    
    # Summary
    ax6 = plt.subplot(2, 3, 6)
    ax6.axis('off')
    summary = (
        "SIMULATION SUMMARY\n\n"
        f"Resolution: 1/4 of paper\n"
        f"Grid: {NX} x {NY}\n"
        f"Domain: {LX*1e9:.0f} x {LY*1e9:.0f} nm\n\n"
        f"FINAL STATE (t={time_ns[-1]:.2f} ns):\n"
        f"  Total Prob:    {total_prob[-1]:.6f}\n"
        f"  Inlet:         {inlet_prob[-1]:.4f}\n"
        f"  Impurity:      {impurity_prob[-1]:.4f}\n"
        f"  Outlet:        {outlet_prob[-1]:.4f}\n\n"
        f"  Transmission:  {transmission[-1]:.4f}\n\n"
        f"CONSERVATION:\n"
        f"  Max |ΔP|:      {np.max(error):.2e}\n"
        f"  Final |ΔP|:    {error[-1]:.2e}\n"
    )
    ax6.text(0.1, 0.5, summary, fontsize=10, family='monospace',
             verticalalignment='center')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'results.png'), dpi=150)
    plt.close()

if __name__ == "__main__":
    # Run for ~1.5× wave transit time (448 cells at c → 700 steps)
    results = run_simulation_gpu(n_steps=700, output_freq=100)
    
    print(f"\nFinal Results:")
    print(f"  Transmission: {results['outlet_prob'][-1]/(results['outlet_prob'][-1] + results['inlet_prob'][-1] + 1e-10):.4f}")
    print(f"  Conservation: {results['total_prob'][-1]:.6f}")
