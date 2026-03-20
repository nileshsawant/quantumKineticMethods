"""
Palpacelli 2012 – Klein Tunneling in Random Impurities (Full Paper Resolution)
===============================================================================
Reproduces Figure 8 of:
  S. Palpacelli, G. Falcucci, S. Succi (2012)
  "Klein Tunneling in the Presence of Random Impurities"
  Int. J. Mod. Phys. C, 23(12), 1250080.

PHYSICS SUMMARY
---------------
Low-energy electrons in graphene obey the 2D massless Dirac equation with
Fermi velocity v_F ≈ 10⁶ m/s acting as the effective speed of light.  Because
they are relativistic Dirac fermions, they exhibit **Klein tunneling**: near-
perfect transmission through potential barriers at normal incidence, regardless
of barrier height.

This simulation sends a Gaussian wave packet (E = 80 meV) through a region
filled with randomly placed square potential barriers (C = 0.5%, V = 50 meV).
The wave packet fragments into a spreading plane front – each impurity
individually scatters and redirects part of the wave, and the combined effect
is measured via the transmission coefficient T(t) = P_outlet / P_total.

NUMERICAL METHOD: Quantum Lattice Boltzmann (QLB)
--------------------------------------------------
The QLB method (Palpacelli 2008; Dellar 2011) is the quantum analogue of the
classical lattice Boltzmann method.  Instead of fluid distribution functions,
we track the 4-component complex Dirac spinor ψ(x, y, t).

Each time step is split into two 1D sweeps via operator splitting (Lie-Trotter):
  1. X-sweep along x (propagation direction):
       (a) Rotate ψ into characteristic frame: ψ̃ = X⁻¹ψ
           X diagonalises -αₓ → components (u₁,u₂) travel right, (d₁,d₂) left
       (b) Collide: apply unitary collision matrix Q (mixes u↔d via mass/potential)
       (c) Stream: shift u components one cell right, d components one cell left
       (d) Rotate back: ψ = X ψ̃
  2. Y-sweep along y (transverse direction):
       Y = Identity (β is already diagonal in the Majorana form), so no rotation.
       Periodic boundary conditions top and bottom.

The collision matrix Q is unitary (|a|²+|b|²=1), so total probability ∫|ψ|²dV
is exactly conserved at every step (no numerical dissipation).

KEY IMPLEMENTATION NOTES
--------------------------
• v_F (not c) is the effective speed: DT = DX/v_F ≈ 9.6×10⁻¹⁶ s
• Energy unit ΔE = ℏ/DT ≈ 686 meV → k₀=0.117 gives E=80 meV
• V_field stores potential ENERGY in Joules (= q × V_volts); g̃₃ = V·DT/(3·ℏ)
  — do NOT multiply by Q_ELECTRON again inside the collision formula
• Inlet: bounce-back (u←d, d=0) in X-characteristic frame
• Outlet: open/absorbing (d=0) + exponential sponge layer
• Y boundaries: periodic wrap-around

References (markdown files in this directory):
  palpacelli2012.md – simulation setup (domain, BCs, wave packet, impurities)
  dellar2011.md     – rotation matrices X,Z; collision matrices Q, Q̂; isotropy
  palpacelli2008.md – original 1D QLB formulation
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import time

# ---------------------------------------------------------------------------
# GPU acceleration via CuPy (optional)
# CuPy mirrors the NumPy API but runs on CUDA GPUs.  If CuPy is not installed
# the code silently falls back to plain NumPy on the CPU; all logic is identical.
# ---------------------------------------------------------------------------
try:
    import cupy as cp
    GPU_AVAILABLE = True
    print("GPU (CuPy) detected and will be used for acceleration!")
except ImportError:
    print("CuPy not available, falling back to NumPy (CPU)")
    cp = np
    GPU_AVAILABLE = False

# ---------------------------------------------------------------------------
# Import Dirac matrices and physical constants from the core solver module.
#
# ALPHA_X, ALPHA_Y, ALPHA_Z  – 4×4 Dirac alpha matrices (streaming operators)
# BETA                        – 4×4 Dirac beta matrix (mass term / y-streaming)
# IDENTITY_4x4                – 4×4 identity
#
# X_MATRIX, X_INV_MATRIX  –  Rotation that diagonalises -αₓ  (Dellar 2011 Eq.22)
# Y_MATRIX, Y_INV_MATRIX  –  Identity (β is already diagonal in Majorana form)
# Z_MATRIX, Z_INV_MATRIX  –  Rotation that diagonalises -αz  (Dellar 2011 Eq.10)
#
# HBAR        = 1.0545718×10⁻³⁴ J·s   (reduced Planck constant)
# C           = 2.99792458×10⁸ m/s   (speed of light, used only in solver module)
# Q_ELECTRON  = 1.60217663×10⁻¹⁹ C   (elementary charge)
# ---------------------------------------------------------------------------
from dirac_qlb_solver import (
    ALPHA_X, ALPHA_Y, ALPHA_Z, BETA, IDENTITY_4x4,
    X_MATRIX, Y_MATRIX, Z_MATRIX,
    X_INV_MATRIX, Y_INV_MATRIX, Z_INV_MATRIX,
    HBAR, C, Q_ELECTRON
)

# =============================================================================
# SIMULATION PARAMETERS  (Paper: Palpacelli 2012, Section 4)
# =============================================================================

# ---------------------------------------------------------------------------
# Grid dimensions
# Full paper resolution for the 2D graphene test: NX=2048, NY=512, NZ=1.
# For a 3D simulation increase NZ and set any required z-extent.
# ---------------------------------------------------------------------------
NX = 2048  # cells along x  (propagation direction)
NY = 512   # cells along y  (transverse, periodic BCs)
NZ = 1     # cells along z  (set to 1 for the 2D Palpacelli test case;
           #                  increase for full 3D runs, periodic BCs)

# ---------------------------------------------------------------------------
# Physical lattice spacing
# Chosen large enough to resolve both σ=48 and d=8 with ~6 cells/impurity
# (Palpacelli 2012 Sec. 4: Δx = 0.96 nm).
# DY = DZ = DX for an isotropic square lattice.
# ---------------------------------------------------------------------------
DX = 0.96e-9  # m  (0.96 nm)
DY = DX
DZ = DX

# ---------------------------------------------------------------------------
# Fermi velocity – the speed of light analogue in graphene
#
# In graphene, the band structure near the Dirac point gives a linear
# dispersion E = ℏ v_F |k|, with v_F ≈ 10⁶ m/s.  The Dirac equation
# for graphene is therefore identical to the relativistic Dirac equation
# but with v_F instead of c.
#
# IMPORTANT: v_F ≈ 10⁶ m/s is ~300× SMALLER than c = 3×10⁸ m/s.
# Using c by mistake makes DT 300× too small → g̃₃ ≈ 0 → impurities invisible.
# ---------------------------------------------------------------------------
V_F = 1.0e6  # m/s  (Fermi velocity in graphene)

# Physical domain extents (for plotting)
LX = NX * DX  # ~1.966 µm
LY = NY * DY  # ~0.492 µm
LZ = NZ * DZ  # DZ when NZ=1; NZ×DZ for a true 3D domain

# ---------------------------------------------------------------------------
# Time step – light-cone (Courant-Friedrichs-Lewy) condition
#
# The QLB streaming step shifts each characteristic component by exactly one
# lattice cell per time step.  This requires DT = DX / v_F so that one cell
# equals one light-cone step.
#
# DT = 0.96×10⁻⁹ / 10⁶ = 9.6×10⁻¹⁶ s
# Energy unit: ΔE = ℏ/DT ≈ 686 meV
# ---------------------------------------------------------------------------
DT = DX / V_F  # ≈ 9.6×10⁻¹⁶ s

# ---------------------------------------------------------------------------
# Particle mass
# For graphene's massless Dirac fermions set M_PARTICLE = 0.
# The mass enters the collision coefficient b̂ = m̃₃ / D; with m=0, b̂=0
# and the collision reduces to a pure phase rotation (no spin-flip scattering).
# For a massive particle (e.g. bilayer graphene gap) set M_PARTICLE > 0.
# ---------------------------------------------------------------------------
M_PARTICLE = 0.0           # kg  (0 = massless graphene)
OMEGA_C = 0.0              # rad/s  (Compton frequency = M*v_F²/ℏ; 0 here)

# ---------------------------------------------------------------------------
# Domain regions  (Palpacelli 2012 Sec. 4)
#
# The domain is divided into three sections along x:
#
#   [0, 512)      INLET      Wave packet is placed here at t=0.
#                            Bounce-back BC at x=0: left-movers reflect as right-movers.
#
#   [512, 1536)   IMPURITIES Random square barriers (concentration C, height V).
#
#   [1536, 2048]  OUTLET     Transmitted wave is measured here.
#                            Open BC at x=N_x-1: left-movers are absorbed.
# ---------------------------------------------------------------------------
INLET_START    = 0
INLET_END      = 512
IMPURITY_START = 512
IMPURITY_END   = 1536
OUTLET_START   = 1536
OUTLET_END     = 2048

# ---------------------------------------------------------------------------
# Wave-packet parameters  (Palpacelli 2012 Sec. 4, Eq. 11)
#
# The initial state is a 3D Gaussian ψ ∝ exp(-r²/4σ²) · exp(ik₀x),
# centred at (3σ, NY/2, NZ/2).  For the 2D case (NZ=1) the z-extent
# collapses to one cell and the packet is effectively 2D.
# σ = 48 lattice cells, k₀ = 0.117/Δx.
#
# Energy:  E = ℏ v_F k₀ = K0_LATTICE × ΔE = 0.117 × 686 meV ≈ 80 meV
#
# Diffraction: a Gaussian beam diffracts at half-angle θ ≈ 1/(σ·k₀).
#   θ = 1/(48 × 0.117) ≈ 0.178 rad ≈ 10° — same as in the paper.
#   Over the 512-cell inlet, the beam spreads ~90 cells transversely (18% of NY).
#
# The ratio σ/d = 48/8 = 6 ensures the packet is wider than each individual
# impurity, so it can flow around impurities like a classical fluid as well as
# tunnel through them quantum mechanically.
# ---------------------------------------------------------------------------
SIGMA_LATTICE = 48    # wave-packet width σ in lattice cells (paper: 48)
K0_LATTICE = 0.117    # k₀·Δx (dimensionless lattice momentum; paper: 0.117 → E=80 meV)

# ---------------------------------------------------------------------------
# Impurity parameters  (Palpacelli 2012 Fig. 8)
#
# Impurities are randomly placed d×d square barriers in the impurity region.
# C = N_imp · d² / A_impurity_region  (fractional area coverage)
#
# V = 50 meV < E = 80 meV → "weak potential" regime.
# The barriers are sub-threshold (classically forbidden for all angles), but
# quantum tunnelling (Klein tunneling) still allows significant transmission.
#
# DEFAULT_BARRIER_HEIGHT is stored as potential ENERGY in Joules = q·V_volts.
# Inside perform_qlb_line_gpu the coupling is g̃₃ = V_energy·DT / (3·ℏ) ≈ 0.024.
# ---------------------------------------------------------------------------
IMPURITY_SIZE_CELLS = 8          # square-barrier side d in cells (paper: 8)
DEFAULT_CONCENTRATION = 0.005    # 0.5% area coverage (paper Fig. 8)
DEFAULT_BARRIER_HEIGHT = 50e-3 * Q_ELECTRON  # 50 meV expressed in Joules

print("\n" + "="*80)
print("PALPACELLI 2012 - FULL RESOLUTION")
print("="*80)
print(f"\nConfiguration:")
print(f"  Grid: {NX} × {NY} × {NZ} (Palpacelli 2012 paper: 2048 × 512 × 1)")
print(f"  Domain: {LX*1e9:.1f} × {LY*1e9:.1f} nm (paper: 1966 × 492 nm)")
print(f"  Lattice spacing: {DX*1e9:.3f} nm")
print(f"  Fermi velocity:  {V_F:.1e} m/s (effective c in graphene Dirac eq.)")
print(f"  Time step:       {DT:.3e} s")
print(f"  Energy unit ΔE:  {(6.582e-16/DT)*1e3:.1f} meV")
print(f"  Wave packet energy: {K0_LATTICE * 6.582e-16/DT * 1e3:.1f} meV (paper: 80 meV)")  # ℏ/DT * K0
print(f"  Wave packet σ: {SIGMA_LATTICE} cells (paper: 48; diffr. angle={1/(SIGMA_LATTICE*K0_LATTICE)*57.3:.0f}°)")
print(f"  Impurity size: {IMPURITY_SIZE_CELLS}×{IMPURITY_SIZE_CELLS} cells (σ/d={SIGMA_LATTICE//IMPURITY_SIZE_CELLS}, paper σ/d=6)")
print(f"  Concentration: {DEFAULT_CONCENTRATION*100:.1f}% (Fig. 8: 0.5%)")
print(f"  Barrier height: {DEFAULT_BARRIER_HEIGHT/Q_ELECTRON*1e3:.0f} meV (Fig. 8: 50 meV)")
# Diagnostic: check coupling strength g_tilde_3 at impurity cells
_g3_check = DEFAULT_BARRIER_HEIGHT * DT / (3 * HBAR)
_energy_ratio = DEFAULT_BARRIER_HEIGHT / (K0_LATTICE * HBAR / DT)
print(f"  Coupling g̃₃ at V=50meV: {_g3_check:.4f} (should be ~0.024)")
print(f"  V/E ratio: {_energy_ratio:.3f} (paper: 50/80 = 0.625, weak potential)")
print(f"  GPU acceleration: {'ENABLED' if GPU_AVAILABLE else 'DISABLED'}")
print("="*80)

# ---------------------------------------------------------------------------
# Copy rotation matrices to GPU memory (if CuPy/CUDA is available).
# All subsequent einsum operations will then run on the GPU.
# On CPU, these are just aliases to the same NumPy arrays.
# ---------------------------------------------------------------------------
if GPU_AVAILABLE:
    X_MATRIX_GPU     = cp.asarray(X_MATRIX)
    Y_MATRIX_GPU     = cp.asarray(Y_MATRIX)
    Z_MATRIX_GPU     = cp.asarray(Z_MATRIX)
    X_INV_MATRIX_GPU = cp.asarray(X_INV_MATRIX)
    Y_INV_MATRIX_GPU = cp.asarray(Y_INV_MATRIX)
    Z_INV_MATRIX_GPU = cp.asarray(Z_INV_MATRIX)
    BETA_GPU         = cp.asarray(BETA)
else:
    X_MATRIX_GPU     = X_MATRIX
    Y_MATRIX_GPU     = Y_MATRIX
    Z_MATRIX_GPU     = Z_MATRIX
    X_INV_MATRIX_GPU = X_INV_MATRIX
    Y_INV_MATRIX_GPU = Y_INV_MATRIX
    Z_INV_MATRIX_GPU = Z_INV_MATRIX
    BETA_GPU         = BETA

def generate_random_impurities_gpu(concentration, barrier_height, seed=42):
    """
    Build the scalar potential field V(x, y) for the impurity region.

    Each impurity is a square barrier of side IMPURITY_SIZE_CELLS cells and
    constant potential energy `barrier_height` (in Joules).  Impurities are
    placed at random positions entirely within the impurity region
    [IMPURITY_START, IMPURITY_END) × [0, NY).

    Concentration is defined as  C = N_imp · d² / A_region  (Palpacelli 2012
    Sec. 4), so the number of impurities is:
        N_imp = C · A_region / d²

    Parameters
    ----------
    concentration   : float  – fractional area coverage  (e.g. 0.005 = 0.5%)
    barrier_height  : float  – potential energy in Joules (= q × V_volts)
    seed            : int    – NumPy random seed for reproducibility

    Returns
    -------
    V_field          : array (NX, NY, NZ) on GPU (or CPU)
    impurity_positions : list of (ix, iy) lower-left corners
    """
    # Allocate potential field on CPU (random number generation is CPU-only)
    V_field = np.zeros((NX, NY, NZ), dtype=np.float64)

    # Number of impurities from the concentration definition
    impurity_region_length = IMPURITY_END - IMPURITY_START
    impurity_region_area   = impurity_region_length * NY
    impurity_area          = IMPURITY_SIZE_CELLS * IMPURITY_SIZE_CELLS
    N_impurities = int(concentration * impurity_region_area / impurity_area)

    print(f"\nGenerating {N_impurities} random impurities...")

    np.random.seed(seed)
    impurity_positions = []

    for _ in range(N_impurities):
        # Random lower-left corner, keeping the full d×d block inside the region
        ix = np.random.randint(IMPURITY_START, IMPURITY_END - IMPURITY_SIZE_CELLS)
        iy = np.random.randint(0, NY - IMPURITY_SIZE_CELLS)

        # Fill the d×d×NZ block with the barrier potential energy.
        # All z-layers get the same value so the potential is z-invariant
        # for the 2D case (NZ=1) and extends uniformly in z for NZ > 1.
        for di in range(IMPURITY_SIZE_CELLS):
            for dj in range(IMPURITY_SIZE_CELLS):
                V_field[ix + di, iy + dj, :] = barrier_height

        impurity_positions.append((ix, iy))

    # Move potential field to GPU for use in the time loop
    if GPU_AVAILABLE:
        V_field = cp.asarray(V_field)

    return V_field, impurity_positions

def initialize_wave_packet_gpu():
    """
    Initialise the 4-component Dirac spinor as a 3D Gaussian wave packet.

    Following Palpacelli 2012 Eq. 11, the initial state is:

        ψ_l(x, y) = A_l / (4πσ²)^(1/2) · exp(-r²/4σ²) · exp(i k₀ x)

    with l=1,2 (bi-spinor components), A₁ = A₂ = 1/√2 for head-on incidence
    (no y-momentum, φ = 0).

    In 4-component form the non-zero entries are components 1 and 3 (0-indexed):
        ψ₁ = ψ₃ = amplitude / √2

    WHY THIS CHOICE?
    The X-rotation X⁻¹ maps the physical spinor to the characteristic frame
    (u₁, u₂, d₁, d₂).  For ψ = (0, A, 0, A)/√2 one can verify analytically:
        u₁ = 0,  u₂ = A,  d₁ = 0,  d₂ = 0
    — a PURE RIGHT-MOVER.  Any left-moving content (non-zero d₁ or d₂) would
    cause the wave packet to immediately split into a forward- and a backward-
    travelling pulse at t=0, which is an initialisation artefact.

    The centre is placed at x₀ = 3σ from the inlet boundary so that the
    Gaussian tail is negligible at x=0 and the inlet bounce-back BC does not
    prematurely reflect the packet.
    """
    # Centre position (in cells)
    x0 = INLET_START + 3 * SIGMA_LATTICE   # 3σ from left edge
    y0 = NY // 2                            # middle of domain height
    z0 = NZ // 2                             # centre in z (= 0 when NZ=1)

    # Physical width σ (in metres)
    D0_initial = SIGMA_LATTICE * DX         # e.g. 48 × 0.96 nm = 46 nm

    # Carrier momentum (SI units) corresponding to K0_LATTICE
    # k₀ = K0_LATTICE / DX  →  p₀ = ℏ k₀
    initial_momentum_x = HBAR * K0_LATTICE / DX

    # Build coordinate arrays on CPU (meshgrid) then compute the envelope
    i_grid, j_grid, k_grid = np.meshgrid(
        np.arange(NX), np.arange(NY), np.arange(NZ), indexing='ij'
    )
    # Displacement from packet centre (metres)
    rx = (i_grid - x0) * DX
    ry = (j_grid - y0) * DY
    rz = (k_grid - z0) * DZ

    # Gaussian envelope: exp(-r² / 4σ²)
    gaussian_envelope = np.exp(-(rx**2 + ry**2 + rz**2) / (4 * D0_initial**2))

    # Plane-wave phase: exp(i k₀ x)  – carrier wave in the +x direction
    phase_factor = np.exp(1j * initial_momentum_x * rx / HBAR)

    # Normalisation prefactor: (2π σ²)^(-3/4) in 3D Gaussian convention
    # (When NZ=1 the z-Gaussian is a single cell and the prefactor still
    # normalises correctly because the wave function is numerically re-normalised
    # below via norm_factor.)
    base_amplitude = (
        (1 / (np.sqrt(2 * np.pi) * D0_initial))**(3/2)
        * gaussian_envelope
        * phase_factor
    )

    # Assemble 4-component spinor: only components 1 and 3 are non-zero
    psi_cpu = np.zeros((NX, NY, NZ, 4), dtype=np.complex128)
    psi_cpu[:, :, :, 1] = base_amplitude / np.sqrt(2)  # up-spin right-mover
    psi_cpu[:, :, :, 3] = base_amplitude / np.sqrt(2)  # down-spin right-mover

    # Normalise so that ∫|ψ|² dV = 1 (total probability = 1)
    norm_factor = np.sum(np.abs(psi_cpu)**2) * DX * DY * DZ
    if norm_factor > 0:
        psi_cpu /= np.sqrt(norm_factor)

    prob_check = np.sum(np.abs(psi_cpu)**2) * DX * DY * DZ
    print(f"\nWave packet initialized:")
    print(f"  Centre cell:  ({x0}, {y0}, {z0})")
    print(f"  Width σ:      {D0_initial*1e9:.2f} nm  ({SIGMA_LATTICE} cells)")
    print(f"  Energy E:     {K0_LATTICE * HBAR/DT / Q_ELECTRON * 1e3:.1f} meV")
    print(f"  Total prob:   {prob_check:.6f}  (should be 1.0)")

    # Transfer to GPU for the simulation loop
    if GPU_AVAILABLE:
        psi = cp.asarray(psi_cpu)
    else:
        psi = psi_cpu

    return psi

def perform_qlb_line_gpu(psi_line, V_line, axis, R_inv, R_op, periodic=False):
    """
    Perform one QLB sub-step along a single 1D line.

    This implements the four-stage QLB update (Dellar 2011, Palpacelli 2012):
        1. Rotate  ψ → ψ̃ = R⁻¹ ψ   (into characteristic frame for this axis)
        2. Collide  ψ̃ → Q ψ̃         (unitary mixing of right/left movers)
        3. Stream   shift u→right, d→left  (exact, one cell per step)
        4. Rotate back  ψ = R ψ̃     (back to physical spinor basis)

    Parameters
    ----------
    psi_line  : array (N, 4)   – 4-component spinor along a 1D line (on GPU/CPU)
    V_line    : array (N,)     – potential ENERGY [J] along the line (q·V_volts)
    axis      : str            – 'x', 'y', or 'z' (determines collision sign pattern)
    R_inv     : array (4,4)    – rotation matrix R⁻¹ = R†  (X_INV or Y_INV=I)
    R_op      : array (4,4)    – rotation matrix R        (X     or Y=I)
    periodic  : bool           – True → periodic BCs (y-direction);
                                 False → open BCs (x-direction, handled externally)

    Returns
    -------
    psi_out : array (N, 4)  – updated spinor (on same device as input)

    COLLISION DETAILS
    -----------------
    The 3D QLB scheme splits the collision term into three equal parts, one per
    directional sweep, each using one-third of the time step.

    Step 1 – compute dimensionless collision parameters at every cell:

        m̃₃ = (1/3) · M · v_F² · DT / ℏ          (rescaled Compton frequency)
        g̃₃ = (1/3) · V_energy · DT / ℏ           (rescaled potential coupling)
                    ^^^^^^^^^^^^
                    V_line is already in Joules = q·V_volts.
                    Do NOT multiply by Q_ELECTRON here.

        Ω₃ = m̃₃² - g̃₃²
        D  = 1 + Ω₃/4 - i·g̃₃
        â  = (1 - Ω₃/4) / D           (diagonal collision factor)
        b̂  = m̃₃ / D                   (off-diagonal, mass-coupling factor; =0 for massless)

    Step 2 – rotate to characteristic frame: ψ̃ = R⁻¹ ψ, giving (ũ₁,ũ₂,d̃₁,d̃₂)

    Step 3 – apply collision.  The sign pattern differs between sweeps because
    X⁻¹Q̂X has the Q-sign pattern while Q̂ itself has a different pattern
    (Dellar 2011: "X⁻¹Q̂X and Z⁻¹Q̂Z have the same sign pattern as Q, but Q̂
    itself does not").

      Y-sweep  (R=I, working directly in Q̂):
        ũ₁' = â·ũ₁ - b̂·d̃₂      ← minus sign on b̂
        ũ₂' = â·ũ₂ + b̂·d̃₁
        d̃₁' = â·d̃₁ - b̂·ũ₂
        d̃₂' = â·d̃₂ + b̂·ũ₁      ← plus sign on b̂

      X/Z-sweep  (already in rotated frame X⁻¹Q̂X = Q-sign pattern):
        ũ₁' = â·ũ₁ + b̂·d̃₂      ← plus sign on b̂
        ũ₂' = â·ũ₂ + b̂·d̃₁
        d̃₁' = â·d̃₁ - b̂·ũ₂
        d̃₂' = â·d̃₂ - b̂·ũ₁      ← minus sign on b̂

      For massless particles (m=0 → b̂=0) both patterns are identical.

    Step 4 – stream: ũ components shift +1, d̃ shift -1 (exactly one cell).

    Step 5 – rotate back: ψ = R ψ̃.
    """
    N_points = psi_line.shape[0]
    xp = cp if GPU_AVAILABLE else np  # use cupy or numpy uniformly

    # ------------------------------------------------------------------
    # STEP 1: compute collision coefficients at each cell along the line
    # ------------------------------------------------------------------
    # Rescaled Compton frequency (one-third of full DT because we split
    # the collision equally among the three sweeps)
    m_tilde_3 = M_PARTICLE * V_F**2 * DT / (3 * HBAR) if HBAR > 0 else 0.0

    # Rescaled potential coupling.
    # V_line is stored as ENERGY in Joules (= q × V_volts), so the formula is
    #   g̃₃ = V_energy · DT / (3 · ℏ)
    # NOT Q_ELECTRON × V_line, which would double-count the charge factor.
    g_tilde_3 = V_line * DT / (3 * HBAR) if HBAR > 0 else xp.zeros_like(V_line)

    Omega_3 = m_tilde_3**2 - g_tilde_3**2

    # Common denominator D = 1 + Ω₃/4 - i·g̃₃
    # Guard against |D| → 0 (cannot happen physically but avoids NaN for V→∞)
    denominator = 1 + Omega_3 / 4 - 1j * g_tilde_3
    denominator = xp.where(xp.abs(denominator) < 1e-15, 1e-15, denominator)

    a_hat = (1 - Omega_3 / 4) / denominator   # |a|² + |b|² = 1  (unitarity)
    b_hat = m_tilde_3 / denominator            # = 0 for massless particles

    # ------------------------------------------------------------------
    # STEP 2: rotate ψ into the characteristic frame for this axis.
    # einsum 'ij,kj->ki':  for each grid cell k, compute  ψ̃[k] = R⁻¹ ψ[k]
    # Result shape: (N, 4) = (N_points, [u₁, u₂, d₁, d₂])
    # ------------------------------------------------------------------
    psi_rotated = xp.einsum('ij,kj->ki', R_inv, psi_line)

    # Extract the four characteristic components as 1D arrays
    u1 = psi_rotated[:, 0]   # right-mover, spin-up
    u2 = psi_rotated[:, 1]   # right-mover, spin-down
    d1 = psi_rotated[:, 2]   # left-mover,  spin-up
    d2 = psi_rotated[:, 3]   # left-mover,  spin-down

    # ------------------------------------------------------------------
    # STEP 3: apply the collision (unitary mixing of right ↔ left movers)
    # ------------------------------------------------------------------
    if axis == 'y':
        # Y-sweep: Y = I, no rotation was applied, so we collide with Q̂ directly.
        # Q̂ has the α_y sign pattern (Palpacelli 2012 Eq. 10):
        #   (u₁,d₂) coupling has a MINUS sign on b̂
        collided_u1 = a_hat * u1 - b_hat * d2
        collided_u2 = a_hat * u2 + b_hat * d1
        collided_d1 = a_hat * d1 - b_hat * u2
        collided_d2 = a_hat * d2 + b_hat * u1
    else:
        # X and Z sweeps: rotated frame — X⁻¹Q̂X and Z⁻¹Q̂Z both give the Q sign
        # pattern (Dellar 2011: all non-identity rotation matrices R satisfy
        # R⁻¹Q̂R = Q-pattern).  (u₁,d₂) coupling has a PLUS sign on b̂.
        collided_u1 = a_hat * u1 + b_hat * d2
        collided_u2 = a_hat * u2 + b_hat * d1
        collided_d1 = a_hat * d1 - b_hat * u2
        collided_d2 = a_hat * d2 - b_hat * u1

    # ------------------------------------------------------------------
    # STEP 4: stream – shift characteristic components by one lattice cell.
    # u components travel in the +direction → each u₁[k] moves to cell k+1.
    # d components travel in the -direction → each d₁[k] moves to cell k-1.
    # This is exact (no interpolation) because DT = DX/v_F was chosen so that
    # exactly one cell is traversed per time step.
    # ------------------------------------------------------------------
    psi_collided_rotated = xp.zeros_like(psi_rotated)

    # Interior streaming (bulk of the domain)
    psi_collided_rotated[1:,  0] = collided_u1[:-1]  # u₁: cell k → k+1
    psi_collided_rotated[1:,  1] = collided_u2[:-1]  # u₂: cell k → k+1
    psi_collided_rotated[:-1, 2] = collided_d1[1:]   # d₁: cell k → k-1
    psi_collided_rotated[:-1, 3] = collided_d2[1:]   # d₂: cell k → k-1

    # Boundary streaming
    if periodic:
        # PERIODIC (Y and Z directions — wrap-around BCs):
        # u components wrap from the last cell to the first
        psi_collided_rotated[0,  0] = collided_u1[-1]  # u₁: last → first
        psi_collided_rotated[0,  1] = collided_u2[-1]  # u₂: last → first
        # d components wrap from the first cell to the last
        psi_collided_rotated[-1, 2] = collided_d1[0]   # d₁: first → last
        psi_collided_rotated[-1, 3] = collided_d2[0]   # d₂: first → last
    # else: OPEN (X-direction) – edge cells receive nothing from outside.
    # u at cell 0 and d at cell N-1 are left as 0, which corresponds to an
    # absorbing/outflow condition.  The actual inlet bounce-back and outlet
    # open BCs are applied externally in the main time loop AFTER each X-sweep.

    # ------------------------------------------------------------------
    # STEP 5: rotate back to the physical spinor basis.
    # einsum 'ij,kj->ki':  for each cell k, compute  ψ[k] = R ψ̃[k]
    # ------------------------------------------------------------------
    psi_out = xp.einsum('ij,kj->ki', R_op, psi_collided_rotated)

    return psi_out

def perform_qlb_sweep_gpu(psi, V_field, axis, R_inv, R_op, periodic=False):
    """
    Perform one QLB sub-step over the FULL (NX, NY, NZ, 4) field for a given axis.

    This is the vectorized version of perform_qlb_line_gpu.  It replaces all
    Python loops over grid lines with a single set of array operations, making
    it efficient both on GPU (CuPy) and CPU (NumPy).

    The algorithm is identical to perform_qlb_line_gpu but operates on the
    entire field array at once:
        1. Rotate   psi → psi_rot = einsum('ij,...j->...i', R_inv, psi)
        2. Collide  psi_rot → Q psi_rot  (axis-dependent sign pattern)
        3. Stream   shift along the sweep axis
        4. Rotate back  psi_out = einsum('ij,...j->...i', R_op, psi_str)

    Parameters
    ----------
    psi      : array (NX, NY, NZ, 4) – full spinor field on GPU/CPU
    V_field  : array (NX, NY, NZ)    – potential energy [J] at every cell
    axis     : str  – 'x' (open BCs), 'y' (periodic), or 'z' (periodic)
    R_inv    : array (4,4) – R⁻¹ rotation matrix for this axis
    R_op     : array (4,4) – R  rotation matrix for this axis
    periodic : bool – True → periodic wrap; False → absorbing open edges

    Returns
    -------
    psi_out : array (NX, NY, NZ, 4) – updated field (same device as input)
    """
    xp = cp if GPU_AVAILABLE else np

    # ------------------------------------------------------------------
    # STEP 1: collision coefficients at every cell  (NX, NY, NZ) arrays
    # ------------------------------------------------------------------
    m_tilde_3 = M_PARTICLE * V_F**2 * DT / (3 * HBAR) if HBAR > 0 else 0.0
    # g_tilde_3[i,j,k] = V_field[i,j,k]*DT / (3*hbar) — already in Joules, no Q_ELECTRON
    g_tilde_3 = V_field * DT / (3 * HBAR) if HBAR > 0 else xp.zeros_like(V_field)

    Omega_3   = m_tilde_3**2 - g_tilde_3**2
    denom     = 1 + Omega_3 / 4 - 1j * g_tilde_3
    denom     = xp.where(xp.abs(denom) < 1e-15, 1e-15, denom)
    a_hat     = (1 - Omega_3 / 4) / denom   # (NX,NY,NZ) – diagonal factor
    b_hat     = m_tilde_3 / denom            # scalar/array – mass coupling (0 if massless)

    # ------------------------------------------------------------------
    # STEP 2: rotate entire field to characteristic frame
    # einsum 'ij,...j->...i' applies the (4,4) matrix to the last axis of
    # psi at every (i,j,k) cell simultaneously.
    # ------------------------------------------------------------------
    psi_rot = xp.einsum('ij,...j->...i', R_inv, psi)   # (NX, NY, NZ, 4)

    u1 = psi_rot[..., 0]   # (NX, NY, NZ)
    u2 = psi_rot[..., 1]
    d1 = psi_rot[..., 2]
    d2 = psi_rot[..., 3]

    # ------------------------------------------------------------------
    # STEP 3: collision — sign pattern depends on axis
    #   Y-sweep: Y=I → apply Q̂ directly (minus on u1–d2 coupling)
    #   X/Z-sweep: rotated frame → X⁻¹Q̂X = Z⁻¹Q̂Z = Q (plus on u1–d2 coupling)
    # ------------------------------------------------------------------
    if axis == 'y':
        cu1 = a_hat * u1 - b_hat * d2
        cu2 = a_hat * u2 + b_hat * d1
        cd1 = a_hat * d1 - b_hat * u2
        cd2 = a_hat * d2 + b_hat * u1
    else:  # 'x' or 'z'
        cu1 = a_hat * u1 + b_hat * d2
        cu2 = a_hat * u2 + b_hat * d1
        cd1 = a_hat * d1 - b_hat * u2
        cd2 = a_hat * d2 - b_hat * u1

    # ------------------------------------------------------------------
    # STEP 4: stream along the sweep axis (array slices, no Python loop)
    #   u components shift +1 along the axis  (right-movers)
    #   d components shift -1 along the axis  (left-movers)
    # ------------------------------------------------------------------
    psi_str = xp.zeros_like(psi_rot)

    if axis == 'x':   # stream along dimension 0
        psi_str[1:,  :, :, 0] = cu1[:-1, :, :]   # u1 at k → k+1
        psi_str[1:,  :, :, 1] = cu2[:-1, :, :]
        psi_str[:-1, :, :, 2] = cd1[1:,  :, :]   # d1 at k → k-1
        psi_str[:-1, :, :, 3] = cd2[1:,  :, :]
        if periodic:
            psi_str[0,  :, :, 0] = cu1[-1, :, :]
            psi_str[0,  :, :, 1] = cu2[-1, :, :]
            psi_str[-1, :, :, 2] = cd1[0,  :, :]
            psi_str[-1, :, :, 3] = cd2[0,  :, :]
        # else open: edge cells left as 0 (inlet bounce-back/outlet absorbing
        # are applied externally after this function returns)

    elif axis == 'y':   # stream along dimension 1
        psi_str[:, 1:,  :, 0] = cu1[:, :-1, :]
        psi_str[:, 1:,  :, 1] = cu2[:, :-1, :]
        psi_str[:, :-1, :, 2] = cd1[:, 1:,  :]
        psi_str[:, :-1, :, 3] = cd2[:, 1:,  :]
        if periodic:
            psi_str[:, 0,  :, 0] = cu1[:, -1, :]
            psi_str[:, 0,  :, 1] = cu2[:, -1, :]
            psi_str[:, -1, :, 2] = cd1[:, 0,  :]
            psi_str[:, -1, :, 3] = cd2[:, 0,  :]

    else:  # 'z' — stream along dimension 2
        psi_str[:, :, 1:,  0] = cu1[:, :, :-1]
        psi_str[:, :, 1:,  1] = cu2[:, :, :-1]
        psi_str[:, :, :-1, 2] = cd1[:, :, 1: ]
        psi_str[:, :, :-1, 3] = cd2[:, :, 1: ]
        if periodic:
            psi_str[:, :, 0,  0] = cu1[:, :, -1]
            psi_str[:, :, 0,  1] = cu2[:, :, -1]
            psi_str[:, :, -1, 2] = cd1[:, :, 0 ]
            psi_str[:, :, -1, 3] = cd2[:, :, 0 ]

    # ------------------------------------------------------------------
    # STEP 5: rotate back to physical spinor basis
    # ------------------------------------------------------------------
    psi_out = xp.einsum('ij,...j->...i', R_op, psi_str)

    return psi_out


def run_simulation_gpu(n_steps=1000, output_freq=100, output_dir="palpacelli_gpu", save_snapshots=True):
    """
    Run the full Palpacelli 2012 simulation.

    The main time loop advances ψ by one operator-splitting step per iteration.
    Each step consists of:

      (A) X-sweep  — QLB sub-step along each x-line  (NX cells, open BCs)
      (B) X boundary conditions  — applied in the X-characteristic frame
      (C) Y-sweep  — QLB sub-step along each y-line  (NY cells, periodic BCs)
      (D) Z-sweep  — QLB sub-step along each z-line  (NZ cells, periodic BCs)
      (E) Sponge damping — exponential attenuation in the last 20 x-cells

    Operator splitting note
    -----------------------
    The 3D QLB scheme uses Lie-Trotter splitting  S_x · S_y · S_z  with the
    full DT and a 1/3 collision weight per sub-step (one third of the collision
    is applied in each of the three sweeps, so the full collision is applied
    once per time step).

    2D reduction (NZ=1, periodic Z)
    --------------------------------
    When NZ=1, the Z-streaming wraps trivially (one cell back to itself), so
    streaming is the identity.  The Z-sweep still applies 1/3 of the collision,
    which is essential — omitting it would reduce the effective barrier coupling
    to 2/3 of the specified value.

    Parameters
    ----------
    n_steps     : int   – number of time steps to run (1800 ≈ 1.73 ps matches paper)
    output_freq : int   – save diagnostics and snapshots every this many steps
    output_dir  : str   – directory for all output files
    save_snapshots : bool – whether to write per-frame PNG images
    """
    os.makedirs(output_dir, exist_ok=True)

    # Create snapshots subdirectory if saving snapshots
    if save_snapshots:
        snapshots_dir = os.path.join(output_dir, 'snapshots')
        os.makedirs(snapshots_dir, exist_ok=True)

    print("\n" + "="*80)
    print("STARTING GPU SIMULATION")
    print("="*80)

    # ------------------------------------------------------------------
    # [1/4] Set up the potential landscape: random square-barrier impurities
    # ------------------------------------------------------------------
    print("\n[1/4] Generating impurities...")
    V_field, impurity_positions = generate_random_impurities_gpu(
        DEFAULT_CONCENTRATION, DEFAULT_BARRIER_HEIGHT, seed=42
    )
    # V_field[i, j, 0] is the potential energy [J] at cell (i,j).
    # The seed is fixed so that runs are reproducible and directly comparable
    # to Fig. 8 of Palpacelli 2012.

    # ------------------------------------------------------------------
    # [2/4] Initialise wave packet to a right-moving Gaussian spinor
    # ------------------------------------------------------------------
    print("\n[2/4] Initializing wave packet...")
    psi = initialize_wave_packet_gpu()
    # psi has shape (NX, NY, NZ=1, 4) and lives on the GPU if available.

    # Record the peak probability density at t=0 so that snapshots use a
    # FIXED colour scale throughout the run.  Without this, the colourmap
    # would rescale whenever the packet spreads or exits, making the animation
    # look like everything is getting brighter even though probability is conserved.
    if GPU_AVAILABLE:
        psi_cpu_temp = cp.asnumpy(psi)
    else:
        psi_cpu_temp = psi
    prob_density_initial = np.sum(np.abs(psi_cpu_temp[:, :, NZ//2, :])**2, axis=-1)
    PROB_DENSITY_MAX = np.max(prob_density_initial) * 1.2  # 20% headroom

    # ------------------------------------------------------------------
    # [3/4] Storage arrays for time-series diagnostics (kept on CPU)
    # ------------------------------------------------------------------
    time_points        = []   # physical time t = step * DT  [s]
    total_probability  = []   # ∫|ψ|² dV  over full domain  (should ≈ 1)
    inlet_probability  = []   # probability in inlet region  (x < IMPURITY_START)
    impurity_probability = [] # probability in impurity region
    outlet_probability = []   # probability in outlet region  (x > IMPURITY_END)

    print("\n[3/4] Running QLB simulation...")
    print(f"  Steps: {n_steps}, Output freq: {output_freq}")

    start_time = time.time()

    for step in range(n_steps):

        xp = cp if GPU_AVAILABLE else np

        # ==============================================================
        # STAGE (A): X-SWEEP  (vectorized over all NY×NZ lines at once)
        # perform_qlb_sweep_gpu replaces the former 'for j in range(NY)' loop.
        # open BCs (periodic=False): edge cells are left as 0 here; the actual
        # inlet bounce-back and outlet absorbing conditions follow in stage (B).
        # ==============================================================
        psi = perform_qlb_sweep_gpu(psi, V_field, 'x',
                                     X_INV_MATRIX_GPU, X_MATRIX_GPU,
                                     periodic=False)

        # ==============================================================
        # STAGE (B): X-DIRECTION BOUNDARY CONDITIONS  (vectorized over NY×NZ)
        #
        # Inlet   (x = 0):  bounce-back — left-movers reflect as right-movers
        #   In characteristic frame: u ← d,  d = 0
        # Outlet  (x = NX-1):  open/absorbing — zero incoming left-movers
        #   In characteristic frame: d = 0
        # ==============================================================
        # Rotate both boundary slices to characteristic frame simultaneously
        # einsum 'ij,...j->...i' broadcasts over all (NY, NZ) transverse cells
        psi_rot_inlet  = xp.einsum('ij,...j->...i', X_INV_MATRIX_GPU, psi[0,  :, :, :])
        psi_rot_outlet = xp.einsum('ij,...j->...i', X_INV_MATRIX_GPU, psi[-1, :, :, :])

        # Inlet: bounce-back — copy d into u, then zero d
        psi_rot_inlet[..., 0] = psi_rot_inlet[..., 2].copy()  # u₁ ← d₁
        psi_rot_inlet[..., 1] = psi_rot_inlet[..., 3].copy()  # u₂ ← d₂
        psi_rot_inlet[..., 2] = 0.0
        psi_rot_inlet[..., 3] = 0.0

        # Outlet: absorbing — zero incoming left-movers
        psi_rot_outlet[..., 2] = 0.0
        psi_rot_outlet[..., 3] = 0.0

        # Rotate back to physical basis
        psi[0,  :, :, :] = xp.einsum('ij,...j->...i', X_MATRIX_GPU, psi_rot_inlet)
        psi[-1, :, :, :] = xp.einsum('ij,...j->...i', X_MATRIX_GPU, psi_rot_outlet)

        # ==============================================================
        # STAGE (C): Y-SWEEP  (vectorized over all NX×NZ lines at once)
        # Y = I so R_inv = R = identity; periodic BCs.
        # ==============================================================
        psi = perform_qlb_sweep_gpu(psi, V_field, 'y',
                                     Y_INV_MATRIX_GPU, Y_MATRIX_GPU,
                                     periodic=True)

        # ==============================================================
        # STAGE (D): Z-SWEEP  (vectorized over all NX×NY lines at once)
        # Periodic BCs.  For NZ=1 the streaming is the identity (one cell
        # wraps to itself), but the collision still applies the remaining
        # 1/3 of the barrier coupling not covered by X and Y sweeps.
        # ==============================================================
        psi = perform_qlb_sweep_gpu(psi, V_field, 'z',
                                     Z_INV_MATRIX_GPU, Z_MATRIX_GPU,
                                     periodic=True)

        # ==============================================================
        # STAGE (E): SPONGE LAYER  (vectorized, no Python loops)
        # Exponential damping in the last sponge_width x-cells.
        # i_local = 0 at inner edge (factor ≈ 1); sponge_width-1 at x=NX-1
        # (factor ≈ e⁻³ ≈ 0.05) — absorbs transmitted wave, minimises reflection.
        # ==============================================================
        sponge_width = 20
        sponge_start = max(0, NX - sponge_width)
        i_local  = xp.arange(sponge_start, NX) - sponge_start   # 0,1,...,sponge_width-1
        damping  = xp.exp(-3.0 * i_local / sponge_width)         # (sponge_width,)
        damping  = damping.reshape(-1, 1, 1, 1)                   # broadcast over NY,NZ,4
        psi[sponge_start:, :, :, :] *= damping

        # ==============================================================
        # DIAGNOSTICS: every output_freq steps, compute regional probabilities
        # ==============================================================
        if step % output_freq == 0:
            # Transfer to CPU for NumPy analysis
            if GPU_AVAILABLE:
                psi_cpu = cp.asnumpy(psi)
            else:
                psi_cpu = psi

            # Probability density [dimensionless per voxel]:
            #   ρ(i,j,k) = (Σ_α |ψ_α|²) · DX · DY · DZ
            # summing over the 4 spinor components α.
            # Integration over a region then gives the probability P = Σ ρ.
            prob_density = np.sum(np.abs(psi_cpu)**2, axis=-1) * DX * DY * DZ

            total_prob    = np.sum(prob_density)
            # Regional sums for the three domain segments (see config block)
            inlet_prob    = np.sum(prob_density[INLET_START:INLET_END,       :, :])
            impurity_prob = np.sum(prob_density[IMPURITY_START:IMPURITY_END, :, :])
            outlet_prob   = np.sum(prob_density[OUTLET_START:OUTLET_END,     :, :])

            time_points.append(step * DT)
            total_probability.append(total_prob)
            inlet_probability.append(inlet_prob)
            impurity_probability.append(impurity_prob)
            outlet_probability.append(outlet_prob)

            # ETA estimate
            elapsed      = time.time() - start_time
            time_per_step = elapsed / (step + 1)
            eta           = time_per_step * (n_steps - step - 1)

            print(f"  Step {step:5d}/{n_steps} | "
                  f"P={total_prob:.4f} | "
                  f"Inlet={inlet_prob:.3f} | "
                  f"Impurity={impurity_prob:.3f} | "
                  f"Outlet={outlet_prob:.3f} | "
                  f"ETA: {eta:.0f}s")

            # Save snapshot image for this frame
            if save_snapshots:
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

    # Transfer final state to CPU for plotting
    if GPU_AVAILABLE:
        psi_final    = cp.asnumpy(psi)
        V_field_cpu  = cp.asnumpy(V_field)
    else:
        psi_final    = psi
        V_field_cpu  = V_field

    # ------------------------------------------------------------------
    # [4/4] Post-processing: analysis plots + optional animation
    # ------------------------------------------------------------------
    print("\n[4/4] Creating plots...")
    plot_results(time_points, total_probability, inlet_probability,
                 impurity_probability, outlet_probability, psi_final, V_field_cpu, output_dir)

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
        'time':           np.array(time_points),
        'total_prob':     np.array(total_probability),
        'inlet_prob':     np.array(inlet_probability),
        'impurity_prob':  np.array(impurity_probability),
        'outlet_prob':    np.array(outlet_probability),
        'psi_final':      psi_final,
        'V_field':        V_field_cpu,
    }

def save_snapshot(psi_cpu, V_field_cpu, step, time_val, total_prob, inlet_prob,
                  impurity_prob, outlet_prob, output_dir, prob_density_max):
    """
    Write a two-panel PNG snapshot of the current simulation state.

    Left panel  — Probability density map:  ρ(x,y) = Σ_α |ψ_α(x,y)|² · DX·DY·DZ
                  Uses a FIXED colour scale anchored at the t=0 peak (prob_density_max).
                  Cyan dashed lines mark the impurity region boundaries.

    Right panel — Potential field [meV] with probability contours overlaid.
                  Provides visual comparison of wave-packet position vs. barriers.

    Parameters
    ----------
    psi_cpu          : ndarray (NX, NY, 1, 4) – spinor field on CPU
    V_field_cpu      : ndarray (NX, NY, 1)    – potential energy field [J]
    step             : int    – current time step number
    time_val         : float  – physical time [s]
    total_prob       : float  – ∫|ψ|² dV  at this step
    inlet_prob       : float  – probability in inlet region
    impurity_prob    : float  – probability in impurity region
    outlet_prob      : float  – probability in outlet region
    output_dir       : str    – directory to save the PNG
    prob_density_max : float  – peak |ψ|² (no vol. factor) from t=0 → sets fixed scale
    """
    
    fig = plt.figure(figsize=(14, 5))

    # ------------------------------------------------------------------
    # LEFT PANEL: probability density ρ(x,y)
    # ρ = Σ_α |ψ_α|² summed over spinor components, × voxel volume DX·DY·DZ
    # so that Σ_{ij} ρ(i,j) ≈ 1.0.  .T transposes so y is the vertical axis.
    # ------------------------------------------------------------------
    ax1 = plt.subplot(1, 2, 1)

    # Sum |ψ|² over 4 spinor components at the midplane z=NZ//2, scale by voxel volume.
    # For NZ=1 this is z=0; for NZ>1 it shows the central z-slice.
    prob_per_pixel = np.sum(np.abs(psi_cpu[:, :, NZ//2, :])**2, axis=-1).T * DX * DY * DZ

    # Fixed colour scale anchored to the t=0 peak: prevents false brightening
    # when the wave packet spreads and the local maximum drops.
    vmax_fixed = prob_density_max * DX * DY * DZ

    im1 = ax1.imshow(prob_per_pixel, origin='lower', cmap='hot', aspect='equal',
                     extent=[0, LX*1e9, 0, LY*1e9], vmin=0, vmax=vmax_fixed)
    # Cyan dashed lines mark the impurity region [IMPURITY_START, IMPURITY_END]
    ax1.axvline(x=IMPURITY_START*DX*1e9, color='cyan', linestyle='--', alpha=0.5, lw=1)
    ax1.axvline(x=IMPURITY_END  *DX*1e9, color='cyan', linestyle='--', alpha=0.5, lw=1)
    ax1.set_xlabel('x (nm)')
    ax1.set_ylabel('y (nm)')
    # Time displayed in picoseconds: DT ≈ 9.6e-16 s → 1000 steps ≈ 0.96 ps
    ax1.set_title(f'Probability Density | step={step} | t={time_val*1e12:.3f} ps')
    plt.colorbar(im1, ax=ax1, label='Probability per Pixel',
                 orientation='horizontal', pad=0.1)

    # ------------------------------------------------------------------
    # RIGHT PANEL: potential landscape V(x,y) [meV]  +  ρ contours overlaid
    # Convert from Joules to meV:  V [meV] = V [J] / Q_ELECTRON × 1e3
    # White contours of ρ show packet position relative to barriers.
    # ------------------------------------------------------------------
    ax2 = plt.subplot(1, 2, 2)
    V_slice = V_field_cpu[:, :, NZ//2].T / Q_ELECTRON * 1e3   # midplane (NY, NX), [meV]
    im2 = ax2.imshow(V_slice, origin='lower', cmap='viridis', aspect='equal',
                     extent=[0, LX*1e9, 0, LY*1e9], alpha=0.7)

    # Overlay ρ(x,y) contours: packet position visible against barrier pattern
    X, Y = np.meshgrid(np.linspace(0, LX*1e9, NX), np.linspace(0, LY*1e9, NY))
    ax2.contour(X, Y, prob_per_pixel, levels=5, colors='white', alpha=0.6, linewidths=1)

    ax2.set_xlabel('x (nm)')
    ax2.set_ylabel('y (nm)')
    ax2.set_title(f'Potential + Wave Packet | t={time_val*1e12:.3f} ps')
    plt.colorbar(im2, ax=ax2, label='Potential (meV)',
                 orientation='horizontal', pad=0.1)

    # Numeric summary text box (upper-right corner)
    stats_text = (f'Step: {step}\nP_total: {total_prob:.4f}\n'
                  f'Inlet: {inlet_prob:.3f}\nImpurity: {impurity_prob:.3f}\n'
                  f'Outlet: {outlet_prob:.3f}')
    ax2.text(0.98, 0.98, stats_text, transform=ax2.transAxes,
             fontsize=9, verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'snapshot_{step:05d}.png'),
                dpi=100, bbox_inches='tight')
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
    """
    Create the six-panel analysis figure saved as results.png.

    Panel layout (2 rows × 3 columns):
      [0,0]  Regional probability vs step   — P_inlet / impurity / outlet vs time
      [0,1]  Transmission coefficient T     — T = P_outlet / P_total  (0→1)
      [0,2]  Conservation error             — |P_total − 1| on log scale
      [1,0]  Final probability density map  — ρ(x,y) heat map at last recorded step
      [1,1]  Impurity potential map         — V(x,y) in meV
      [1,2]  Text summary                   — key parameters + final state numbers

    Transmission definition
    -----------------------
    T = P_outlet / max(P_total, ε)  (ε = 1e-10 avoids division by zero)
    Converges from 0→1 as the packet flows inlet→outlet; equals the fraction
    of the transmitted probability at any time.  For massless graphene carriers
    at normal incidence, T = 1 exactly (perfect Klein tunneling).

    Parameters
    ----------
    time_pts      : list[float] – physical times [s] at each diagnostic step
    total_prob    : list[float] – ∫|ψ|² dV at each recorded step
    inlet_prob    : list[float] – P in inlet region     (x < IMPURITY_START)
    impurity_prob : list[float] – P in impurity region  (IMPURITY_START ≤ x < IMPURITY_END)
    outlet_prob   : list[float] – P in outlet region    (x ≥ IMPURITY_END)
    psi           : ndarray (NX, NY, 1, 4) – final spinor field on CPU
    V_field       : ndarray (NX, NY, 1)    – potential energy field [J]
    output_dir    : str   – directory to write results.png
    """
    
    time_ps = np.array(time_pts) * 1e12  # s → ps for readable axis labels
    n_steps_total = len(time_pts)
    # Recover integer step numbers from physical times (DT = DX/V_F)
    step_numbers = [int(round(t / DT)) for t in time_pts]

    fig = plt.figure(figsize=(16, 10))

    # ------------------------------------------------------------------
    # PANEL [0,0]: regional probability time series
    # P_total ≈ 1 always (conservation check).
    # P_inlet decreases as the packet travels right.
    # P_outlet increases when the packet crosses the impurity region.
    # ------------------------------------------------------------------
    ax1 = plt.subplot(2, 3, 1)
    ax1.plot(step_numbers, total_prob,    'k-', lw=2, label='Total')
    ax1.plot(step_numbers, inlet_prob,    'b-',       label='Inlet')
    ax1.plot(step_numbers, impurity_prob, 'r-',       label='Impurity')
    ax1.plot(step_numbers, outlet_prob,   'g-',       label='Outlet')
    ax1.set_xlabel('Time step')
    ax1.set_ylabel('Probability')
    ax1.set_title('Regional Probability Evolution')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # ------------------------------------------------------------------
    # PANEL [0,1]: transmission coefficient T = P_outlet / P_total
    # For massless Dirac carriers at normal incidence, T = 1 exactly
    # (perfect Klein tunneling).  At finite impurity concentration and
    # mixed incident angles, T < 1.  Compare with paper Fig. 8.
    # ------------------------------------------------------------------
    transmission = np.array(outlet_prob) / (np.maximum(np.array(total_prob), 1e-10))
    ax2 = plt.subplot(2, 3, 2)
    ax2.plot(step_numbers, transmission, 'g-', lw=2)
    ax2.set_xlabel('Time step')
    ax2.set_ylabel('T = P_outlet / P_total')
    ax2.set_title('Transmission Coefficient')
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim([0, 1.1])
    
    # ------------------------------------------------------------------
    # PANEL [0,2]: probability conservation error |P_total − 1| (log scale)
    # QLB collision is exactly unitary → error should be near floating-point
    # precision in the lossless interior.  The error grows once the packet
    # starts leaving through the outlet (absorbed by the sponge layer).
    # ------------------------------------------------------------------
    ax3 = plt.subplot(2, 3, 3)
    error = np.abs(np.array(total_prob) - 1.0)
    ax3.semilogy(step_numbers, error, 'k-', lw=2)
    ax3.set_xlabel('Time step')
    ax3.set_ylabel('|P - 1|')
    ax3.set_title('Conservation Error')
    ax3.grid(True, alpha=0.3)

    # ------------------------------------------------------------------
    # PANEL [1,0]: final probability density ρ(x,y)
    # ρ = Σ_α |ψ_α|² · DX·DY·DZ  (summed over 4 spinor components α)
    # Cyan dashed lines mark the impurity region boundaries.
    # ------------------------------------------------------------------
    ax4 = plt.subplot(2, 3, 4)
    prob_final = np.sum(np.abs(psi[:, :, NZ//2, :])**2, axis=-1).T * DX * DY * DZ
    im = ax4.imshow(prob_final, origin='lower', cmap='hot', aspect='equal',
                    extent=[0, LX*1e9, 0, LY*1e9])
    ax4.axvline(x=IMPURITY_START*DX*1e9, color='cyan', linestyle='--', alpha=0.7)
    ax4.axvline(x=IMPURITY_END  *DX*1e9, color='cyan', linestyle='--', alpha=0.7)
    ax4.set_xlabel('x (nm)')
    ax4.set_ylabel('y (nm)')
    ax4.set_title(f'Final Probability Density (step {step_numbers[-1]})')
    plt.colorbar(im, ax=ax4)

    # ------------------------------------------------------------------
    # PANEL [1,1]: impurity potential field V(x,y) [meV]
    # Displayed in meV: V [meV] = V [J] / Q_ELECTRON × 1e3.
    # Each non-zero square is one impurity of side d_cells × d_cells cells.
    # ------------------------------------------------------------------
    ax5 = plt.subplot(2, 3, 5)
    V_slice = V_field[:, :, NZ//2].T / Q_ELECTRON * 1e3   # midplane [meV], shape (NY, NX)
    im = ax5.imshow(V_slice, origin='lower', cmap='viridis', aspect='equal',
                    extent=[0, LX*1e9, 0, LY*1e9])
    ax5.set_xlabel('x (nm)')
    ax5.set_ylabel('y (nm)')
    ax5.set_title('Impurity Potential Field (meV)')
    plt.colorbar(im, ax=ax5)

    # ------------------------------------------------------------------
    # PANEL [1,2]: plain-text simulation summary
    # ------------------------------------------------------------------
    ax6 = plt.subplot(2, 3, 6)
    ax6.axis('off')
    summary = (
        "SIMULATION SUMMARY\n\n"
        f"Grid:   {NX} × {NY}  (full paper res.)\n"
        f"Domain: {LX*1e9:.0f} × {LY*1e9:.0f} nm\n"
        f"σ = {SIGMA_LATTICE} cells  |  d = {IMPURITY_SIZE_CELLS} cells\n"
        f"C = {DEFAULT_CONCENTRATION*100:.1f}%  |  E = 80 meV  |  m = 0\n"
        f"V_barrier = {DEFAULT_BARRIER_HEIGHT/Q_ELECTRON*1e3:.0f} meV\n\n"
        f"FINAL STATE  (step {step_numbers[-1]},  t = {time_ps[-1]:.3f} ps)\n"
        f"  Total Prob:     {total_prob[-1]:.6f}\n"
        f"  Inlet:          {inlet_prob[-1]:.4f}\n"
        f"  Impurity:       {impurity_prob[-1]:.4f}\n"
        f"  Outlet:         {outlet_prob[-1]:.4f}\n\n"
        f"  T = P_out/P_tot = {transmission[-1]:.4f}\n\n"
        f"CONSERVATION\n"
        f"  Max  |ΔP|: {np.max(error):.2e}\n"
        f"  Final|ΔP|: {error[-1]:.2e}\n"
    )
    ax6.text(0.05, 0.5, summary, fontsize=10, family='monospace',
             verticalalignment='center')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'results.png'), dpi=150)
    plt.close()

if __name__ == "__main__":
    # ------------------------------------------------------------------
    # Entry point: run the Palpacelli 2012 Klein-tunneling simulation.
    #
    # n_steps = 1800  →  t_final = 1800 × DT ≈ 1800 × 9.6e-16 s ≈ 1.73 ps
    # This matches the simulation time used in Fig. 8 of the paper, where
    # the wave packet has fully traversed the 200-nm impurity region.
    #
    # output_freq = 50  →  36 diagnostic / snapshot frames over the run.
    # Increase to reduce disk usage; decrease for smoother animations.
    #
    # Output directory structure:
    #   palpacelli_gpu/
    #     snapshots/   — one PNG per output step (used for animation)
    #     animation.gif — animated wave-packet evolution
    #     results.png   — six-panel analysis figure
    # ------------------------------------------------------------------
    results = run_simulation_gpu(n_steps=1800, output_freq=50)

    # Final one-line summary printed to stdout
    final_T = results['outlet_prob'][-1] / max(results['total_prob'][-1], 1e-10)
    print(f"\nFinal Results:")
    print(f"  Transmission T = P_outlet / P_total = {final_T:.4f}")
    print(f"  Probability conservation (should be ≈ 1.0): {results['total_prob'][-1]:.6f}")
