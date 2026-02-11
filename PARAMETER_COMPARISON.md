# Parameter Comparison: Palpacelli 2012 vs Current Implementation

## From Paper: "Klein Tunneling in the Presence of Random Impurities" (Palpacelli et al. 2012)

### Section 4: Klein Paradox in Random Media - Key Parameters

#### Physical Parameters:
- **Lattice spacing**: `Δx = 0.96 nm`
- **Energy**: `E = 0.117` (80 meV in physical units)
- **Wave packet spread**: `σ = 48` lattice spacings = `48 × 0.96 nm = 46.08 nm`
- **Particle mass**: 
  - Massless case: `m = 0`
  - Massive case: `m = 0.1` (mc² ≈ 0.07 eV = 70 meV)

#### Grid Configuration:
- **Total domain**: `2048 × 512` cells
- **Inlet region**: `[0, 512) × 512`
- **Impurity region**: `[512, 1536) × 512` (1024 cells long)
- **Outlet region**: `[1536, 2048] × 512`

#### Impurity Configuration:
- **Impurity size**: `d = 8` lattice spacings (square shape)
- **Impurity concentration**: `C = N·d²/A` where:
  - `N` = number of impurities
  - `A = Ly × Lz` = total area
- **Tested concentrations**: `C = 0.001, 0.005, 0.01, 0.05` (0.1%, 0.5%, 1%, 5%)
- **Barrier heights tested**: `V = 25, 50, 100, 200, 285 meV`

#### Boundary Conditions:
- **Top/Bottom**: Periodic
- **Inlet**: Bounce-back
- **Outlet**: Open (to avoid reflections)

---

## Current Implementation Parameters

### From `dirac_qlb_solver.py`:

```python
# Spatial dimensions
NX = 64   # Number of grid points in x-direction
NY = 64   # Number of grid points in y-direction
NZ = 1    # For 2D graphene example

# Physical domain size
LX = 1e-9  # 1 nm
LY = 1e-9  # 1 nm
LZ = 1e-9  # 1 nm

# Discretization
DX = LX / NX  # = 1.5625e-11 m = 0.015625 nm
DY = LY / NY
DZ = LZ / NZ if NZ > 1 else DX

# Particle mass
M_PARTICLE = 9.1093837e-31  # kg (electron mass, modified to 0 for graphene)

# Time step
DT = DX / C  # = 5.21e-20 s
```

### From `setup_graphene_simulation()`:

```python
# Massless Dirac fermions
M_PARTICLE = 0.0  # Massless for ideal graphene

# Impurity Configuration
impurity_concentration = 0.05  # 5% impurities
impurity_potential_strength = 0.1 * Q_ELECTRON  # 0.1 eV barrier

# Define impurity region
impurity_region_start_x = NX // 3  # = 21
impurity_region_end_x = 2 * NX // 3  # = 42
# Impurity region length = 21 cells

# Wave packet initialization
D0_graphene = 8 * DX  # Initial spread = 0.125 nm
initial_momentum_x_graphene = 1.0 * HBAR / D0_graphene
```

---

## Comparison Summary

| Parameter | Palpacelli 2012 | Current Implementation | Match? |
|-----------|-----------------|------------------------|---------|
| **Grid spacing** | 0.96 nm | 0.015625 nm | ❌ 61× smaller |
| **Domain size** | 2048 × 512 | 64 × 64 | ❌ Much smaller |
| **Wave packet σ** | 46.08 nm | 0.125 nm | ❌ 369× smaller |
| **Impurity size** | 7.68 nm (8 cells) | ~0.125 nm (8 cells) | ✓ Same cell count |
| **Mass (graphene)** | m = 0 | m = 0 | ✓ Correct |
| **Impurity conc.** | 0.1% - 5% | 5% | ✓ Within range |
| **Potential** | 25-285 meV | 100 meV | ✓ Within range |
| **Impurity region** | 50% of domain | 33% of domain | ~ Similar |
| **Boundary conditions** | Periodic (y), bounce-back (inlet), open (outlet) | Periodic (all) | ❌ Different |

---

## Key Discrepancies

### 1. **Scale Mismatch (CRITICAL)**
The current implementation operates at **~60× smaller spatial scale**:
- **Paper**: dx = 0.96 nm, domain = ~2 μm × 0.5 μm
- **Current**: dx = 0.016 nm, domain = 1 nm × 1 nm

This means:
- Wave packet is ~370× smaller relative to the domain
- Time scale is also ~60× faster
- Physical phenomena may be under-resolved

### 2. **Grid Resolution**
- **Paper**: 2048 × 512 = 1,048,576 cells
- **Current**: 64 × 64 = 4,096 cells (256× fewer cells)

### 3. **Boundary Conditions**
- **Paper**: Open outlet boundary (realistic for transmission measurement)
- **Current**: Periodic boundaries (unphysical for transmission studies)

### 4. **Wave Packet Size**
- **Paper**: σ = 48 cells = 46 nm (extends over multiple impurities)
- **Current**: σ = 8 cells = 0.125 nm (comparable to single impurity)

This is crucial because:
- In the paper, wave packets can "split and flow around" obstacles
- In current implementation, wave packet is too small for this behavior

---

## Recommendations

### To Match Palpacelli 2012 Setup:

1. **Increase spatial scale**:
   ```python
   DX = 0.96e-9  # 0.96 nm (instead of 0.015625 nm)
   NX = 2048
   NY = 512
   LX = NX * DX  # ~2 μm
   LY = NY * DX  # ~0.5 μm
   ```

2. **Adjust wave packet spread**:
   ```python
   D0_graphene = 48 * DX  # 46.08 nm
   ```

3. **Set impurity size**:
   ```python
   impurity_size = 8  # cells (7.68 nm)
   ```

4. **Implement proper boundary conditions**:
   - Top/Bottom: Periodic
   - Left (inlet): Bounce-back or absorbing
   - Right (outlet): Open/transmitting boundary

5. **Update impurity region**:
   ```python
   impurity_region_start = NX // 4  # ~512 cells
   impurity_region_end = 3 * NX // 4  # ~1536 cells
   ```

6. **Energy considerations**:
   - Paper uses E = 80 meV for massless
   - This corresponds to k₀ ≈ 1 in lattice units
   - Adjust initial momentum accordingly

---

## Current Status

The current implementation is a **scaled-down demonstration** suitable for:
- ✅ Testing the QLB algorithm
- ✅ Verifying probability conservation
- ✅ Educational purposes
- ✅ Quick validation runs

But it **does not match** the Palpacelli 2012 Klein tunneling study in terms of:
- ❌ Physical scale (60× too small)
- ❌ Grid resolution (256× fewer cells)
- ❌ Wave packet to impurity size ratio
- ❌ Boundary conditions for transmission studies

To replicate the paper's results on Klein tunneling and conductivity loss, the scale must be increased significantly.
