"""Test which spinor combinations create right-moving eigenstates."""
import numpy as np
from dirac_qlb_solver import ALPHA_X, X_MATRIX, X_INV_MATRIX, C, M_PARTICLE, HBAR

# For massless Dirac: H = c * α_x * p_x
# Eigenvalue equation: H |ψ⟩ = E |ψ⟩
# c * α_x * p_x |ψ⟩ = E |ψ⟩

# For momentum p_x > 0 (right-moving):
# E = ±c * |p_x| (positive/negative energy branches)

print("Testing Dirac eigenstates for right-moving wave packet")
print("=" * 60)

# Test different spinor combinations
test_states = {
    "(1,0,0,0)": np.array([1, 0, 0, 0]),
    "(0,1,0,0)": np.array([0, 1, 0, 0]),
    "(0,0,1,0)": np.array([0, 0, 1, 0]),
    "(0,0,0,1)": np.array([0, 0, 0, 1]),
    "(1,0,1,0)/√2": np.array([1, 0, 1, 0]) / np.sqrt(2),
    "(1,0,0,1)/√2": np.array([1, 0, 0, 1]) / np.sqrt(2),
    "(0,1,1,0)/√2": np.array([0, 1, 1, 0]) / np.sqrt(2),
    "(0,1,0,1)/√2": np.array([0, 1, 0, 1]) / np.sqrt(2),
}

print("\nChecking which states are eigenstates of α_x:")
for name, psi in test_states.items():
    result = ALPHA_X @ psi
    # Check if result is proportional to psi (eigenstate)
    ratio = result / (psi + 1e-10)
    is_eigenstate = np.allclose(ratio[np.abs(psi) > 0.01], ratio[np.abs(psi) > 0.01][0])
    
    if is_eigenstate:
        eigenvalue = ratio[np.abs(psi) > 0.01][0].real
        print(f"  {name:20s} → EIGENSTATE with λ = {eigenvalue:+.2f}")
    else:
        print(f"  {name:20s} → not an eigenstate")

print("\n" + "=" * 60)
print("Characteristic basis decomposition (X_INV @ ψ):")
print("Components: (u1, u2, d1, d2)")
print("  u1, u2 → move RIGHT")
print("  d1, d2 → move LEFT")
print("=" * 60)

for name, psi in test_states.items():
    char_basis = X_INV_MATRIX @ psi
    u_components = np.abs(char_basis[:2])**2
    d_components = np.abs(char_basis[2:])**2
    right_prob = np.sum(u_components)
    left_prob = np.sum(d_components)
    
    print(f"\n{name:20s}")
    print(f"  Right (u1,u2): {right_prob:.3f}")
    print(f"  Left  (d1,d2): {left_prob:.3f}")
    if right_prob > 0.95:
        print(f"  → PURE RIGHT-MOVING ✓")
    elif left_prob > 0.95:
        print(f"  → PURE LEFT-MOVING")
    else:
        print(f"  → MIXED (will split!)")
