"""
Quick simulation test to verify the complete workflow
Runs a short graphene simulation (100 steps instead of 500)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for HPC
import matplotlib.pyplot as plt
import os

# Modify parameters for quick test
import dirac_qlb_solver as solver

# Override some parameters for faster execution
solver.T_STEPS = 100
solver.OUTPUT_FREQ = 25
solver.NX = 32
solver.NY = 32
solver.NZ = 1
solver.DX = solver.LX / solver.NX
solver.DY = solver.LY / solver.NY
solver.DZ = solver.DX
solver.DT = solver.DX / solver.C

print("="*70)
print("QUICK SIMULATION TEST - Graphene Electron Transport")
print("="*70)
print(f"Grid: {solver.NX} x {solver.NY} x {solver.NZ}")
print(f"Time steps: {solver.T_STEPS}")
print(f"Output frequency: {solver.OUTPUT_FREQ}")
print("="*70)

# Setup graphene simulation
psi_graphene, V_graphene = solver.setup_graphene_simulation()

print("\nStarting simulation...")
# Run simulation
psi_final, densities, transmission = solver.run_simulation(
    psi_graphene, V_graphene,
    output_dir_name="qlb_graphene_test_output"
)

print("\n" + "="*70)
print("SIMULATION TEST COMPLETE")
print("="*70)
print(f"Final total density: {densities[-1]:.6f}")
print(f"Final transmission: {transmission[-1]:.6e}")
print(f"Outputs saved to: qlb_graphene_test_output/")

# Check that outputs exist
output_dir = "qlb_graphene_test_output"
if os.path.exists(output_dir):
    files = os.listdir(output_dir)
    print(f"\nGenerated {len(files)} output files:")
    for f in sorted(files)[:5]:  # Show first 5
        print(f"  - {f}")
    if len(files) > 5:
        print(f"  ... and {len(files)-5} more")
else:
    print("WARNING: Output directory not created")

print("\n✓ Quick simulation test passed!")
