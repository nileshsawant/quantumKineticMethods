# Klein Tunneling Simulation - Palpacelli 2012

## Quick Start

Run the Klein tunneling simulation (1/4 resolution, ~90 seconds):
```bash
module load anaconda3/2024.06.1
conda activate ./dirac_qlb_env
python run_palpacelli_gpu.py
```

Results saved to: `palpacelli_gpu/results.png`

## What This Simulates

Klein tunneling of massless Dirac fermions (like electrons in graphene) through a random distribution of potential barriers. This replicates:

**Palpacelli, S., et al. (2012). "Klein Tunneling in the Presence of Random Impurities"**

## Configuration

| Parameter | Value | Paper Value |
|-----------|-------|-------------|
| Grid | 512 × 128 | 2048 × 512 (1/4 resolution) |
| Lattice spacing | 0.96 nm | 0.96 nm ✓ |
| Domain | 0.49 × 0.12 μm | 1.97 × 0.49 μm |
| Wave packet σ | 12 cells | 48 cells (scaled) |
| Impurities | 409 random, 2×2 cells | 409 random, 8×8 cells |
| Concentration | 5% | 5% ✓ |
| Barrier height | 100 meV | 100 meV ✓ |

## Performance

- **CPU**: ~23.8 steps/s (2000 steps in 84s)
- **GPU (with CuPy)**: ~100-200 steps/s estimated (10× faster)

To enable GPU:
```bash
pip install cupy-cuda12x  # Match your CUDA version
```

## Files

- `palpacelli2012_config.py` - Full paper parameters
- `run_palpacelli_gpu.py` - Main simulation script (CPU/GPU)
- `palpacelli_gpu/results.png` - Simulation results (6 plots)

## Key Results

✅ **Probability conserved**: P = 1.0000 (to 4 decimal places)  
✅ **Klein tunneling observed**: Wave packet propagates through barriers  
✅ **Scattering visible**: Probability reflects/transmits at impurities  
✅ **Periodic dynamics**: Wave bounces between boundaries  

## Physics

The simulation solves the 3D Dirac equation using the Quantum Lattice Boltzmann (QLB) method:

```
iℏ ∂ψ/∂t = (c α·∇ + mc² β + V(r) β) ψ
```

For graphene: m = 0 (massless Dirac fermions)

## Next Steps

1. **Parameter sweep**: Test different concentrations (0.1%, 0.5%, 1%, 5%)
2. **Barrier heights**: Test 25, 50, 100, 200, 285 meV
3. **Full resolution**: Run 2048×512 grid with GPU (5-10 min)
4. **Compare with paper**: Validate transmission vs. time (Figure 5)

## References

- Palpacelli, S., Mendoza, M., Herrmann, H. J., & Succi, S. (2012). Klein Tunneling in the Presence of Random Impurities. Int. J. Mod. Phys. C, 23(12).
- Paper PDF: `kleinTunnelingPresenceRandomImpurities_palpacelli_klein_2012.pdf`
