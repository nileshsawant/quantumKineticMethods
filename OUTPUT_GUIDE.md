# Simulation Output Guide

## Updated Features

The simulation now saves **time-resolved snapshots** during execution, not just final results!

## Output Structure

```
palpacelli_gpu/
├── results.png          # Final comprehensive analysis (6 plots)
├── animation.gif        # Animated GIF showing evolution (created automatically)
└── snapshots/           # Individual frames at each output step
    ├── snapshot_00000.png  # Initial state (t=0)
    ├── snapshot_00100.png  # Step 100
    ├── snapshot_00200.png  # Step 200
    ├── snapshot_00300.png  # Step 300
    └── ...                 # Up to snapshot_02000.png
```

## Snapshot Content

Each snapshot PNG (90-110 KB) contains:

**Left panel:** Probability density heatmap
- Shows wave packet location and intensity
- Cyan dashed lines mark impurity region boundaries

**Right panel:** Potential field + wave packet overlay
- Colored background: Potential barriers (meV)
- White contours: Wave packet probability density
- Stats box: Current step, probabilities in each region

## Final Results (results.png)

Comprehensive 6-panel analysis:
1. **Regional Probability vs Time** - Inlet/Impurity/Outlet evolution
2. **Transmission Coefficient** - T(t) from 0 to 1
3. **Conservation Error** - |P - 1| (should be < 10⁻¹⁵)
4. **Final Probability Density** - 2D heatmap
5. **Potential Field** - Random impurity distribution
6. **Summary Statistics** - Final state numbers

## Animation (animation.gif)

Automatically generated from all snapshots:
- Duration: ~200ms per frame
- Loops continuously
- Shows complete wave packet evolution through random impurities
- File size: ~2-5 MB for 20 frames

## Customization

Edit `run_palpacelli_gpu.py`:

```python
# At the end of file, modify:
results = run_simulation_gpu(
    n_steps=2000,        # Total simulation steps
    output_freq=100,     # Save snapshot every 100 steps
    save_snapshots=True  # Set False to disable snapshots
)
```

**Trade-offs:**
- `output_freq=50` → More snapshots, slower simulation, larger storage
- `output_freq=200` → Fewer snapshots, faster simulation, smaller storage
- `save_snapshots=False` → Only final results.png, fastest

## File Sizes

For 2000 steps with `output_freq=100` (20 snapshots):

| Item | Size | Count |
|------|------|-------|
| Individual snapshot | ~100 KB | 20 files |
| Total snapshots | ~2 MB | - |
| Animation GIF | ~3-5 MB | 1 file |
| Final results | ~373 KB | 1 file |
| **Total** | **~5-7 MB** | - |

## Viewing Results

### On local machine (after scp):
```bash
# View snapshots
eog palpacelli_gpu/snapshots/snapshot_*.png

# View animation
firefox palpacelli_gpu/animation.gif
# or
eog palpacelli_gpu/animation.gif

# View final results
eog palpacelli_gpu/results.png
```

### On remote server (with X11 forwarding):
```bash
ssh -X user@server
display palpacelli_gpu/animation.gif
```

### Download all results:
```bash
scp -r user@server:/path/to/palpacelli_gpu/ ./
```

## Performance Impact

Saving snapshots adds ~5-10% overhead:
- Without snapshots: ~23.8 steps/s → 84s total
- With snapshots (freq=100): ~22 steps/s → 91s total
- Minimal impact, highly worthwhile for visualization!

## Klein Tunneling Observations

Watch for in snapshots:
1. **Wave packet initialization** in inlet region (step 0)
2. **Propagation** toward impurity region (steps 100-500)
3. **Scattering** at random barriers (steps 500-1000)
4. **Transmission/reflection** dynamics (steps 1000-1500)
5. **Periodic revivals** from boundary reflections (steps 1500-2000)

The probability numbers in each snapshot show how the wave packet distributes across regions over time, revealing Klein tunneling and quantum scattering behavior!
