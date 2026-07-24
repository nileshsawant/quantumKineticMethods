"""
2D classical-vs-circuit overlay: oblique Klein tunneling through a barrier.
==========================================================================

A massless Dirac wave packet is launched at an angle toward a vertical potential
barrier and evolved two ways -- the classical 2D QLB algorithm and the ported
quantum circuit (qiskit-aer, GPU) -- and the probability density |psi(x,y)|^2 is
shown as heatmaps side by side.  Unlike normal incidence (perfect Klein
transmission), oblique incidence partially reflects, so the packet visibly splits
into transmitted and reflected lobes; the circuit reproduces the classical result
to machine precision.

Run:
    module load qiskit/aer-gpu
    python -m qlb_port.plot_validation_2d
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from . import operators as ops
from . import backend as bk
from . import twod

NX = NY = 5                    # 32 x 32 lattice
Nx = Ny = 2 ** NX
XB0, XBW = 18, 3               # barrier x-range [18, 21)
G = 0.9                        # barrier coupling


def initial_packet(x0=6, y0=16, sigma=3.0, kx=0.6, ky=0.5):
    sp = ops.X_ROTATION @ (np.array([0, 0, 1, 1], dtype=complex) / np.sqrt(2))
    X, Y = np.meshgrid(np.arange(Nx), np.arange(Ny), indexing="xy")   # (Ny, Nx)
    env = np.exp(-(((X - x0) ** 2 + (Y - y0) ** 2) / (2 * sigma ** 2))) * np.exp(1j * (kx * X + ky * Y))
    psi = np.zeros((Ny, Nx, 4), dtype=complex)
    for c in range(4):
        psi[:, :, c] = env * sp[c]
    return (psi / np.linalg.norm(psi)).reshape(-1)


def density(sv):
    return (np.abs(sv.reshape(Ny, Nx, 4)) ** 2).sum(axis=2)   # (Ny, Nx)


def main():
    print("Device:", bk.device_report())
    snaps = [0, 13, 26]
    V = twod.barrier_field_2d(NX, NY, range(XB0, XB0 + XBW), G)
    psi0 = initial_packet()

    # classical snapshots
    classical = {}
    psi = psi0.copy(); tprev = 0
    for t in snaps:
        for _ in range(t - tprev):
            psi = twod.classical_step_2d_potential(psi, NX, NY, V)
        tprev = t
        classical[t] = density(psi)
    # circuit snapshots (one transpile, one run)
    step = twod.sweep2d_circuit_potential(NX, NY, V)
    svs = bk.evolve_snapshots(step, psi0, snaps)
    quantum = {t: density(svs[t]) for t in snaps}

    vmax = max(classical[0].max(), classical[snaps[-1]].max())
    fig, axes = plt.subplots(len(snaps), 3, figsize=(12, 11))
    for r, t in enumerate(snaps):
        d_c, d_q = classical[t], quantum[t]
        diff = np.abs(d_c - d_q)
        for col, (dat, title) in enumerate(
                [(d_c, f"classical  t={t}"),
                 (d_q, f"circuit  t={t}"),
                 (diff, f"|difference|  (max {diff.max():.1e})")]):
            ax = axes[r, col]
            cmap = "inferno" if col < 2 else "viridis"
            vm = vmax if col < 2 else max(diff.max(), 1e-16)
            im = ax.imshow(dat, origin="lower", aspect="auto", cmap=cmap,
                           vmin=0, vmax=vm, extent=[0, Nx, 0, Ny])
            ax.axvspan(XB0, XB0 + XBW, color="cyan", alpha=0.25)
            ax.set_ylabel("y"); ax.set_xlabel("x")
            ax.set_title(title, fontsize=10)
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    # transmission/reflection at final time (from the circuit)
    dfin = quantum[snaps[-1]]
    refl = dfin[:, :XB0].sum(); trans = dfin[:, XB0 + XBW:].sum()
    maxdev = max(np.max(np.abs(classical[t] - quantum[t])) for t in snaps)
    fig.suptitle(f"2D oblique Klein tunneling: classical QLB vs quantum circuit (GPU)\n"
                 f"reflected R≈{refl:.2f}, transmitted T≈{trans:.2f}; "
                 f"classical vs circuit max |Δ| = {maxdev:.1e}", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = "qlb_port/validation_2d.png"
    fig.savefig(out, dpi=140)
    print(f"Saved {out}   R≈{refl:.2f} T≈{trans:.2f}  max|Δ|={maxdev:.2e}")


if __name__ == "__main__":
    main()
