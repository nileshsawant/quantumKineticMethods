"""
3D classical-vs-circuit overlay: Dirac wave packet through a planar barrier.
===========================================================================

A massless 3D Dirac wave packet is launched (with a transverse component) toward a
planar potential barrier and evolved with both the classical 3D QLB algorithm and
the ported quantum circuit (qiskit-aer, GPU).  Because a full 3D density is hard to
draw, we show two projections of |psi(x,y,z)|^2:

  * top row: the x-profile (summed over y and z) at several times -- classical line
    vs circuit markers -- with the barrier shaded; and
  * bottom row: the xy-projection (summed over z) at the final time, classical vs
    circuit vs |difference|.

Run:
    module load qiskit/aer-gpu
    python -m qlb_port.plot_validation_3d
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from . import operators as ops
from . import backend as bk
from . import threed

NX, NY, NZ = 4, 3, 3           # 16 x 8 x 8 lattice (12 qubits)
Nx, Ny, Nz = 2 ** NX, 2 ** NY, 2 ** NZ
XB0, XBW, G = 9, 2, 0.9        # planar barrier at x in [9, 11)


def initial_packet(x0=2, sigma=1.6, kx=0.6, ky=0.4):
    sp = ops.X_ROTATION @ (np.array([0, 0, 1, 1], dtype=complex) / np.sqrt(2))
    Z, Y, X = np.meshgrid(np.arange(Nz), np.arange(Ny), np.arange(Nx), indexing="ij")
    env = np.exp(-(((X - x0) ** 2 + (Y - Ny // 2) ** 2 + (Z - Nz // 2) ** 2) / (2 * sigma ** 2)))
    env = env * np.exp(1j * (kx * X + ky * Y))
    psi = np.zeros((Nz, Ny, Nx, 4), dtype=complex)
    for c in range(4):
        psi[:, :, :, c] = env * sp[c]
    return (psi / np.linalg.norm(psi)).reshape(-1)


def dens_xyz(sv):
    return (np.abs(sv.reshape(Nz, Ny, Nx, 4)) ** 2).sum(axis=3)   # (Nz, Ny, Nx)


def main():
    print("Device:", bk.device_report())
    snaps = [0, 7, 14]
    V = threed.planar_barrier_field_3d(NX, NY, NZ, range(XB0, XB0 + XBW), G)
    psi0 = initial_packet()

    classical = {}
    psi = psi0.copy(); tprev = 0
    for t in snaps:
        for _ in range(t - tprev):
            psi = threed.classical_step_3d_potential(psi, NX, NY, NZ, V)
        tprev = t
        classical[t] = dens_xyz(psi)
    step = threed.sweep3d_circuit_potential(NX, NY, NZ, V)
    svs = bk.evolve_snapshots(step, psi0, snaps)
    quantum = {t: dens_xyz(svs[t]) for t in snaps}

    maxdev = max(np.max(np.abs(classical[t] - quantum[t])) for t in snaps)
    x = np.arange(Nx)
    colors = plt.cm.viridis(np.linspace(0, 0.8, len(snaps)))
    fig = plt.figure(figsize=(13, 8))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.1])

    # top: x-profile (sum over y,z), classical line + circuit markers
    axp = fig.add_subplot(gs[0, :])
    for t, col in zip(snaps, colors):
        axp.plot(x, classical[t].sum(axis=(0, 1)), "-", color=col, lw=2, label=f"classical t={t}")
        axp.plot(x, quantum[t].sum(axis=(0, 1)), "o", color=col, mfc="none", ms=5, mew=1.3,
                 label=f"circuit t={t}")
    axp.axvspan(XB0, XB0 + XBW, color="red", alpha=0.12, label="barrier")
    axp.set_xlabel("x"); axp.set_ylabel(r"$\sum_{y,z}|\psi|^2$")
    axp.set_title(f"3D Dirac through a planar barrier: x-profile (classical vs circuit)   "
                  f"max |Δ| = {maxdev:.1e}", fontsize=11)
    axp.legend(fontsize=7, ncol=3, loc="upper right"); axp.grid(alpha=0.3)

    # bottom: xy-projection (sum over z) at final time
    tf = snaps[-1]
    cxy, qxy = classical[tf].sum(0), quantum[tf].sum(0)   # (Ny, Nx)
    diff = np.abs(cxy - qxy)
    vmax = max(cxy.max(), qxy.max())
    for col, (dat, title, cmap, vm) in enumerate([
            (cxy, f"classical  t={tf}", "inferno", vmax),
            (qxy, f"circuit  t={tf}", "inferno", vmax),
            (diff, f"|difference| (max {diff.max():.0e})", "viridis", max(diff.max(), 1e-16))]):
        ax = fig.add_subplot(gs[1, col])
        im = ax.imshow(dat, origin="lower", aspect="auto", cmap=cmap, vmin=0, vmax=vm,
                       extent=[0, Nx, 0, Ny])
        ax.axvspan(XB0, XB0 + XBW, color="cyan", alpha=0.25)
        ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_title(title, fontsize=10)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    fig.suptitle("3D QLB: classical algorithm vs ported quantum circuit (GPU)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = "qlb_port/validation_3d.png"
    fig.savefig(out, dpi=140)
    print(f"Saved {out}   max|Δ|={maxdev:.2e}")


if __name__ == "__main__":
    main()
