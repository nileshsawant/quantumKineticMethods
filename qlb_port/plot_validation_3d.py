"""
3D classical-vs-circuit overlay: a wave packet sent along the body diagonal.
===========================================================================

A massless 3D Dirac wave packet is launched along the (1,1,1) body diagonal of a
cubic lattice and evolved with both the classical 3D QLB algorithm and the ported
quantum circuit (qiskit-aer, GPU).  The "diagonal mover" is the +1 eigenstate of the
velocity operator (ALPHA_X + BETA + ALPHA_Z)/sqrt(3) (in Dellar's representation the
per-axis velocity operators are ALPHA_X, BETA, ALPHA_Z), so it has equal group
velocity along x, y and z and all three sweeps' streaming are genuinely exercised.

The full 3D density is shown as three orthogonal projections (xy, xz, yz) at the
start and end times; the packet moves from one corner toward the opposite, and the
circuit reproduces the classical result to machine precision.

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

NX = NY = NZ = 4                # 16^3 cubic lattice (14 qubits)
Nx = Ny = Nz = 2 ** NX
T = 8
K0 = 0.5


def diagonal_spinor():
    """+1 eigenstate of (ALPHA_X + BETA + ALPHA_Z)/sqrt(3): equal velocity on x,y,z."""
    A = ops.ALPHA_X + ops.BETA + ops.ALPHA_Z
    w, V = np.linalg.eigh(A)
    return V[:, int(np.argmax(w))]


def initial_packet(x0=4, sigma=2.0):
    sp = diagonal_spinor()
    Z, Y, X = np.meshgrid(np.arange(Nz), np.arange(Ny), np.arange(Nx), indexing="ij")
    env = np.exp(-(((X - x0) ** 2 + (Y - x0) ** 2 + (Z - x0) ** 2) / (2 * sigma ** 2)))
    env = env * np.exp(1j * K0 * (X + Y + Z))
    psi = np.zeros((Nz, Ny, Nx, 4), dtype=complex)
    for c in range(4):
        psi[:, :, :, c] = env * sp[c]
    return (psi / np.linalg.norm(psi)).reshape(-1)


def dens(sv):
    return (np.abs(sv.reshape(Nz, Ny, Nx, 4)) ** 2).sum(axis=3)   # (Nz, Ny, Nx)


def com(d):
    tot = d.sum()
    return (float((d.sum((0, 1)) * np.arange(Nx)).sum() / tot),
            float((d.sum((0, 2)) * np.arange(Ny)).sum() / tot),
            float((d.sum((1, 2)) * np.arange(Nz)).sum() / tot))


def main():
    print("Device:", bk.device_report())
    psi0 = initial_packet()

    # classical densities at t=0 and t=T
    pc = psi0.copy()
    for _ in range(T):
        pc = threed.classical_step_3d(pc, NX, NY, NZ, 0.0)
    d_c0, d_cT = dens(psi0), dens(pc)
    # circuit densities at t=0 and t=T (one transpile, one run)
    svs = bk.evolve_snapshots(threed.sweep3d_circuit(NX, NY, NZ), psi0, [0, T])
    d_q0, d_qT = dens(svs[0]), dens(svs[T])

    maxdev = max(np.max(np.abs(d_c0 - d_q0)), np.max(np.abs(d_cT - d_qT)))
    c0, cT = com(d_q0), com(d_qT)
    print(f"COM start {tuple(f'{v:.2f}' for v in c0)} -> end {tuple(f'{v:.2f}' for v in cT)}"
          f"   max|Δ|={maxdev:.2e}")

    fig, axes = plt.subplots(2, 3, figsize=(13, 8.5))

    # ---- Row 1: orthogonal projections at t=T -------------------------------
    # filled heatmap = classical, overlaid contours = quantum circuit.
    def proj(d):
        return [d.sum(0), d.sum(1), d.sum(2)]      # xy (Ny,Nx), xz (Nz,Nx), yz (Nz,Ny)
    labels = [("x", "y"), ("x", "z"), ("y", "z")]
    pc_xy, qc_xy = proj(d_cT), proj(d_qT)
    vmax = max(p.max() for p in pc_xy)
    for j in range(3):
        ax = axes[0, j]
        im = ax.imshow(pc_xy[j], origin="lower", aspect="auto", cmap="inferno",
                       vmin=0, vmax=vmax, extent=[0, pc_xy[j].shape[1], 0, pc_xy[j].shape[0]])
        lv = np.linspace(0.15, 0.95, 5) * vmax
        ax.contour(np.linspace(0.5, pc_xy[j].shape[1] - 0.5, pc_xy[j].shape[1]),
                   np.linspace(0.5, pc_xy[j].shape[0] - 0.5, pc_xy[j].shape[0]),
                   qc_xy[j], levels=lv, colors="cyan", linewidths=1.0)
        ax.set_xlabel(labels[j][0]); ax.set_ylabel(labels[j][1])
        ax.set_title(f"{labels[j][0]}{labels[j][1]}-projection, t={T}\n"
                     f"fill = classical, cyan lines = circuit", fontsize=9)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    # ---- Row 2: 1D marginals along x, y, z (classical line + circuit markers)
    axis_len = [Nx, Ny, Nz]
    axis_name = ["x", "y", "z"]
    marg = lambda d, k: d.sum(axis=tuple(a for a in (0, 1, 2) if a != (2 - k)))  # k=0->x,1->y,2->z
    col0, colT = "0.55", "C3"
    for k in range(3):
        ax = axes[1, k]
        r = np.arange(axis_len[k])
        ax.plot(r, marg(d_c0, k), "-", color=col0, lw=2, label="classical t=0")
        ax.plot(r, marg(d_q0, k), "o", color=col0, mfc="none", ms=5, mew=1.2, label="circuit t=0")
        ax.plot(r, marg(d_cT, k), "-", color=colT, lw=2, label=f"classical t={T}")
        ax.plot(r, marg(d_qT, k), "o", color=colT, mfc="none", ms=5, mew=1.2, label=f"circuit t={T}")
        ax.axvline(c0[k], color=col0, ls=":", lw=1)
        ax.axvline(cT[k], color=colT, ls=":", lw=1)
        ax.set_xlabel(f"{axis_name[k]}")
        ax.set_ylabel(rf"$\sum_{{\neq {axis_name[k]}}}|\psi|^2$")
        ax.set_title(f"{axis_name[k]}-marginal: COM {c0[k]:.1f} → {cT[k]:.1f}", fontsize=9)
        ax.grid(alpha=0.3)
        if k == 0:
            ax.legend(fontsize=7, loc="upper right")

    fig.suptitle(
        "3D Dirac packet along the (1,1,1) diagonal: classical QLB vs quantum circuit (GPU)\n"
        f"COM {tuple(round(v, 1) for v in c0)} → {tuple(round(v, 1) for v in cT)} "
        f"(equal motion in x, y, z);  classical vs circuit max |Δ| = {maxdev:.1e}",
        fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = "qlb_port/validation_3d.png"
    fig.savefig(out, dpi=140)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
