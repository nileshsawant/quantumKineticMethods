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

    # classical at t=0 and t=T
    pc = psi0.copy()
    for _ in range(T):
        pc = threed.classical_step_3d(pc, NX, NY, NZ, 0.0)
    d_c0, d_cT = dens(psi0), dens(pc)
    # circuit at t=0 and t=T (one transpile, one run)
    svs = bk.evolve_snapshots(threed.sweep3d_circuit(NX, NY, NZ), psi0, [0, T])
    d_q0, d_qT = dens(svs[0]), dens(svs[T])

    maxdev = max(np.max(np.abs(d_c0 - d_q0)), np.max(np.abs(d_cT - d_qT)))
    c0, cT = com(d_q0), com(d_qT)
    print(f"COM start {tuple(f'{v:.2f}' for v in c0)} -> end {tuple(f'{v:.2f}' for v in cT)}"
          f"   max|Δ|={maxdev:.2e}")

    # three orthogonal projections at start (row 0) and end (row 1), circuit density
    def proj(d):
        return [d.sum(0), d.sum(1), d.sum(2)]   # xy (Ny,Nx), xz (Nz,Nx), yz (Nz,Ny)
    labels = [("x", "y"), ("x", "z"), ("y", "z")]
    rows = [("t = 0", d_q0, c0), (f"t = {T}", d_qT, cT)]
    vmax = max(proj(d_q0)[0].max(), proj(d_qT)[0].max())

    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    for r, (tlab, d, c) in enumerate(rows):
        pr = proj(d)
        cc = [(c[0], c[1]), (c[0], c[2]), (c[1], c[2])]   # COM in each projection
        for j in range(3):
            ax = axes[r, j]
            im = ax.imshow(pr[j], origin="lower", aspect="auto", cmap="inferno",
                           vmin=0, vmax=vmax,
                           extent=[0, pr[j].shape[1], 0, pr[j].shape[0]])
            ax.plot(cc[j][0], cc[j][1], "c+", ms=12, mew=2)   # mark center of mass
            ax.set_xlabel(labels[j][0]); ax.set_ylabel(labels[j][1])
            ax.set_title(f"{tlab}   {labels[j][0]}{labels[j][1]}-projection", fontsize=10)
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    fig.suptitle(
        "3D Dirac packet sent along the (1,1,1) diagonal: classical QLB vs quantum circuit (GPU)\n"
        f"COM {tuple(round(v, 1) for v in c0)} -> {tuple(round(v, 1) for v in cT)} "
        f"(equal motion in x, y, z);  classical vs circuit max |Δ| = {maxdev:.1e}",
        fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = "qlb_port/validation_3d.png"
    fig.savefig(out, dpi=140)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
