"""
Reflecting-wall boundary: a wave packet bouncing in a box (classical vs circuit).
=================================================================================

Boundary conditions are implemented in the 3D framework with a per-axis choice of
'periodic' or 'reflecting' (bounce-back hard walls); a 1D or 2D problem is just a 3D
lattice with thin periodic transverse dimensions (as in classical LB codes).

Here we run a 1D-in-3D box: x has reflecting walls, y and z are thin and periodic,
and the packet is uniform in y, z.  A massless Dirac packet is launched toward the
+x wall; it reflects (bounce-back is unitary) and returns.  We show the space-time
density |psi(x, t)|^2 for the classical algorithm and the quantum circuit side by
side -- the reflecting trajectory (a 'bounce') is identical in both.

Run:
    module load qiskit/aer-gpu
    python -m qlb_port.plot_validation_bc
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from . import operators as ops
from . import backend as bk
from . import threed

NX, NY, NZ = 5, 1, 1           # 32 sites in x, thin (2-site) periodic y, z
Nx = 2 ** NX
T = 50
BC = ("reflecting", "periodic", "periodic")


def initial_packet(x0=6, sigma=1.8, k0=0.6):
    sp = ops.X_ROTATION @ (np.array([0, 0, 1, 1], dtype=complex) / np.sqrt(2))
    xg = np.arange(Nx)
    env = np.exp(-((xg - x0) ** 2) / (2 * sigma ** 2)) * np.exp(1j * k0 * xg)
    psi = np.zeros((2, 2, Nx, 4), dtype=complex)          # (Nz=2, Ny=2, Nx, 4), uniform in y,z
    for c in range(4):
        psi[:, :, :, c] = env[None, None, :] * sp[c]
    return (psi / np.linalg.norm(psi)).reshape(-1)


def xprofile(flat):
    return (np.abs(flat.reshape(2, 2, Nx, 4)) ** 2).sum((0, 1, 3))   # (Nx,)


def main():
    print("Device:", bk.device_report())
    psi0 = initial_packet()

    # classical space-time
    prof_c = [xprofile(psi0)]
    pc = psi0.copy()
    for _ in range(T):
        pc = threed.classical_step_3d(pc, NX, NY, NZ, 0.0, bc=BC)
        prof_c.append(xprofile(pc))
    prof_c = np.array(prof_c)
    # circuit space-time (one transpile, saves at every step)
    svs = bk.evolve_snapshots(threed.sweep3d_circuit(NX, NY, NZ, bc=BC), psi0, list(range(T + 1)))
    prof_q = np.array([xprofile(svs[t]) for t in range(T + 1)])

    diff = np.abs(prof_c - prof_q)
    maxdev = diff.max()
    print(f"space-time max |Δ| = {maxdev:.2e}")

    fig, axes = plt.subplots(1, 3, figsize=(13, 5.5), sharey=True)
    vmax = prof_c.max()
    for ax, dat, title, cmap, vm in [
            (axes[0], prof_c, "classical QLB", "inferno", vmax),
            (axes[1], prof_q, "quantum circuit", "inferno", vmax),
            (axes[2], diff, f"|difference| (max {maxdev:.0e})", "viridis", max(maxdev, 1e-16))]:
        im = ax.imshow(dat, origin="lower", aspect="auto", cmap=cmap, vmin=0, vmax=vm,
                       extent=[0, Nx, 0, T])
        ax.axvline(0.5, color="cyan", lw=2)          # reflecting walls
        ax.axvline(Nx - 0.5, color="cyan", lw=2)
        ax.set_xlabel("x"); ax.set_title(title, fontsize=11)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    axes[0].set_ylabel("time step  t")
    fig.suptitle("Reflecting-wall box (bounce-back): classical QLB vs quantum circuit (GPU)\n"
                 f"packet reflects off the +x wall and returns; probability conserved; "
                 f"max |Δ| = {maxdev:.1e}", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = "qlb_port/validation_bc.png"
    fig.savefig(out, dpi=140)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
