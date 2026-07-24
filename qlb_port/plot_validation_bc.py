"""
Reflecting-wall boundary in 3D: an oblique wave packet bouncing in a box.
========================================================================

Boundary conditions are implemented in the 3D framework with a per-axis choice of
'periodic' or 'reflecting' (bounce-back hard walls).  Here we stress-test the
reflecting walls on *all three* axes at once, and against the quantum circuit.

A massless Dirac wave packet is launched from a low corner along an **oblique**
direction (v_x > v_y > v_z, none equal), so it approaches the +x, +y and +z walls
at three different angles and three different times.  Because the walls are hit in
sequence, one run exercises every reflecting plane (yz, xz, xy walls) independently
and shows the bounce-back BC works for each angle of incidence.

The oblique mover is built as an eigenstate u(k) of the one-step massless Fourier
operator M(k) = Mz(kz) My(ky) Mx(kx); the (+,+,+) branch (all Dirac velocities
<alpha_a> > 0) has measured group velocity ~(0.50, 0.36, 0.26).

The figure shows, for several time steps (columns), the density projected onto the
xy, xz and yz planes (rows): classical QLB as filled color, the quantum circuit as
overlaid contours -- they coincide to ~1e-14.  The bottom row traces the packet's
center-of-mass trajectory (the oblique bounce) in each plane, classical line vs
circuit markers.

Run (on a GPU compute node):
    srun ... python -m qlb_port.plot_validation_bc
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import PowerNorm
from matplotlib.collections import LineCollection

from . import operators as ops
from . import backend as bk
from . import threed

NX = NY = NZ = 5                    # 32 sites per axis  (2 + 5+5+5 = 17 qubits)
Nx = Ny = Nz = 2 ** NX
T = 60
TIMES = [0, 12, 24, 36, 48, 60]     # snapshot columns
K = (1.00, 0.68, 0.15)              # carrier momentum -> oblique v, 3 distinct wall angles
R0 = (7.0, 7.0, 7.0)                # launch from the low corner
SIGMA = 3.2
BC = ("reflecting", "reflecting", "reflecting")


def select_oblique_mover(k):
    """Eigenstate of the one-step Fourier operator M(k) that moves into the (+,+,+)
    octant: the branch whose Dirac velocities <alpha_x>,<alpha_y>,<alpha_z> are all
    positive.  Returns the 4-spinor u and its (approximate) velocity vector."""
    M = np.eye(4, dtype=complex)
    for axis, ka in zip(("x", "y", "z"), k):
        R = ops.ROTATIONS[axis]
        S = np.diag(np.exp(-1j * ka * ops.streaming_signs(axis)))
        M = (R @ S @ R.conj().T) @ M
    _, V = np.linalg.eig(M)
    A = {"x": ops.ALPHA_X, "y": ops.BETA, "z": ops.ALPHA_Z}
    best, best_u, best_v = -np.inf, None, None
    for j in range(4):
        u = V[:, j] / np.linalg.norm(V[:, j])
        v = np.array([np.real(u.conj() @ A[a] @ u) for a in ("x", "y", "z")])
        if v.min() > 0 and v[0] > best:      # all-positive octant, fastest in x
            best, best_u, best_v = v[0], u, v
    return best_u, best_v


def initial_packet(u):
    """3D Gaussian envelope * plane-wave carrier * spinor u, flat (i=((z*Ny+y)*Nx+x)*4+c)."""
    kx, ky, kz = K
    ex = np.exp(-((np.arange(Nx) - R0[0]) ** 2) / (2 * SIGMA ** 2)) * np.exp(1j * kx * np.arange(Nx))
    ey = np.exp(-((np.arange(Ny) - R0[1]) ** 2) / (2 * SIGMA ** 2)) * np.exp(1j * ky * np.arange(Ny))
    ez = np.exp(-((np.arange(Nz) - R0[2]) ** 2) / (2 * SIGMA ** 2)) * np.exp(1j * kz * np.arange(Nz))
    psi = np.zeros((Nz, Ny, Nx, 4), dtype=complex)
    for c in range(4):
        psi[..., c] = ez[:, None, None] * ey[None, :, None] * ex[None, None, :] * u[c]
    return (psi / np.linalg.norm(psi)).reshape(-1)


def density(flat):
    return (np.abs(flat.reshape(Nz, Ny, Nx, 4)) ** 2).sum(3)   # (Nz, Ny, Nx)


def projections(d):
    """Return xy (sum over z), xz (sum over y), yz (sum over x) planes."""
    return d.sum(0), d.sum(1), d.sum(2)          # (Ny,Nx), (Nz,Nx), (Nz,Ny)


def com(d):
    tot = d.sum()
    cx = (d.sum((0, 1)) * np.arange(Nx)).sum() / tot
    cy = (d.sum((0, 2)) * np.arange(Ny)).sum() / tot
    cz = (d.sum((1, 2)) * np.arange(Nz)).sum() / tot
    return np.array([cx, cy, cz])


def main():
    print("Device:", bk.device_report())
    u, v = select_oblique_mover(K)
    print(f"oblique mover: v=({v[0]:+.3f},{v[1]:+.3f},{v[2]:+.3f})  |v|={np.linalg.norm(v):.3f}")
    psi0 = initial_packet(u)

    # --- classical: COM every step (smooth trajectory) + density at snapshot times ---
    dens_c = {0: density(psi0)}
    traj_c = [com(dens_c[0])]
    p = psi0.copy()
    for t in range(1, T + 1):
        p = threed.classical_step_3d(p, NX, NY, NZ, 0.0, bc=BC)
        d = density(p)
        traj_c.append(com(d))
        if t in TIMES:
            dens_c[t] = d
    traj_c = np.array(traj_c)

    # --- circuit: density + COM at snapshot times (one transpile, one run) ---
    step = threed.sweep3d_circuit(NX, NY, NZ, bc=BC)
    svs = bk.evolve_snapshots(step, psi0, TIMES)
    dens_q = {t: density(svs[t]) for t in TIMES}
    traj_q = np.array([com(dens_q[t]) for t in TIMES])

    v_meas = (traj_c[4] - traj_c[0]) / 4.0            # true lattice group velocity
    print(f"measured group velocity v=({v_meas[0]:+.2f},{v_meas[1]:+.2f},{v_meas[2]:+.2f}) "
          f"|v|={np.linalg.norm(v_meas):.2f}")

    maxdev = max(np.abs(dens_c[t] - dens_q[t]).max() for t in TIMES)
    print(f"snapshot density  max |Δ| = {maxdev:.2e}")

    # -------------------------------------------------------------------------
    # figure: rows = {xy, xz, yz} planes, cols = snapshot times; bottom = COM path
    # -------------------------------------------------------------------------
    nt = len(TIMES)
    fig = plt.figure(figsize=(2.15 * nt + 0.6, 10.6))
    gs = GridSpec(4, nt, figure=fig, height_ratios=[1, 1, 1, 1.6], hspace=0.32, wspace=0.12)

    plane_names = ["xy  (sum over z)", "xz  (sum over y)", "yz  (sum over x)"]
    plane_axes = [("x", "y", Nx, Ny), ("x", "z", Nx, Nz), ("y", "z", Ny, Nz)]
    # per-plane vmax for a common color scale across time
    vmax = []
    for pi in range(3):
        vmax.append(max(projections(dens_c[t])[pi].max() for t in TIMES))

    for row in range(3):                         # plane
        hlab, vlab, Nh, Nv = plane_axes[row]
        for col, t in enumerate(TIMES):          # time
            ax = fig.add_subplot(gs[row, col])
            pc = projections(dens_c[t])[row]
            pq = projections(dens_q[t])[row]
            ax.imshow(pc, origin="lower", aspect="equal", cmap="inferno",
                      norm=PowerNorm(0.5, vmin=0, vmax=vmax[row]),
                      extent=[-0.5, Nh - 0.5, -0.5, Nv - 0.5])
            lv = np.array([0.25, 0.55, 0.85]) * vmax[row]
            ax.contour(np.arange(Nh), np.arange(Nv), pq, levels=lv,
                       colors="cyan", linewidths=0.7, alpha=0.9)
            for w in (-0.5, Nh - 0.5):
                ax.axvline(w, color="deepskyblue", lw=1.5)
            for w in (-0.5, Nv - 0.5):
                ax.axhline(w, color="deepskyblue", lw=1.5)
            ax.set_xticks([]); ax.set_yticks([])
            if col == 0:
                ax.set_ylabel(f"{plane_names[row]}\n{vlab}", fontsize=9)
            if row == 0:
                ax.set_title(f"t = {t}", fontsize=10)
            if row == 2:
                ax.set_xlabel(hlab, fontsize=8)

    # bottom row: COM trajectory in each plane, colored by time (shows the bounce)
    tt = np.arange(T + 1)
    spans = [slice(0, 2), slice(2, 4), slice(4, 6)]
    idx = [(0, 1, "x", "y", Nx, Ny), (0, 2, "x", "z", Nx, Nz), (1, 2, "y", "z", Ny, Nz)]
    traj_axes, lc_last = [], None
    for (a, b, hl, vl, Nh, Nv), sp in zip(idx, spans):
        ax = fig.add_subplot(gs[3, sp])
        pts = np.array([traj_c[:, a], traj_c[:, b]]).T.reshape(-1, 1, 2)
        segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
        lc = LineCollection(segs, cmap="viridis", norm=plt.Normalize(0, T))
        lc.set_array(tt[:-1])
        lc.set_linewidth(2.5)
        ax.add_collection(lc)
        lc_last = lc
        ax.plot(traj_q[:, a], traj_q[:, b], "o", mfc="none", mec="crimson", mew=1.6,
                ms=8, label="circuit COM")
        ax.plot([traj_c[0, a]], [traj_c[0, b]], "*", color="lime", ms=15, mec="k",
                label="launch")
        ax.add_patch(plt.Rectangle((0, 0), Nh - 1, Nv - 1, fill=False,
                                   ec="deepskyblue", lw=1.5))
        ax.set_xlim(-1, Nh)
        ax.set_ylim(-1, Nv)
        ax.set_xlabel(hl)
        ax.set_ylabel(vl)
        ax.set_aspect("equal")
        ax.set_title(f"COM path ({hl}{vl}): oblique bounce", fontsize=9)
        if sp == spans[0]:
            ax.legend(fontsize=7, loc="upper left", framealpha=0.9)
        traj_axes.append(ax)
    cb = fig.colorbar(lc_last, ax=traj_axes, orientation="horizontal",
                      fraction=0.06, pad=0.26, aspect=60)
    cb.set_label("time step t  (COM path color: launch -> return)", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    fig.suptitle(
        "Reflecting-wall box in 3D: oblique packet bounces off all three wall pairs "
        "(classical QLB vs quantum circuit, GPU)\n"
        f"group velocity v = ({v_meas[0]:.2f}, {v_meas[1]:.2f}, {v_meas[2]:.2f}), "
        f"|v| = {np.linalg.norm(v_meas):.2f} (all components distinct -> 3 angles of incidence);  "
        f"filled = classical (√-scaled), cyan contours = circuit;  max |Δρ| = {maxdev:.1e}",
        fontsize=11)
    out = "qlb_port/validation_bc.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
