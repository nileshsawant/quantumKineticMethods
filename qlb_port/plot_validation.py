"""
Overlay the classical QLB algorithm against the ported quantum circuit.
=======================================================================

Evolves a 1D Dirac wave packet two ways and plots the probability density on top of
each other so agreement is visible by eye:

  * classical: the QLB sub-step loop (rotate -> collide -> stream -> rotate) on a
    (N, 4) spinor field, and
  * quantum:   the ported sweep circuit run on the qiskit-aer (GPU) statevector.

Top row: massless free particle (rigid +x propagation).
Bottom row: massive free particle (Zitterbewegung / dispersion — exercises the
collision gate).  Classical = solid line, quantum circuit = open circles.

Run:
    module load qiskit/aer-gpu
    python -m qlb_port.plot_validation
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from . import operators as ops
from . import backend as bk
from . import sweep

AXIS = "x"
N_POS = 6                 # lattice N = 64
N = 2 ** N_POS


def classical_step(psi, m_tilde, g_tilde=0.0):
    """One classical QLB sub-step on a (N, 4) spinor field (periodic)."""
    R = ops.ROTATIONS[AXIS]
    R_inv = R.conj().T
    Q_char = ops.collision_operator_char(AXIS, m_tilde, g_tilde)
    signs = ops.streaming_signs(AXIS)
    pr = psi @ R_inv.T                      # rotate into char frame
    pr = pr @ Q_char.T                      # collide
    out = np.empty_like(pr)                 # stream +/-1 per component (periodic)
    for c in range(4):
        out[:, c] = np.roll(pr[:, c], int(signs[c]))
    return out @ R.T                        # rotate back


def initial_packet(x0=16, sigma=4.0, k0=0.5):
    """Physical +x-moving Gaussian Dirac spinor field, shape (N, 4), normalized."""
    sp = ops.X_ROTATION @ (np.array([0, 0, 1, 1], dtype=complex) / np.sqrt(2))
    x = np.arange(N)
    env = np.exp(-((x - x0) ** 2) / (2 * sigma ** 2)) * np.exp(1j * k0 * x)
    psi = env[:, None] * sp[None, :]
    psi /= np.linalg.norm(psi)
    return psi


def field_to_statevector(psi):
    """(N,4) field -> statevector with index i = x*4 + c."""
    return psi.reshape(-1)


def statevector_to_density(sv):
    """|psi|^2 summed over the 4 spinor components, per position."""
    return (np.abs(sv.reshape(N, 4)) ** 2).sum(axis=1)


def run(m_tilde, snapshots):
    psi0 = initial_packet()
    sv0 = field_to_statevector(psi0)
    x = np.arange(N)

    classical, quantum = {}, {}
    # classical loop, capturing snapshots
    psi = psi0.copy()
    tprev = 0
    for t in snapshots:
        for _ in range(t - tprev):
            psi = classical_step(psi, m_tilde)
        tprev = t
        classical[t] = (np.abs(psi) ** 2).sum(axis=1)
    # quantum circuit, one run per snapshot
    for t in snapshots:
        if t == 0:
            sv = sv0
        else:
            sv = bk.apply_to_statevector(sweep.evolution_circuit(AXIS, N_POS, t, m_tilde=m_tilde), sv0)
        quantum[t] = statevector_to_density(sv)
    # quantitative agreement
    max_dev = max(np.max(np.abs(classical[t] - quantum[t])) for t in snapshots)
    return x, classical, quantum, max_dev


def main():
    print("Device:", bk.device_report())
    snapshots = [0, 12, 24, 36]
    cases = [("Massless free particle (m=0)", 0.0),
             ("Massive free particle (m=0.35): Zitterbewegung", 0.35)]

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    colors = plt.cm.viridis(np.linspace(0, 0.85, len(snapshots)))

    for ax, (title, m) in zip(axes, cases):
        x, classical, quantum, max_dev = run(m, snapshots)
        for t, col in zip(snapshots, colors):
            ax.plot(x, classical[t], "-", color=col, lw=2,
                    label=f"classical  t={t}")
            ax.plot(x, quantum[t], "o", color=col, mfc="none", ms=5, mew=1.3,
                    label=f"circuit    t={t}")
        ax.set_title(f"{title}    (max |Δ| = {max_dev:.2e})", fontsize=11)
        ax.set_ylabel(r"$|\psi(x)|^2$")
        ax.grid(alpha=0.3)
        print(f"{title}: max |classical - circuit| = {max_dev:.3e}")

    axes[-1].set_xlabel("lattice site  x")
    # compact legend (classical vs circuit only, using the first two handles per axis)
    handles, labels = axes[0].get_legend_handles_labels()
    axes[0].legend(handles, labels, ncol=2, fontsize=8, loc="upper left")
    fig.suptitle("Classical QLB algorithm vs ported quantum circuit  (1D Dirac, GPU)",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = "qlb_port/validation_overlay.png"
    fig.savefig(out, dpi=150)
    print(f"\nSaved {out}")


if __name__ == "__main__":
    main()
