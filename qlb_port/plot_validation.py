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
from . import potential

AXIS = "x"
N_POS = 6                 # lattice N = 64
N = 2 ** N_POS


def classical_step(psi, m_tilde=0.0, V_tilde=None):
    """One classical QLB sub-step on a (N, 4) spinor field (periodic).

    V_tilde : None for free particle, or a per-site potential coupling array.
    """
    R = ops.ROTATIONS[AXIS]
    R_inv = R.conj().T
    signs = ops.streaming_signs(AXIS)
    pr = psi @ R_inv.T                          # rotate into char frame
    if V_tilde is None:                         # uniform collide
        pr = pr @ ops.collision_operator_char(AXIS, m_tilde, 0.0).T
    elif m_tilde == 0.0:                        # massless: position-dependent phase
        pr = pr * potential.collision_phases(V_tilde)[:, None]
    else:                                       # massive: per-site 2-qubit collision
        for x in range(N):
            pr[x] = ops.collision_operator_char(AXIS, m_tilde, float(V_tilde[x])) @ pr[x]
    out = np.empty_like(pr)                     # stream +/-1 per component (periodic)
    for c in range(4):
        out[:, c] = np.roll(pr[:, c], int(signs[c]))
    return out @ R.T                            # rotate back


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


def run(m_tilde, snapshots, x0=16, sigma=4.0, k0=0.5, V_tilde=None):
    """Evolve the packet classically and via circuit; return spinor fields (N,4) per snapshot.

    Returns (x, classical_fields, circuit_fields) where each *_fields[t] is a (N,4) array.
    """
    psi0 = initial_packet(x0=x0, sigma=sigma, k0=k0)
    sv0 = field_to_statevector(psi0)

    classical = {}
    psi = psi0.copy()
    tprev = 0
    for t in snapshots:
        for _ in range(t - tprev):
            psi = classical_step(psi, m_tilde, V_tilde)
        tprev = t
        classical[t] = psi.copy()

    if V_tilde is None:
        step = sweep.sweep_circuit(AXIS, N_POS, m_tilde=m_tilde)
    else:
        step = potential.sweep_circuit_potential(AXIS, N_POS, V_tilde, m_tilde=m_tilde)
    svs = bk.evolve_snapshots(step, sv0, snapshots)
    quantum = {t: svs[t].reshape(N, 4) for t in snapshots}
    return np.arange(N), classical, quantum


def density(field):
    """|psi(x)|^2 summed over spinor components, for a (N,4) field."""
    return (np.abs(field) ** 2).sum(axis=1)


def phase_shift(field_bar, field_free, thresh=1e-3):
    """
    Amplitude-weighted local phase shift arg(sum_c psi_bar,c conj(psi_free,c)) vs x,
    masked (NaN) where the free density is below `thresh` (phase undefined there).
    """
    z = np.sum(field_bar * np.conj(field_free), axis=1)
    dphi = np.angle(z)
    dphi[density(field_free) < thresh] = np.nan
    return dphi


def main():
    print("Device:", bk.device_report())
    snaps = [0, 18, 36]
    tf = snaps[-1]
    colors = plt.cm.viridis(np.linspace(0, 0.8, len(snaps)))

    # Massless runs: free and barrier (reused for the density and phase panels)
    Vk = potential.impurity_field(N_POS, range(40, 44), g_value=0.9)
    x, cf_free, qf_free = run(0.0, snaps, x0=20, sigma=4.0, k0=0.6, V_tilde=None)
    _, cf_bar,  qf_bar = run(0.0, snaps, x0=20, sigma=4.0, k0=0.6, V_tilde=Vk)
    # Massive runs: free (Zitterbewegung) and barrier (reflection)
    _, cf_mfree, qf_mfree = run(0.35, snaps, x0=20, sigma=4.0, k0=0.6, V_tilde=None)
    Vm = potential.impurity_field(N_POS, range(40, 44), g_value=2.0)
    _, cf_mbar, qf_mbar = run(0.6, snaps, x0=20, sigma=4.0, k0=0.6, V_tilde=Vm)

    fig, axes = plt.subplots(2, 2, figsize=(13, 8.5))
    barspan = (39.5, 43.5)

    def dev(cf, qf):
        return max(np.max(np.abs(density(cf[t]) - density(qf[t]))) for t in snaps)

    # (0,0) Massless: |psi|^2 unchanged by the barrier (perfect transmission, T=1)
    ax = axes[0, 0]
    for t, col in zip(snaps, colors):
        ax.plot(x, density(cf_bar[t]), "-", color=col, lw=2, label=f"barrier  t={t}")
        ax.plot(x, density(qf_bar[t]), "o", color=col, mfc="none", ms=4, mew=1.1)
        ax.plot(x, density(cf_free[t]), ":", color=col, lw=1.2, label=f"free      t={t}")
    ax.axvspan(*barspan, color="red", alpha=0.12, label="barrier")
    ax.set_title(f"Massless: $|\\psi|^2$ identical with/without barrier (Klein, T=1)\n"
                 f"classical vs circuit max |Δ| = {dev(cf_bar, qf_bar):.1e}", fontsize=10)
    ax.set_ylabel(r"$|\psi(x)|^2$")
    ax.legend(fontsize=7, ncol=2, loc="upper left")

    # (0,1) Massless: the phase the barrier imprints (this is what Klein tunneling affects)
    ax = axes[0, 1]
    dphi_c = phase_shift(cf_bar[tf], cf_free[tf], thresh=2e-2)
    dphi_q = phase_shift(qf_bar[tf], qf_free[tf], thresh=2e-2)
    ax.axhline(0.0, color="gray", ls=":", lw=1.2, label="free (no barrier)")
    ax.plot(x, dphi_c, "-", color="C3", lw=2, label="barrier classical")
    ax.plot(x, dphi_q, "o", color="C3", mfc="none", ms=5, mew=1.3, label="barrier circuit")
    ax.axvspan(*barspan, color="red", alpha=0.12)
    ax.ticklabel_format(axis="y", useOffset=False)
    ax.set_ylim(-np.pi, 0.4)
    mean_shift = float(np.nanmean(dphi_c))
    pdev = np.nanmax(np.abs(dphi_c - dphi_q))
    ax.set_title(f"Massless: barrier imprints a phase on the transmitted packet\n"
                 f"mean ΔΦ ≈ {mean_shift:.2f} rad ≈ {np.degrees(mean_shift):.0f}° "
                 f"(|ψ|² unchanged); circuit vs classical max |Δ| = {pdev:.0e}", fontsize=10)
    ax.set_ylabel(r"$\Delta\Phi(x)$  [rad]")
    ax.legend(fontsize=7, loc="lower left")

    # (1,0) Massive free: Zitterbewegung / dispersion
    ax = axes[1, 0]
    for t, col in zip(snaps, colors):
        ax.plot(x, density(cf_mfree[t]), "-", color=col, lw=2, label=f"classical t={t}")
        ax.plot(x, density(qf_mfree[t]), "o", color=col, mfc="none", ms=4, mew=1.1)
    ax.set_title(f"Massive free (m=0.35): Zitterbewegung / dispersion\n"
                 f"classical vs circuit max |Δ| = {dev(cf_mfree, qf_mfree):.1e}", fontsize=10)
    ax.set_ylabel(r"$|\psi(x)|^2$")

    # (1,1) Massive barrier: visible reflection (mass gap breaks Klein protection)
    ax = axes[1, 1]
    for t, col in zip(snaps, colors):
        ax.plot(x, density(cf_mbar[t]), "-", color=col, lw=2, label=f"classical t={t}")
        ax.plot(x, density(qf_mbar[t]), "o", color=col, mfc="none", ms=4, mew=1.1)
    ax.axvspan(*barspan, color="red", alpha=0.12)
    refl = (np.abs(qf_mbar[tf][:38]) ** 2).sum()
    ax.set_title(f"Massive (m=0.6): barrier reflects (Klein broken), R≈{refl:.2f}\n"
                 f"classical vs circuit max |Δ| = {dev(cf_mbar, qf_mbar):.1e}", fontsize=10)
    ax.set_ylabel(r"$|\psi(x)|^2$")

    for ax in axes.flat:
        ax.set_xlabel("lattice site  x")
        ax.grid(alpha=0.3)
    fig.suptitle("Classical QLB algorithm vs ported quantum circuit  (1D Dirac, GPU)",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = "qlb_port/validation_overlay.png"
    fig.savefig(out, dpi=150)
    print(f"Saved {out}   massless dphi max|Δ|={pdev:.2e}   massive R≈{refl:.2f}")


if __name__ == "__main__":
    main()
