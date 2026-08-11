"""
Replicate the Gerritsma et al. (Nature 463, 68 (2010)) trapped-ion experiment.
==============================================================================

That experiment analog-simulated the 1+1D Dirac equation with a single trapped ion and
observed **Zitterbewegung** (the trembling motion of a relativistic wave packet) together
with the relativistic-to-non-relativistic crossover.  Here we reproduce the same physics
digitally, with the QLB Dirac scheme, running **both**

  * the classical QLB solver (the exact one-step operator ``sweep.sweep_operator``), and
  * the ported quantum circuit (``sweep.sweep_circuit`` on the qiskit-aer statevector
    emulator),

from an identical initial state, and confirm that the two agree to machine precision.

Physics (1+1D, lattice units, hbar = c = dx = dt = 1).  The one-step Bloch symbol along x
is ``U(k) = R . diag(e^{-i k s_c}) . Q_char(m) . R^{-1}`` with per-component streaming
signs ``s_c``; its eigenphases are ``+/- E(k)`` with the Dirac dispersion

    E(k) = sqrt( sin^2 k + m~^2 ) .

A packet that mixes the +E and -E eigenmodes at a carrier k0 shows Zitterbewegung with

    omega_ZB = 2 E(k0)        (trembling frequency)
    R_ZB     = m~ / (2 E^2)   (trembling amplitude in position; = A_v / omega_ZB,
                               with velocity amplitude A_v = |<+|alpha_x|->| = m~/E)

so a massless packet ([alpha_x, H] = 0) does not tremble, and increasing the mass at fixed
momentum raises omega_ZB and moves R_ZB through its maximum -- the crossover.

Run:
    python -m qlb_port.plot_validation_zitterbewegung
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit, transpile

from . import operators as ops
from . import backend as bk
from . import sweep

# ---- fixed experiment parameters ------------------------------------------------
N_POS = 8                       # lattice N = 2**8 = 256 sites (2 + 8 = 10 qubits)
N = 2 ** N_POS
K0 = 0.6                        # carrier wavenumber (pc = sin k0)
SIGMA = 14.0                    # Gaussian width (sites); wide => clean single-frequency ZB
X0 = N // 2

_R = ops.ROTATIONS["x"]
_RINV = _R.conj().T
_SIGNS = ops.streaming_signs("x").astype(float)
ALPHA_X = ops.ALPHA_X


def bloch_symbol(k, m):
    """4x4 one-step operator at wavenumber k: U(k) = R . S(k) . Q_char(m) . R^-1."""
    S = np.diag(np.exp(-1j * k * _SIGNS))
    return _R @ S @ ops.collision_operator_char("x", m, 0.0) @ _RINV


def energy_eigenmodes(k0, m):
    """Return (u_plus, u_minus, E): the +E and -E spinor eigenmodes at k0 and E=E(k0)."""
    w, V = np.linalg.eig(bloch_symbol(k0, m))
    ph = np.angle(w)
    return V[:, int(np.argmax(ph))], V[:, int(np.argmin(ph))], 0.5 * (ph.max() - ph.min())


def initial_state(m, mix=(1.0, 1.0)):
    """Gaussian * plane-wave * (mix of +/-E eigenmodes); flat statevector (i = x*4 + c)."""
    up, un, E = energy_eigenmodes(K0, m)
    u0 = mix[0] * up + mix[1] * un
    u0 = u0 / np.linalg.norm(u0)
    x = np.arange(N)
    env = np.exp(-((x - X0) ** 2) / (2 * SIGMA ** 2)) * np.exp(1j * K0 * x)
    psi = (env[:, None] * u0[None, :]).reshape(-1)
    return psi / np.linalg.norm(psi), E


def observables(sv):
    """Return (<x>, <alpha_x>) for a flat statevector sv (layout i = x*4 + c)."""
    field = sv.reshape(N, 4)
    rho = (np.abs(field) ** 2).sum(1)
    x_mean = (np.arange(N) * rho).sum() / rho.sum()
    v_mean = np.real(np.einsum("xi,ij,xj->", field.conj(), ALPHA_X, field))
    return x_mean, v_mean


def run_classical(m, T, mix=(1.0, 1.0)):
    """Evolve with the exact classical one-step operator; return <x>, <alpha_x>, E, states."""
    psi, E = initial_state(m, mix)
    U = sweep.sweep_operator("x", N_POS, m_tilde=m)     # (4N x 4N) exact solver step
    xs, vs, states = [], [], []
    for _ in range(T + 1):
        states.append(psi)
        x_mean, v_mean = observables(psi)
        xs.append(x_mean); vs.append(v_mean)
        psi = U @ psi
    return np.array(xs), np.array(vs), E, states


def run_circuit(m, T, mix=(1.0, 1.0)):
    """Evolve with the ported circuit on the emulator; return <x>, <alpha_x>, E, states.

    One QLB step is transpiled once and reused for all T steps (with a statevector save
    after each), so the whole evolution is a single simulator run and no large circuit is
    re-transpiled.
    """
    psi, E = initial_state(m, mix)
    sim = bk.statevector_simulator()
    tstep = transpile(sweep.sweep_circuit("x", N_POS, m_tilde=m), sim, optimization_level=0)
    n = tstep.num_qubits
    qc = QuantumCircuit(n)
    qc.set_statevector(psi)
    qc.save_statevector(label="t0")
    for t in range(1, T + 1):
        qc.compose(tstep, inplace=True)
        qc.save_statevector(label=f"t{t}")
    data = sim.run(qc).result().data()
    states = [np.asarray(data[f"t{t}"]) for t in range(T + 1)]
    xs = np.array([observables(s)[0] for s in states])
    vs = np.array([observables(s)[1] for s in states])
    return xs, vs, E, states


def zb_frequency(v, E):
    """Dominant angular frequency (rad/step) of the detrended <alpha_x>(t) signal."""
    t = np.arange(len(v), dtype=float)
    vd = (v - np.polyval(np.polyfit(t, v, 1), t)) * np.hanning(len(v))
    n = 8 * len(vd)
    sp = np.abs(np.fft.rfft(vd, n=n))
    f = np.fft.rfftfreq(n, 1.0) * 2 * np.pi
    band = (f > 0.4 * 2 * E) & (f < 1.6 * 2 * E)
    return f[band][int(np.argmax(sp[band]))]


def main():
    print("Device:", bk.device_report())
    print(f"1+1D Dirac Zitterbewegung; lattice N={N} ({2 + N_POS} qubits), k0={K0}, "
          f"pc=sin k0={np.sin(K0):.3f}\n")

    masses = [0.0, 0.1, 0.2, 0.35, 0.6, 1.0, 2.0]
    T = 90
    print(" m~    E_lat  omegaZB_circ  omegaZB_clas   2E    2sqrt(sin^2k+m^2)  max|drho|  fidelity")
    rows = []
    for m in masses:
        xc, vc, E, cstates = run_classical(m, T)
        xq, vq, E2, qstates = run_circuit(m, T)
        wq = zb_frequency(vq, E) if m > 0 else 0.0
        wc = zb_frequency(vc, E) if m > 0 else 0.0
        maxdev = max(
            float(np.max(np.abs((np.abs(cstates[t].reshape(N, 4)) ** 2).sum(1)
                                - (np.abs(qstates[t].reshape(N, 4)) ** 2).sum(1))))
            for t in range(T + 1))
        fid = bk.state_fidelity(cstates[T], qstates[T])
        wth = 2 * np.sqrt(np.sin(K0) ** 2 + m ** 2)
        print(f"{m:5.2f} {E:6.3f}   {wq:8.3f}    {wc:8.3f}   {2*E:6.3f}   {wth:6.3f}       "
              f"{maxdev:.1e}   {fid:.6f}", flush=True)
        rows.append((m, E, wq))

    # ---- figure: classical lines + circuit markers ------------------------------
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 4.6))
    Tf = 40                              # display window (Zitterbewegung fades as the +/-E parts separate)
    for m, col in [(0.0, "0.5"), (0.2, "C0"), (0.35, "C1"), (0.6, "C3")]:
        xc, vc, E, _ = run_classical(m, Tf)
        xq, vq, E2, _ = run_circuit(m, Tf)
        tt = np.arange(len(xc))
        axA.plot(tt, xc - X0, color=col, lw=2,
                 label=("massless (no ZB)" if m == 0 else fr"$\tilde m$={m}"))
        axA.plot(tt[::2], xq[::2] - X0, "o", color=col, mfc="none", ms=5, mew=1.1)
    axA.axhline(0, color="k", lw=0.6, ls=":")
    axA.set_xlabel("time step  t"); axA.set_ylabel(r"$\langle x(t)\rangle - x_0$  (sites)")
    axA.set_title("Zitterbewegung: classical (lines) vs circuit (markers)")
    axA.legend(fontsize=8); axA.grid(alpha=0.3)

    mm = np.array([r[0] for r in rows[1:]]); wq = np.array([r[2] for r in rows[1:]])
    mg = np.linspace(0.05, 2.0, 60)
    Eg = np.array([energy_eigenmodes(K0, mmg)[2] for mmg in mg])   # exact lattice dispersion
    axB.plot(mg, 2 * Eg, "b-", label=r"$\omega_{ZB}=2E$ (lattice Dirac)")
    axB.plot(mm, wq, "bo", ms=7, label=r"$\omega_{ZB}$ (circuit)")
    axB.set_xlabel(r"mass  $\tilde m$"); axB.set_ylabel(r"$\omega_{ZB}$  (rad/step)", color="b")
    axB.tick_params(axis="y", labelcolor="b")
    axB.set_title(r"Relativistic $\to$ non-relativistic crossover")
    ax2 = axB.twinx()
    ax2.plot(mg, mg / (2 * Eg ** 2), "r--", label=r"$R_{ZB}=\tilde m/2E^2$ (Dirac)")
    ax2.set_ylabel(r"$R_{ZB}$  (sites)", color="r"); ax2.tick_params(axis="y", labelcolor="r")
    h1, l1 = axB.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    axB.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper left"); axB.grid(alpha=0.3)

    fig.tight_layout()
    out = "qlb_port/validation_zitterbewegung.png"
    fig.savefig(out, dpi=140)
    print(f"\nSaved {out}")


if __name__ == "__main__":
    main()
