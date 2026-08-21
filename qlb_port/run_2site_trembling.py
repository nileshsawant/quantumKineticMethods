"""
2-site (n_pos=1) Zitterbewegung trembling on hardware: <alpha_x>(t), t=0..tmax.
==============================================================================

On the two-site lattice the only carrier is k0=0, at which the velocity
<alpha_x> starts at zero and oscillates with period 2*pi/2E ~= 4.1 steps,
swinging almost the full [-1,+1]. Crucially, each time point is an *independent*
shallow circuit (prepare -> t steps -> rotate into the alpha_x basis -> measure
one spinor bit), only one or two two-qubit gates after transpilation, so unlike
the four-site demo every step survives on hardware. Scanning one full period
therefore traces the trembling itself, not just its first sign flip.

Each time t is submitted as its own circuit (measurement ends the run, so a time
series cannot be read from one circuit). Raw device counts go to a SEPARATE file
(default hw_2site_trembling_raw.txt) so earlier runs are not overwritten.

Usage (after `module load qiskit/aer-gpu` and exporting the credentials):
    # emulator only (no cloud, no credits):
    PYTHONPATH=. python3 -m qlb_port.run_2site_trembling
    # on hardware (tmax+1 jobs):
    PYTHONPATH=. python3 -m qlb_port.run_2site_trembling --submit --shots 2000
"""

import argparse
import datetime

import numpy as np
from qiskit import transpile

from . import run_hardware_zitterbewegung as H


def _cz(qc):
    t = transpile(qc, basis_gates=["cz", "rz", "sx", "x"], optimization_level=3)
    return sum(v for g, v in t.count_ops().items() if g == "cz")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--npos", type=int, default=1, help="position qubits (1 => 2 sites).")
    ap.add_argument("--mass", type=float, default=0.8)
    ap.add_argument("--k0", type=float, default=np.pi / 2,
                    help="target carrier (snapped to nearest lattice momentum).")
    ap.add_argument("--tmax", type=int, default=4, help="scan t=0..tmax (one period ~4).")
    ap.add_argument("--shots", type=int, default=2000)
    ap.add_argument("--submit", action="store_true",
                    help="ACTUALLY submit to the cloud backend (spends credits).")
    ap.add_argument("--backend-name", default=H.DEFAULT_BACKEND)
    ap.add_argument("--raw", default="qlb_port/hw_2site_trembling_raw.txt",
                    help="raw-counts file (kept separate from earlier runs).")
    ap.add_argument("--out", default="qlb_port/hw_2site_trembling.png")
    ap.add_argument("--data", default="qlb_port/hw_2site_trembling.npz")
    args = ap.parse_args()

    npos, m, shots, tmax = args.npos, args.mass, args.shots, args.tmax
    N = 2 ** npos
    k0 = H.nearest_lattice_k(npos, args.k0)
    psi_m, E = H.mode_state(npos, k0, m)                 # +/-E mode mix for the velocity
    ts = list(range(tmax + 1))

    # exact velocity trajectory, and a smooth analytic curve from a 3-parameter
    # fit {1, cos 2Et, sin 2Et} over a few extra steps (the trembling is a single
    # tone at 2E, so the fit is exact up to round-off).
    tf = np.arange(max(tmax, 8) + 1)
    traj = H.classical_evolution(psi_m, npos, m, int(tf[-1]))
    avf = np.array([H.alpha_x_expectation(traj[t], npos) for t in tf])
    av_exact = avf[:tmax + 1]
    Amat = np.column_stack([np.ones_like(tf, float), np.cos(2 * E * tf), np.sin(2 * E * tf)])
    coef, *_ = np.linalg.lstsq(Amat, avf, rcond=None)
    t_smooth = np.linspace(0, tmax, 400)
    av_smooth = coef[0] + coef[1] * np.cos(2 * E * t_smooth) + coef[2] * np.sin(2 * E * t_smooth)

    # one velocity circuit per time; emulator (shot-based, like hardware)
    circs = {t: H.alpha_x_circuit(npos, psi_m, m, t) for t in ts}
    av_aer = np.array([H.alpha_x_from_counts(H.run_aer(circs[t], shots), shots) for t in ts])

    period = 2 * np.pi / (2 * E)
    print(f"2-site trembling: npos={npos} (N={N} sites, {2 + npos} qubits)  "
          f"k0={k0:.4f}  m~={m}  E={E:.4f}  period={period:.2f} steps  shots={shots}")
    print("   t   exact<ax>   emul<ax>   cz")
    for i, t in enumerate(ts):
        print(f"  {t:2d}   {av_exact[i]:+8.4f}   {av_aer[i]:+8.4f}   {_cz(circs[t]):2d}")

    av_hw = None
    if args.submit:
        service = H.get_service()
        av_hw = np.full(len(ts), np.nan)
        try:
            with open(args.raw, "w") as fh:
                fh.write(f"# 2-site trembling raw hardware counts  backend={args.backend_name}  "
                         f"{datetime.datetime.now().isoformat(timespec='seconds')}\n")
                fh.write(f"# npos={npos} N={N} mass={m} k0={k0:.6f} E={E:.6f} "
                         f"shots={shots} tmax={tmax}\n")
                fh.write("# <alpha_x> = P(q1=1) - P(q1=0) after rotating into the alpha_x basis\n")
                print(f"\nSubmitting {len(ts)} velocity circuits to {args.backend_name} "
                      f"({shots} shots each) ...")
                for i, t in enumerate(ts):
                    print(f"  velocity t={t} circuit:")
                    cv = H.run_hardware(service, args.backend_name, circs[t], shots)
                    H._dump_counts(fh, f"alpha_x_t={t}_measure_q1", cv, 1, shots)
                    av_hw[i] = H.alpha_x_from_counts(cv, shots)
                    fh.write(f"  <alpha_x>_hw(t={t}) = {av_hw[i]:+.4f}\n")
                    fh.flush()
            print(f"  raw hardware counts written to {args.raw}")
            print("   t   exact<ax>   hw<ax>")
            for i, t in enumerate(ts):
                print(f"  {t:2d}   {av_exact[i]:+8.4f}   {av_hw[i]:+8.4f}")
        finally:
            service.close()

    save = dict(t=np.array(ts), av_exact=av_exact, av_aer=av_aer, E=float(E),
                k0=k0, npos=npos, mass=m, shots=shots)
    if av_hw is not None:
        save["av_hw"] = av_hw
    np.savez(args.data, **save)

    # figure: the trembling curve, same colour code as the other hardware figures
    # (exact grey, emulator blue C0, hardware red C3)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 4.6))
    ax.plot(t_smooth, av_smooth, color="0.6", lw=1.6, zorder=1,
            label=r"exact $\langle\alpha_x\rangle(t)$")
    ax.plot(ts, av_exact, "o", color="0.4", ms=5, zorder=2, label="exact (steps)")
    ax.plot(ts, av_aer, "s", color="C0", ms=8, mfc="none", mew=1.7, zorder=3,
            label=f"emulator ({shots} shots)")
    if av_hw is not None:
        ax.plot(ts, av_hw, "D", color="C3", ms=8, zorder=4, label="hardware")
    ax.axhline(0, color="k", lw=0.6, ls=":")
    ax.set_xlabel("step  $t$")
    ax.set_ylabel(r"$\langle\alpha_x\rangle(t)$")
    ax.set_ylim(-1.18, 1.18)
    ax.set_xticks(ts)
    ax.set_title(rf"Two-site Zitterbewegung on hardware "
                 rf"($N=2$, $\tilde m={m}$, period $\approx{period:.1f}$ steps)")
    ax.legend(frameon=False, ncol=2, loc="lower right")
    fig.tight_layout()
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"\nwrote {args.out} and {args.data}")


if __name__ == "__main__":
    main()
