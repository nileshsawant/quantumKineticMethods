"""
2-site (n_pos=1) position sloshing on hardware: rho(x,t) and <x>(t), t=0..tmax.
==============================================================================

Companion to run_2site_trembling.py. Where that script follows the *velocity*
<alpha_x>(t) of the mode-mix state, this one follows the *position* of a
localized packet: on two sites the streaming increment swaps the two sites each
step, so a packet sloshes back and forth with period two, and its mean position
<x> = rho(1) oscillates between its two extremes -- the position side of the
Zitterbewegung. Each time t is a separate density circuit (prepare -> t sweeps
-> measure all qubits), shallow enough to survive on hardware.

Raw device counts go to a SEPARATE file (default hw_2site_density_raw.txt) so no
earlier run is overwritten.

Usage (after `module load qiskit/aer-gpu` and exporting the credentials):
    # emulator only (no cloud, no credits):
    PYTHONPATH=. python3 -m qlb_port.run_2site_density
    # on hardware (tmax+1 jobs):
    PYTHONPATH=. python3 -m qlb_port.run_2site_density --submit --shots 2000
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
    ap.add_argument("--sigma", type=float, default=1.0)
    ap.add_argument("--tmax", type=int, default=4, help="scan t=0..tmax.")
    ap.add_argument("--shots", type=int, default=2000)
    ap.add_argument("--submit", action="store_true",
                    help="ACTUALLY submit to the cloud backend (spends credits).")
    ap.add_argument("--backend-name", default=H.DEFAULT_BACKEND)
    ap.add_argument("--raw", default="qlb_port/hw_2site_density_raw.txt",
                    help="raw-counts file (kept separate from earlier runs).")
    ap.add_argument("--out", default="qlb_port/hw_2site_density.png")
    ap.add_argument("--data", default="qlb_port/hw_2site_density.npz")
    args = ap.parse_args()

    npos, m, shots, tmax = args.npos, args.mass, args.shots, args.tmax
    N = 2 ** npos
    xgrid = np.arange(N)
    k0 = H.nearest_lattice_k(npos, args.k0)
    psi_d, _ = H.packet_state(npos, k0, m, args.sigma)   # localized packet
    ts = list(range(tmax + 1))

    # exact density trajectory and mean position <x> = sum_x x rho(x)
    traj = H.classical_evolution(psi_d, npos, m, tmax)
    rho_exact = np.array([H.density(traj[t], npos) for t in ts])
    x_exact = rho_exact @ xgrid

    # one density circuit per time; emulator (shot-based, like hardware)
    circs = {t: H.density_circuit(npos, psi_d, m, t) for t in ts}
    rho_aer = np.array([H.density_from_counts(H.run_aer(circs[t], shots), npos, shots) for t in ts])
    x_aer = rho_aer @ xgrid

    print(f"2-site density sloshing: npos={npos} (N={N} sites, {2 + npos} qubits)  "
          f"k0={k0:.4f}  m~={m}  shots={shots}")
    print("   t   rho_exact            <x>ex   <x>em   cz")
    for i, t in enumerate(ts):
        print(f"  {t:2d}   {np.round(rho_exact[i], 3).tolist()!s:20s} "
              f"{x_exact[i]:+.3f}  {x_aer[i]:+.3f}   {_cz(circs[t]):2d}")

    rho_hw = x_hw = None
    if args.submit:
        service = H.get_service()
        rho_hw = np.full((len(ts), N), np.nan)
        try:
            with open(args.raw, "w") as fh:
                fh.write(f"# 2-site density raw hardware counts  backend={args.backend_name}  "
                         f"{datetime.datetime.now().isoformat(timespec='seconds')}\n")
                fh.write(f"# npos={npos} N={N} mass={m} k0={k0:.6f} sigma={args.sigma} "
                         f"shots={shots} tmax={tmax}\n")
                fh.write("# density layout: site x = high bits, spinor c = low 2 bits\n")
                print(f"\nSubmitting {len(ts)} density circuits to {args.backend_name} "
                      f"({shots} shots each) ...")
                for i, t in enumerate(ts):
                    print(f"  density t={t} circuit:")
                    cd = H.run_hardware(service, args.backend_name, circs[t], shots)
                    H._dump_counts(fh, f"density_t={t}", cd, 2 + npos, shots)
                    rho_hw[i] = H.density_from_counts(cd, npos, shots)
                    fh.write("  rho_hw(x) = " + " ".join(f"{v:.4f}" for v in rho_hw[i]) + "\n")
                    fh.flush()
            x_hw = rho_hw @ xgrid
            print(f"  raw hardware counts written to {args.raw}")
            print("   t   rho_hw               <x>hw")
            for i, t in enumerate(ts):
                print(f"  {t:2d}   {np.round(rho_hw[i], 3).tolist()!s:20s} {x_hw[i]:+.3f}")
        finally:
            service.close()

    save = dict(t=np.array(ts), rho_exact=rho_exact, rho_aer=rho_aer,
                x_exact=x_exact, x_aer=x_aer, k0=k0, npos=npos, mass=m, shots=shots)
    if rho_hw is not None:
        save["rho_hw"] = rho_hw
        save["x_hw"] = x_hw
    np.savez(args.data, **save)

    # standalone figure: mean position <x>(t) sloshing (exact grey / emu C0 / hw C3)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 4.6))
    ax.plot(ts, x_exact, "-o", color="0.5", lw=1.3, ms=5, zorder=2, label=r"exact  $\langle x\rangle$")
    ax.plot(ts, x_aer, "s", color="C0", ms=8, mfc="none", mew=1.7, zorder=3,
            label=f"emulator ({shots} shots)")
    if x_hw is not None:
        ax.plot(ts, x_hw, "D", color="C3", ms=8, zorder=4, label="hardware")
    ax.axhline(0.5, color="k", lw=0.6, ls=":")
    ax.set_xlabel("step  $t$")
    ax.set_ylabel(r"mean position  $\langle x\rangle=\rho(1)$")
    ax.set_ylim(-0.02, 1.02)
    ax.set_xticks(ts)
    ax.set_title(rf"Two-site position sloshing ($N=2$, $\tilde m={m}$, period $2$ steps)")
    ax.legend(frameon=False, ncol=2, loc="upper right")
    fig.tight_layout()
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"\nwrote {args.out} and {args.data}")


if __name__ == "__main__":
    main()
