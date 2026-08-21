"""
2-site (n_pos=1) t=0 state check: hardware vs emulator vs exact.
================================================================

Measures the *prepared* initial state at t=0 (no QLB step applied) two ways --
the position density rho(x) and the mean velocity <alpha_x> -- on the Aer
emulator and, with --submit, on the cloud backend, and overlays them in the same
figure style as run_hardware_zitterbewegung (exact grey, emulator blue, hardware
red).

At n_pos=1 the state preparation is only about one two-qubit gate, so t=0 isolates
the preparation-plus-readout floor with the smallest possible circuit.

Raw device counts are written to a SEPARATE file (default hw_2site_t0_raw.txt) so
that the earlier run's hw_demo_raw.txt is not overwritten.

Usage (after `module load qiskit/aer-gpu` and exporting the credentials):
    # emulator only (no cloud, no credits):
    PYTHONPATH=. python3 -m qlb_port.run_2site_t0
    # on hardware:
    PYTHONPATH=. python3 -m qlb_port.run_2site_t0 --submit --shots 2000
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
    ap.add_argument("--shots", type=int, default=2000)
    ap.add_argument("--submit", action="store_true",
                    help="ACTUALLY submit to the cloud backend (spends credits).")
    ap.add_argument("--backend-name", default=H.DEFAULT_BACKEND)
    ap.add_argument("--raw", default="qlb_port/hw_2site_t0_raw.txt",
                    help="raw-counts file (kept separate from earlier runs).")
    ap.add_argument("--out", default="qlb_port/hw_2site_t0.png")
    ap.add_argument("--data", default="qlb_port/hw_2site_t0.npz")
    args = ap.parse_args()

    npos, m, shots = args.npos, args.mass, args.shots
    N = 2 ** npos
    k0 = H.nearest_lattice_k(npos, args.k0)
    psi_d, _ = H.packet_state(npos, k0, m, args.sigma)   # localized state for the density
    psi_m, E = H.mode_state(npos, k0, m)                 # +/-E mode mix for the velocity

    # exact t=0 (no evolution)
    rho_exact = H.density(psi_d, npos)
    av_exact = H.alpha_x_expectation(psi_m, npos)

    # t=0 circuits: prepare, apply 0 steps, measure
    qc_d = H.density_circuit(npos, psi_d, m, 0)
    qc_v = H.alpha_x_circuit(npos, psi_m, m, 0)

    # emulator (Aer), shot-based like hardware
    rho_aer = H.density_from_counts(H.run_aer(qc_d, shots), npos, shots)
    av_aer = H.alpha_x_from_counts(H.run_aer(qc_v, shots), shots)

    print(f"2-site t=0 check: npos={npos} (N={N} sites, {2 + npos} qubits)  "
          f"k0={k0:.4f}  m~={m}  shots={shots}")
    print(f"  t=0 gate cost (cz, all-to-all): density={_cz(qc_d)}  velocity={_cz(qc_v)}")
    print(f"  exact    : rho(x)={np.round(rho_exact, 3).tolist()}  <alpha_x>={av_exact:+.3f}")
    print(f"  emulator : rho(x)={np.round(rho_aer, 3).tolist()}  <alpha_x>={av_aer:+.3f}")

    rho_hw = av_hw = None
    if args.submit:
        service = H.get_service()
        try:
            with open(args.raw, "w") as fh:
                fh.write(f"# 2-site t=0 raw hardware counts  backend={args.backend_name}  "
                         f"{datetime.datetime.now().isoformat(timespec='seconds')}\n")
                fh.write(f"# npos={npos} N={N} mass={m} k0={k0:.6f} shots={shots} "
                         f"(t=0, no QLB step)\n")
                fh.write("# density layout: site x = high bits, spinor c = low 2 bits\n")
                print(f"\nSubmitting 2 circuits to {args.backend_name} ({shots} shots each) ...")
                print("  density t=0 circuit:")
                cd = H.run_hardware(service, args.backend_name, qc_d, shots)
                H._dump_counts(fh, "density_t0", cd, 2 + npos, shots)
                rho_hw = H.density_from_counts(cd, npos, shots)
                fh.write("  rho_hw(x) = " + " ".join(f"{v:.4f}" for v in rho_hw) + "\n")
                fh.flush()
                print("  velocity t=0 circuit:")
                cv = H.run_hardware(service, args.backend_name, qc_v, shots)
                H._dump_counts(fh, "alpha_x_t0_measure_q1", cv, 1, shots)
                av_hw = H.alpha_x_from_counts(cv, shots)
                fh.write(f"  <alpha_x>_hw(t0) = {av_hw:+.4f}\n")
                fh.flush()
            print(f"  raw hardware counts written to {args.raw}")
            print(f"  hardware : rho(x)={np.round(rho_hw, 3).tolist()}  <alpha_x>={av_hw:+.3f}")
            print(f"  density total-variation distance  exact-vs-emulator: "
                  f"{H._tvd(rho_exact, rho_aer):.3f}   exact-vs-hardware: "
                  f"{H._tvd(rho_exact, rho_hw):.3f}")
        finally:
            service.close()

    save = dict(x=np.arange(N), rho_exact=rho_exact, rho_aer=rho_aer,
                av_exact=float(av_exact), av_aer=float(av_aer), shots=shots,
                npos=npos, mass=m, k0=k0)
    if rho_hw is not None:
        save["rho_hw"] = rho_hw
        save["av_hw"] = float(av_hw)
    np.savez(args.data, **save)

    # figure: same style as run_hardware_zitterbewegung (exact grey / emulator blue / hardware red)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.2))
    xs = np.arange(N)
    if rho_hw is None:
        w = 0.38
        axA.bar(xs - w / 2, rho_exact, w, color="0.65", label="exact")
        axA.bar(xs + w / 2, rho_aer, w, color="C0", label=f"emulator ({shots} shots)")
    else:
        w = 0.28
        axA.bar(xs - w, rho_exact, w, color="0.65", label="exact")
        axA.bar(xs, rho_aer, w, color="C0", label=f"emulator ({shots} shots)")
        axA.bar(xs + w, rho_hw, w, color="C3", label="hardware")
    axA.set_xticks(xs)
    axA.set_xlabel("site  $x$")
    axA.set_ylabel(r"$\rho(x)$ at $t=0$")
    axA.set_title(f"$t=0$ density, {N}-site line ({2 + npos} qubits)")
    axA.legend(frameon=False)

    labels, vals, cols = ["exact", "emulator"], [av_exact, av_aer], ["0.65", "C0"]
    if av_hw is not None:
        labels.append("hardware"); vals.append(av_hw); cols.append("C3")
    axB.bar(range(len(vals)), vals, color=cols, width=0.6)
    axB.axhline(0, color="k", lw=0.6, ls=":")
    axB.axhline(av_exact, color="0.65", lw=1.0, ls="--")
    axB.set_xticks(range(len(vals)))
    axB.set_xticklabels(labels)
    axB.set_ylim(-1.05, 1.05)
    axB.set_ylabel(r"$\langle\alpha_x\rangle$ at $t=0$")
    axB.set_title(rf"$t=0$ velocity (exact $={av_exact:+.2f}$)")
    tag = "vs hardware" if rho_hw is not None else "(emulator only)"
    fig.suptitle(f"2-site $t=0$ state: exact vs emulator {tag}", y=1.02)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"\nwrote {args.out} and {args.data}")


if __name__ == "__main__":
    main()
