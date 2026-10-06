"""
Parity symmetry verification and run-to-run variation of the four-site (n_pos=2) instance on QCS.

Reads (no cloud, no credits), all in qlb_port/:
    qcs_quil_n2_parity/qcs_results_parity_pinned.json   even-parity states, all qubits read, + 16 calibration jobs
    qcs_quil_n2_parity/qcs_results_repeat_pinned.json   t=3,4 parity circuits twice more, step-1 velocity t=3,4
    qcs_quil_n2_t4check/qcs_results_t4check.json        t=4 velocity: even/odd sectors, step-1 state (all / q1 read)
    qcs_quil_n2_fourier/qcs_results_n2f_pinned.json     step-1 run (readout mitigation), for comparison

Mitigation: the 16-outcome distribution is corrected with the confusion matrix of the same batch
(non-negative least squares), then shots with odd spinor parity q0 xor q1 are discarded.

Also simulates (density matrix, no shots) the even-state velocity with 2% two-qubit depolarizing
and a coherent RZ(phi) on q1 after every two-qubit gate, to show the sensitivity of t=4.

Writes <outdir>/hw_4site_qcs_parity.png and prints the numbers used in the paper.

Usage:
    PYTHONPATH=. python3 -m qlb_port.plot_4site_qcs_parity
"""

import argparse
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import nnls

from .plot_4site_qcs_readout import confusion, tvd

P = "qlb_port/"
S = 2000
IDX = np.arange(16)
C, X = IDX & 3, IDX >> 2
EVEN = (C == 0) | (C == 3)


def dist(job):
    p = np.zeros(2 ** job["ro_bits"])
    for k, v in job["counts"].items():
        p[int(k)] = v / S
    return p


def observables(p, A=None, post=False, sector=EVEN):
    """(<alpha_x>, rho, kept fraction) from a 16-outcome distribution."""
    if A is not None:
        q, _ = nnls(A, p)
        p = q / q.sum()
    kept = p[sector].sum()
    q = np.where(sector, p, 0) / kept if post else p
    return q[(C >> 1) == 1].sum() - q[(C >> 1) == 0].sum(), np.array([q[X == k].sum() for k in range(4)]), kept


def simulate_phase_sensitivity(phis=(-0.15, -0.08, 0.0, 0.08, 0.15)):
    from qiskit import QuantumCircuit, transpile
    from qiskit.circuit.library import RZGate, StatePreparation, UnitaryGate
    from qiskit.quantum_info import DensityMatrix, Operator
    from qiskit_aer.noise import depolarizing_error
    from . import operators as ops, run_hardware_zitterbewegung as H, sweep

    dep = depolarizing_error(0.02, 2).to_quantumchannel()
    psi0 = H.even_sector(H.mode_state(2, 0.0, 0.8)[0], 2)
    out = {}
    for t in range(5):
        qc = QuantumCircuit(4)
        qc.append(StatePreparation(psi0), range(4))
        for _ in range(t):
            qc.compose(sweep.sweep_circuit("x", 2, m_tilde=0.8, streaming_method="fourier"), inplace=True)
        qc.append(UnitaryGate(ops.ROTATIONS["x"].conj().T), [0, 1])
        qc = transpile(qc, basis_gates=["rx", "ry", "rz", "cx"], optimization_level=3, seed_transpiler=0)
        for phi in phis:
            rho = DensityMatrix.from_label("0000")
            for inst in qc.data:
                q = [qc.find_bit(b).index for b in inst.qubits]
                rho = rho.evolve(Operator(inst.operation), q)
                if len(q) == 2:
                    if 1 in q:
                        rho = rho.evolve(Operator(RZGate(phi)), [1])
                    rho = rho.evolve(dep, q)
            out[(t, phi)] = observables(np.real(np.diag(rho.data)))[0]
    return phis, out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default="proposal/quairk/ancillaStreaming/figures")
    ap.add_argument("--no-sim", action="store_true", help="skip the phase-sensitivity simulation")
    args = ap.parse_args()

    R1 = json.load(open(P + "qcs_quil_n2_parity/qcs_results_parity_pinned.json"))
    RR = json.load(open(P + "qcs_quil_n2_parity/qcs_results_repeat_pinned.json"))
    T4 = json.load(open(P + "qcs_quil_n2_t4check/qcs_results_t4check.json"))
    N1 = json.load(open(P + "qcs_quil_n2_fourier/qcs_results_n2f_pinned.json"))
    TW = [{j["name"]: j for j in json.load(open(P + f"qcs_quil_n2_twirl/qcs_results_twirl_{r}.json"))["jobs"]}
          for r in "AB"]
    J1 = {j["name"]: j for j in R1["jobs"]}
    JR = {j["name"]: j for j in RR["jobs"]}
    JT = {j["name"]: j for j in T4["jobs"]}
    JN = {j["name"]: j for j in N1["jobs"]}
    A = confusion(J1, S)
    ts = np.arange(5)
    print(f"parity run: {R1['total_execution_us'] / 1e6:.1f} s, repeat: {RR['total_execution_us'] / 1e6:.1f} s, "
          f"t=4 check: {T4['total_execution_us'] / 1e6:.1f} s QPU execution")

    # ---- run 1: velocity state and packet ----
    V = {k: [] for k in ("exact", "raw", "ro", "post", "ropost", "kept")}
    print("\nvelocity state:  t  exact  raw  readout  post  readout+post  kept")
    for t in ts:
        j = J1[f"palpha_t{t}"]; p = dist(j)
        row = (j["exact_alpha"], observables(p)[0], observables(p, A)[0], observables(p, post=True)[0],
               observables(p, A, post=True)[0], observables(p, post=True)[2])
        for k, v in zip(V, row):
            V[k].append(v)
        print(f"  {t}  " + "  ".join(f"{v:+.3f}" for v in row[:-1]) + f"  {row[-1]:.3f}")
    print("packet velocity:  t  exact  raw  post  readout+post")
    for t in ts:
        j = J1[f"pdens_t{t}"]; p = dist(j)
        print(f"  {t}  {j['exact_alpha']:+.3f}  {observables(p)[0]:+.3f}  {observables(p, post=True)[0]:+.3f}"
              f"  {observables(p, A, post=True)[0]:+.3f}")
    D = {k: [] for k in ("raw", "ro", "post", "ropost", "kept")}
    print("packet density TVD:  t  raw  readout  post  readout+post  kept   | readout+post rho")
    for t in ts:
        j = J1[f"pdens_t{t}"]; p = dist(j); ex = j["exact_rho"]
        r = [observables(p)[1], observables(p, A)[1], observables(p, post=True)[1], observables(p, A, post=True)[1]]
        for k, v in zip(("raw", "ro", "post", "ropost"), r):
            D[k].append(tvd(ex, v))
        D["kept"].append(observables(p, post=True)[2])
        print(f"  {t}  " + "  ".join(f"{tvd(ex, v):.3f}" for v in r) + f"  {D['kept'][-1]:.3f}   {np.round(r[3], 3).tolist()}")

    # ---- repeats at t=3,4 (three runs of each parity circuit) ----
    rep = {}
    for nm in ("palpha_t3", "palpha_t4", "pdens_t3", "pdens_t4"):
        runs = [J1[nm], JR[nm + "_r1"], JR[nm + "_r2"]]
        rep[nm] = np.array([(observables(dist(j))[0], observables(dist(j), A, post=True)[0],
                             tvd(j["exact_rho"], observables(dist(j))[1]),
                             tvd(j["exact_rho"], observables(dist(j), A, post=True)[1])) for j in runs])
        r = rep[nm]
        print(f"{nm}: v raw {np.round(r[:, 0], 3)}  v ro+post {np.round(r[:, 1], 3)}  "
              f"TVD raw {np.round(r[:, 2], 3)}  TVD ro+post {np.round(r[:, 3], 3)}")
        print(f"      mean+-std: v raw {r[:, 0].mean():+.3f}+-{r[:, 0].std(ddof=1):.3f}  v ro+post "
              f"{r[:, 1].mean():+.3f}+-{r[:, 1].std(ddof=1):.3f}  TVD raw {r[:, 2].mean():.3f}+-{r[:, 2].std(ddof=1):.3f}"
              f"  TVD ro+post {r[:, 3].mean():.3f}+-{r[:, 3].std(ddof=1):.3f}")

    # ---- every raw t=4 velocity measurement of today (incl. untwirled controls of the twirl runs) ----
    odd = ~EVEN
    t4 = {
        "even parity": [J1["palpha_t4"]["qpu"], JR["palpha_t4_r1"]["qpu"], JR["palpha_t4_r2"]["qpu"],
                        JT["even_all_t4_r1"]["qpu"], JT["even_all_t4_r2"]["qpu"]]
                       + [w["even_all_t4_native"]["qpu"] for w in TW],
        "odd parity": [JT["odd_all_t4_r1"]["qpu"], JT["odd_all_t4_r2"]["qpu"]],
        "original state,\nall qubits read": [JT["step1_all_t4_r1"]["qpu"], JT["step1_all_t4_r2"]["qpu"]],
        "original state,\n$q_1$ read": [JN["alpha_t4"]["qpu"], JR["alpha_t4"]["qpu"],
                                      JT["step1_q1_t4_r1"]["qpu"], JT["step1_q1_t4_r2"]["qpu"]]
                                     + [w["step1_q1_t4_native"]["qpu"] for w in TW],
    }
    exact4 = J1["palpha_t4"]["exact_alpha"]
    print(f"\nt=4 raw velocity (exact {exact4:+.3f}; shot noise ~{np.sqrt(1 / S):.3f}):")
    for k, v in t4.items():
        print(f"  {k.replace(chr(10), ' '):28s} {np.round(v, 3).tolist()}  mean {np.mean(v):+.3f}  std {np.std(v, ddof=1):.3f}")
    for nm, sec in (("even_all_t4_r1", EVEN), ("even_all_t4_r2", EVEN), ("odd_all_t4_r1", odd), ("odd_all_t4_r2", odd)):
        v, _, kept = observables(dist(JT[nm]), post=True, sector=sec)
        print(f"  {nm}: post-selected {v:+.3f} (kept {kept:.3f})")
    print(f"step-1 t=3 velocity: {JN['alpha_t3']['qpu']:+.3f} (first run), {JR['alpha_t3']['qpu']:+.3f} (repeat)")

    if not args.no_sim:
        phis, sim = simulate_phase_sensitivity()
        print("\nsimulated even-state velocity, 2% depolarizing + coherent RZ(phi) on q1 per 2q gate:")
        for phi in phis:
            print(f"  phi={phi:+.2f}: " + "  ".join(f"{sim[(t, phi)]:+.3f}" for t in ts))

    # ---- figure ----
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(16, 4.4),
                                        gridspec_kw={"width_ratios": [1.05, 1.05, 1.0]})
    E = float(R1["E"])
    tsm = np.linspace(0, 4, 400)
    axA.plot(tsm, np.sin(2 * E * tsm), color="0.6", lw=1.6, label=r"exact  $\sin 2Et$")
    axA.plot(ts, V["exact"], "o", color="0.4", ms=5)
    axA.plot(ts, V["raw"], "o", color="C4", ms=9, mfc="none", mew=1.7, label="raw")
    axA.plot(ts, V["ropost"], "o", color="C4", ms=8, label="readout-mitigated + parity post-selected")
    for t, nm in ((3, "palpha_t3"), (4, "palpha_t4")):
        axA.plot([t + 0.12] * 2, rep[nm][1:, 1], "o", color="C4", ms=4, alpha=0.6)
        axA.plot([t + 0.12] * 2, rep[nm][1:, 0], "o", color="C4", ms=4, mfc="none", alpha=0.6)
    axA.axhline(0, color="k", lw=0.6, ls=":")
    axA.set_xticks(ts)
    axA.set_xlabel("step  $t$")
    axA.set_ylabel(r"$\langle\alpha_x\rangle(t)$")
    axA.set_ylim(-1.18, 1.18)
    axA.set_title("(a) velocity, even-parity state")
    axA.legend(frameon=False, loc="lower left", fontsize=8.5)

    w = 0.2
    cols = (("raw", "C4", 0.3), ("ro", "C2", 0.5), ("post", "C4", 0.6), ("ropost", "C4", 1.0))
    labels = {"raw": "raw", "ro": "readout-mitigated", "post": "parity post-selected",
              "ropost": "both"}
    for k, (key, col, al) in enumerate(cols):
        axB.bar(ts + (k - 1.5) * w, D[key], w, color=col, alpha=al, label=labels[key] + " (run 1)")
    for t, nm in ((3, "pdens_t3"), (4, "pdens_t4")):
        for k, key in ((0, 2), (3, 3)):
            r = rep[nm][:, key]
            x = t + (k - 1.5) * w
            axB.plot([x + 0.07] * len(r), r, "o", color="k", ms=3.5, mfc="none", mew=0.9,
                     label="runs 1-3" if (t, k) == (3, 0) else None)
            axB.errorbar(x, r.mean(), yerr=r.std(ddof=1), fmt="_", color="k", ms=8, lw=1,
                         label="mean $\\pm$ std, runs 1-3" if (t, k) == (3, 0) else None)
    axB.set_xticks(ts)
    axB.set_xlabel("step  $t$")
    axB.set_ylabel("TVD from exact density")
    axB.set_title("(b) packet density error")
    axB.legend(frameon=False, loc="upper left", fontsize=8.5)

    for k, (nm, v) in enumerate(t4.items()):
        axC.plot(np.full(len(v), k) + np.linspace(-0.12, 0.12, len(v)), v, "o", color="C4", ms=7)
    axC.axhline(exact4, color="0.4", lw=1.2, label=f"exact  {exact4:+.3f}")
    axC.axhline(0, color="k", lw=0.6, ls=":")
    axC.set_xticks(range(len(t4)))
    axC.set_xticklabels(list(t4), fontsize=8.5)
    axC.set_ylabel(r"$\langle\alpha_x\rangle(4)$, raw")
    axC.set_xlim(-0.5, len(t4) - 0.5)
    axC.set_ylim(-0.45, 0.35)
    axC.set_title("(c) repeated runs at $t=4$")
    axC.legend(frameon=False, loc="lower left", fontsize=8.5)
    fig.tight_layout()
    out = f"{args.outdir}/hw_4site_qcs_parity.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
