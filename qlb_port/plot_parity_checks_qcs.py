"""
Analysis of the in-circuit parity checks of the QLB velocity circuits on Rigetti QCS
(see export_parity_checks.py).

For each step t and variant (no checks / checks after every step but the last) the velocity
<alpha_x> = P(q1=1) - P(q1=0) is computed
  * raw, from all shots;
  * after symmetry verification, keeping the shots with even final parity q0 xor q1;
  * with the checks as well, keeping only shots whose check ancillas also all read 0,
together with the kept fractions.  Rounds are pooled; the per-round values give the run-to-run spread.
With --noise-check the programs are compiled with quilc (as for the device) and evolved as density
matrices with 0.8% two-qubit depolarizing noise after every CZ, 0.1% after every RX, and 1.5% readout
flips, instead of using device data (needs a quilc server).

Usage:
    python plot_parity_checks_qcs.py [--results file] [--noise-check] [--outdir dir]
"""

import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "qcs_quil"))

RES = "qlb_port/qcs_quil_n2_pcheck/qcs_results_pcheck.json"
LAYOUT = [65, 74, 92, 83, 56, 64, 66]


def metrics(prob, n, n_anc):
    """prob over integers with bit k = qubit k (0 q0, 1 q1, 2-3 position, 4.. ancillas)."""
    idx = np.arange(len(prob))
    b = lambda k: (idx >> k) & 1
    v = lambda m: float((prob[m] * (2 * b(1)[m] - 1)).sum() / prob[m].sum())
    even = (b(0) ^ b(1)) == 0
    ok = even.copy()
    for k in range(4, 4 + n_anc):
        ok &= b(k) == 0
    allm = np.ones(len(prob), bool)
    return {"raw": v(allm), "par": v(even), "chk": v(ok),
            "k_par": float(prob[even].sum()), "k_chk": float(prob[ok].sum())}


def noisy_distribution(quil_text, n, p2=0.008, p1=0.001, pm=0.015):
    from pyquil import Program, get_qc
    from run_qcs import _pin_text
    from cdr_qcs import parse
    qc = get_qc("Cepheus-1-108Q", compiler_timeout=300)
    nat = str(qc.compiler.quil_to_native_quil(Program(_pin_text(quil_text, dict(enumerate(LAYOUT)))),
                                              protoquil=True))
    gates, _ = parse(nat)
    pos = {q: k for k, q in enumerate(LAYOUT[:n])}
    I2 = np.eye(2)
    P = [I2, np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])]
    rho = np.zeros([2] * (2 * n), complex)
    rho[(0,) * (2 * n)] = 1                       # axes: ket qubits n-1..0 then bra qubits n-1..0

    def ax(q):
        return n - 1 - pos[q]

    def apply1(r, U, q):
        a = ax(q)
        r = np.moveaxis(np.tensordot(U, r, axes=([1], [a])), 0, a)
        return np.moveaxis(np.tensordot(U.conj(), r, axes=([1], [n + a])), 0, n + a)

    def dep1(r, q, p):
        return (1 - p) * r + p / 4 * sum(apply1(r, S, q) for S in P)

    for g, ang, qs in gates:
        if g == "CZ":
            d = np.ones([2] * (2 * n))
            ka, kb = ax(qs[0]), ax(qs[1])
            sl = [slice(None)] * (2 * n)
            sl[ka], sl[kb] = 1, 1
            d[tuple(sl)] *= -1
            sl = [slice(None)] * (2 * n)
            sl[n + ka], sl[n + kb] = 1, 1
            d[tuple(sl)] *= -1
            rho = rho * d
            new = (1 - p2) * rho
            for A in P:
                for B in P:
                    new = new + p2 / 16 * apply1(apply1(rho, A, qs[0]), B, qs[1])
            rho = new
        else:
            c, s = math.cos(ang / 2), math.sin(ang / 2)
            U = (np.array([[c, -1j * s], [-1j * s, c]]) if g == "RX"
                 else np.diag([np.exp(-1j * ang / 2), np.exp(1j * ang / 2)]))
            rho = apply1(rho, U, qs[0])
            if g == "RX":
                rho = dep1(rho, qs[0], p1)
    prob = np.real(np.einsum("".join(chr(97 + i) for i in range(n)) * 2 + "->"
                             + "".join(chr(97 + i) for i in range(n)), rho))   # axis m <-> qubit n-1-m
    flip = np.array([[1 - pm, pm], [pm, 1 - pm]])
    for a in range(n):
        prob = np.moveaxis(np.tensordot(flip, prob, axes=([1], [a])), 0, a)
    return prob.reshape(-1)                          # flat index: bit k <-> qubit k


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default=RES)
    ap.add_argument("--noise-check", action="store_true")
    ap.add_argument("--diag", action="store_true",
                    help="summarize the diagnostic runs (export_check_diagnostics.py) instead")
    ap.add_argument("--outdir", default="proposal/quairk/ancillaStreaming/figures")
    args = ap.parse_args()
    if args.diag:
        for f, cols in (("qlb_port/qcs_quil_n2_pdiag/qcs_results_pdiag.json",
                         "q0 q1 p0 p1 | a56 a64 a66"),
                        ("qlb_port/qcs_quil_n2_czrep/qcs_results_czrep.json", "q0 | a56 a64")):
            R_ = json.load(open(f))
            print(f"{f}: {R_.get('total_execution_us', 0) / 1e6:.1f} s QPU; P(1) of {cols}")
            for j in R_["jobs"]:
                b = np.array([[int(c) for c in s] for s in j["bitstrings"]])
                extra = (f"   parity odd {float((b[:, 0] != b[:, 1]).mean()):.3f}" if b.shape[1] == 7 else "")
                print(f"   {j['name']:12s} " + " ".join(f"{x:.3f}" for x in b.mean(axis=0)) + extra)
        return
    R_ = json.load(open(args.results))
    jobs = [j for j in R_["jobs"] if j.get("role") == "pcheck"]
    base = os.path.dirname(args.results)

    if args.noise_check:
        print("noise model (quilc-compiled programs, density matrix)")
        print(" t  variant   exact    raw      parity   checks   kept(parity) kept(checks)")
        for j in [j for j in jobs if j["round"] == 0]:
            n = 4 + j["n_anc"]
            m = metrics(noisy_distribution(open(os.path.join(base, j["quil"])).read(), n), n, j["n_anc"])
            print(f" {j['t']}  {'checks' if j['checks'] else 'none  '}  {j['exact']:+.3f}   {m['raw']:+.3f}"
                  f"   {m['par']:+.3f}   {m['chk']:+.3f}    {m['k_par']:.3f}        {m['k_chk']:.3f}")
        return

    print(f"{len(jobs)} jobs, {R_.get('total_execution_us', 0) / 1e6:.1f} s QPU execution")
    groups = {}
    for j in jobs:
        if "bitstrings" not in j:
            continue
        n = 4 + j["n_anc"]
        ints = np.array([sum(int(c) << k for k, c in enumerate(s)) for s in j["bitstrings"]])
        prob = np.bincount(ints, minlength=2 ** n) / len(ints)
        groups.setdefault((j["t"], j["checks"]), []).append((metrics(prob, n, j["n_anc"]), prob, j))

    print("\n t  variant   exact    raw      parity   checks   kept(parity) kept(checks)   per-round checks")
    rows = {}
    for (t, c), lst in sorted(groups.items()):
        j0 = lst[0][2]
        n = 4 + j0["n_anc"]
        m = metrics(np.mean([p for _, p, _ in lst], axis=0), n, j0["n_anc"])
        rows[(t, c)] = (m, j0["exact"], [x["chk"] for x, _, _ in lst], [x["par"] for x, _, _ in lst])
        print(f" {t}  {'checks' if c else 'none  '}  {j0['exact']:+.3f}   {m['raw']:+.3f}   {m['par']:+.3f}"
              f"   {m['chk']:+.3f}    {m['k_par']:.3f}        {m['k_chk']:.3f}         "
              + " ".join(f"{x:+.3f}" for x in rows[(t, c)][2]))
    print("\n errors caught only by the checks (fraction of shots with even final parity but a check = 1):")
    for (t, c), (m, _, _, _) in sorted(rows.items()):
        if c:
            print(f"   t = {t}: {m['k_par'] - m['k_chk']:.3f}")
    print("\n per job: P(check ancilla = 1) for checks after steps 1, 2, ..., and P(final parity odd)")
    for j in jobs:
        if "bitstrings" in j:
            b = np.array([[int(c) for c in s] for s in j["bitstrings"]])
            print(f"   {j['name']:14s} " + " ".join(f"{x:.3f}" for x in b[:, 4:].mean(axis=0))
                  + f"   odd {float((b[:, 0] != b[:, 1]).mean()):.3f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    E = float(R_["E"])
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.5, 4.2), gridspec_kw={"width_ratios": [1.3, 1]})
    ts = np.linspace(0, 4, 300)
    axA.plot(ts, np.sin(2 * E * ts), color="0.6", lw=1.5, label=r"exact $\sin 2Et$")
    tt = sorted({t for t, _ in rows})
    axA.plot(np.array(tt) - 0.1, [rows[(t, 0)][0]["par"] for t in tt], "o", color="C0", mfc="none", ms=8,
             label="no checks, final parity")
    tc = sorted(t for t, c in rows if c)
    axA.plot(np.array(tc) + 0.1, [rows[(t, 1)][0]["par"] for t in tc], "s", color="C1", mfc="none", ms=8,
             label="checks present, final parity only")
    axA.plot(np.array(tc) + 0.1, [rows[(t, 1)][0]["chk"] for t in tc], "s", color="C1", ms=7,
             label="checks present, all checks")
    axA.axhline(0, color="k", lw=0.6, ls=":")
    axA.set_xlabel("step $t$")
    axA.set_ylabel(r"$\langle\alpha_x\rangle(t)$")
    axA.set_title("(a) velocity after post-selection")
    axA.legend(frameon=False, fontsize=8, loc="lower left")
    axB.plot(tt, [rows[(t, 0)][0]["k_par"] for t in tt], "o-", color="C0", label="no checks, final parity")
    axB.plot(tc, [rows[(t, 1)][0]["k_par"] for t in tc], "s--", color="C1", mfc="none",
             label="checks present, final parity only")
    axB.plot(tc, [rows[(t, 1)][0]["k_chk"] for t in tc], "s-", color="C1", label="checks present, all checks")
    axB.set_xlabel("step $t$")
    axB.set_ylabel("fraction of shots kept")
    axB.set_title("(b) shots kept")
    axB.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    out = os.path.join(args.outdir, "hw_4site_qcs_pcheck.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
