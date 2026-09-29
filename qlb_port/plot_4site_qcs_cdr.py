"""
Clifford data regression of the twirled four-site velocity circuits (t=1..4) on Rigetti QCS.

Reads qlb_port/qcs_quil_n2_cdr/qcs_results_cdr.json: per t, 16 twirled variants of the target
(125 shots each) and 32 twirled near-Clifford training circuits (250 shots each) whose exact values are
stored in the manifest.  For each t a straight line  ideal = a + b * measured  is fitted to the training
circuits (least squares) and applied to the pooled target value.  Uncertainties: bootstrap over the
training circuits and the binomial shot counts (2000 resamples).

Writes <outdir>/hw_4site_qcs_cdr.png and prints the numbers used in the paper.  With --noise-check it
instead applies the same regression to an exact density-matrix simulation of the programs with a
simple noise model (2% two-qubit depolarizing and amplitude damping gamma=0.01 on every qubit after
each CZ, readout errors 0.3% / 2.2%).

Usage:
    PYTHONPATH=. python3 -m qlb_port.plot_4site_qcs_cdr [--results file] [--noise-check]
"""

import argparse
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = "qlb_port/qcs_quil_n2_cdr/qcs_results_cdr.json"


def counts(j):
    n1 = j["counts"].get("1", 0)
    return n1, sum(j["counts"].values())


def noise_check(d, qubits=(100, 101, 102, 103), meas=102, p2=0.02, gamma=0.01, e01=0.003, e10=0.022):
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "qcs_quil"))
    from cdr_qcs import parse
    idx = {q: k for k, q in enumerate(qubits)}
    I2, X, Y, Z = np.eye(2), np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])

    def op(mats):
        out = np.array([[1]])
        for q in qubits:
            out = np.kron(out, mats.get(q, I2))
        return out

    def unitary(name, ang, q):
        if name == "CZ":
            return np.diag([-1 if (i >> (3 - idx[q[0]])) & 1 and (i >> (3 - idx[q[1]])) & 1 else 1
                            for i in range(16)]).astype(complex)
        c, s = np.cos(ang / 2), np.sin(ang / 2)
        u = (np.array([[c, -1j * s], [-1j * s, c]]) if name == "RX"
             else np.diag([np.exp(-1j * ang / 2), np.exp(1j * ang / 2)]))
        return op({q[0]: u})

    damp = [(op({q: np.array([[1, 0], [0, np.sqrt(1 - gamma)]])}), op({q: np.array([[0, np.sqrt(gamma)], [0, 0]])}))
            for q in qubits]

    def run(gates):
        rho = np.zeros((16, 16), complex)
        rho[0, 0] = 1
        for g in gates:
            U = unitary(*g)
            rho = U @ rho @ U.conj().T
            if g[0] == "CZ":
                new = (1 - p2) * rho
                for P in (I2, X, Y, Z):
                    for Q in (I2, X, Y, Z):
                        K = op({g[2][0]: P, g[2][1]: Q})
                        new += p2 / 16 * K @ rho @ K.conj().T
                rho = new
                for K0, K1 in damp:
                    rho = K0 @ rho @ K0.conj().T + K1 @ rho @ K1.conj().T
        p = np.real(np.diag(rho))
        p1 = sum(p[i] for i in range(16) if (i >> (3 - idx[meas])) & 1)
        return 2 * (p1 * (1 - e10) + (1 - p1) * e01) - 1

    M = json.load(open(os.path.join(d, "manifest.json")))["jobs"]
    print(" t  exact   noisy   CDR     fit a    fit b")
    for t in range(1, 5):
        tgt = [j for j in M if j["group"] == f"alpha_t{t}"][0]
        vt = run(parse(open(os.path.join(d, tgt["quil"])).read())[0])
        tr = [j for j in M if j["group"] == f"alpha_t{t}_cdr"]
        xi = np.array([run(parse(open(os.path.join(d, j["quil"])).read())[0]) for j in tr])
        yi = np.array([j["ideal"] for j in tr])
        b, a = np.polyfit(xi, yi, 1)
        print(f" {t}  {tgt['exact']:+.3f}  {vt:+.3f}  {a + b * vt:+.3f}  {a:+.3f}   {b:.3f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default="proposal/quairk/ancillaStreaming/figures")
    ap.add_argument("--results", default=RES)
    ap.add_argument("--noise-check", action="store_true")
    args = ap.parse_args()
    if args.noise_check:
        noise_check(os.path.dirname(RES))
        return

    R = json.load(open(args.results))
    J = R["jobs"]
    print(f"CDR batch: {len(J)} jobs, {R.get('total_execution_us', 0) / 1e6:.1f} s QPU execution")
    rng = np.random.default_rng(11)
    rows = {}
    print("\n t  exact   target(raw)   fit a      fit b     CDR estimate   (bootstrap sd)   train rms resid")
    for t in range(1, 5):
        tgt = [j for j in J if j["group"] == f"alpha_t{t}"]
        trn = [j for j in J if j["group"] == f"alpha_t{t}_cdr"]
        n1t = sum(counts(j)[0] for j in tgt); nt = sum(counts(j)[1] for j in tgt)
        vt = 2 * n1t / nt - 1
        xi = np.array([2 * counts(j)[0] / counts(j)[1] - 1 for j in trn])
        yi = np.array([j["ideal"] for j in trn])
        ni = np.array([counts(j)[1] for j in trn])
        b, a = np.polyfit(xi, yi, 1)
        est = a + b * vt
        resid = yi - (a + b * xi)
        boot = []
        for _ in range(2000):
            pick = rng.integers(0, len(trn), len(trn))
            xb = 2 * rng.binomial(ni[pick], (xi[pick] + 1) / 2) / ni[pick] - 1
            bb, ab = np.polyfit(xb, yi[pick], 1)
            vtb = 2 * rng.binomial(nt, n1t / nt) / nt - 1
            boot.append(ab + bb * vtb)
        sd = float(np.std(boot))
        exact = tgt[0]["exact"]
        rows[t] = dict(exact=exact, vt=vt, a=a, b=b, est=est, sd=sd, xi=xi, yi=yi)
        print(f" {t}  {exact:+.3f}   {vt:+.3f}        {a:+.3f}    {b:.3f}     {est:+.3f}         ({sd:.3f})"
              f"          {np.sqrt(np.mean(resid ** 2)):.3f}")
        vv = np.array([2 * counts(j)[0] / counts(j)[1] - 1 for j in tgt])
        print(f"    target variants: sd {vv.std(ddof=1):.3f} (shot noise {np.sqrt(np.mean((1 - vv ** 2) / 125)):.3f});"
              f" training ideals {yi.min():+.2f}..{yi.max():+.2f}; |CDR - exact| {abs(est - exact):.3f}"
              f" = {abs(est - exact) / sd:.1f} sd; raw error {abs(vt - exact):.3f};"
              f" shot-noise part of train resid {b * np.sqrt(np.mean((1 - xi ** 2) / ni)):.3f}")
        # q1 readout errors (0->1 0.3%, 1->0 2.2%) measured in the ZNE batch, for comparison only
        print(f"    target with ZNE-batch readout correction: {(vt - (0.003 - 0.022)) / (1 - 0.003 - 0.022):+.3f}")
    print(f"  mean |error| over t=1..4: raw {np.mean([abs(r['vt'] - r['exact']) for r in rows.values()]):.3f},"
          f" CDR {np.mean([abs(r['est'] - r['exact']) for r in rows.values()]):.3f}")

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.5, 4.4), gridspec_kw={"width_ratios": [1.15, 1.0]})
    E = float(R["E"])
    tsm = np.linspace(0, 4, 400)
    ts = np.array(sorted(rows))
    axA.plot(tsm, np.sin(2 * E * tsm), color="0.6", lw=1.6, label=r"exact  $\sin 2Et$")
    axA.plot(ts, [rows[t]["exact"] for t in ts], "o", color="0.4", ms=5)
    axA.plot(ts - 0.08, [rows[t]["vt"] for t in ts], "o", color="C4", mfc="none", ms=8, mew=1.6,
             label="twirled (raw)")
    axA.errorbar(ts + 0.08, [rows[t]["est"] for t in ts], yerr=[rows[t]["sd"] for t in ts], fmt="s",
                 color="C1", ms=7, capsize=3, label="Clifford data regression")
    axA.axhline(0, color="k", lw=0.6, ls=":")
    axA.set_xticks(range(5))
    axA.set_xlabel("step  $t$")
    axA.set_ylabel(r"$\langle\alpha_x\rangle(t)$")
    axA.set_ylim(-1.35, 1.35)
    axA.set_title("(a) velocity, twirled and corrected")
    axA.legend(frameon=False, loc="lower left", fontsize=8.5)

    xs = np.linspace(-1, 1, 10)
    for t, col in zip(ts, ("C0", "C1", "C2", "C3")):
        r = rows[t]
        axB.plot(r["xi"], r["yi"], "o", color=col, ms=4, alpha=0.6)
        axB.plot(xs, r["a"] + r["b"] * xs, "-", color=col, lw=1.2,
                 label=f"$t={t}$: $a={r['a']:+.2f}$, $b={r['b']:.2f}$")
    axB.plot([-1, 1], [-1, 1], ":", color="k", lw=0.8)
    axB.set_xlabel("measured on the device (twirled training circuits)")
    axB.set_ylabel("ideal value")
    axB.set_xlim(-1.05, 1.05)
    axB.set_ylim(-1.1, 1.1)
    axB.set_title("(b) training circuits and linear fits")
    axB.legend(frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout()
    out = f"{args.outdir}/hw_4site_qcs_cdr.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
