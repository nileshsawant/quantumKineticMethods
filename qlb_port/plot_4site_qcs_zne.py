"""
Zero-noise extrapolation and mirror rescaling of the twirled four-site velocity circuits (t=1..4) on QCS.

Reads qlb_port/qcs_quil_n2_zne/qcs_results_zne.json: for each t, 16 Pauli-twirled variants (125 shots
each) at CZ fold factors 1, 3, 5, 16 twirled variants of the mirror circuit (circuit then inverse,
ideal <alpha_x> = -1), and a q1 readout calibration.  Every pooled value is readout-corrected with the
2x2 confusion matrix of q1.  Estimators:
    ZNE, linear in the fold factor through s = 1, 3;
    ZNE, exponential  v(s) = v0 exp(-k s)  fitted to s = 1, 3, 5;
    rescaling  v1 / sqrt(m),  m = -(mirror value), the decay of a circuit of twice the depth.
Error bars: parametric bootstrap of the binomial shot counts (2000 resamples).

Writes <outdir>/hw_4site_qcs_zne.png and prints the numbers used in the paper.

Usage:
    PYTHONPATH=. python3 -m qlb_port.plot_4site_qcs_zne
"""

import argparse
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = "qlb_port/qcs_quil_n2_zne/qcs_results_zne.json"
SCALES = (1, 3, 5)


def pooled_p1(jobs):
    """Pooled (shots reading 1, total shots) of the single measured bit."""
    n1 = sum(j["counts"].get("1", 0) for j in jobs)
    n = sum(sum(j["counts"].values()) for j in jobs)
    return n1, n


def value(p1, M):
    q = np.clip(np.linalg.solve(M, [1 - p1, p1]), 0, None)
    q = q / q.sum()
    return q[1] - q[0]


def estimators(v):
    """v: dict scale -> value, plus 'mirror'.  Returns (lin, exp, rescaled)."""
    lin = v[1] - (v[3] - v[1]) / 2
    y = np.array([v[s] for s in SCALES])
    if np.all(np.sign(y) == np.sign(y[0])) and np.all(y != 0):
        k, lna = np.polyfit(np.array(SCALES, float), np.log(np.abs(y)), 1)
        ex = np.sign(y[0]) * np.exp(lna)
    else:
        ex = np.nan
    m = -v["mirror"]
    resc = v[1] / np.sqrt(m) if m > 0 else np.nan
    return lin, ex, resc


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default="proposal/quairk/ancillaStreaming/figures")
    ap.add_argument("--results", default=RES)
    args = ap.parse_args()

    R = json.load(open(args.results))
    J = R["jobs"]
    print(f"QCS ZNE batch: {len(J)} jobs, {R['total_execution_us'] / 1e6:.1f} s QPU execution")
    cal = {j["prepared"]: j for j in J if j["group"] == "cal"}
    c0, c1 = cal[0]["counts"], cal[1]["counts"]
    e01 = c0.get("1", 0) / sum(c0.values())
    e10 = c1.get("0", 0) / sum(c1.values())
    M = np.array([[1 - e01, e10], [e01, 1 - e10]])
    print(f"q1 readout: P(1|0)={e01:.3f}  P(0|1)={e10:.3f}")

    rng = np.random.default_rng(7)
    B = 2000
    rows = {}
    print("\n t  exact   s=1     s=3     s=5    mirror | ZNE lin   ZNE exp   rescaled   (bootstrap sd)")
    for t in range(1, 5):
        grp = {s: [j for j in J if j["group"] == (f"alpha_t{t}" if s == 1 else f"alpha_t{t}_x{s}")] for s in SCALES}
        grp["mirror"] = [j for j in J if j["group"] == f"alpha_t{t}_mirror"]
        counts = {k: pooled_p1(v) for k, v in grp.items()}
        v = {k: value(n1 / n, M) for k, (n1, n) in counts.items()}
        est = estimators(v)
        boot = []
        for _ in range(B):
            vb = {k: value(rng.binomial(n, n1 / n) / n, M) for k, (n1, n) in counts.items()}
            boot.append(estimators(vb))
        boot = np.array(boot)
        sd = np.nanstd(boot, axis=0)
        exact = grp[1][0]["exact"]
        sdv = {k: np.sqrt(max(1 - v[k] ** 2, 1e-4) / counts[k][1]) for k in v}
        rows[t] = dict(exact=exact, v=v, sdv=sdv, est=est, sd=sd)
        print(f" {t}  {exact:+.3f}  {v[1]:+.3f}  {v[3]:+.3f}  {v[5]:+.3f}  {v['mirror']:+.3f} | "
              f"{est[0]:+.3f}    {est[1]:+.3f}    {est[2]:+.3f}    ({sd[0]:.3f}, {sd[1]:.3f}, {sd[2]:.3f})")
    for t, r in rows.items():
        parts = []
        for s in SCALES:
            vals = [j["qpu"] for j in J if j["group"] == (f"alpha_t{t}" if s == 1 else f"alpha_t{t}_x{s}")]
            parts.append(f"s={s}: {np.std(vals, ddof=1):.3f} (shot {np.sqrt(max(1 - np.mean(vals) ** 2, 0) / 125):.3f})")
        print(f" t={t} variant spread, raw: " + "  ".join(parts) + f"   mirror decay m={-r['v']['mirror']:.3f}")

    # ---- figure: (a) velocity vs t, (b) decay with fold factor ----
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.5, 4.4), gridspec_kw={"width_ratios": [1.15, 1.0]})
    E = float(R["E"])
    tsm = np.linspace(0, 4, 400)
    ts = np.array(sorted(rows))
    axA.plot(tsm, np.sin(2 * E * tsm), color="0.6", lw=1.6, label=r"exact  $\sin 2Et$")
    axA.plot(ts, [rows[t]["exact"] for t in ts], "o", color="0.4", ms=5)
    axA.errorbar(ts - 0.12, [rows[t]["v"][1] for t in ts], yerr=[rows[t]["sdv"][1] for t in ts], fmt="o",
                 color="C4", mfc="none", ms=8, mew=1.6, capsize=3, label="twirled, readout-corrected")
    axA.errorbar(ts, [rows[t]["est"][1] for t in ts], yerr=[rows[t]["sd"][1] for t in ts], fmt="s",
                 color="C1", ms=7, capsize=3, label="zero-noise extrapolation (exponential)")
    axA.errorbar(ts + 0.12, [rows[t]["est"][2] for t in ts], yerr=[rows[t]["sd"][2] for t in ts], fmt="^",
                 color="C2", ms=8, capsize=3, label="mirror rescaling")
    axA.axhline(0, color="k", lw=0.6, ls=":")
    axA.set_xticks(range(5))
    axA.set_xlabel("step  $t$")
    axA.set_ylabel(r"$\langle\alpha_x\rangle(t)$")
    axA.set_ylim(-1.35, 1.35)
    axA.set_title("(a) velocity, twirled and mitigated")
    axA.legend(frameon=False, loc="lower left", fontsize=8.5)

    ss = np.linspace(0, 5.3, 100)
    for t, col in zip(ts, ("C0", "C1", "C2", "C3")):
        r = rows[t]
        y = np.array([r["v"][s] for s in SCALES])
        axB.errorbar(SCALES, y, yerr=[r["sdv"][s] for s in SCALES], fmt="o", color=col, capsize=3, label=f"$t={t}$")
        if np.isfinite(r["est"][1]):
            k, lna = np.polyfit(np.array(SCALES, float), np.log(np.abs(y)), 1)
            axB.plot(ss, np.sign(y[0]) * np.exp(lna + k * ss), "-", color=col, lw=1)
        axB.plot(0, r["exact"], "*", color=col, ms=11, mec="k", mew=0.5)
    axB.axhline(0, color="k", lw=0.6, ls=":")
    axB.set_xlabel("CZ fold factor  $s$  (noise scale)")
    axB.set_ylabel(r"$\langle\alpha_x\rangle$, twirled")
    axB.set_xlim(-0.3, 5.5)
    axB.set_title("(b) extrapolation to zero noise (stars: exact)")
    axB.legend(frameon=False, fontsize=8.5, loc="center right")
    fig.tight_layout()
    out = f"{args.outdir}/hw_4site_qcs_zne.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
