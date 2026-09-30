"""
Phase-error model for the four-site velocity measurements: does one coherent error explain the
deviations at t = 2 and t = 4 together?

For every batch that measured the velocity at t = 1..4 (and t = 0 where available) the raw values are
fitted, weighted by shot noise, with
  M0: v = lam * sin(2E t) + delta                  (damping and offset)
  M1: v = lam * sin(2E t + phi) + delta            (constant phase shift)
  M2: v = lam * sin((2E + eps) t) + delta          (phase error eps per step)
A phase error shifts the steep points t = 2 and t = 4 in opposite directions (Eq. eq:phase of the paper),
while damping and offset cannot; the residuals at t = 2 and 4 under M0 and M2 test this.

Usage:
    PYTHONPATH=. python -m qlb_port.fit_phase_model
"""

import json

import numpy as np
from scipy.optimize import least_squares

B = "qlb_port/"


def pooled(jobs):
    n1 = sum(j["counts"].get("1", 0) for j in jobs)
    n = sum(sum(j["counts"].values()) for j in jobs)
    return 2 * n1 / n - 1, n


def batches():
    out = {}
    d = np.load(B + "hw_4site_fourier_trembling.npz")
    out["Open Quantum (Sec. 3)"] = [(int(t), float(v), int(d["shots"])) for t, v in zip(d["t"], d["av_hw"])]
    J = json.load(open(B + "qcs_quil_n2_fourier/qcs_results_n2f_pinned.json"))["jobs"]
    out["QCS raw (Sec. 4.3)"] = [(j["t"], j["qpu"], 2000) for j in J if j["kind"] == "alpha"]
    J = json.load(open(B + "qcs_quil_n2_parity/qcs_results_parity_pinned.json"))["jobs"]
    out["QCS even parity, raw (Sec. 4.5)"] = [(j["t"], j["qpu"], 2000) for j in J if j["kind"] == "palpha"]
    J = json.load(open(B + "qcs_quil_n2_zne/qcs_results_zne.json"))["jobs"]
    out["twirled, ZNE batch (Sec. 4.8)"] = [
        (t, *pooled([j for j in J if j.get("group") == f"alpha_t{t}" and j.get("scale") == 1]))
        for t in range(1, 5)]
    J = json.load(open(B + "qcs_quil_n2_cdr/qcs_results_cdr.json"))["jobs"]
    out["twirled, CDR batch (Sec. 4.10)"] = [
        (t, *pooled([j for j in J if j.get("group") == f"alpha_t{t}"])) for t in range(1, 5)]
    J = json.load(open(B + "qcs_quil_n2_pcheck/qcs_results_pcheck.json"))["jobs"]
    rows = []
    for t in range(1, 5):
        bits = [s for j in J if j["t"] == t and not j["checks"] for s in j["bitstrings"]]
        rows.append((t, float(np.mean([2 * int(s[1]) - 1 for s in bits])), len(bits)))
    out["even parity, raw, new qubits (Sec. 5.5)"] = rows
    return out


def fit(t, v, s, E, model):
    def f(p):
        if model == 0:
            return p[0] * np.sin(2 * E * t) + p[1]
        if model == 1:
            return p[0] * np.sin(2 * E * t + p[2]) + p[1]
        return p[0] * np.sin((2 * E + p[2]) * t) + p[1]
    p0 = [0.8, 0.0] + ([] if model == 0 else [0.0])
    r = least_squares(lambda p: (f(p) - v) / s, p0)
    return r.x, float((r.fun ** 2).sum()), f(r.x)


def main():
    E = float(json.load(open(B + "qcs_quil_n2_cdr/qcs_results_cdr.json"))["E"])
    print(f"2E = {2 * E:.4f}; exact v(t) = sin(2Et)\n")
    for name, rows in batches().items():
        t = np.array([r[0] for r in rows], float)
        v = np.array([r[1] for r in rows])
        n = np.array([r[2] for r in rows])
        s = np.sqrt(np.maximum(1 - v ** 2, 0.05) / n)
        print(f"{name}: t = {t.astype(int).tolist()}, v = {np.round(v, 3).tolist()}")
        for model, label in ((0, "M0 damping+offset"), (1, "M1 + constant phase"), (2, "M2 + phase per step")):
            p, chi2, fv = fit(t, v, s, E, model)
            res = v - fv
            extra = "" if model == 0 else (f", phi = {p[2]:+.3f}" if model == 1 else f", eps = {p[2]:+.3f} per step")
            r2 = res[t == 2][0] if 2 in t else float("nan")
            r4 = res[t == 4][0] if 4 in t else float("nan")
            print(f"   {label:22s} lam = {p[0]:.3f}, delta = {p[1]:+.3f}{extra};  chi2 = {chi2:7.1f}"
                  f" (dof {len(t) - len(p)});  residual t=2 {r2:+.3f}, t=4 {r4:+.3f}")
        print()


if __name__ == "__main__":
    main()
