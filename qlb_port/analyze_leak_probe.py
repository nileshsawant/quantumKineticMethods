"""Analyze leak_probe results: three-state IQ classification against the |0>, |1>, |2> reference
programs of the same run, target-partner correlations, and an IQ figure.
Usage: python analyze_leak_probe.py [results dir]"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

D = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "qcs_repmid_results")
PAIRS = {"A": [(86, 77), (104, 95)], "C": [(86, 95), (68, 77)]}

fig, ax = plt.subplots(2, 3, figsize=(13, 8.4))
for ri, ro in enumerate(("noshelve", "default")):
    res = json.load(open(os.path.join(D, f"leak_probe_{ro}.json")))
    jobs = {" ".join(map(str, j["case"])): j for j in res["jobs"]}
    iq = lambda j, q: np.array(j["iq"][f"q{q}_unclassified"])
    means, sig = {}, {}
    for q in res["qubits"]:
        refs = [iq(jobs[f"ref {s}"], q) for s in (0, 1, 2)]
        means[q] = np.array([r.mean(0) for r in refs])
        sig[q] = np.sqrt(np.mean([r.var(0).sum() / 2 for r in refs[:2]]))

    def cls(j, q):
        """0, 1, 2 by the nearest reference mean."""
        return np.argmin(np.linalg.norm(iq(j, q)[:, None] - means[q][None], axis=2), axis=1)

    print(f"== readout {ro}: reference separations in units of sigma")
    for q in res["qubits"]:
        m = means[q]
        print(f"  q{q}: |0>-|1> {np.linalg.norm(m[1] - m[0]) / sig[q]:.1f}  "
              f"|1>-|2> {np.linalg.norm(m[2] - m[1]) / sig[q]:.1f}  "
              f"|0>-|2> {np.linalg.norm(m[2] - m[0]) / sig[q]:.1f}")
    print("  ref state -> fraction classified 0/1/2")
    for s in (0, 1, 2):
        print(f"  ref {s}: " + "  ".join(
            f"q{q} " + "/".join(f"{np.mean(cls(jobs[f'ref {s}'], q) == c):.2f}" for c in range(3))
            for q in res["qubits"]))
    print("  case            target: P(0) P(2)   partner: P(flip)  P(target 0 and partner flipped)")
    for name, j in jobs.items():
        if name.startswith("ref"):
            continue
        _, g, p, n = name.split()
        row = []
        for t, b in PAIRS[g]:
            ct, cb = cls(j, t), cls(j, b)
            flip = (np.minimum(cb, 1) != j["ideal"][str(b)])
            row.append(f"{t}<-{b}: {np.mean(ct == 0):.2f} {np.mean(ct == 2):.2f}  "
                       f"{np.mean(flip):.2f} {np.mean((ct == 0) & flip):.2f}")
        print(f"  {name:14s} " + "   ".join(row))
    for qi, q in enumerate((104, 86, 68)):
        a = ax[ri, qi]
        for s, c in zip((0, 1, 2), ("0.75", "0.45", "tab:green")):
            a.scatter(*iq(jobs[f"ref {s}"], q).T, s=4, c=c, label=f"ref |{s}>")
        g = "A" if q != 68 else "C"
        j = jobs[f"cz {'C' if q == 86 else g} 0 8"]
        a.scatter(*iq(j, q).T, s=5, c="tab:red", alpha=0.6, label="|1>, 8 CZ")
        a.set_title(f"{ro} readout, q{q}")
        a.legend(fontsize=7, markerscale=3)
plt.tight_layout()
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "leak_probe_iq.png")
plt.savefig(out, dpi=70)
print("wrote", out)
