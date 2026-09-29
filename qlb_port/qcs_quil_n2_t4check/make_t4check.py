"""Export the t=4 diagnostic circuits: even/odd parity sectors and the step-1 state, read on all qubits,
plus the step-1 q1-only circuit.  Usage: PYTHONPATH=. python /tmp/make_t4.py <outdir>"""
import json, os, sys
import numpy as np
from qlb_port import operators as ops, run_hardware_zitterbewegung as H
from qlb_port.export_quil import circuit_to_quil

out = sys.argv[1]; os.makedirs(out, exist_ok=True)
R = ops.ROTATIONS["x"]; Ri = R.conj().T; npos = 2; m = 0.8; t = 4


def sector(psi, keep):
    g = (Ri @ psi.reshape(4, 4).T).T
    g[:, [k for k in range(4) if k not in keep]] = 0
    o = (R @ g.T).T.reshape(-1)
    return o / np.linalg.norm(o)


base, E = H.mode_state(npos, 0.0, m)
exact = float(H.alpha_x_expectation(H.classical_evolution(base, npos, m, t)[t], npos))
specs = [("even_all", sector(base, (0, 3)), H.parity_circuit, "palpha"),
         ("odd_all", sector(base, (1, 2)), H.parity_circuit, "palpha"),
         ("step1_all", base, H.parity_circuit, "palpha"),
         ("step1_q1", base, H.alpha_x_circuit, "alpha")]
jobs = []
for nm, psi, build, kind in specs:
    quil, measured = circuit_to_quil(build(npos, psi, m, t, streaming_method="fourier"))
    open(os.path.join(out, f"{nm}_t4.quil"), "w").write(quil)
    for r in (1, 2):
        jobs.append({"name": f"{nm}_t4_r{r}", "kind": kind, "t": t, "quil": f"{nm}_t4.quil",
                     "npos": npos, "ro_bits": len(measured), "exact": exact, "exact_alpha": exact,
                     "exact_rho": [0.25] * 4})
json.dump({"description": "t=4 velocity diagnostic: parity sectors and readout width",
           "npos": npos, "n_qubits": 4, "mass": m, "k0": 0.0, "E": float(E), "streaming": "fourier",
           "jobs": jobs}, open(os.path.join(out, "manifest.json"), "w"), indent=2)
print(f"exact <alpha_x>(4) = {exact:+.4f};", [j["name"] for j in jobs])
