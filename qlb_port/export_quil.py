"""
Export the 2-site Dirac-QLB circuits to plain Quil, for running on the Rigetti
QCS Cepheus-1 QPU with pyQuil (>= 4.7, on-demand access).
=============================================================================

Why this split.  The circuits are built and *validated* here in qiskit; QCS
speaks Quil.  Rather than install a whole qiskit->QCS bridge (the old
`qiskit-rigetti` is abandoned and pins ancient qiskit), we translate each
already-validated circuit 1:1 into generic Quil text.  quilc (on the QCS side)
then compiles that generic Quil to the Cepheus-1 native ISA {RX, RZ, iSWAP} and
routes it onto the chip -- exactly the job quilc exists to do.  Nothing about the
physics is re-derived on the Quil side.

The translation is faithful because qiskit and Quil share rotation/entangler
conventions:  qiskit RX/RY/RZ(theta) = exp(-i theta P / 2) = Quil RX/RY/RZ(theta),
and qiskit cx(c, t) = Quil CNOT c t.  Global phase is unobservable and dropped.
We transpile to the universal basis {rx, ry, rz, cx} first so every instruction
has a direct Quil equivalent; quilc re-optimizes into the native set anyway, so
this intermediate basis does not affect the circuit that finally runs.

Layout reminder (from run_hardware_zitterbewegung):  flat index i = x*4 + c,
so qubit 0,1 = spinor c, qubit 2.. = position x.  The alpha_x circuit measures
qubit 1 (its sign is the alpha_x eigenvalue after the R^-1 rotation); the density
circuit measures every qubit and rho(x) = sum_c |psi(x,c)|^2 with x = i >> 2.

This script talks to NO cloud service and needs NO credentials.  It writes:
    <outdir>/<name>.quil      one generic-Quil program per (observable, t)
    <outdir>/manifest.json    job list + exact and emulator reference values

Run it here (qiskit env), copy <outdir> to the machine that has pyQuil + your
QCS credentials (Rigetti's hosted QCS JupyterLab is easiest -- pyQuil, quilc and
your creds are already there), then run run_qcs.py against the manifest.

Usage:
    module load qiskit/aer-gpu
    PYTHONPATH=. python3 -m qlb_port.export_quil                 # + emulator refs
    PYTHONPATH=. python3 -m qlb_port.export_quil --no-emulator   # exact refs only
"""

import argparse
import json
import os

import numpy as np
from qiskit import transpile

from . import run_hardware_zitterbewegung as H

# universal basis whose every gate has a 1:1 Quil equivalent; quilc re-optimizes
# to the Cepheus-1 native set {rx, rz, iswap}, so this choice is cosmetic.
_BRIDGE_BASIS = ["rx", "ry", "rz", "cx"]


def circuit_to_quil(qc):
    """Translate a qiskit circuit to a generic-Quil program string.

    Returns (quil_text, measured_qubits) where measured_qubits[k] is the qiskit
    qubit read into ro[k].  Raises if an unexpected gate survives transpilation.
    """
    t = transpile(qc, basis_gates=_BRIDGE_BASIS, optimization_level=3)
    gate_lines, meas = [], []                       # meas: (clbit, qubit) pairs
    for instr in t.data:
        name = instr.operation.name
        qubits = [t.find_bit(q).index for q in instr.qubits]
        if name in ("barrier", "id", "delay"):
            continue
        if name == "measure":
            clbit = t.find_bit(instr.clbits[0]).index
            meas.append((clbit, qubits[0]))
            continue
        if name in ("rx", "ry", "rz"):
            theta = float(instr.operation.params[0])
            gate_lines.append(f"{name.upper()}({theta:.17g}) {qubits[0]}")
        elif name == "cx":
            gate_lines.append(f"CNOT {qubits[0]} {qubits[1]}")
        else:
            raise ValueError(f"unexpected gate after transpile to {_BRIDGE_BASIS}: {name!r}")

    meas.sort()                                     # order ro[] by classical bit
    lines = [f"DECLARE ro BIT[{len(meas)}]"]
    lines += gate_lines
    lines += [f"MEASURE {q} ro[{clbit}]" for clbit, q in meas]
    return "\n".join(lines) + "\n", [q for _, q in meas]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--npos", type=int, default=1, help="position qubits (1 => 2 sites).")
    ap.add_argument("--mass", type=float, default=0.8)
    ap.add_argument("--k0", type=float, default=np.pi / 2,
                    help="target carrier (snapped to nearest lattice momentum).")
    ap.add_argument("--sigma", type=float, default=1.0, help="density-packet width.")
    ap.add_argument("--tmax", type=int, default=4, help="scan t=0..tmax.")
    ap.add_argument("--shots", type=int, default=2000,
                    help="emulator shots for the reference column (QPU shots set in run_qcs.py).")
    ap.add_argument("--outdir", default="qlb_port/qcs_quil")
    ap.add_argument("--no-emulator", action="store_true",
                    help="skip the local Aer reference column (exact only).")
    args = ap.parse_args()

    npos, m, tmax, shots = args.npos, args.mass, args.tmax, args.shots
    os.makedirs(args.outdir, exist_ok=True)
    k0 = H.nearest_lattice_k(npos, args.k0)
    xgrid = np.arange(2 ** npos)
    ts = list(range(tmax + 1))

    # velocity (trembling) state and density (sloshing) packet -- identical to the
    # run_2site_* drivers, so the QPU reproduces those exact/emulator numbers.
    psi_m, E = H.mode_state(npos, k0, m)
    psi_d, _ = H.packet_state(npos, k0, m, args.sigma)

    traj_v = H.classical_evolution(psi_m, npos, m, tmax)
    traj_d = H.classical_evolution(psi_d, npos, m, tmax)

    jobs = []
    for t in ts:
        # -- velocity: <alpha_x>(t), measures qubit 1 --
        qc = H.alpha_x_circuit(npos, psi_m, m, t)
        quil, measured = circuit_to_quil(qc)
        assert measured == [1], f"alpha circuit measured {measured}, expected [1]"
        fname = f"alpha_t{t}.quil"
        with open(os.path.join(args.outdir, fname), "w") as fh:
            fh.write(quil)
        rec = {"name": f"alpha_t{t}", "kind": "alpha", "t": t, "quil": fname,
               "npos": npos, "ro_bits": 1,
               "exact": float(H.alpha_x_expectation(traj_v[t], npos))}
        if not args.no_emulator:
            rec["emulator"] = float(H.alpha_x_from_counts(H.run_aer(qc, shots), shots))
        jobs.append(rec)

        # -- density: rho(x,t) and <x>(t), measures every qubit --
        qc = H.density_circuit(npos, psi_d, m, t)
        quil, measured = circuit_to_quil(qc)
        assert measured == list(range(2 + npos)), \
            f"density circuit measured {measured}, expected {list(range(2 + npos))}"
        fname = f"dens_t{t}.quil"
        with open(os.path.join(args.outdir, fname), "w") as fh:
            fh.write(quil)
        rho = H.density(traj_d[t], npos)
        rec = {"name": f"dens_t{t}", "kind": "density", "t": t, "quil": fname,
               "npos": npos, "ro_bits": 2 + npos,
               "exact_rho": [float(v) for v in rho],
               "exact": float(rho @ xgrid)}
        if not args.no_emulator:
            rho_e = H.density_from_counts(H.run_aer(qc, shots), npos, shots)
            rec["emulator_rho"] = [float(v) for v in rho_e]
            rec["emulator"] = float(rho_e @ xgrid)
        jobs.append(rec)

    manifest = {
        "description": "2-site Dirac-QLB circuits as generic Quil for Rigetti QCS Cepheus-1.",
        "npos": npos, "n_qubits": 2 + npos, "mass": m, "k0": k0, "sigma": args.sigma,
        "E": float(E), "tmax": tmax, "emulator_shots": shots,
        "alpha": "<alpha_x> = P(ro[0]=1) - P(ro[0]=0)  (single spinor bit)",
        "density": "rho(x) from i = sum_k ro[k]*2^k, x = i >> 2;  <x> = sum_x x*rho(x)",
        "jobs": jobs,
    }
    with open(os.path.join(args.outdir, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)

    # console summary
    print(f"exported {len(jobs)} Quil programs to {args.outdir}/  "
          f"(npos={npos}, N={2 ** npos} sites, {2 + npos} qubits, k0={k0:.4f}, "
          f"m~={m}, E={E:.4f})")
    hdr = "  job          t   exact" + ("" if args.no_emulator else "    emulator")
    print(hdr)
    for j in jobs:
        line = f"  {j['name']:11s} {j['t']:2d}   {j['exact']:+7.4f}"
        if not args.no_emulator:
            line += f"   {j['emulator']:+7.4f}"
        print(line)
    print(f"\nmanifest: {os.path.join(args.outdir, 'manifest.json')}")
    print("next: copy this folder to a pyQuil + QCS environment and run run_qcs.py")


if __name__ == "__main__":
    main()
