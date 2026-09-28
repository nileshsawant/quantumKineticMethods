"""
Gate-level drawings of the four-site (n_pos=2) t=1 density circuit with Fourier
streaming, as in Fig. hw-native of the two-site paper.

Writes to <outdir>:
    hw_4site_t1_generic.png   compiled to {Rz,Ry,Rx,CX}, all-to-all connectivity
    hw_4site_t1_native.png    compiled against the device target on the pinned qubits
                              (0,1,11,10); only with --device, which needs
                              OPENQUANTUM_CLIENT_ID / OPENQUANTUM_CLIENT_SECRET in the env

The submitted circuit was compiled without a fixed seed, and the device compilation varies
between 7 and 8 two-qubit gates with the seed, so the drawing uses the first seed that
reproduces the submitted statistics (8 two-qubit gates, depth 18; see
hw_4site_fourier_submit_log.txt).

Usage:
    PYTHONPATH=. python3 -m qlb_port.draw_4site_circuit [--device]
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister, transpile

from . import run_hardware_zitterbewegung as H

N_POS, MASS, K0, SIGMA, T = 2, 0.8, 0.0, 1.0, 1
LAYOUT = [0, 1, 11, 10]
SUBMITTED_2Q, SUBMITTED_DEPTH = 8, 18


def labeled(qc):
    """Same circuit on registers named q (spinor) and p (position)."""
    out = QuantumCircuit(QuantumRegister(2, "q"), QuantumRegister(N_POS, "p"),
                         ClassicalRegister(qc.num_clbits, "c"))
    return out.compose(qc, qubits=range(qc.num_qubits), clbits=range(qc.num_clbits))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default="proposal/quairk/ancillaStreaming/figures")
    ap.add_argument("--device", action="store_true")
    args = ap.parse_args()

    psi0, _ = H.packet_state(N_POS, K0, MASS, SIGMA)
    qc = labeled(H.density_circuit(N_POS, psi0, MASS, T, streaming_method="fourier"))

    tq = transpile(qc, basis_gates=["rz", "ry", "rx", "cx"], optimization_level=3,
                   seed_transpiler=0)
    print(f"generic: {dict(tq.count_ops())} depth={tq.depth()}")
    tq.draw("mpl", fold=-1).savefig(os.path.join(args.outdir, "hw_4site_t1_generic.png"),
                                    dpi=150, bbox_inches="tight")

    if args.device:
        backend = H.get_service().return_backend(H.DEFAULT_BACKEND)
        for seed in range(200):
            tn = transpile(qc, backend, optimization_level=3, initial_layout=LAYOUT,
                           seed_transpiler=seed)
            n2 = sum(1 for inst in tn.data if inst.operation.num_qubits == 2)
            if (n2, tn.depth()) == (SUBMITTED_2Q, SUBMITTED_DEPTH):
                break
        else:
            raise RuntimeError("no seed reproduces the submitted gate count and depth")
        print(f"device (seed {seed}): {dict(tn.count_ops())} depth={tn.depth()}")
        tn.draw("mpl", fold=-1, idle_wires=False).savefig(
            os.path.join(args.outdir, "hw_4site_t1_native.png"), dpi=150, bbox_inches="tight")


if __name__ == "__main__":
    main()
