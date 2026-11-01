"""
Port verified QLB circuits from Qiskit to Guppy and run them on Quantinuum's Selene emulator.

A Qiskit circuit is transpiled to {Rz, Ry, Rx, CX}, rebuilt gate by gate as a pytket circuit
(pytket angles are in half-turns), loaded into Guppy with `guppy.load_pytket`, compiled for the
Helios platform and sampled on Selene's statevector simulator. Outcome frequencies are compared with
the exact probabilities of the same Qiskit circuit.

Run on a compute node with the qiskit module (guppylang installed in its venv):
    module load qiskit/aer-gpu
    PYTHONPATH=. python3 -m qlb_port.guppy_port [--shots 4000]
"""
import argparse
import time

import numpy as np
from pytket import Circuit
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import StatePreparation, UnitaryGate
from qiskit.quantum_info import Statevector

from guppylang import guppy
from guppylang.std.builtins import comptime, result
from guppylang.std.quantum import array, collect_measurements, measure_array, qubit

from . import operators as ops
from . import sweep, threed

BASIS = ["rz", "ry", "rx", "cx"]
_TK = {"rz": "Rz", "ry": "Ry", "rx": "Rx"}


def to_pytket(qc):
    """Unitary Qiskit circuit (no measurements) -> pytket circuit on the same qubit indices."""
    t = transpile(qc, basis_gates=BASIS, optimization_level=3, seed_transpiler=7)
    tk = Circuit(t.num_qubits)
    for inst in t.data:
        name = inst.operation.name
        qs = [t.find_bit(q).index for q in inst.qubits]
        if name == "cx":
            tk.CX(*qs)
        elif name in _TK:
            getattr(tk, _TK[name])(float(inst.operation.params[0]) / np.pi, qs[0])
        elif name != "barrier":
            raise ValueError(f"unexpected gate {name}")
    return tk


def run_selene(tk, shots, seed=1):
    """Measure every qubit after `tk`; return outcome counts indexed little-endian like Qiskit."""
    n = tk.n_qubits
    body = guppy.load_pytket("qlb_body", tk)

    @guppy
    def main() -> None:
        q = array(qubit() for _ in range(comptime(n)))
        body(q)
        result("m", collect_measurements(measure_array(q)))

    res = main.emulator(n_qubits=n).statevector_sim().with_shots(shots).with_seed(seed).run()
    counts = np.zeros(2 ** n)
    for bits, c in res.register_counts()["m"].items():
        counts[sum(1 << i for i, b in enumerate(bits) if b == "1")] += c
    return counts


def gaussian(n_pos, k0=0.5):
    N = 2 ** n_pos
    x = np.arange(N)
    env = np.exp(-(x - N / 4) ** 2 / (2 * (N / 8) ** 2) + 1j * k0 * x)
    return env / np.linalg.norm(env)


def instance(dim, n_pos, steps, m_tilde, velocity):
    """Packet preparation + `steps` QLB steps (+ velocity readout rotation) as one Qiskit circuit."""
    w, v = np.linalg.eigh(ops.ALPHA_X)
    spinor = v[:, np.argmax(w)]  # <alpha_x> = +1 at t=0, so a massive packet trembles
    step = (sweep.sweep_circuit("x", n_pos, m_tilde=m_tilde, streaming_method="fourier") if dim == 1
            else threed.sweep3d_circuit(n_pos, n_pos, n_pos, m_tilde=m_tilde, streaming_method="fourier"))
    qc = QuantumCircuit(step.num_qubits)
    qc.append(StatePreparation(spinor), [0, 1])
    for a in range(dim):
        qc.append(StatePreparation(gaussian(n_pos)), list(range(2 + a * n_pos, 2 + (a + 1) * n_pos)))
    for _ in range(steps):
        qc.compose(step, inplace=True)
    if velocity:
        # rotates alpha_x to diag(-1,-1,1,1): <alpha_x> = P(q1=1) - P(q1=0)
        qc.append(UnitaryGate(np.linalg.inv(ops.X_ROTATION)), [0, 1])
    return qc


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--shots", type=int, default=4000)
    args = ap.parse_args()
    cases = [  # (label, dim, n_pos, steps, m_tilde, velocity readout)
        ("1D 4 sites, velocity, t=1", 1, 2, 1, 0.3, True),
        ("1D 16 sites, density, t=4", 1, 4, 4, 0.3, False),
        ("1D 16 sites, velocity, t=2", 1, 4, 2, 0.3, True),
        ("1D 16 sites, velocity, t=4", 1, 4, 4, 0.3, True),
        ("1D 16 sites, velocity, t=8", 1, 4, 8, 0.3, True),
        ("3D 4^3 sites, density, t=2", 3, 2, 2, 0.3, False),
        ("3D 8^3 sites, density, t=2", 3, 3, 2, 0.3, False),
    ]
    print(f"{'case':30s} {'qubits':>6s} {'CX':>5s} {'TVD':>7s} {'TVD noise':>9s} "
          f"{'<a_x> exact':>11s} {'<a_x> Selene':>12s} {'time s':>7s}")
    for label, dim, n_pos, steps, m_tilde, vel in cases:
        qc = instance(dim, n_pos, steps, m_tilde, vel)
        p = Statevector(qc).probabilities()
        tk = to_pytket(qc)
        t0 = time.time()
        counts = run_selene(tk, args.shots)
        dt = time.time() - t0
        f = counts / counts.sum()
        tvd = 0.5 * np.abs(f - p).sum()
        # mean TVD of `shots` exact samples, as the shot-noise reference
        rng = np.random.default_rng(0)
        noise = np.mean([0.5 * np.abs(rng.multinomial(args.shots, p) / args.shots - p).sum() for _ in range(50)])
        ncx = sum(1 for c in tk.get_commands() if c.op.type.name == "CX")
        if vel:
            q1 = (np.arange(len(p)) >> 1) & 1
            ax_exact, ax_sel = (p * (2 * q1 - 1)).sum(), (f * (2 * q1 - 1)).sum()
            vtxt = f"{ax_exact:+11.3f} {ax_sel:+12.3f}"
        else:
            vtxt = f"{'':11s} {'':12s}"
        print(f"{label:30s} {qc.num_qubits:6d} {ncx:5d} {tvd:7.3f} {noise:9.3f} {vtxt} {dt:7.1f}", flush=True)


if __name__ == "__main__":
    main()
