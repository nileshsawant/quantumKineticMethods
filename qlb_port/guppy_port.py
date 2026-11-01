"""
Port verified QLB circuits from Qiskit to Guppy and run them on Quantinuum's Selene emulator.

A Qiskit circuit is transpiled to {Rz, Ry, Rx, CX}, rebuilt gate by gate as a pytket circuit
(pytket angles are in half-turns), loaded into Guppy with `guppy.load_pytket`, compiled for the
Helios platform and sampled on Selene's statevector simulator. Outcome frequencies are compared with
the exact probabilities of the same Qiskit circuit.

Run on a compute node with the qiskit module (guppylang installed in its venv):
    module load qiskit/aer-gpu
    PYTHONPATH=. python3 -m qlb_port.guppy_port [--shots 4000]
    PYTHONPATH=. python3 -m qlb_port.guppy_port --noise [--shots 2000]   # Helios-like depolarizing noise
    PYTHONPATH=. python3 -m qlb_port.guppy_port --counts                 # two-qubit counts of whole runs

The --noise study converts the published Helios average infidelities (arXiv:2511.05465: two-qubit
7.9e-4, single-qubit 2.5e-5, state preparation and measurement 4.8e-4) to depolarizing probabilities
p = r (d+1)/d, assuming Selene draws a uniformly random non-identity Pauli with probability p. It has
no memory, leakage, crosstalk or coherent errors, so it is a lower bound on the noise of the device.
"""
import argparse
import json
import os
import time

import numpy as np
from pytket import Circuit
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import StatePreparation, UnitaryGate
from qiskit.quantum_info import Statevector

from guppylang import guppy
from guppylang.std.builtins import comptime, result
from guppylang.std.quantum import array, collect_measurements, measure_array, qubit
from selene_sim import DepolarizingErrorModel

from . import operators as ops
from . import sweep, threed

_TK = {"rz": "Rz", "ry": "Ry", "rx": "Rx"}
HELIOS = dict(p_2q=7.9e-4 * 5 / 4, p_1q=2.5e-5 * 3 / 2, p_meas=4.8e-4, p_init=4.8e-4)


def to_pytket(qc, two_qubit="cx"):
    """Unitary Qiskit circuit (no measurements) -> pytket circuit on the same qubit indices.
    two_qubit='cx' or 'rzz' (arbitrary-angle ZZ, the Helios native entangler)."""
    t = transpile(qc, basis_gates=["rz", "ry", "rx", two_qubit], optimization_level=3, seed_transpiler=7)
    tk = Circuit(t.num_qubits)
    for inst in t.data:
        name = inst.operation.name
        qs = [t.find_bit(q).index for q in inst.qubits]
        if name == "cx":
            tk.CX(*qs)
        elif name == "rzz":
            tk.ZZPhase(float(inst.operation.params[0]) / np.pi, *qs)
        elif name in _TK:
            getattr(tk, _TK[name])(float(inst.operation.params[0]) / np.pi, qs[0])
        elif name != "barrier":
            raise ValueError(f"unexpected gate {name}")
    return tk


def run_selene(tk, shots, seed=1, error_model=None):
    """Measure every qubit after `tk`; return outcome counts indexed little-endian like Qiskit and
    the number of native two-qubit (rzz) operations per shot after compilation for Helios."""
    n = tk.n_qubits
    body = guppy.load_pytket("qlb_body", tk)

    @guppy
    def main() -> None:
        q = array(qubit() for _ in range(comptime(n)))
        body(q)
        result("m", collect_measurements(measure_array(q)))

    em = main.emulator(n_qubits=n).statevector_sim().with_shots(shots).with_seed(seed).with_metrics()
    em = em.with_n_processes(min(16, os.cpu_count() or 1))
    if error_model is not None:
        em = em.with_error_model(error_model)
    res = em.run()
    counts = np.zeros(2 ** n)
    for bits, c in res.register_counts()["m"].items():
        counts[sum(1 << i for i, b in enumerate(bits) if b == "1")] += c
    return counts, res.metrics()[0]["user_program"]["rzz_count"]


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


def two_qubit_count(qc, two_qubit):
    t = transpile(qc, basis_gates=["rz", "ry", "rx", two_qubit], optimization_level=3, seed_transpiler=7)
    return int(t.count_ops().get(two_qubit, 0))


def counts_study(out="qlb_port/guppy_counts.json"):
    """Two-qubit counts (CX and arbitrary-angle ZZ) of one step on its own, of the packet preparation, and
    of whole runs (preparation + T steps, density readout) transpiled at once. In a whole run the
    inverse Fourier transform that ends one streaming step cancels the transform that starts the next,
    so the cost is an overhead N0 plus a marginal cost k per step, fitted from T=4 and T=8."""
    cases = [("1D, 4 sites", 1, 2, 4), ("1D, 16 sites", 1, 4, 20), ("1D, 64 sites", 1, 6, 20),
             ("3D, 4^3 sites", 3, 2, 8), ("3D, 8^3 sites", 3, 3, 8), ("3D, 16^3 sites", 3, 4, 8)]
    rows = []
    for label, dim, n_pos, t_plan in cases:
        step = (sweep.sweep_circuit("x", n_pos, m_tilde=0.3, streaming_method="fourier") if dim == 1
                else threed.sweep3d_circuit(n_pos, n_pos, n_pos, m_tilde=0.3, streaming_method="fourier"))
        r = {"label": label, "dim": dim, "n_pos": n_pos, "qubits": step.num_qubits, "steps_planned": t_plan}
        for g in ("cx", "rzz"):
            n = {T: two_qubit_count(instance(dim, n_pos, T, 0.3, False), g) for T in (0, 4, 8, t_plan)}
            k = (n[8] - n[4]) / 4
            r.update({f"{g}_single_step": two_qubit_count(step, g), f"{g}_prep": n[0], f"{g}_T4": n[4],
                      f"{g}_T8": n[8], f"{g}_per_step": k, f"{g}_overhead": n[4] - 4 * k,
                      f"{g}_total": n[t_plan], f"{g}_total_linear": n[4] - 4 * k + k * t_plan})
        rows.append(r)
        print(f"{label:16s} q={r['qubits']:2d} | single step CX {r['cx_single_step']:3d} ZZ {r['rzz_single_step']:3d} | "
              f"prep CX {r['cx_prep']:3d} ZZ {r['rzz_prep']:3d} | per step CX {r['cx_per_step']:5.1f} "
              f"ZZ {r['rzz_per_step']:5.1f} | overhead CX {r['cx_overhead']:5.0f} ZZ {r['rzz_overhead']:5.0f} | "
              f"T={t_plan}: CX {r['cx_total']:4d} ({r['cx_total_linear']:.0f}) ZZ {r['rzz_total']:4d} "
              f"({r['rzz_total_linear']:.0f}) | steps x single step ZZ {r['rzz_prep'] + t_plan * r['rzz_single_step']}",
              flush=True)
    json.dump(rows, open(out, "w"), indent=1)
    print(f"wrote {out}")


def x_marginal(prob, n_pos):
    """Marginal probability along x (register on qubits 2 .. 1+n_pos)."""
    x = (np.arange(len(prob)) >> 2) & (2 ** n_pos - 1)
    return np.bincount(x, weights=prob, minlength=2 ** n_pos)


def noise_study(shots, out="qlb_port/guppy_noise"):
    em = lambda s: DepolarizingErrorModel(random_seed=s, **HELIOS)  # noqa: E731
    print("depolarizing probabilities:", HELIOS)
    rng = np.random.default_rng(0)

    # 1D Zitterbewegung, 16 sites (6 qubits)
    ts = np.arange(0, 11)
    zb = {"t": ts, "exact": [], "noisy": [], "sd": [], "rzz": []}
    print(f"\n1D 16 sites, m=0.3: {'t':>3s} {'rzz':>5s} {'<a_x> exact':>11s} {'Selene+noise':>12s} {'sd':>6s}")
    for t in ts:
        qc = instance(1, 4, int(t), 0.3, True)
        p = Statevector(qc).probabilities()
        counts, nzz = run_selene(to_pytket(qc, "rzz"), shots, seed=int(t) + 1, error_model=em(int(t) + 11))
        q1 = 2 * ((np.arange(len(p)) >> 1) & 1) - 1
        ex, no = (p * q1).sum(), (counts * q1).sum() / counts.sum()
        sd = np.sqrt(max(1 - no ** 2, 0) / shots)
        for k, v in zip(("exact", "noisy", "sd", "rzz"), (ex, no, sd, nzz)):
            zb[k].append(v)
        print(f"{'':27s}{t:3d} {nzz:5d} {ex:+11.3f} {no:+12.3f} {sd:6.3f}", flush=True)

    # 3D packets: x-marginal of the density
    print(f"\n3D density x-marginal: {'lattice':>8s} {'t':>3s} {'qubits':>6s} {'rzz':>5s} {'TVD':>6s} "
          f"{'TVD shots':>9s} {'<x> exact':>9s} {'<x> noisy':>9s}")
    d3 = {}
    for n_pos, steps in ((3, (2, 4, 8)), (4, (2, 4, 8))):
        for t in steps:
            qc = instance(3, n_pos, t, 0.3, False)
            p = Statevector(qc).probabilities()
            counts, nzz = run_selene(to_pytket(qc, "rzz"), shots, seed=t + 1, error_model=em(t + 21))
            mx_e, mx_n = x_marginal(p, n_pos), x_marginal(counts / counts.sum(), n_pos)
            ref = np.mean([0.5 * np.abs(rng.multinomial(shots, mx_e) / shots - mx_e).sum() for _ in range(50)])
            xs = np.arange(2 ** n_pos)
            d3[f"{2 ** n_pos}_{t}"] = (mx_e, mx_n)
            print(f"{'':22s}{2 ** n_pos:5d}^3 {t:3d} {qc.num_qubits:6d} {nzz:5d} "
                  f"{0.5 * np.abs(mx_n - mx_e).sum():6.3f} {ref:9.3f} {(xs * mx_e).sum():9.2f} "
                  f"{(xs * mx_n).sum():9.2f}", flush=True)

    np.savez(out + ".npz", shots=shots, **{k: np.array(v) for k, v in zb.items()},
             **{f"m3d_{k}": np.array(v) for k, v in d3.items()})
    plot_noise(out)


def plot_noise(out="qlb_port/guppy_noise"):
    """Draw the noise-study figure from the saved .npz."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    d = np.load(out + ".npz")
    shots = int(d["shots"])
    fig, (a, b) = plt.subplots(1, 2, figsize=(9, 3.2))
    a.plot(d["t"], d["exact"], "o-", color="0.5", label="exact")
    a.errorbar(d["t"], d["noisy"], yerr=d["sd"], fmt="s", color="C3", ms=4, capsize=2,
               label=f"Selene, Helios-like noise ({shots} shots)")
    a.set_xlabel("time step $t$")
    a.set_ylabel(r"$\langle\alpha_x\rangle$")
    a.set_title("(a) 1D, 16 sites, 6 qubits", loc="left", fontsize=10)
    a.legend(frameon=False, fontsize=8)
    mx_e, mx_n = d["m3d_16_8"]
    xs = np.arange(16)
    b.plot(xs, mx_e, "-", color="0.5", label="exact")
    b.plot(xs, mx_n, "o", color="C3", mfc="none", label="Selene, Helios-like noise")
    b.set_xlabel("$x$")
    b.set_ylabel("marginal density")
    b.set_title("(b) 3D, $16^3$ sites, 14 qubits, $t=8$", loc="left", fontsize=10)
    b.set_ylim(0, 1.35 * max(mx_e.max(), mx_n.max()))
    b.legend(frameon=False, fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(out + ".png", dpi=200)
    print(f"wrote {out}.png")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--shots", type=int, default=4000)
    ap.add_argument("--noise", action="store_true", help="run the Helios-like noise study instead")
    ap.add_argument("--counts", action="store_true", help="two-qubit gate counts of whole runs")
    ap.add_argument("--plot-noise", action="store_true", help="redraw the noise figure from guppy_noise.npz")
    args = ap.parse_args()
    if args.plot_noise:
        plot_noise()
        return
    if args.counts:
        counts_study()
        return
    if args.noise:
        noise_study(args.shots)
        return
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
        counts, _ = run_selene(tk, args.shots)
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
