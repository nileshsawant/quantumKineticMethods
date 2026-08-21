"""
Run a (reduced-size) QLB Dirac / Zitterbewegung circuit on real quantum hardware.
================================================================================

The full Zitterbewegung example (``plot_validation_zitterbewegung.py``) is a
10-qubit, 90-step evolution: well over a million two-qubit gates.  That is far
beyond the coherence of any current device and is an emulator-only result.  This
script instead runs a *drastically reduced* circuit whose only purpose is to show
that the QLB Dirac step itself *executes on real hardware* and that the measured
output tracks the noiseless emulator.  It talks to the Open Quantum cloud
(Rigetti Cepheus-1) through the qiskit plugin, and by default runs nothing on the
cloud at all -- it stays local on the Aer emulator until you pass ``--submit``.

Credentials are read ONLY from environment variables and are never stored here::

    export OPENQUANTUM_CLIENT_ID="s_..."        # in YOUR shell, not in this file
    export OPENQUANTUM_CLIENT_SECRET="..."

See https://docs.openquantum.com/core-sdk/authentication/ .

Differences from the emulator script (why that one cannot run on hardware):
  * it uses ``set_statevector`` / ``save_statevector`` (simulator-only); here the
    state is prepared with a real ``StatePreparation`` gate and read out with real
    ``measure`` operations and finite shots;
  * ``<alpha_x>`` is off-diagonal, so it is measured by rotating the spinor into
    the alpha_x eigenbasis (append R^-1 on the spinor qubits) and reading qubit 1:
    ``<alpha_x> = P(q1=1) - P(q1=0)``.

Usage (load the module first: ``module load qiskit/aer-gpu``)::

    # local emulator only, no cloud, no credits spent (start here):
    python -m qlb_port.run_hardware_zitterbewegung --observable alphax --npos 2 --tmax 6
    python -m qlb_port.run_hardware_zitterbewegung --observable density --npos 2 --steps 1

    # once credentials are exported, a cheap auth + plumbing check (Bell state):
    python -m qlb_port.run_hardware_zitterbewegung --bell --submit --shots 200

    # list the backends your account can see:
    python -m qlb_port.run_hardware_zitterbewegung --list-backends

    # the reduced QLB demo on hardware (spends credits; keep it small):
    python -m qlb_port.run_hardware_zitterbewegung --observable density --npos 2 --steps 1 \
        --submit --shots 2000
"""

import argparse
import os
import sys

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import StatePreparation, UnitaryGate
from qiskit_aer import AerSimulator

from . import operators as ops
from . import sweep

# ---- Dirac-frame constants (x-sweep), shared with the emulator script ----------
_R = ops.ROTATIONS["x"]
_RINV = _R.conj().T
_SIGNS = ops.streaming_signs("x").astype(float)
ALPHA_X = ops.ALPHA_X

DEFAULT_BACKEND = "rigetti:cepheus-1-108q"
# native two-qubit gate names we might see after transpiling to a real target
_TWO_Q = ("cx", "cz", "ecr", "iswap", "xx_plus_yy", "rzz", "rxx")


# ================================================================================
# physics helpers  (small-lattice versions of the emulator script's builders)
# ================================================================================
def bloch_symbol(k, m):
    """4x4 one-step operator at wavenumber k: U(k) = R . diag(e^{-i k s_c}) . Qchar . R^-1."""
    S = np.diag(np.exp(-1j * k * _SIGNS))
    return _R @ S @ ops.collision_operator_char("x", m, 0.0) @ _RINV


def energy_eigenmodes(k0, m):
    """The +E / -E eigenmodes of U(k0) that are partners in one collision block
    (largest velocity coupling |<+|alpha_x|->|), plus the eigenphase E."""
    w, V = np.linalg.eig(bloch_symbol(k0, m))
    ph = np.angle(w)
    order = np.argsort(ph)
    E = 0.5 * (ph[order[-1]] - ph[order[0]])
    lo, hi = order[:2], order[2:]
    best = (-1.0, int(lo[0]), int(hi[0]))
    for i in lo:
        for j in hi:
            c = abs(np.vdot(V[:, i], ALPHA_X @ V[:, j]))
            if c > best[0]:
                best = (c, int(i), int(j))
    return V[:, best[1]], V[:, best[2]], E


def nearest_lattice_k(n_pos, k_target):
    """Nearest lattice momentum 2*pi*j/N to k_target (so a plane wave is periodic)."""
    N = 2 ** n_pos
    ks = 2 * np.pi * np.arange(N) / N
    d = np.abs(((ks - k_target + np.pi) % (2 * np.pi)) - np.pi)
    return float(ks[int(np.argmin(d))])


def mode_state(n_pos, k0, m, mix=(1.0, 1.0)):
    """Plane wave (at a lattice momentum) x (+/-E eigenmode mix): the clean ZB state.
    Under one QLB step this stays a single Bloch mode, so <alpha_x>(t) trembles at 2E
    independent of lattice size (the size only limits which k0 are available)."""
    up, un, E = energy_eigenmodes(k0, m)
    u0 = mix[0] * up + mix[1] * un
    u0 = u0 / np.linalg.norm(u0)
    N = 2 ** n_pos
    env = np.exp(1j * k0 * np.arange(N))
    psi = (env[:, None] * u0[None, :]).reshape(-1)
    return psi / np.linalg.norm(psi), E


def packet_state(n_pos, k0, m, sigma, mix=(1.0, 1.0)):
    """Localized Gaussian x plane wave x (+/-E eigenmode mix): a moving/spreading blob,
    used for the position-density demo (mirrors the emulator script, tiny lattice)."""
    up, un, E = energy_eigenmodes(k0, m)
    u0 = mix[0] * up + mix[1] * un
    u0 = u0 / np.linalg.norm(u0)
    N = 2 ** n_pos
    x = np.arange(N)
    env = np.exp(-((x - N / 2) ** 2) / (2 * sigma ** 2)) * np.exp(1j * k0 * x)
    psi = (env[:, None] * u0[None, :]).reshape(-1)
    return psi / np.linalg.norm(psi), E


def classical_evolution(psi0, n_pos, m, tmax):
    """Exact reference: {t: U^t psi0} using the classical one-step operator."""
    U = sweep.sweep_operator("x", n_pos, m_tilde=m)
    out, p = {}, psi0.copy()
    for t in range(tmax + 1):
        out[t] = p
        p = U @ p
    return out


def alpha_x_expectation(psi, n_pos):
    """<alpha_x> for a flat statevector (layout i = x*4 + c)."""
    field = psi.reshape(2 ** n_pos, 4)
    return float(np.real(np.einsum("xi,ij,xj->", field.conj(), ALPHA_X, field)))


def density(psi, n_pos):
    """Position density rho(x) = sum_c |psi(x,c)|^2 (layout i = x*4 + c)."""
    field = psi.reshape(2 ** n_pos, 4)
    return (np.abs(field) ** 2).sum(1)


# ================================================================================
# circuit builders  (real state-prep + real measurement, hardware-ready)
# ================================================================================
def _prep(n_pos, psi0):
    """Circuit that prepares psi0 with a StatePreparation gate (no reset -> hardware-ok)."""
    n = 2 + n_pos
    qc = QuantumCircuit(n, name="prep")
    qc.append(StatePreparation(psi0), range(n))
    return qc


def alpha_x_circuit(n_pos, psi0, m, t):
    """Prepare psi0, apply t QLB steps, rotate spinor into the alpha_x eigenbasis,
    and measure qubit 1 (its sign = the alpha_x eigenvalue)."""
    qc = _prep(n_pos, psi0)
    step = sweep.sweep_circuit("x", n_pos, m_tilde=m)
    for _ in range(t):
        qc.compose(step, inplace=True)
    qc.append(UnitaryGate(_RINV, label="Rinv_meas"), [0, 1])
    meas = QuantumCircuit(2 + n_pos, 1)
    meas.compose(qc, inplace=True)
    meas.measure(1, 0)
    return meas


def density_circuit(n_pos, psi0, m, t):
    """Prepare psi0, apply t QLB steps, and measure every qubit (position + spinor)."""
    qc = _prep(n_pos, psi0)
    step = sweep.sweep_circuit("x", n_pos, m_tilde=m)
    for _ in range(t):
        qc.compose(step, inplace=True)
    n = 2 + n_pos
    meas = QuantumCircuit(n, n)
    meas.compose(qc, inplace=True)
    meas.measure(range(n), range(n))
    return meas


def bell_circuit():
    qc = QuantumCircuit(2, 2)
    qc.h(0)
    qc.cx(0, 1)
    qc.measure([0, 1], [0, 1])
    return qc


# ================================================================================
# execution helpers
# ================================================================================
def _counts_to_ints(counts):
    """Normalize qiskit counts keys (hex '0x3', binary '11', or ints) to integer -> count."""
    out = {}
    for k, v in counts.items():
        if isinstance(k, str):
            i = int(k, 16) if k.startswith("0x") else int(k.replace(" ", ""), 2)
        else:
            i = int(k)
        out[i] = out.get(i, 0) + int(v)
    return out


def run_aer(qc, shots):
    """Sample the circuit on the local Aer emulator (shot-based, like hardware)."""
    sim = AerSimulator()
    tqc = transpile(qc, sim)
    return _counts_to_ints(sim.run(tqc, shots=shots).result().get_counts())


def get_service():
    """Build an OpenQuantumService from env-var credentials (never hard-coded)."""
    cid = os.environ.get("OPENQUANTUM_CLIENT_ID")
    sec = os.environ.get("OPENQUANTUM_CLIENT_SECRET")
    if not cid or not sec:
        sys.exit("ERROR: set OPENQUANTUM_CLIENT_ID and OPENQUANTUM_CLIENT_SECRET in your "
                 "shell (export ...); never put secrets in this file.")
    from openquantum_sdk.auth import ClientCredentials
    from openquantum_sdk_qiskit import OpenQuantumService
    return OpenQuantumService(creds=ClientCredentials(client_id=cid, client_secret=sec))


def run_hardware(service, backend_name, qc, shots):
    """Transpile to the real backend target, report cost, submit, return integer counts."""
    backend = service.return_backend(backend_name)
    tqc = transpile(qc, backend, optimization_level=3)
    n2 = sum(v for g, v in tqc.count_ops().items() if g in _TWO_Q)
    print(f"    [{backend_name}] transpiled: qubits={tqc.num_qubits} depth={tqc.depth()} "
          f"2q-gates={n2}")
    job = backend.run(tqc, shots=shots)
    return _counts_to_ints(job.result().get_counts())


# ================================================================================
# observable post-processing
# ================================================================================
def alpha_x_from_counts(counts_int, shots):
    n1 = counts_int.get(1, 0)
    n0 = counts_int.get(0, 0)
    tot = n0 + n1
    return (n1 - n0) / tot if tot else float("nan")


def density_from_counts(counts_int, n_pos, shots):
    N = 2 ** n_pos
    rho = np.zeros(N)
    for i, c in counts_int.items():
        rho[i >> 2] += c            # x = i >> 2 (spinor = 2 low bits)
    return rho / max(rho.sum(), 1)


def _tvd(p, q):
    """Total-variation distance between two probability vectors."""
    return 0.5 * float(np.abs(np.asarray(p) - np.asarray(q)).sum())


# ================================================================================
# drivers
# ================================================================================
def run_alphax(args, service):
    k0 = nearest_lattice_k(args.npos, args.k0)
    psi0, E = mode_state(args.npos, k0, args.mass)
    exact = classical_evolution(psi0, args.npos, args.mass, args.tmax)
    period = (2 * np.pi / (2 * E)) if E > 1e-9 else float("inf")
    print(f"\nZitterbewegung <alpha_x>(t):  npos={args.npos} (N={2**args.npos} sites, "
          f"{2+args.npos} qubits)  k0={k0:.4f}  m~={args.mass}")
    print(f"  E={E:.4f}  ->  omega_ZB=2E={2*E:.4f} rad/step  (trembling period "
          f"{period:.2f} steps)")
    if args.submit and period < 2.5:
        print("  NOTE: period <= ~2 steps at this lattice size => the trembling is "
              "time-aliased; expect only a sign alternation, not a smooth curve.")
    header = "   t   exact     aer" + ("     hardware" if args.submit else "")
    print(header)
    for t in range(args.tmax + 1):
        ex = alpha_x_expectation(exact[t], args.npos)
        qc = alpha_x_circuit(args.npos, psi0, args.mass, t)
        av = alpha_x_from_counts(run_aer(qc, args.shots), args.shots)
        line = f"  {t:2d}  {ex:+.3f}   {av:+.3f}"
        if args.submit:
            hv = alpha_x_from_counts(
                run_hardware(service, args.backend_name, qc, args.shots), args.shots)
            line += f"   {hv:+.3f}"
        print(line, flush=True)


def run_density(args, service):
    k0 = nearest_lattice_k(args.npos, args.k0)
    psi0, E = packet_state(args.npos, k0, args.mass, args.sigma)
    exact = classical_evolution(psi0, args.npos, args.mass, args.steps)[args.steps]
    rho_exact = density(exact, args.npos)
    qc = density_circuit(args.npos, psi0, args.mass, args.steps)
    rho_aer = density_from_counts(run_aer(qc, args.shots), args.npos, args.shots)
    print(f"\nPosition density rho(x) after {args.steps} QLB step(s):  npos={args.npos} "
          f"(N={2**args.npos} sites, {2+args.npos} qubits)  k0={k0:.4f}  m~={args.mass}")
    rho_hw = None
    if args.submit:
        rho_hw = density_from_counts(
            run_hardware(service, args.backend_name, qc, args.shots), args.npos, args.shots)
    print("   x    exact     aer" + ("    hardware" if args.submit else ""))
    for x in range(2 ** args.npos):
        line = f"  {x:2d}   {rho_exact[x]:.3f}   {rho_aer[x]:.3f}"
        if rho_hw is not None:
            line += f"    {rho_hw[x]:.3f}"
        print(line)
    print(f"  total-variation distance  exact vs aer: {_tvd(rho_exact, rho_aer):.3f}"
          + (f"   exact vs hardware: {_tvd(rho_exact, rho_hw):.3f}" if rho_hw is not None
             else ""))


def run_bell(args, service):
    qc = bell_circuit()
    print("\nBell-state plumbing check (should be ~50/50 on '00' and '11'):")
    aer = run_aer(qc, args.shots)
    print("  aer     :", {format(i, "02b"): c for i, c in sorted(aer.items())})
    if args.submit:
        hw = run_hardware(service, args.backend_name, qc, args.shots)
        print("  hardware:", {format(i, "02b"): c for i, c in sorted(hw.items())})


# ================================================================================
# figures  (labeled circuit drawing + emulator results plot, hardware-overlay ready)
# ================================================================================
def draw_circuit_figure(n_pos, mass, k0, out_stage, out_native, steps=1):
    """Save two drawings of the reduced density circuit: a labeled stage-level view
    (state prep | one QLB step as Rinv, Q, stream, R | measure) and the decomposed
    native-gate view that actually runs on the device."""
    import matplotlib
    matplotlib.use("Agg")
    from . import streaming as st
    psi0, _ = packet_state(n_pos, k0, mass, 1.0)
    n = 2 + n_pos
    R = ops.ROTATIONS["x"]
    Rinv = R.conj().T
    Qc = ops.collision_operator_char("x", mass, 0.0)

    sp = StatePreparation(psi0)
    sp.label = "state prep"
    disp = QuantumCircuit(n, n)
    disp.append(sp, range(n))
    for _ in range(steps):
        disp.append(UnitaryGate(Rinv, label="Rinv"), [0, 1])
        disp.append(UnitaryGate(Qc, label="Q"), [0, 1])
        stream = st.streaming_circuit("x", n_pos).to_gate()
        stream.label = "stream"
        disp.append(stream, range(n))
        disp.append(UnitaryGate(R, label="R"), [0, 1])
    disp.measure(range(n), range(n))
    fig = disp.draw("mpl", fold=-1)
    fig.savefig(out_stage, dpi=150, bbox_inches="tight")

    tqc = transpile(density_circuit(n_pos, psi0, mass, steps),
                    basis_gates=["cz", "rz", "sx", "x"], optimization_level=3)
    text = str(tqc.draw(output="text", fold=120))
    with open(out_native, "w") as fh:
        fh.write(text)
    return disp, tqc, text


def make_results_figure(n_pos, mass, k0, shots, tmax, out_png, out_npz,
                        rho_hw=None, av_hw=None, t_hw=None):
    """Results for the reduced instance, saved as a figure and as data. Hardware points
    (rho_hw; av_hw at steps t_hw) are overlaid on the same axes when provided."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    N = 2 ** n_pos
    psi_d, E = packet_state(n_pos, k0, mass, 1.0)
    rho_exact = density(classical_evolution(psi_d, n_pos, mass, 1)[1], n_pos)
    rho_aer = density_from_counts(run_aer(density_circuit(n_pos, psi_d, mass, 1), shots),
                                  n_pos, shots)
    psi_m, E = mode_state(n_pos, k0, mass)
    ev = classical_evolution(psi_m, n_pos, mass, tmax)
    ts = np.arange(tmax + 1)
    av_exact = np.array([alpha_x_expectation(ev[t], n_pos) for t in ts])
    av_aer = np.array([alpha_x_from_counts(
        run_aer(alpha_x_circuit(n_pos, psi_m, mass, int(t)), shots), shots) for t in ts])
    save = dict(x=np.arange(N), rho_exact=rho_exact, rho_aer=rho_aer, t=ts,
                av_exact=av_exact, av_aer=av_aer, shots=shots, n_pos=n_pos, mass=mass,
                k0=k0, E=float(E))
    if rho_hw is not None:
        save["rho_hw"] = rho_hw
    if av_hw is not None:
        save["av_hw"] = av_hw
        save["t_hw"] = np.asarray(t_hw)
    np.savez(out_npz, **save)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.2))
    xs = np.arange(N)
    if rho_hw is None:
        w = 0.38
        axA.bar(xs - w / 2, rho_exact, w, color="0.65", label="exact")
        axA.bar(xs + w / 2, rho_aer, w, color="C0", label=f"emulator ({shots} shots)")
    else:
        w = 0.28
        axA.bar(xs - w, rho_exact, w, color="0.65", label="exact")
        axA.bar(xs, rho_aer, w, color="C0", label=f"emulator ({shots} shots)")
        axA.bar(xs + w, rho_hw, w, color="C3", label="hardware")
    axA.set_xticks(xs)
    axA.set_xlabel("site  $x$")
    axA.set_ylabel(r"$\rho(x)$ after one step")
    axA.set_title(f"Density, {N}-site line ({2 + n_pos} qubits)")
    axA.legend(frameon=False)

    axB.plot(ts, av_exact, "-", color="0.4", lw=1.5, label="exact")
    axB.plot(ts, av_aer, "o", color="C0", mfc="none", ms=7,
             label=f"emulator ({shots} shots)")
    if av_hw is not None:
        axB.plot(np.asarray(t_hw), av_hw, "s", color="C3", ms=7, label="hardware")
    axB.axhline(0, color="k", lw=0.6, ls=":")
    axB.set_xlabel("step  $t$")
    axB.set_ylabel(r"$\langle\alpha_x(t)\rangle$")
    axB.set_ylim(-1.05, 1.05)
    period = 2 * np.pi / (2 * E) if E > 1e-9 else float("inf")
    axB.set_title(rf"Velocity trembling  ($2E={2 * E:.2f}$, period ${period:.1f}$)")
    axB.legend(frameon=False)
    tag = "with hardware" if rho_hw is not None else "hardware points to be overlaid"
    fig.suptitle(f"Reduced QLB Dirac instance: exact vs emulator ({tag})", y=1.02)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    if rho_hw is not None:
        print(f"\n  density total-variation distance   exact vs emulator: "
              f"{_tvd(rho_exact, rho_aer):.3f}   exact vs hardware: "
              f"{_tvd(rho_exact, rho_hw):.3f}")
    return rho_exact, rho_aer, av_exact, av_aer


def _dump_counts(fh, label, counts_int, nbits, shots):
    """Append raw measurement counts (as bitstrings) to an open text file, flush to disk."""
    fh.write(f"\n[{label}] shots={shots} nbits={nbits}\n")
    for i in sorted(counts_int):
        fh.write(f"  {format(i, '0%db' % nbits)}  {counts_int[i]}\n")
    fh.flush()
    os.fsync(fh.fileno())


def collect_hardware(service, backend_name, n_pos, mass, k0, shots, steps, hw_tmax=None,
                     raw_path="qlb_port/hw_demo_raw.txt"):
    """Run the reduced circuits on the cloud backend and write the raw counts to a text
    file immediately (before any post-processing), so the numbers survive even if the
    later plotting fails. Returns (rho_hw, av_hw, t_hw)."""
    import datetime
    nqd = 2 + n_pos
    with open(raw_path, "w") as fh:
        fh.write(f"# raw hardware counts  backend={backend_name}  "
                 f"{datetime.datetime.now().isoformat(timespec='seconds')}\n")
        fh.write(f"# n_pos={n_pos} mass={mass} k0={k0:.6f} shots={shots} steps={steps}\n")
        fh.write("# density qubit layout: x = high 2 bits, spinor c = low 2 bits\n")

        psi_d, _ = packet_state(n_pos, k0, mass, 1.0)
        print("  density circuit:")
        counts_d = run_hardware(service, backend_name,
                                density_circuit(n_pos, psi_d, mass, steps), shots)
        _dump_counts(fh, f"density_after_{steps}_step", counts_d, nqd, shots)
        rho_hw = density_from_counts(counts_d, n_pos, shots)
        fh.write("  rho_hw(x) = " + " ".join(f"{v:.4f}" for v in rho_hw) + "\n")
        fh.flush()

        av_hw = t_hw = None
        if hw_tmax is not None:
            psi_m, _ = mode_state(n_pos, k0, mass)
            t_hw = np.arange(hw_tmax + 1)
            av_hw = np.empty(len(t_hw))
            for i, t in enumerate(t_hw):
                print(f"  alpha_x circuit t={int(t)}:")
                counts_a = run_hardware(service, backend_name,
                                        alpha_x_circuit(n_pos, psi_m, mass, int(t)), shots)
                _dump_counts(fh, f"alpha_x_t={int(t)}_measure_q1", counts_a, 1, shots)
                av_hw[i] = alpha_x_from_counts(counts_a, shots)
                fh.write(f"  <alpha_x>(t={int(t)}) = {av_hw[i]:+.4f}\n")
                fh.flush()
    print(f"  raw hardware counts written to {raw_path}")
    return rho_hw, av_hw, t_hw


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--observable", choices=["alphax", "density"], default="alphax",
                    help="alphax: <alpha_x>(t) trembling; density: rho(x) after --steps.")
    ap.add_argument("--npos", type=int, default=2,
                    help="position qubits (lattice N=2**npos). Keep <=2 for hardware.")
    ap.add_argument("--k0", type=float, default=np.pi / 2,
                    help="target carrier wavenumber (snapped to nearest lattice momentum).")
    ap.add_argument("--mass", type=float, default=0.8, help="dimensionless mass m~.")
    ap.add_argument("--sigma", type=float, default=1.0,
                    help="Gaussian width in sites (density observable only).")
    ap.add_argument("--tmax", type=int, default=6, help="max step for the alphax curve.")
    ap.add_argument("--steps", type=int, default=1, help="number of steps for density.")
    ap.add_argument("--shots", type=int, default=2000)
    ap.add_argument("--bell", action="store_true", help="run only the Bell plumbing check.")
    ap.add_argument("--submit", action="store_true",
                    help="ACTUALLY submit to the cloud backend (spends credits). Off by "
                         "default: everything runs locally on Aer.")
    ap.add_argument("--list-backends", action="store_true",
                    help="list backends visible to your account, then exit.")
    ap.add_argument("--backend-name", default=DEFAULT_BACKEND)
    ap.add_argument("--figures", action="store_true",
                    help="draw the reduced circuit and make the emulator results plot "
                         "(local only), then exit.")
    ap.add_argument("--hw-alphax", action="store_true",
                    help="with --figures --submit, also run the <alpha_x>(t) curve on "
                         "hardware (more circuits, more credits).")
    ap.add_argument("--hw-tmax", type=int, default=3,
                    help="max step for the hardware <alpha_x> curve (with --hw-alphax).")
    args = ap.parse_args()

    service = None
    if args.list_backends:
        service = get_service()
        print("Backends visible to your account:")
        for b in service.backends():
            print("  -", getattr(b, "name", b))
        service.close()
        return

    if args.figures:
        k0 = nearest_lattice_k(args.npos, args.k0)
        _, tqc, text = draw_circuit_figure(
            args.npos, args.mass, k0, "qlb_port/hw_demo_circuit.png",
            "qlb_port/hw_demo_circuit_native.txt", steps=args.steps)
        print("Reduced density circuit, decomposed to native gates "
              f"(cz/rz/sx/x), depth {tqc.depth()}:\n")
        print(text)
        rho_hw = av_hw = t_hw = None
        if args.submit:
            service = get_service()
            try:
                print(f"\nSubmitting reduced circuits to {args.backend_name} "
                      f"({args.shots} shots each) ...")
                rho_hw, av_hw, t_hw = collect_hardware(
                    service, args.backend_name, args.npos, args.mass, k0, args.shots,
                    args.steps, hw_tmax=(args.hw_tmax if args.hw_alphax else None))
            finally:
                service.close()
        make_results_figure(args.npos, args.mass, k0, args.shots, args.tmax,
                            "qlb_port/hw_demo_results.png", "qlb_port/hw_demo_data.npz",
                            rho_hw=rho_hw, av_hw=av_hw, t_hw=t_hw)
        print("\nwrote qlb_port/hw_demo_circuit.png (labeled), "
              "hw_demo_circuit_native.txt, hw_demo_results.png, hw_demo_data.npz")
        return

    if args.submit:
        service = get_service()

    try:
        if args.bell:
            run_bell(args, service)
        elif args.observable == "alphax":
            run_alphax(args, service)
        else:
            run_density(args, service)
    finally:
        if service is not None:
            service.close()

    if not args.submit:
        print("\n(local emulator only; add --submit to run on the cloud backend once "
              "OPENQUANTUM_CLIENT_ID / OPENQUANTUM_CLIENT_SECRET are exported.)")


if __name__ == "__main__":
    main()
