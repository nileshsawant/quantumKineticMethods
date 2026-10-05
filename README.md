# quantumKineticMethods

**Exact quantum circuits for the quantum lattice Boltzmann (QLB) method for the Dirac equation.**

[![Documentation](https://img.shields.io/badge/docs-online-blue)](https://nileshsawant.github.io/quantumKineticMethods/)

This repository ports the three-dimensional Succi–Dellar Dirac QLB scheme, operation by
operation, to **exact quantum circuits** on qubits, and verifies on a state-vector emulator
that the circuits reproduce the classical QLB solver to machine precision. It accompanies the
paper *Exact quantum circuits for lattice Boltzmann realization of the Dirac
equation* ([arXiv:2608.06570](https://doi.org/10.48550/arXiv.2608.06570)).

It also contains the code and raw data of a follow-up study (manuscript in preparation):
cheaper **arithmetic streaming** circuits, runs of a four-site instance on a superconducting
processor (Rigetti Cepheus-1), a systematic test of seven **error-mitigation** techniques, and
small **error-detection and error-correction** experiments on the same device
([below](#running-on-quantum-hardware)).

A QLB time step is a fixed sequence of **unitary** operations — a basis rotation, a collision,
a streaming shift, and the inverse rotation — so it maps naturally onto a quantum circuit. Here
that mapping is made explicit and checked, layer by layer, against a classical reference.

![3D reflecting box](qlb_port/validation_bc.png)

*An oblique massless Dirac packet launched from a corner of a 32³ box bounces off all three
reflecting walls and returns. The ported circuit (cyan contours) matches the classical solver
(filled colour) to max|Δρ| = 6.7×10⁻¹⁷, with probability conserved.*

## What it does

- **Amplitude encoding.** The four Dirac spinor components live in 2 shared "spinor" qubits;
  each lattice axis of $N$ sites is addressed by $\log_2 N$ "position" qubits. A $32^3$ lattice
  ($32{,}768$ sites) is 17 qubits.
- **Every operation as a gate circuit.** The fixed rotations and the two-qubit collision (via
  the KAK decomposition); streaming as a **controlled increment** of the position register; the
  position-dependent potential as a **phase oracle** (massless) or a **position-multiplexed
  collision** (massive); and periodic and reflecting (bounce-back) boundaries as unitary
  circuits.
- **Composition.** Single-axis sweeps are tensored — one position register per axis on a shared
  spinor register — into 1D, 2D and 3D time steps.
- **A general harness.** A small "unitary → verified circuit" routine compiles any target
  unitary to a chosen gate set and independently checks it against its classical matrix.

No claim of computational advantage is made. The contribution is exact representability,
together with measured gate counts and a validated, reusable set of circuit primitives.

## Two ideas worth a look

- **Streaming is $+1$ on the position register.** Moving every amplitude one site is the map
  $|x\rangle \mapsto |x+1\rangle$, i.e. binary "add one" — a ripple of multi-controlled-X gates
  on the address bits, applied to all sites at once. Two adder-based alternatives cut its
  cost: a ripple-carry increment on $n_{\mathrm{pos}}-1$ clean ancilla
  (`ancilla_streaming_circuit(axis, n_pos)`) and an ancilla-free Fourier (Draper) increment
  (`streaming_circuit(..., method="fourier")`). At $n_{\mathrm{pos}}=6$ they reduce one streaming
  step from 2203 two-qubit gates to 154 and 69
  (`python -m qlb_port.validate_ancilla_streaming`).
- **Bounce-back is $+1$ on a *folded* ring.** Folding the direction qubit into the position
  register as its most significant bit turns a hard-wall reflection into a single plain
  increment: interior movers advance one site, and a mover that reaches a wall crosses the fold
  and returns with its direction reversed. One increment does advection and reflection at once.

## Validation

Circuit vs. classical solver, maximum density deviation over all sites and recorded times
(state fidelity is 1 to twelve digits throughout):

| Test | Lattice (qubits) | Physics | max\|Δρ\| |
|---|---|---|---|
| 1D free / barrier / phase | $2^6$ line (8) | massless & massive | $3.7\times10^{-12}$ |
| 2D oblique Klein | $2^5\times2^5$ (12) | massless barrier | $4.7\times10^{-16}$ |
| 3D diagonal mover | $2^4\times2^4\times2^4$ (14) | massless free | $1.0\times10^{-17}$ |
| 3D reflecting box | $2^5\times2^5\times2^5$ (17) | massless, bounce-back | $6.7\times10^{-17}$ |

<p align="center">
  <img src="qlb_port/validation_overlay.png" width="49%" />
  <img src="qlb_port/validation_2d.png" width="49%" />
</p>

*Left: 1D — Klein transmission of a massless packet through a barrier, the phase the barrier
imprints, massive-particle Zitterbewegung, and massive-barrier reflection. Right: 2D — a
massless packet at oblique incidence on a barrier splits into reflected and transmitted lobes.
Classical solver and circuit coincide in every panel.*

## Zitterbewegung and the trapped-ion comparison

Gerritsma *et al.* ([*Nature* **463**, 68 (2010)](https://doi.org/10.1038/nature08688)) simulated
the 1+1D Dirac equation on a single trapped ion and observed **Zitterbewegung** — the trembling
motion of a relativistic wave packet — together with the crossover from relativistic to
non-relativistic behaviour. That **analog** experiment is reproduced here **digitally**: the same
circuit construction advances the 1+1D Dirac dynamics, and because the ported circuit equals the
classical solver exactly, the two are run side by side.

On a 256-site line (10 qubits) at fixed carrier momentum, the mass is swept so that $mc^2$ crosses
$pc$ through the crossover. Two observables are read off the evolving packet:

- the **trembling frequency** $\omega_{\mathrm{ZB}} = 2E$, which follows the exact lattice
  dispersion and rises with the mass — a massless packet does not tremble at all, since $\alpha_x$
  is then conserved;
- the **position amplitude** $R_{\mathrm{ZB}} = \hat{b}/(2E\sin E)$, which grows from zero, peaks
  near $\tilde{m} \approx 0.7$, and dies away in both the ultrarelativistic and non-relativistic
  limits — the same crossover Gerritsma *et al.* report.

Across the whole sweep the ported circuit reproduces the classical solver to
$\max|\Delta\rho| \le 2.3\times10^{-13}$ at unit state fidelity.

<p align="center">
  <img src="qlb_port/validation_zitterbewegung.png" width="85%" />
</p>

*Left: the mean position of a massless packet (flat) and three massive packets, which tremble and
then relax to a drift as their ±E branches separate. Right: the measured trembling frequency
(filled points) tracking 2E, and the position amplitude (open squares) tracking the exact
b̂/(2E·sin E) curve, rising to a maximum near m̃ ≈ 0.7 and falling off in both limits.*

Regenerate it with `python -m qlb_port.plot_validation_zitterbewegung`.

## Running on quantum hardware

The smallest instances were run on the Rigetti Cepheus-1-108Q superconducting processor,
through the Open Quantum cloud platform (Qiskit) and directly through Rigetti QCS (pyQuil).
Every circuit was first checked on a noiseless emulator, and the raw counts of every job are
stored in `qlb_port/`.

- **Two- and four-site Zitterbewegung.** A two-site instance (3 qubits) and a four-site ring
  with Fourier streaming (4 qubits) recover the trembling velocity over a full period; on four
  sites the raw hardware keeps 70 to 80% of the ideal amplitude.
- **Qubit selection.** Pinning the circuit to qubits chosen from the live calibration removes
  the readout failures of automatic placement. `qlb_port/select_qubits.py` ranks all SWAP-free
  placements of any circuit by the calibrated fidelities they use.
- **Error mitigation.** Readout calibration, symmetry verification, Pauli twirling, zero-noise
  extrapolation, mirror rescaling, dynamical decoupling, and Clifford data regression, each
  applied on its own to the same four-site circuits. Twirling followed by Clifford data
  regression brings all four steps of the period to within about one standard deviation of the
  ideal values, reducing the mean error from 0.19 to 0.05.
- **Error detection and correction.** Bit- and phase-flip repetition codes on up to 17 qubits,
  error detection with the [[4,2,2]] code (logical Bell-state fidelity 0.90 to 0.98), and
  parity checks inside the QLB circuits.

Step-by-step commands for both platforms, all workflows, and lessons from the runs are in
[`qlb_port/HARDWARE.md`](qlb_port/HARDWARE.md). Vendor-neutral guidance for AI coding agents
(qubit pinning, pre-flight checks, mitigation, pitfalls) is packaged as the agent skill
[`.github/skills/quantum-hardware-execution`](.github/skills/quantum-hardware-execution/SKILL.md).

## Repository layout

```
qlb_port/                 circuit library
  operators.py            Dirac matrices, rotations, collision (Majorana form)
  streaming.py            streaming as a controlled increment (MCX, ripple-carry with ancilla,
                          Fourier); bounce-back (folded ring)
  potential.py            phase oracle / position-multiplexed collision
  sweep.py                one single-axis QLB substep (rotate – collide – stream – rotate)
  twod.py, threed.py      2D / 3D time steps from tensored position registers
  port.py                 "unitary -> verified circuit" harness (compile + check)
  backend.py              Qiskit Aer state-vector emulator interface
  test_port.py            every operator and composition checked vs its classical matrix
  validate_ancilla_streaming.py  streaming constructions: correctness and gate counts
  plot_validation*.py     regenerate the validation figures
  validation_*.png        the figures shown above

  HARDWARE.md             how to run on hardware (Open Quantum and Rigetti QCS)
  run_2site_*.py, run_hardware_zitterbewegung.py   Qiskit drivers for the Open Quantum runs
  export_quil.py          circuits -> Quil programs with exact reference values
  select_qubits.py        rank qubit placements from live calibration (any vendor)
  qcs_quil/               QCS runner (run_qcs.py) and generators for twirling, ZNE, mirrors,
                          dynamical decoupling, CDR, repetition and [[4,2,2]] codes
  export_parity_checks.py, export_check_diagnostics.py   in-circuit parity checks
  qcs_quil_*/             programs, manifests and raw results of every QCS batch
  hw_*                    raw counts, data and figures of the Open Quantum runs
  plot_4site_*.py, plot_*_qcs.py, fit_phase_model.py   analyses of the hardware data
  mitigation_examples.py  numbers of the worked examples in the paper
.github/skills/           agent skill: running circuits on noisy quantum hardware
dirac_qlb_solver.py       classical QLB Dirac solver (the reference the circuits are checked against)
test_qlb_validation.py    physics validation of the solver against the exact Dirac equation
test_dirac_qlb_solver.py  solver unit tests
```

## Quick start

```bash
pip install numpy scipy matplotlib qiskit qiskit-aer   # use qiskit-aer-gpu for the large circuits
```

```bash
# circuits reproduce the classical scheme (fast, runs on CPU)
python -m qlb_port.test_port

# physics validation of the classical solver vs the exact Dirac equation
python test_qlb_validation.py

# regenerate the 1D / 2D / 3D validation figures
python -m qlb_port.plot_validation
python -m qlb_port.plot_validation_2d
python -m qlb_port.plot_validation_3d

# reproduce the Zitterbewegung / trapped-ion comparison
python -m qlb_port.plot_validation_zitterbewegung

# ripple-carry and Fourier streaming: correctness and gate counts
python -m qlb_port.validate_ancilla_streaming

# qubit-selection search checked against brute force
python -m qlb_port.test_select_qubits
```

The larger circuits (for example the 17-qubit reflecting box) are emulated fastest on a GPU
through `qiskit-aer-gpu`. The hardware workflows additionally need `pyquil` (with `quilc`),
`stim`, and `pymatching`; see [`qlb_port/HARDWARE.md`](qlb_port/HARDWARE.md).

## The scheme, in one paragraph

The state is a four-component Dirac spinor on a lattice. Following Dellar (2011), the scheme
uses the **Majorana form** of the Dirac equation, in which the three matrices multiplying the
spatial derivatives are real, so each spatial gradient becomes a $\pm1$ lattice shift after a
fixed spinor rotation. A time step is split into per-axis substeps; each rotates into the
characteristic frame, applies the $\mathrm{SU}(2)$ collision that carries the mass and the
potential, streams by one site, and rotates back. Every factor is unitary, so the discrete
$\ell_2$ norm of the field is conserved exactly — which is precisely what a quantum circuit on
$n$ qubits implements.

## Reference

This code accompanies the paper below. Please cite it if you use this repository.

Sawant, N. et al. (2026) “Exact quantum circuits for lattice Boltzmann realization of the Dirac equation.” arXiv. Available at: https://doi.org/10.48550/ARXIV.2608.06570.

The hardware, error-mitigation, and error-correction results are described in: Sawant, N.,
Young, E., Griffin, K. P., and Martin, M., “Arithmetic streaming and error mitigation for
simulating the Dirac equation through a quantum lattice Boltzmann realization” (in preparation).

## License

This project is released under the [Apache License 2.0](LICENSE).
