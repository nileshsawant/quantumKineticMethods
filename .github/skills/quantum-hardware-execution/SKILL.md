---
name: quantum-hardware-execution
description: "Use when running quantum circuits on real quantum hardware (QPU) of any vendor (Rigetti, IBM, IQM, IonQ, Quantinuum, ...) or planning such runs: qubit selection and pinning from live calibration data, layout/routing, transpilation and two-qubit gate-count optimization, pre-flight checks (noiseless simulation, compile-only, noise models), cost-safe job submission, readout error mitigation (confusion matrix), symmetry verification / post-selection, Pauli twirling / randomized compiling, zero-noise extrapolation, mirror-circuit rescaling, dynamical decoupling, Clifford data regression, shot noise and bootstrap uncertainties, run-to-run drift, and small error-detection/correction experiments (repetition code, [[4,2,2]] code, parity checks)."
---

# Running circuits on noisy quantum hardware

A vendor-neutral procedure for getting trustworthy expectation values out of a noisy quantum processor.
It is distilled from a complete study (circuit optimization, qubit pinning, seven mitigation techniques,
and small error-correction experiments) on a superconducting device. Links point to working examples in
the [quantumKineticMethods](https://github.com/nileshsawant/quantumKineticMethods) repository (pinned
to commit `ee60b02`); they use Qiskit for circuit construction and pyQuil/Quil for one vendor, but every
step has an equivalent in other SDKs.

## Ground rules

1. **Know the exact answer first.** For every circuit, compute the ideal observable by state-vector
   simulation before touching hardware. Without it you cannot tell a good run from a bad one.
2. **Change one thing at a time** and compare against a **raw baseline taken in the same session, on the
   same qubits**. Device noise drifts on the scale of minutes to hours.
3. **Never choose the estimator after looking at the exact result.** Fix the procedure, then run it.
4. **Every number comes from a script** that reads the saved raw data. Save raw counts (and job IDs)
   to disk the moment each job returns, before any post-processing.
5. **Spend hardware time last and carefully.** Validate everything that is free (simulation,
   compilation) first; submit a single job, inspect it, then the batch.
6. **Credentials never pass through the agent.** The user places tokens/keys in their own files or
   environment; never echo, read, or commit them. Do not tunnel around network blocks.

## Step 1: Minimize the circuit

Two-qubit gates dominate the error budget (typically 5-10x worse than one-qubit gates), so the
two-qubit count and depth are the primary cost metrics.

- Transpile at the highest optimization level to a universal basis (e.g. `{rz, ry, rx, cx}`) and record
  two-qubit count and depth. Compare alternative constructions of the same operation (e.g. arithmetic
  via ripple-carry vs Fourier-basis adders); asymptotics can mislead at small sizes.
- After every transformation, verify the unitary (or output state) against the reference to machine
  precision. Example: [`port_unitary`](https://github.com/nileshsawant/quantumKineticMethods/blob/ee60b02183c77d3cb9e1faaa31fe66484390e7c0/qlb_port/port.py#L40),
  [`two_qubit_count`](https://github.com/nileshsawant/quantumKineticMethods/blob/ee60b02183c77d3cb9e1faaa31fe66484390e7c0/qlb_port/port.py#L61),
  [`verify`](https://github.com/nileshsawant/quantumKineticMethods/blob/ee60b02183c77d3cb9e1faaa31fe66484390e7c0/qlb_port/port.py#L67).
- Translating between SDKs: go through a universal basis with identical rotation conventions
  ($R_P(\theta)=e^{-i\theta P/2}$), drop only the global phase, and re-simulate the translated program.
  Example: [`circuit_to_quil`](https://github.com/nileshsawant/quantumKineticMethods/blob/ee60b02183c77d3cb9e1faaa31fe66484390e7c0/qlb_port/export_quil.py#L54).
- Merge cancelling layers across logical blocks (e.g. a basis rotation at the end of one step and its
  inverse at the start of the next).

## Step 2: Choose and pin the physical qubits

Automatic placement is usually blind to readout quality and can land on a faulty qubit. Choose the
qubits yourself from the **live** calibration.

1. Fetch the calibration: per-qubit readout fidelity, per-edge two-qubit fidelity, per-qubit one-qubit
   (randomized-benchmarking) fidelity. (Qiskit: `backend.target` / `backend.properties()`; Rigetti QCS:
   instruction-set architecture, see [HARDWARE.md](https://github.com/nileshsawant/quantumKineticMethods/blob/ee60b02183c77d3cb9e1faaa31fe66484390e7c0/qlb_port/HARDWARE.md).)
2. Build the circuit's **interaction graph** (which qubit pairs share a two-qubit gate). If it is a
   path or tree, it may embed into the coupling map without any SWAP.
3. Enumerate embeddings of that graph into the coupling map (subgraph isomorphism; for chains, walk the
   graph) and score each by the product of the fidelities it actually uses, weighted by gate counts:
   two-qubit edges used, readout of measured qubits, one-qubit fidelity of all qubits. The assignment
   of logical to physical qubits matters, so different orientations of the same chain score differently.
   Ready-made tool: [`select_qubits.py`](https://github.com/nileshsawant/quantumKineticMethods/blob/9f92b5bd032e134376fd47d56c0f7cd45931c202/qlb_port/select_qubits.py)
   ([`select_layouts`](https://github.com/nileshsawant/quantumKineticMethods/blob/9f92b5bd032e134376fd47d56c0f7cd45931c202/qlb_port/select_qubits.py#L224),
   branch-and-bound over SWAP-free embeddings, checked against brute force in
   [`test_select_qubits.py`](https://github.com/nileshsawant/quantumKineticMethods/blob/9f92b5bd032e134376fd47d56c0f7cd45931c202/qlb_port/test_select_qubits.py)).
   It reads a generic calibration JSON, a Rigetti QCS ISA
   ([`calibration_from_qcs_isa`](https://github.com/nileshsawant/quantumKineticMethods/blob/9f92b5bd032e134376fd47d56c0f7cd45931c202/qlb_port/select_qubits.py#L82)),
   or a Qiskit `Target`
   ([`calibration_from_qiskit_target`](https://github.com/nileshsawant/quantumKineticMethods/blob/9f92b5bd032e134376fd47d56c0f7cd45931c202/qlb_port/select_qubits.py#L114)),
   and takes the circuit as an edge list, Quil, OpenQASM, or a Qiskit circuit
   ([`interaction_from_qiskit`](https://github.com/nileshsawant/quantumKineticMethods/blob/9f92b5bd032e134376fd47d56c0f7cd45931c202/qlb_port/select_qubits.py#L174)).
   Example: `python3 qlb_port/select_qubits.py --qcs-isa isa.json --edges 0-1,1-3,2-3 --top 5`.
   If no SWAP-free layout exists, reduce the interaction graph (different construction) or accept routing.
4. Pin the best embedding (Qiskit: `initial_layout=`; Quil: `PRAGMA INITIAL_REWIRING "NAIVE"` with
   remapped qubit indices, see [`_pin_text`](https://github.com/nileshsawant/quantumKineticMethods/blob/ee60b02183c77d3cb9e1faaa31fe66484390e7c0/qlb_port/qcs_quil/run_qcs.py#L85)).
5. Check the compiled program: no SWAPs, two-qubit count equal to the all-to-all count.
6. Optional cheap pre-check: a state-preparation-only circuit that reads every qubit of the layout.
   A qubit collapsing to $0$ (energy relaxation before readout) shows up immediately.

## Step 3: Pre-flight checks (free) and submission

- Run the **compiled native** program on a noiseless simulator with the same pinning; it must
  reproduce the exact values within shot noise. This catches translation, pinning, and readout-mapping
  bugs.
- Compile for the real device (compile-only) to catch ISA/connectivity errors.
- Simulate with a simple noise model (depolarizing per gate, thermal relaxation, readout error) to
  predict whether a technique can help at all. Real devices are usually worse and more correlated than
  such models; use them for direction, not magnitude.
- **Test vendor features with one real job.** Local simulators often accept features the hardware
  path does not: mid-circuit measurement, delays, barriers/fences, verbatim blocks, compiler bypass.
  (Observed: programs with mid-circuit measurement compiled and ran but executed only one shot; native
  programs submitted without the vendor compiler also returned one shot.)
- Submit one job, verify the number of shots and bits per shot, then submit the batch. Record job IDs
  at submission ([example](https://github.com/nileshsawant/quantumKineticMethods/blob/ee60b02183c77d3cb9e1faaa31fe66484390e7c0/qlb_port/qcs_quil/run_qcs.py#L221)).

## Step 4: Statistics and reproducibility

- **Shot noise** of a $\pm1$ observable $v=P(1)-P(0)$ from $N$ shots: $\sigma=\sqrt{(1-v^2)/N}$
  (binomial). $N=2000$ gives $\sigma\approx0.022$ at $v=0$.
- For derived estimates, use a **bootstrap**: redraw counts (parametric: binomial with the measured
  probability) and recompute; when a model is fitted on training data, resample the training set too.
- Report deviations in units of the standard deviation; use $\chi^2$ against the number of degrees
  of freedom when fitting models.
- **Measure run-to-run variation.** Repeat key circuits, interleaved, minutes and hours apart. If the
  spread exceeds shot noise, single-run differences of that size must not be interpreted.
- **Sensitivity analysis.** Where the ideal signal is steep (e.g. near zero crossings of an
  oscillation), a small coherent phase error shifts the value strongly; at extrema it barely matters.
  Fit damping/offset/phase models across batches to identify the error type:
  [`fit_phase_model.fit`](https://github.com/nileshsawant/quantumKineticMethods/blob/ee60b02183c77d3cb9e1faaa31fe66484390e7c0/qlb_port/fit_phase_model.py#L55).
  Opposite-sign residuals at points of opposite slope indicate a phase error; same-sign residuals an
  offset.

## Step 5: Error mitigation ladder

Apply in this order; each later stage assumes the earlier ones. All techniques improve *estimates of
expectation values*, not the quantum state, and all cost extra shots.

| # | Technique | Removes | Assumes / fails when | Example |
|---|---|---|---|---|
| 1 | Readout calibration | readout errors, incl. correlated ones | stable readout; full matrix costs $2^n$ circuits | [`confusion`](https://github.com/nileshsawant/quantumKineticMethods/blob/ee60b02183c77d3cb9e1faaa31fe66484390e7c0/qlb_port/plot_4site_qcs_readout.py#L37), [`mitigate`](https://github.com/nileshsawant/quantumKineticMethods/blob/ee60b02183c77d3cb9e1faaa31fe66484390e7c0/qlb_port/plot_4site_qcs_readout.py#L56) |
| 2 | Symmetry verification | errors that change a conserved quantity | a conserved quantity exists; misses errors that commute with it | [`observables`](https://github.com/nileshsawant/quantumKineticMethods/blob/ee60b02183c77d3cb9e1faaa31fe66484390e7c0/qlb_port/plot_4site_qcs_parity.py#L47) |
| 3 | Pauli twirling (randomized compiling) | coherence of errors (turns drift/sign flips into reproducible damping) | Clifford two-qubit gates; Paulis inserted after compilation | [`twirl`](https://github.com/nileshsawant/quantumKineticMethods/blob/ee60b02183c77d3cb9e1faaa31fe66484390e7c0/qlb_port/qcs_quil/twirl_qcs.py#L45), [`preserve`](https://github.com/nileshsawant/quantumKineticMethods/blob/ee60b02183c77d3cb9e1faaa31fe66484390e7c0/qlb_port/qcs_quil/twirl_qcs.py#L61) |
| 4 | Zero-noise extrapolation | damping | pure, smooth damping; fails with offsets or sign changes | [`fold`](https://github.com/nileshsawant/quantumKineticMethods/blob/ee60b02183c77d3cb9e1faaa31fe66484390e7c0/qlb_port/qcs_quil/twirl_qcs.py#L78), [`estimators`](https://github.com/nileshsawant/quantumKineticMethods/blob/ee60b02183c77d3cb9e1faaa31fe66484390e7c0/qlb_port/plot_4site_qcs_zne.py#L44) |
| 5 | Mirror rescaling | damping | global depolarizing noise; overshoot past the physical range shows it is not | [`mirror`](https://github.com/nileshsawant/quantumKineticMethods/blob/ee60b02183c77d3cb9e1faaa31fe66484390e7c0/qlb_port/qcs_quil/twirl_qcs.py#L100) |
| 6 | Dynamical decoupling | slow phase drift on idle qubits | timing control; idle errors dominate | [`build`](https://github.com/nileshsawant/quantumKineticMethods/blob/ee60b02183c77d3cb9e1faaa31fe66484390e7c0/qlb_port/qcs_quil/dd_qcs.py#L56) |
| 7 | Clifford data regression | damping **and** offset | training circuits see the same noise as the target | [`cdr_qcs`](https://github.com/nileshsawant/quantumKineticMethods/blob/ee60b02183c77d3cb9e1faaa31fe66484390e7c0/qlb_port/qcs_quil/cdr_qcs.py#L87), [`analysis`](https://github.com/nileshsawant/quantumKineticMethods/blob/ee60b02183c77d3cb9e1faaa31fe66484390e7c0/qlb_port/plot_4site_qcs_cdr.py#L92) |

Technique notes:

1. **Readout.** Prepare each basis state ($X$ on the 1-bits), measure, histogram: column $b$ of the
   confusion matrix $A$. Solve $p=Aq$ by non-negative least squares and renormalize (plain inversion can
   give negative probabilities). Per-qubit (tensored) calibration misses correlated readout errors
   (e.g. one qubit misread when its neighbor is excited); inspect the full matrix for $n\lesssim6$.
   For a single qubit: $v=(v_{\rm meas}-e_{01}+e_{10})/(1-e_{01}-e_{10})$. Run calibration in the same batch.
2. **Symmetry verification.** Find an observable $S$ conserved by every step. An error $E$ is detected
   only if it anticommutes with $S$. If the natural initial state mixes symmetry sectors, prepare a
   state in one sector with the same observable dynamics (verify by exact simulation). Discarding shots
   inflates shot noise by $1/\sqrt{\text{kept fraction}}$. Checking the symmetry mid-circuit through
   ancillas (deferred measurement) only helps if the check gates are much better than the circuit;
   diagnose check-gate errors in isolation
   ([`cz_repetition_programs`](https://github.com/nileshsawant/quantumKineticMethods/blob/ee60b02183c77d3cb9e1faaa31fe66484390e7c0/qlb_port/export_check_diagnostics.py#L39)).
3. **Twirling.** Around each Clifford two-qubit gate $G$ insert random Paulis $P$ before and
   $GPG^\dagger$ after (for CZ: $X\otimes I\to X\otimes Z$, $Y\otimes I\to Y\otimes Z$, $Z$ passes). Compile
   first, insert into the native program, and protect it from the compiler (verbatim/preserve block,
   or disable optimization), then confirm the executed gate sequence. Use 16+ variants at equal total
   shots. **Diagnostic:** if the spread between variants exceeds their shot noise, coherent errors are
   present. Twirled values are reproducible but damped toward zero.
4. **ZNE.** Fold self-inverse gates ($G\to G^{s}$, $s=1,3,5$) or $G\to GG^\dagger G$; twirl each copy.
   Linear: $v_0=v_1-(v_3-v_1)/2$. Exponential: $v_s=v_0e^{-\kappa s}$, undefined if the $v_s$ change
   sign and meaningless if the last point is near zero. Dividing by a damping factor multiplies any
   offset by the same factor.
5. **Rescaling.** Run the circuit followed by its inverse (mirror); ideal result known. Under global
   depolarizing noise the mirror decays as $\lambda^2$, so $v\approx v_{\rm meas}/\sqrt{m}$. Estimates
   beyond the physical range mean the mirror is more sensitive (it involves all qubits) than the
   observable.
6. **Dynamical decoupling.** Pairs of $X$ pulses on idle qubits cancel quasi-static $Z$ rotations. If
   the SDK only offers barriers/fences, build CZ-by-CZ layers and run three versions (as compiled,
   layered without pulses, layered with pulses) to separate the effect of scheduling from that of the
   echoes. A large layering effect means the error depends on the gate schedule (crosstalk), not idling.
7. **CDR.** Replace non-Clifford gates of the target by Clifford ones with probability ~0.7 to make
   training circuits whose ideal values are computable (exactly for few qubits, by stabilizer
   simulation for many). Fit $v_{\rm ideal}=a+b\,v_{\rm meas}$ on raw values (absorbs readout too).
   On hardware where $R_Z$ is a virtual frame change, changing only $R_Z$ angles keeps the physical pulse
   sequence identical, which makes the training noise match the target
   ([`is_clifford`](https://github.com/nileshsawant/quantumKineticMethods/blob/ee60b02183c77d3cb9e1faaa31fe66484390e7c0/qlb_port/qcs_quil/cdr_qcs.py#L83),
   [`ideal_value`](https://github.com/nileshsawant/quantumKineticMethods/blob/ee60b02183c77d3cb9e1faaa31fe66484390e7c0/qlb_port/qcs_quil/cdr_qcs.py#L59)).
   Twirl training circuits and target alike. Check that training ideal values span the target range,
   and test the whole procedure on a density-matrix noise model first
   ([`noise_check`](https://github.com/nileshsawant/quantumKineticMethods/blob/ee60b02183c77d3cb9e1faaa31fe66484390e7c0/qlb_port/plot_4site_qcs_cdr.py#L37)).
   Training-circuit scatter beyond shot noise reveals circuit-dependent (coherent) noise.

Recommended sequence for a new circuit family: **pin from calibration -> readout calibration -> twirl
-> CDR** (or rescaling if offsets are negligible), with symmetry post-selection where a symmetry exists.

## Step 6: Small error-detection/correction experiments

- **Capability check first:** repeated syndrome extraction needs mid-circuit measurement and reset.
  Without it, measure everything at the end: the final data qubits already contain the last-round
  syndrome, so decoding reduces to a majority vote (confirm with a matching decoder).
- **Repetition codes** (bit-flip and phase-flip, both logical states):
  [`repcode ops`](https://github.com/nileshsawant/quantumKineticMethods/blob/ee60b02183c77d3cb9e1faaa31fe66484390e7c0/qlb_port/qcs_quil/repcode_qcs.py#L36),
  with a Stim twin for noise-model predictions
  ([`to_stim`](https://github.com/nileshsawant/quantumKineticMethods/blob/ee60b02183c77d3cb9e1faaa31fe66484390e7c0/qlb_port/qcs_quil/repcode_qcs.py#L88)).
  Compare the logical error with the independent-error prediction from measured per-qubit error rates
  ([`majority_fail_independent`](https://github.com/nileshsawant/quantumKineticMethods/blob/ee60b02183c77d3cb9e1faaa31fe66484390e7c0/qlb_port/plot_repcode_qcs.py#L37));
  an excess means correlated errors. Energy relaxation makes errors asymmetric: the logical state whose
  qubits sit in $\ket{1}$ fails. The compiler's frame choice decides which state that is; compute the
  excited-state fraction of the compiled program
  ([`excited_fraction`](https://github.com/nileshsawant/quantumKineticMethods/blob/ee60b02183c77d3cb9e1faaa31fe66484390e7c0/qlb_port/plot_repcode_qcs.py#L64)).
- **[[4,2,2]] error detection** on a logical Bell state (two physical Bell pairs, fault-tolerant
  preparation); fidelity $F=\tfrac14(1+\langle XX\rangle-\langle YY\rangle+\langle ZZ\rangle)$ from three
  bases, post-selecting on the four-qubit parity:
  [`c422 ops`](https://github.com/nileshsawant/quantumKineticMethods/blob/ee60b02183c77d3cb9e1faaa31fe66484390e7c0/qlb_port/qcs_quil/c422_qcs.py#L30),
  [`corr`](https://github.com/nileshsawant/quantumKineticMethods/blob/ee60b02183c77d3cb9e1faaa31fe66484390e7c0/qlb_port/plot_c422_qcs.py#L30).
  Add controlled noise with gate pairs that cancel ($G G^\dagger$) inside a verbatim block.
- Compare like with like: detection gain on the same qubits and data, not encoded vs a different pair.

## Pitfalls observed in practice

- Automatic placement read a register bit on a qubit that always returned $0$; pinning fixed it.
- Calibrated two-qubit fidelities (~99%) understated in-context errors: a single CNOT used as a parity
  check failed 8-10% of the time. Measure the building block you actually use.
- The compiler removes gate pairs that cancel unless they are protected.
- Untwirled deep circuits gave values of opposite sign in consecutive jobs; twirling made them agree.
- Coherent phase errors on the spinor/data qubits are invisible to parity checks and post-selection
  (they commute with $Z$-type symmetries).
- Mirror rescaling overshot to $>1$ where the noise was not global depolarizing.
- Single-batch success is not reproducibility; state which results come from a single batch.

## Reporting checklist

- Device, access route, date/time of each batch, processor time used, physical qubits and their
  calibration values.
- Compiled two-qubit counts and depth; whether routing added SWAPs.
- Shots per circuit, number of variants/training circuits, uncertainties and how they were obtained.
- Raw and mitigated values side by side with the exact value; failures reported as well as successes.
- Scripts that regenerate every number from the saved raw counts, e.g.
  [`mitigation_examples.py`](https://github.com/nileshsawant/quantumKineticMethods/blob/ee60b02183c77d3cb9e1faaa31fe66484390e7c0/qlb_port/mitigation_examples.py)
  for worked examples.
