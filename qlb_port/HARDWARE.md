# Running the Dirac QLB circuits on quantum hardware

This guide covers everything in `qlb_port/` that touches real hardware: the two-site and
four-site Zitterbewegung runs, qubit selection, the error-mitigation study, and the
error-detection and error-correction experiments. All of it ran on the Rigetti
Cepheus-1-108Q processor through two routes:

| Route | SDK | Used for |
|---|---|---|
| **Open Quantum** (cloud broker) | Qiskit, OpenQASM 3 | two-site runs, first four-site run |
| **Rigetti QCS** (direct) | pyQuil, Quil | pinned four-site runs, all mitigation and error-correction experiments |

Vendor-neutral guidance (qubit pinning, pre-flight checks, mitigation, pitfalls) is in the
agent skill [`.github/skills/quantum-hardware-execution/SKILL.md`](../.github/skills/quantum-hardware-execution/SKILL.md).

Everything is built and checked for free first (state-vector emulator, noiseless QVM,
compile-only), and raw counts are written to disk as soon as each job returns, before any
post-processing. Credentials are never hard-coded or committed: Open Quantum reads
environment variables, QCS reads `~/.qcs/{settings,secrets}.toml`.

---

## 0. Emulator only (no cloud, no credits)

```bash
module load qiskit/aer-gpu
PYTHONPATH=. python3 -m qlb_port.run_2site_trembling      # velocity <alpha_x>(t), two sites
PYTHONPATH=. python3 -m qlb_port.run_2site_density        # position <x>(t), two sites

# four sites (npos=2) with Fourier-basis streaming; the run_2site_* drivers take any --npos
PYTHONPATH=. python3 -m qlb_port.run_2site_trembling --npos 2 --k0 0 --streaming fourier
PYTHONPATH=. python3 -m qlb_port.run_2site_density   --npos 2 --k0 0 --streaming fourier

# validate the ripple-carry and Fourier streaming circuits and count their gates
PYTHONPATH=. python3 -u -m qlb_port.validate_ancilla_streaming
```

Each driver prints an exact-vs-emulator table and writes an `.npz`.

---

## Path A -- Open Quantum (qiskit)

The `run_2site_*` and `run_hardware_zitterbewegung` drivers submit qiskit circuits
straight to the Open Quantum backend.

```bash
# credentials (set in your shell; never commit these):
export OPENQUANTUM_CLIENT_ID=...
export OPENQUANTUM_CLIENT_SECRET=...

module load qiskit/aer-gpu

# two sites (spends credits); backend = rigetti:cepheus-1-108q
PYTHONPATH=. python3 -m qlb_port.run_2site_trembling --submit --shots 2000
PYTHONPATH=. python3 -m qlb_port.run_2site_density   --submit --shots 2000

# four sites, Fourier streaming, pinned to the chain 0-1-10-11 (q0,q1,p0,p1 -> 0,1,11,10)
PYTHONPATH=. python3 -m qlb_port.run_2site_t0        --npos 2 --k0 0 --submit
PYTHONPATH=. python3 -m qlb_port.run_2site_trembling --npos 2 --k0 0 --streaming fourier --layout 0,1,11,10 --submit
PYTHONPATH=. python3 -m qlb_port.run_2site_density   --npos 2 --k0 0 --streaming fourier --layout 0,1,11,10 --submit

# generic driver: Bell warm-up, backend list, single observable
PYTHONPATH=. python3 -m qlb_port.run_hardware_zitterbewegung --bell
PYTHONPATH=. python3 -m qlb_port.run_hardware_zitterbewegung --list-backends
```

Raw device counts are written to `hw_2site_*_raw.txt` and `hw_4site_fourier_*_raw.txt`
(the `hw_4site_fourier_q19_*` files are the first, automatically placed run, which hit a
faulty qubit and is kept for reference). If a client secret ever appears in a terminal
banner, rotate it in the Open Quantum portal.

---

## Path B -- Rigetti QCS directly (pyQuil)

The circuits are built in qiskit but QCS speaks Quil, so they are exported to
generic Quil first (`quilc` on the QCS side compiles that to the native
`{RX, RZ, CZ}` set and routes it).

### B.1  Export circuits to Quil (qiskit env, no credentials)

```bash
module load qiskit/aer-gpu
PYTHONPATH=. python3 -m qlb_port.export_quil --npos 1 --outdir qlb_port/qcs_quil
```

Writes `qcs_quil/<name>.quil` plus `manifest.json` (carrying the exact and Aer
reference values). The four-site circuits of the mitigation study were exported with

```bash
PYTHONPATH=. python3 -m qlb_port.export_quil --npos 2 --k0 0 --streaming fourier --readout-cal \
    --outdir qlb_port/qcs_quil_n2_fourier
```

`--readout-cal` adds the $2^n$ readout-calibration circuits and `--parity` writes the
even-parity states used for symmetry verification.

### B.2  Easiest: run in Rigetti's hosted QCS JupyterLab

pyQuil, `quilc`, the QVM, and your credentials are all preinstalled there.
Upload the `qcs_quil` folder, open a Terminal, then:

```bash
conda activate python3          # pyQuil is in the 'python3' env, NOT base
cd qcs_quil
python run_qcs.py --bell                       # access smoke test (tiny)
python run_qcs.py --qvm                         # free noiseless rehearsal
python run_qcs.py --compile-only --qubits 103,102,100,101   # free: compile for the device
python run_qcs.py --only alpha_t4 --qubits 103,102,100,101  # ONE real job first, inspect it
python run_qcs.py --shots 2000 --qubits 103,102,100,101 --out qcs_results.json   # the batch
```

Device defaults to `Cepheus-1-108Q`; `--dir` selects the program folder. Only `qc.run` on
the real device spends credits (billed by execution time); `--qvm` and `--compile-only` are
free. `run_qcs.py` compiles every program with `quilc` as ProtoQuil (all measurements at the
end), pins it with `PRAGMA INITIAL_REWIRING "NAIVE"` when `--qubits` is given, records the job
IDs at submission, and stores the raw bit strings and the processor time of every job in
`--out`.

The hub can also be driven from another machine through the JupyterHub REST API (upload
files, run code in a kernel, download the results) with a personal API token kept in a
private file; never commit or print the token.

### B.3  Run locally (pyQuil venv + a quilc server)

```bash
python3 -m venv qcs-venv && qcs-venv/bin/pip install 'pyquil>=4.7'
mkdir -p ~/.qcs   # then place your settings.toml + secrets.toml here (chmod 600)
```

`quilc` is needed for compilation. Without Docker, run it from the official image
through Apptainer and point pyQuil at it:

```bash
module load apptainer
apptainer pull quilc.sif docker://rigetti/quilc:latest
nohup apptainer exec quilc.sif quilc -S -p 5599 &        # RPCQ server
export QCS_SETTINGS_APPLICATIONS_QUILC_URL=tcp://127.0.0.1:5599
qcs-venv/bin/python run_qcs.py --qvm                     # local validation
```

NOTE: compilation and the QVM work from anywhere, but the QPU *execution* endpoint
(Rigetti netblock `64.4.164.0/22`) may be blocked by an HPC egress firewall; if
`qc.run` on the real device times out, submit from the JupyterLab environment
(B.2) instead. If the QCS token refresh fails with "Invalid cross-device link", run
with `TMPDIR=$HOME`.

### B.4  Choosing the qubits from the live calibration

Automatic placement is blind to readout quality; in the first four-site run it read a
register bit on a qubit that always returned 0. `select_qubits.py` finds every SWAP-free
placement of a circuit's interaction graph on the device and ranks them by the calibrated
fidelities they use (two-qubit gates weighted by their counts, readout of measured qubits,
one-qubit gates). Fetching the calibration is read-only and free; it needs the pyQuil venv.

```bash
# live calibration, saved for later reuse, circuit read from a Quil program
python3 qlb_port/select_qubits.py --qcs Cepheus-1-108Q --save-calibration isa.json \
    --quil qlb_port/qcs_quil_n2_fourier/alpha_t4.quil --top 5
# saved calibration, circuit given as its interaction graph
python3 qlb_port/select_qubits.py --qcs-isa isa.json --edges 0-1,1-3,2-3 --top 5
# checks against brute force
PYTHONPATH=. python3 -m qlb_port.test_select_qubits
```

The printed `--qubits` list is in program-qubit order and can be passed to `run_qcs.py`.
The module also reads a generic calibration JSON and a Qiskit `Target`
(`calibration_from_qiskit_target`) for other vendors. The raw calibration is available as

```python
from qcs_sdk import QCSClient
from qcs_sdk.qpu.isa import get_instruction_set_architecture as gisa
isa = gisa(client=QCSClient.load(), quantum_processor_id="Cepheus-1-108Q")
# MEASURE sites carry 'fRO', CZ sites 'fCZ', benchmarks 'fRB' (1q), 'T1', 'T2'.
```

---

## C. Error-mitigation workflows (four-site circuits, QCS)

All runs used the pinned chain `--qubits 103,102,100,101` (q0, q1, p0, p1). Each generator
writes Quil programs plus a `manifest.json` into a folder; the folder is run with
`run_qcs.py --dir <folder>` as in B.2, and an analysis script reads the saved results.
The twirling, DD and CDR generators compile with `quilc` to native programs and wrap them
in `PRAGMA PRESERVE_BLOCK`, so they need the pyQuil venv and a `quilc` server (B.3).

| Technique | Generate | Data folder | Analyse |
|---|---|---|---|
| readout calibration | `export_quil.py --readout-cal` | `qcs_quil_n2_fourier/` | `plot_4site_qcs_readout.py` |
| symmetry verification | `export_quil.py --parity --readout-cal` | `qcs_quil_n2_parity/` | `plot_4site_qcs_parity.py` |
| run-to-run repeats at t=4 | `qcs_quil_n2_t4check/make_t4check.py <outdir>` | `qcs_quil_n2_t4check/` | `plot_4site_qcs_parity.py` |
| Pauli twirling | `qcs_quil/twirl_qcs.py` | `qcs_quil_n2_twirl/` | `plot_4site_qcs_twirl.py` |
| zero-noise extrapolation, mirror rescaling | `twirl_qcs.py --scales 1,3,5 --mirror --cal` | `qcs_quil_n2_zne/` | `plot_4site_qcs_zne.py` |
| dynamical decoupling | `qcs_quil/dd_qcs.py` | `qcs_quil_n2_dd/` | `plot_4site_qcs_dd.py` |
| Clifford data regression | `qcs_quil/cdr_qcs.py` | `qcs_quil_n2_cdr/` | `plot_4site_qcs_cdr.py [--noise-check]` |
| phase-error model fits | (results above) | | `fit_phase_model.py` |
| worked-example numbers | | | `mitigation_examples.py` |

Examples:

```bash
Q=103,102,100,101; SRC=qlb_port/qcs_quil_n2_fourier
python3 qlb_port/qcs_quil/twirl_qcs.py --src $SRC --jobs alpha_t4 --qubits $Q --out qlb_port/qcs_quil_n2_twirl
python3 qlb_port/qcs_quil/twirl_qcs.py --src $SRC --jobs alpha_t1,alpha_t2,alpha_t3,alpha_t4 --qubits $Q \
    --scales 1,3,5 --mirror --cal --out qlb_port/qcs_quil_n2_zne
python3 qlb_port/qcs_quil/dd_qcs.py   --src $SRC --jobs alpha_t4 --qubits $Q --out qlb_port/qcs_quil_n2_dd
python3 qlb_port/qcs_quil/cdr_qcs.py  --src $SRC --jobs alpha_t1,alpha_t2,alpha_t3,alpha_t4 --qubits $Q \
    --out qlb_port/qcs_quil_n2_cdr
PYTHONPATH=. python3 qlb_port/plot_4site_qcs_cdr.py --outdir /tmp/figs --noise-check
```

The analysis scripts write figures into the paper folder by default; pass `--outdir` to
write elsewhere.

---

## D. Error-detection and error-correction experiments (QCS)

These need `stim` and `pymatching` in the pyQuil venv for the noise-model twins and the
decoder.

| Experiment | Generate | Data folder | Analyse |
|---|---|---|---|
| bit- and phase-flip repetition codes, d = 3-9 | `qcs_quil/repcode_qcs.py --out <dir>` | `qcs_quil_repcode/` | `plot_repcode_qcs.py [--sim] [--populations]` |
| [[4,2,2]] logical Bell state | `qcs_quil/c422_qcs.py --out <dir>` | `qcs_quil_c422/` | `plot_c422_qcs.py [--sim]` |
| parity checks inside the QLB circuits | `export_parity_checks.py --outdir <dir>` | `qcs_quil_n2_pcheck/` | `plot_parity_checks_qcs.py [--noise-check]` |
| diagnostics of the checks | `export_check_diagnostics.py [--cz-repetition]` | `qcs_quil_n2_pdiag/`, `qcs_quil_n2_czrep/` | `plot_parity_checks_qcs.py --diag` |

The repetition codes ran on the 17-qubit chain
`100,101,102,93,84,83,74,65,56,57,58,67,68,77,86,95,104`, the [[4,2,2]] code on
`100,101,102,93`, and the parity checks on `65,74,83,92` with ancillas `56,64,66`.

---

## E. Lessons from the QCS runs

- **Measure only at the end.** Programs with gates after a `MEASURE` compile (non-ProtoQuil)
  but execute a single shot on the on-demand path; local compile and QVM do not catch this.
- **Always compile with `quilc`.** Native programs submitted without it also returned one
  shot. Protect hand-placed gates (twirling Paulis, echo pulses, cancelling CZ pairs) with
  `PRAGMA PRESERVE_BLOCK`, otherwise the compiler removes or merges them.
- `DELAY` and an argument-less `FENCE` are rejected; `FENCE` on an explicit qubit list works.
- **Submit one job first** (`--only`), check the shot count and bits per shot, then the batch.
- **Pin the qubits** (B.4) and re-select them from the calibration of the day.
- Results drift between jobs; repeat key circuits and compare against baselines from the same
  session.

---

## Plotting

```bash
# two-site figures, from the .npz files (either platform):
PYTHONPATH=. python3 qlb_port/plot_2site_results.py
# or straight from a QCS results file:
python3 qlb_port/plot_qcs_results.py qcs_results.json
# four-site Open Quantum figures and table values, and the gate-level circuit drawings:
PYTHONPATH=. python3 qlb_port/plot_4site_results.py --outdir /tmp/figs
PYTHONPATH=. python3 qlb_port/draw_4site_circuit.py --outdir /tmp/figs
```
