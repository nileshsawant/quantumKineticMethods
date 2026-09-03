# Running the Dirac QLB circuits on quantum hardware

Two cloud paths are supported for the reduced two--site demo (and larger `npos`):
**Open Quantum** (a broker fronting Rigetti Cepheus-1) and **Rigetti Quantum Cloud
Services (QCS)** directly. Both execute the same circuits; only the submission layer
differs. Everything can be built and cross-checked on the local Aer emulator first,
with no cloud account and no credits.

Credentials are never hard-coded or committed: Open Quantum reads environment
variables, QCS reads `~/.qcs/{settings,secrets}.toml`.

---

## 0. Emulator only (no cloud, no credits)

```bash
module load qiskit/aer-gpu
PYTHONPATH=. python3 -m qlb_port.run_2site_trembling      # velocity <alpha_x>(t)
PYTHONPATH=. python3 -m qlb_port.run_2site_density        # position <x>(t)
```

Each prints an exact-vs-emulator table and writes an `.npz`. The paper figure
(Figure 13 style) is assembled from those two `.npz` files by
`qlb_port/plot_2site_results.py`.

---

## Path A -- Open Quantum (qiskit)

The `run_2site_*` and `run_hardware_zitterbewegung` drivers submit qiskit circuits
straight to the Open Quantum backend.

```bash
# credentials (set in your shell; never commit these):
export OPENQUANTUM_CLIENT_ID=...
export OPENQUANTUM_CLIENT_SECRET=...

module load qiskit/aer-gpu

# submit to hardware (spends credits); backend = rigetti:cepheus-1-108q
PYTHONPATH=. python3 -m qlb_port.run_2site_trembling --submit --shots 2000
PYTHONPATH=. python3 -m qlb_port.run_2site_density   --submit --shots 2000

# generic driver: Bell warm-up, backend list, single observable
PYTHONPATH=. python3 -m qlb_port.run_hardware_zitterbewegung --bell
PYTHONPATH=. python3 -m qlb_port.run_hardware_zitterbewegung --list-backends
```

Raw device counts are written to `hw_2site_*_raw.txt` (kept separate from earlier
runs). If a client secret ever appears in a terminal banner, rotate it in the
Open Quantum portal.

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
reference values). Use `--npos 2` for the four-site (4-qubit) case, etc.

### B.2  Easiest: run in Rigetti's hosted QCS JupyterLab

pyQuil, `quilc`, the QVM, and your credentials are all preinstalled there.
Upload the `qcs_quil` folder, open a Terminal, then:

```bash
conda activate python3          # pyQuil is in the 'python3' env, NOT base
cd qcs_quil
python run_qcs.py --bell                       # access smoke test (tiny)
python run_qcs.py --qvm                         # free noiseless rehearsal
python run_qcs.py --shots 2000                  # real run -> qcs_results.json
python run_qcs.py --shots 2000 --qubits 90,99,101   # pin high-fidelity qubits
```

Device defaults to `Cepheus-1-108Q`. Only `qc.run` on the real device spends
credits; `--qvm` and `--compile-only` are free.

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
(B.2) instead.

### Picking pinning qubits from the live calibration

The auto compiler can place the circuit on weak qubits (readout especially). Pull
the live characterization and pin the good ones:

```python
from qcs_sdk import QCSClient
from qcs_sdk.qpu.isa import get_instruction_set_architecture as gisa
isa = gisa(client=QCSClient.load(), quantum_processor_id="Cepheus-1-108Q")
# MEASURE sites carry 'fRO' (readout fidelity); CZ sites carry 'fCZ' (2q fidelity).
```

Give the position register the best-readout qubit and the spinor pair a
high-fidelity CZ edge, then pass them as `--qubits <spinor0>,<spinor1>,<position>`
(the run pins logical -> physical with NAIVE rewiring).

---

## Plotting

```bash
# paper Figure 13 style, from the .npz files (either platform):
PYTHONPATH=. python3 qlb_port/plot_2site_results.py
# or straight from a QCS results file:
python3 qlb_port/plot_qcs_results.py qcs_results.json
```
