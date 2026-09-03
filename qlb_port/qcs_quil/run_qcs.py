"""
Run the exported 2-site Dirac-QLB Quil programs on Rigetti QCS (Cepheus-1).
===========================================================================

This is the QCS-side companion to export_quil.py.  Run it *where pyQuil and your
QCS credentials live* -- easiest is Rigetti's hosted QCS JupyterLab, where
pyQuil, quilc, the QVM and your ~/.qcs credentials are all preconfigured.  (On a
local machine you instead need pyQuil >= 4.7, ~/.qcs/{settings,secrets}.toml, and
a running quilc server.)  This script imports NO qiskit and needs only the
`qcs_quil/` folder produced by export_quil.py.

On-demand access (no reservation): a job submitted without an active reservation
is automatically routed to the on-demand queue; you are billed only for QPU
execution time (queue time, failed and pre-empted runs are free), reported below
as execution_duration_microseconds.

Cost discipline: everything except `qc.run` is free.  Validate first with
`--qvm` (noiseless simulator, same ISA/topology) and pre-stage with
`--compile-only` (quilc only, no QPU) BEFORE running on the real device.

Workflow:
    python run_qcs.py --list                              # discover the exact QPU name
    python run_qcs.py --bell --device Cepheus-1-108Q      # 1-job access smoke test
    python run_qcs.py --qvm                               # free noiseless dry run
    python run_qcs.py --compile-only --device Cepheus-1-108Q  # pre-compile, no QPU time
    python run_qcs.py --device Cepheus-1-108Q --shots 2000    # the real QPU run
"""

import argparse
import json
import os

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))


def _load_ro(result):
    """Return the ro readout array (shots x bits); get_register_map is the current API."""
    getreg = getattr(result, "get_register_map", None)
    if getreg is not None:
        ro = getreg().get("ro")
        if ro is not None:
            return np.asarray(ro)
    data = getattr(result, "readout_data", None)   # deprecated fallback
    if data is not None and data.get("ro") is not None:
        return np.asarray(data["ro"])
    raise RuntimeError("no 'ro' register found in execution result")


def _observable(job, bits):
    """(value, rho) for a job given its ro bit array.  rho is None for alpha jobs."""
    bits = np.asarray(bits).astype(int)
    if job["kind"] == "alpha":
        b = bits[:, 0]                          # single spinor bit after R^-1
        return float((2 * b - 1).mean()), None  # P(1) - P(0)
    shots, nb = bits.shape                       # density: rebuild flat index i
    ints = (bits * (1 << np.arange(nb))).sum(axis=1)
    N = 2 ** job["npos"]
    rho = np.bincount(ints >> 2, minlength=N).astype(float) / shots   # x = i >> 2
    return float(rho @ np.arange(N)), rho


def _bell_program(Program):
    p = Program()
    ro = p.declare("ro", "BIT", 2)
    from pyquil.gates import H, CNOT, MEASURE
    p += H(0)
    p += CNOT(0, 1)
    p += MEASURE(0, ro[0])
    p += MEASURE(1, ro[1])
    return p


_TWO_Q_NAMES = {"CNOT", "CZ", "ISWAP", "SWAP", "XY", "CPHASE", "PSWAP"}


def _pin_text(quil_text, mapping):
    """Rewrite a Quil program onto specific physical qubits (logical->physical via
    `mapping`) and force NAIVE rewiring, so quilc uses exactly those qubits."""
    def remap(tok):
        return str(mapping.get(int(tok), int(tok)))
    out = ['PRAGMA INITIAL_REWIRING "NAIVE"']
    for line in quil_text.splitlines():
        s = line.strip()
        if not s or s.startswith(("DECLARE", "PRAGMA", "HALT", "RESET", "#")):
            out.append(line)
            continue
        toks = s.split()
        base = toks[0].split("(")[0]
        if base == "MEASURE":
            toks[1] = remap(toks[1])
        elif base in _TWO_Q_NAMES:
            toks[1], toks[2] = remap(toks[1]), remap(toks[2])
        else:
            toks[-1] = remap(toks[-1])
        out.append(" ".join(toks))
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default=_HERE, help="folder with manifest.json + *.quil.")
    ap.add_argument("--device", default="Cepheus-1-108Q", help="QCS QPU name (see --list).")
    ap.add_argument("--shots", type=int, default=2000)
    ap.add_argument("--qvm", action="store_true",
                    help="run on the noiseless QVM (same ISA/topology) -- free, for validation.")
    ap.add_argument("--compile-only", action="store_true",
                    help="quilc-compile every program but do not run (no QPU time).")
    ap.add_argument("--bell", action="store_true",
                    help="run a single Bell circuit as an access smoke test, then exit.")
    ap.add_argument("--list", action="store_true",
                    help="list available quantum computers and exit.")
    ap.add_argument("--timeout", type=float, default=120.0,
                    help="compiler/execution timeout seconds (on-demand may queue).")
    ap.add_argument("--qubits", default=None,
                    help="pin logical qubits 0,1,2,... to these physical qubits, "
                         "comma-separated (e.g. 90,99,101); forces NAIVE rewiring.")
    ap.add_argument("--out", default=os.path.join(_HERE, "qcs_results.json"))
    args = ap.parse_args()

    qmap = None
    if args.qubits:
        qmap = {i: int(p) for i, p in enumerate(args.qubits.split(","))}

    try:
        from pyquil import Program, get_qc
    except ImportError:
        raise SystemExit(
            "pyQuil not found.  Run this in Rigetti's QCS JupyterLab (preinstalled), "
            "or `pip install 'pyquil>=4.7'` in a venv with ~/.qcs credentials + a quilc server.")

    if args.list:
        from pyquil import list_quantum_computers
        print("available quantum computers:")
        for name in list_quantum_computers():
            print("  ", name)
        print("(or use the CLI: `qcs list-quantum-processors`)")
        return

    qc = get_qc(args.device, as_qvm=args.qvm,
                compiler_timeout=args.timeout, execution_timeout=args.timeout)
    where = f"{args.device} (QVM)" if args.qvm else args.device
    print(f"target: {where}   shots={args.shots}")
    if qmap:
        print(f"  pinned logical->physical qubits: {qmap}  (INITIAL_REWIRING NAIVE)")

    if args.bell:
        prog = _bell_program(Program)
        if qmap:
            prog = Program(_pin_text(str(prog), qmap))
        prog.wrap_in_numshots_loop(args.shots)
        res = qc.run(qc.compile(prog))
        bits = _load_ro(res)
        p00 = float(((bits[:, 0] == 0) & (bits[:, 1] == 0)).mean())
        p11 = float(((bits[:, 0] == 1) & (bits[:, 1] == 1)).mean())
        dur = getattr(res, "execution_duration_microseconds", None)
        print(f"  Bell: P(00)={p00:.3f}  P(11)={p11:.3f}  (ideal 0.5/0.5)"
              + (f"   {dur} us" if dur else ""))
        print("access OK." if p00 + p11 > 0.7 else "check output: correlations look weak.")
        return

    with open(os.path.join(args.dir, "manifest.json")) as fh:
        manifest = json.load(fh)
    jobs = manifest["jobs"]
    print(f"loaded {len(jobs)} jobs from {args.dir}/manifest.json  "
          f"(npos={manifest['npos']}, m~={manifest['mass']}, E={manifest['E']:.4f})")

    # compile every program first (free); only then spend QPU time running them.
    executables = {}
    for job in jobs:
        with open(os.path.join(args.dir, job["quil"])) as fh:
            raw = fh.read()
        if qmap:
            raw = _pin_text(raw, qmap)
        prog = Program(raw)
        prog.wrap_in_numshots_loop(args.shots)
        executables[job["name"]] = qc.compile(prog)
        print(f"  compiled {job['name']}")
    if args.compile_only:
        print("compile-only: all programs compiled, no QPU time used.")
        return

    total_us = 0.0
    ref_key = "emulator" if "emulator" in jobs[0] else "exact"
    print(f"\n  job          t   kind      QPU        exact    {ref_key}")
    for job in jobs:
        res = qc.run(executables[job["name"]])
        val, rho = _observable(job, _load_ro(res))
        dur = getattr(res, "execution_duration_microseconds", None)
        total_us += float(dur) if dur else 0.0
        job["qpu"] = val
        if rho is not None:
            job["qpu_rho"] = [float(v) for v in rho]
        print(f"  {job['name']:11s} {job['t']:2d}   {job['kind']:8s} "
              f"{val:+7.4f}   {job['exact']:+7.4f}   {job.get(ref_key, float('nan')):+7.4f}")

    manifest["device"] = args.device
    manifest["qvm"] = args.qvm
    manifest["shots"] = args.shots
    manifest["total_execution_us"] = total_us
    with open(args.out, "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"\ntotal QPU execution time: {total_us:.0f} us  ({total_us / 1e6:.4f} s billed)")
    print(f"results written to {args.out}")


if __name__ == "__main__":
    main()
