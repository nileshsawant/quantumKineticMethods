"""
Compile, check and run the mid-circuit repetition-code programs of repcode_mid_qcs.py.

  --compile   quilc-compile every program for the device on the pinned chain (control flow kept,
              protoquil off) and save the native programs to <dir>/native/.  Reports the number of
              two-qubit gates against the abstract count; more means quilc inserted SWAPs, which
              is not safe across JUMP/LABEL blocks (see below).  Needs QCS credentials, no QPU time.
  --qvm-check run the saved native programs, relabelled onto qubits 0..n-1, on the noiseless QVM
              and check that every ancilla reads 0 and every data qubit reads the logical state.
              The QVM returns only the last pass of an in-program loop, so loop programs are
              checked for one pass.  Needs a local qvm server, no credentials.
  --qpu       run the saved native programs on the device (spends QPU time; needs --confirm).
              Loop programs are submitted once with shot_count = shots; unrolled programs (K shots
              each) are submitted ceil(shots / K) times.  Shots are rebuilt from the raw readout
              stream.

On the 9q-square QVM, where quilc must route, the conditional-reset programs came out wrong
(noiseless syndromes nonzero), so programs with control flow must compile without SWAPs.

On Cepheus-1 a gate that follows a mid-circuit MEASURE on the same qubit too closely mostly fails
(qubit 84, 50 x [RX(pi); MEASURE]: the readout changed after only 10% of the flips); with a 2 us
DELAY after each MEASURE it changed after 63% (q84) and 84% (q100) (mcm_probe2.py).  --qpu therefore
adds DELAY q <--measure-delay> after every MEASURE.  Each MEASURE calibration also plays a 1->2
pulse before readout; repeating it afterwards (--unshelve) made results worse (55%).

Usage:
    python run_repmid_qcs.py --dir ../qcs_quil_repmid --compile
    python run_repmid_qcs.py --dir ../qcs_quil_repmid --qvm-check
    python run_repmid_qcs.py --dir ../qcs_quil_repmid --qpu --only name1,name2 --confirm
"""

import argparse
import json
import os
import re

import numpy as np

from repcode_mid_qcs import readout_stream_to_bits, unchain, unroll_native
from run_qcs import _pin_text

CHAIN = [100, 101, 102, 93, 84, 83, 74, 65, 56, 57, 58, 67, 68, 77, 86, 95, 104]
TWOQ = re.compile(r"^(CZ|ISWAP|XY|CPHASE|CNOT|SWAP)\b")


GATE = re.compile(r"^([A-Z][A-Z\-]*(?:\(.*\))?)\s+([\d\s]+)$")


def relabel(native):
    """Map the physical qubits of a native program onto 0..n-1 (sorted)."""
    lines = native.splitlines()
    phys = set()
    for line in lines:
        if line.startswith("MEASURE"):
            phys.add(int(line.split()[1]))
        elif GATE.match(line):
            phys |= {int(q) for q in GATE.match(line).group(2).split()}
    idx = {p: k for k, p in enumerate(sorted(phys))}
    out = []
    for line in lines:
        if line.startswith("MEASURE"):
            t = line.split()
            out.append(f"MEASURE {idx[int(t[1])]} {t[2]}")
        elif GATE.match(line):
            g, qs = GATE.match(line).groups()
            out.append(g + " " + " ".join(str(idx[int(q)]) for q in qs.split()))
        else:
            out.append(line)
    return "\n".join(out) + "\n", len(phys)


def shelving_pulses(qc):
    """{qubit: PULSE line} of the 1->2 pulse that each MEASURE calibration plays before readout."""
    lines = str(qc.compiler.get_calibration_program()).splitlines()
    out = {}
    for i, ln in enumerate(lines):
        m = re.match(r"DEFCAL MEASURE (\d+) ", ln)
        if m:
            for x in lines[i + 1:i + 4]:
                if "charge_tx_f12" in x and "PULSE" in x:
                    out[int(m.group(1))] = x.strip().replace("NONBLOCKING ", "")
    return out


def measure_defcals_without_shelving(qc, qubits, duration=None):
    """Quil-T lines redefining MEASURE on these qubits as the device calibration minus its 1->2
    pulse, so a qubit read as 1 stays in |1>.  duration (s) shortens the readout pulse and
    capture; the device's 0/1 threshold is then off, so classify from the raw IQ points."""
    lines = str(qc.compiler.get_calibration_program()).splitlines()
    out = []
    for q in qubits:
        for i, ln in enumerate(lines):
            if re.match(rf"DEFCAL MEASURE {q} ", ln):
                out.append(ln)
                for x in lines[i + 1:]:
                    if not x.startswith((" ", "\t")) or not x.strip():
                        break
                    if "charge_tx_f12" in x:
                        continue
                    if duration and ("readout_tx" in x or "readout_rx" in x) and "duration:" in x:
                        x = re.sub(r"duration: [0-9.eE+-]+", f"duration: {duration:g}", x)
                    out.append(x)
                break
    return out


def add_dd(native, data, window, n_slots, na, R):
    """Echo the data qubits (RX(pi), RX(pi) at window/4 and 3*window/4) while the ancillas are read
    out: inserted after the MEASURE of the last ancilla of each round (ro slots are round-major)."""
    dd = [f"DELAY {d} {window / 4:g}\nRX(pi) {d}\nDELAY {d} {window / 2:g}\nRX(pi) {d}\nDELAY {d} {window / 4:g}"
          for d in sorted(data)]
    out = []
    for ln in native.splitlines():
        out.append(ln)
        m = re.match(r"MEASURE \d+ ro\[(\d+)\]", ln)
        if m:
            k = int(m.group(1)) % n_slots
            if k < R * na and k % na == na - 1:
                out += dd
    return "\n".join(out) + "\n"


def add_unshelve(native, f12):
    """Repeat the 1->2 pulse after every MEASURE so a qubit read as 1 is back in |1>, not |2>."""
    out = []
    for ln in native.splitlines():
        out.append(ln)
        if ln.startswith("MEASURE") and int(ln.split()[1]) in f12:
            out.append(f12[int(ln.split()[1])])
    return "\n".join(out) + "\n"


def add_delay(native, seconds):
    """DELAY the measured qubit after every MEASURE, so its readout resonator can empty."""
    out = []
    for ln in native.splitlines():
        out.append(ln)
        if ln.startswith("MEASURE") and seconds > 0:
            out.append(f"DELAY {ln.split()[1]} {seconds:g}")
    return "\n".join(out) + "\n"


def check(bits, job):
    na, R, s = job["d"] - 1, job["rounds"], job["state"]
    return float(bits[:, :R * na].any(axis=1).mean()), float((bits[:, R * na:] != s).any(axis=1).mean())


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True)
    ap.add_argument("--device", default="Cepheus-1-108Q")
    ap.add_argument("--compile", action="store_true")
    ap.add_argument("--qubits", default=None,
                    help="comma-separated physical chain for --compile (default CHAIN)")
    ap.add_argument("--gap", type=float, default=0.0,
                    help="--compile: DELAY (s) between unrolled shots so the qubits relax to |0>")
    ap.add_argument("--noshelve", action="store_true",
                    help="--qpu: redefine MEASURE without its 1->2 pulse")
    ap.add_argument("--ro-duration", type=float, default=None,
                    help="--qpu with --noshelve: readout pulse/capture length (s)")
    ap.add_argument("--dd", type=float, default=0.0,
                    help="--qpu: echo the data qubits over this window (s) during ancilla readout")
    ap.add_argument("--qvm-check", action="store_true")
    ap.add_argument("--qpu", action="store_true")
    ap.add_argument("--confirm", action="store_true", help="required with --qpu")
    ap.add_argument("--dry", action="store_true", help="--qpu: translate only, run nothing")
    ap.add_argument("--unshelve", action="store_true",
                    help="repeat the 1->2 readout pulse after each MEASURE")
    ap.add_argument("--measure-delay", type=float, default=2e-6,
                    help="DELAY (s) on a qubit after each of its MEASUREs")
    ap.add_argument("--tag", default="", help="suffix of the results file")
    ap.add_argument("--only", default=None)
    ap.add_argument("--shots", type=int, default=None, help="override the manifest shot count")
    ap.add_argument("--timeout", type=float, default=300.0)
    args = ap.parse_args()

    from pyquil import Program, get_qc
    man = json.load(open(os.path.join(args.dir, "manifest.json")))
    jobs = man["jobs"]
    if args.only:
        keep = set(args.only.split(","))
        jobs = [j for j in jobs if j["name"] in keep]
    ndir = os.path.join(args.dir, "native")

    if args.compile:
        os.makedirs(ndir, exist_ok=True)
        qc = get_qc(args.device, compiler_timeout=args.timeout)
        chain = [int(x) for x in args.qubits.split(",")] if args.qubits else CHAIN
        for j in jobs:
            raw = open(os.path.join(args.dir, j["quil"])).read()
            nat = str(qc.compiler.quil_to_native_quil(
                Program(_pin_text(raw, dict(enumerate(chain)))), protoquil=False))
            n_abs = sum(1 for line in raw.splitlines() if line.startswith("CNOT"))
            n_nat = sum(1 for line in nat.splitlines() if TWOQ.match(line))
            if j.get("unroll", 1) > 1:
                nat = unroll_native(nat, j["n_slots"], j["unroll"], args.gap)
            j["gap"] = args.gap
            j["qubits"] = chain[:j["n_qubits"]]
            open(os.path.join(ndir, j["quil"]), "w").write(nat)
            flag = "" if n_nat == n_abs else "   <-- extra two-qubit gates (routing)"
            print(f"  {j['name']:30s} CNOT {n_abs:3d}  native 2q {n_nat:3d}{flag}")
        json.dump(man, open(os.path.join(args.dir, "manifest.json"), "w"), indent=2)

    if args.qvm_check:
        from repcode_mid_qcs import to_stim
        for j in jobs:
            nat, n = relabel("\n".join(ln for ln in open(os.path.join(ndir, j["quil"])).read().splitlines()
                                       if not ln.startswith("DELAY")))
            qc = get_qc(f"{n}q-qvm", execution_timeout=args.timeout)
            p = Program(nat)
            p.wrap_in_numshots_loop(1 if j["loop"] else 200)
            r = qc.run(p, memory_map={"shot_count": [3]} if j["loop"] else None)
            bits = np.asarray(r.get_register_map()["ro"]).astype(int)
            if j.get("chain"):
                bits = np.concatenate([unchain(row.reshape(-1, j["n_slots"]), j["code"], j["d"],
                                               j["rounds"], j["state"], j.get("inject", False))
                                       for row in bits])
            bits = bits.reshape(-1, j["n_slots"])
            ref = to_stim(j["code"], j["d"], j["rounds"], j["state"], j["reset"],
                          inject=j.get("inject", False)).compile_sampler().sample(1)[0].astype(int)
            bad = float((bits != ref).any(axis=1).mean())
            print(f"  {j['name']:34s} {n:2d} qubits  shots {len(bits):5d}  reference "
                  f"{''.join(map(str, ref))}  P(shot differs) {bad:.3f}  "
                  f"{'ok' if bad == 0 else 'WRONG'}")

    if args.qpu:
        if not args.confirm:
            raise SystemExit("--qpu spends QPU time; add --confirm")
        qc = get_qc(args.device, compiler_timeout=args.timeout, execution_timeout=args.timeout)
        f12 = shelving_pulses(qc) if args.unshelve else {}
        out_path = os.path.join(args.dir, f"qcs_results_repmid{'_' + args.tag if args.tag else ''}.json")
        done = {}
        if os.path.exists(out_path):
            done = {j["name"]: j for j in json.load(open(out_path))["jobs"]}
        for j in jobs:
            shots = args.shots or j["shots"]
            K = j.get("unroll", 1)
            text = add_delay(add_unshelve(open(os.path.join(ndir, j["quil"])).read(), f12),
                             args.measure_delay)
            if args.dd:
                text = add_dd(text, set(j["qubits"][0::2]), args.dd, j["n_slots"], j["d"] - 1,
                              j["rounds"])
            if args.noshelve:
                mq = sorted({int(ln.split()[1]) for ln in text.splitlines() if ln.startswith("MEASURE")})
                text = "\n".join(measure_defcals_without_shelving(qc, mq, args.ro_duration)) + "\n" + text
            nat = Program(text)
            nat.wrap_in_numshots_loop(1)
            exe = qc.compiler.native_quil_to_executable(nat)
            if args.dry:
                print(f"  translated {j['name']}", flush=True)
                continue
            reps = 1 if j["loop"] else -(-shots // K)    # one job per K unrolled shots
            parts, us, raw = [], 0.0, []
            for _ in range(reps):
                r = qc.run(exe, memory_map={"shot_count": [shots]} if j["loop"] else None)
                us += float(getattr(r, "execution_duration_microseconds", 0) or 0)
                q = r.data.result_data.to_qpu()
                raw.append({"mappings": dict(q.mappings),
                            "streams": {k: [int(x) for x in np.asarray(v.inner()).ravel()]
                                        for k, v in q.readout_values.items()
                                        if k.endswith("_classified")},
                            "iq": {k: [[round(float(z.real), 3), round(float(z.imag), 3)]
                                       for z in np.asarray(v.inner()).ravel()]
                                   for k, v in q.readout_values.items()
                                   if k.endswith("_unclassified")}})
                try:
                    b = readout_stream_to_bits(q, j["n_slots"] * K)
                except ValueError as e:                  # keep the raw data for diagnosis
                    print(f"  {j['name']}: {e}; stream lengths "
                          f"{ {k: len(v) for k, v in raw[-1]['streams'].items()} }")
                    continue
                b = b.reshape(-1, j["n_slots"])
                raw[-1]["bits"] = ["".join(map(str, row)) for row in b]
                if j.get("chain") and not j.get("gap"):
                    try:
                        b = unchain(b, j["code"], j["d"], j["rounds"], j["state"],
                                    j.get("inject", False))
                    except ImportError:                  # no stim here: unchain raw[...]["bits"] later
                        pass
                parts.append(b)
            bits = np.concatenate(parts) if parts else np.zeros((0, j["n_slots"]), int)
            pa, pd = check(bits, j) if len(bits) else (float("nan"), float("nan"))
            j = dict(j, shots_returned=len(bits), jobs_run=reps, execution_us=us, raw=raw,
                     bitstrings=["".join(map(str, row)) for row in bits])
            done[j["name"]] = j
            json.dump({"device": args.device, "jobs": list(done.values())},
                      open(out_path, "w"), indent=1)
            print(f"  {j['name']:30s} shots {len(bits):5d}  P(any ancilla=1) {pa:.3f}  "
                  f"P(data wrong) {pd:.3f}  {j['execution_us']:.0f} us")


if __name__ == "__main__":
    main()
