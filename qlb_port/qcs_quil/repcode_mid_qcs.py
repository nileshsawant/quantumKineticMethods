"""
Repetition-code memory experiments with mid-circuit measurement of the ancillas (Rigetti QCS).

Same codes, layout and check gates as repcode_qcs.py, but every round ends by measuring the d-1
ancillas into their own readout slots, so a fault is located in time as well as in space.
Two ways of handling an ancilla after it is measured:
  * reset:    conditional X on the ancilla if it read 1 (JUMP-UNLESS / X / LABEL), the active reset
              of the Rigetti dynamic-control-flow guide; detector of round r is m_r xor m_(r-1);
  * no reset: the ancilla keeps its value and accumulates the parities, so m_r is the running sum;
              detector of round r is m_r xor m_(r-2).  This needs no control flow inside a shot.

The QPU ignores the shot count for programs that are not ProtoQuil, so with --loop the program
repeats itself `shot_count` times (the guide's template).  A loop of this kind returned a single
pass in the September test, so --unroll K instead puts K shots one after another into one program
with their own readout slots (the one-shot program is compiled once and its native code copied K
times, see unroll_native; quilc is slow on long programs with many branches).  In both cases every
qubit is returned to |0> at the end of a shot by a conditional X on its last readout, and the raw
readout stream has one value per MEASURE execution; group it with readout_stream_to_bits().

Readout record of one shot: R rounds x (d-1) ancillas (round-major), then the d data qubits.

The Cepheus-1 translator accepts a branch on a readout only when a single qubit has been measured
(two or more measured qubits fail to schedule whatever DELAY is added, and RESET has no calibration
mid-program).  --chain therefore writes branch-free programs: no-reset scheme, no return to |0>
between the unrolled shots, so each shot starts in the basis state its qubits were last read in.
The circuits are Clifford, so that start state flips a fixed set of records (start_response), and
unchain() removes it shot by shot using the previous shot's raw readout.

Writes <out>/<name>.quil and <out>/manifest.json.  Free (no QPU time).

Usage:
    python repcode_mid_qcs.py --distances 3,5 --rounds 1,3 --out ../qcs_quil_repmid [--loop]
"""

import argparse
import json
import os
import re

import numpy as np

from repcode_qcs import CODES, layout


def ops_mid(code, d, R, state, reset, inject=False):
    """Abstract op list: ('X',q) ('H',q) ('Z',q) ('CX',c,t) ('IDLE',[q]) ('MM',q) ('CRESET',q)
    ('M',q).  inject puts a known error on data qubit 0 after preparation (X for the bit-flip
    code, Z for the phase-flip code) so that the noiseless checks exercise every branch."""
    data, anc = layout(d)
    out = []
    for q in data:
        if state:
            out.append(("X", q))
        if code == "phase":
            out.append(("H", q))
    if inject:
        out.append(("X" if code == "bit" else "Z", data[0]))
    for _ in range(R):
        if code == "bit":
            for k, a in enumerate(anc):
                out.append(("CX", data[k], a))
            for k, a in enumerate(anc):
                out.append(("CX", data[k + 1], a))
        else:
            for a in anc:
                out.append(("H", a))
            for k, a in enumerate(anc):
                out.append(("CX", a, data[k]))
            for k, a in enumerate(anc):
                out.append(("CX", a, data[k + 1]))
            for a in anc:
                out.append(("H", a))
        out.append(("IDLE", list(data)))
        for a in anc:
            out.append(("MM", a))
        out.append(("IDLE_M", list(data)))            # data wait while the ancillas are read out
        if reset:
            for a in anc:
                out.append(("CRESET", a))
    if code == "phase":
        for q in data:
            out.append(("H", q))
    for q in data:
        out.append(("M", q))
    return out


def to_quil(op_list, n_qubits, loop=False, end_reset=False):
    nm = sum(1 for o in op_list if o[0] in ("MM", "M"))
    lines = [f"DECLARE ro BIT[{nm}]"]
    if loop:
        lines += ["DECLARE shot_count INTEGER[1]", "LABEL @START"]
    k, last, nlab, fresh = 0, {}, 0, set()
    for o in op_list:
        if o[0] in ("X", "H", "Z"):
            lines.append(f"{o[0]} {o[1]}")
        elif o[0] == "CX":
            lines.append(f"CNOT {o[1]} {o[2]}")
        elif o[0] in ("MM", "M"):
            lines.append(f"MEASURE {o[1]} ro[{k}]")
            last[o[1]] = k
            fresh.discard(o[1])
            k += 1
        elif o[0] == "CRESET":
            lines += [f"JUMP-UNLESS @L{nlab} ro[{last[o[1]]}]", f"X {o[1]}", f"LABEL @L{nlab}"]
            fresh.add(o[1])
            nlab += 1
    if loop or end_reset:
        for q in range(n_qubits):                     # return to |0> for the next shot
            if q in last and q not in fresh:
                lines += [f"JUMP-UNLESS @L{nlab} ro[{last[q]}]", f"X {q}", f"LABEL @L{nlab}"]
                nlab += 1
    if loop:
        lines += ["SUB shot_count[0] 1", "JUMP-UNLESS @END shot_count[0]", "JUMP @START",
                  "LABEL @END"]
    return "\n".join(lines) + "\n"


def unroll_native(native, n_slots, K, gap=0.0):
    """K copies of a compiled one-shot program (ending with the return to |0>) in one program:
    labels and ro slots renamed per copy, and HALT replaced by a jump to the next copy.  gap > 0
    adds a DELAY of that many seconds on every measured qubit between copies (passive reset)."""
    body = [ln for ln in native.splitlines() if ln.strip() and not ln.startswith(("DECLARE", "PRAGMA"))]
    while body and body[-1].strip() == "HALT":
        body.pop()
    head = [ln for ln in native.splitlines() if ln.startswith("PRAGMA")]
    qubits = sorted({int(ln.split()[1]) for ln in body if ln.startswith("MEASURE")})
    out = [f"DECLARE ro BIT[{n_slots * K}]"] + head
    for c in range(K):
        for ln in body:
            ln = re.sub(r"@([\w\-]+)", rf"@\1-c{c}", ln)
            ln = re.sub(r"ro\[(\d+)\]", lambda m: f"ro[{int(m.group(1)) + c * n_slots}]", ln)
            out.append(f"JUMP @NEXT-c{c}" if ln.strip() == "HALT" else ln)
        if any(ln.strip() == "HALT" for ln in body):
            out.append(f"LABEL @NEXT-c{c}")
        if gap > 0 and c < K - 1:
            out += [f"DELAY {q} {gap:g}" for q in qubits]
    return "\n".join(out) + "\n"


def to_stim(code, d, R, state, reset, noise=None, inject=False, start=()):
    """Stim circuit with one detector per ancilla per round, final data detectors, and the final
    value of data qubit 0 as the observable.  noise: p1, p2 (depolarizing after gates), pm
    (readout record flip; a misread ancilla is also mis-reset), pz, px (data errors per round),
    km (data errors during each mid-circuit readout, as a multiple of pz, px).  start: qubits that
    begin in |1>."""
    import stim
    nz = noise or {}
    c = stim.Circuit()
    data, anc = layout(d)
    na = len(anc)
    for q in start:
        c.append("X", [q])
    for o in ops_mid(code, d, R, state, reset, inject):
        if o[0] == "Z":
            c.append("Z", [o[1]])
        elif o[0] in ("X", "H"):
            c.append(o[0], [o[1]])
            if nz.get("p1"):
                c.append("DEPOLARIZE1", [o[1]], nz["p1"])
        elif o[0] == "CX":
            c.append("CX", [o[1], o[2]])
            if nz.get("p2"):
                c.append("DEPOLARIZE2", [o[1], o[2]], nz["p2"])
        elif o[0] in ("IDLE", "IDLE_M"):
            f = nz.get("km", 1.0) if o[0] == "IDLE_M" else 1.0
            if nz.get("pz") and f:
                c.append("Z_ERROR", o[1], nz["pz"] * f)
            if nz.get("px") and f:
                c.append("X_ERROR", o[1], nz["px"] * f)
        elif o[0] in ("MM", "M"):
            c.append("M", [o[1]], nz.get("pm", 0.0))
        elif o[0] == "CRESET":                        # the round's ancilla records are the last na
            c.append("CX", [stim.target_rec(-na + anc.index(o[1])), o[1]])
    total = R * na + d

    def rec(idx):
        return stim.target_rec(idx - total)

    def m(r, j):
        return [rec(r * na + j)] if 0 <= r < R else []

    lag = 1 if reset else 2
    for r in range(R):
        for j in range(na):
            c.append("DETECTOR", m(r, j) + m(r - lag, j), [2 * j + 1, r])
    for j in range(na):
        prev = m(R - 1, j) + (m(R - 2, j) if not reset else [])
        c.append("DETECTOR", [rec(R * na + j), rec(R * na + j + 1)] + prev, [2 * j + 1, R])
    c.append("OBSERVABLE_INCLUDE", [rec(R * na)], 0)
    return c


def last_slot(d, R):
    """Record index of the last readout of every qubit in one shot (no-reset scheme)."""
    data, anc = layout(d)
    out = {q: R * (d - 1) + i for i, q in enumerate(data)}
    out.update({a: (R - 1) * (d - 1) + j for j, a in enumerate(anc)})
    return out


def start_response(code, d, R, state, inject=False):
    """(reference record, {qubit: records flipped when that qubit starts in |1>}), noiseless."""
    def rec(start):
        s = to_stim(code, d, R, state, False, inject=inject, start=start).compile_sampler()
        a = s.sample(2).astype(int)
        assert (a[0] == a[1]).all(), "record is not deterministic"
        return a[0]
    ref = rec(())
    return ref, {q: rec((q,)) ^ ref for q in range(2 * d - 1)}


def unchain(bits, code, d, R, state, inject=False):
    """Remove the start-state flips from consecutive chained shots of one job (shot 0 starts in
    |0...0>; shot k starts in the state shot k-1 read out)."""
    _, resp = start_response(code, d, R, state, inject)
    last = last_slot(d, R)
    out = bits.copy()
    for k in range(1, len(bits)):
        for q, idx in last.items():
            if bits[k - 1, idx]:
                out[k] ^= resp[q]
    return out


def readout_stream_to_bits(qpu_data, n_slots):
    """Rebuild a shots x n_slots bit array from QPUResultData.  mappings sends ro[k] to a readout
    stream; all measurements of one qubit share its stream (e.g. q101_classified), which holds one
    value per MEASURE execution in time order.  Every MEASURE here runs once per shot and ro indices
    follow program order, so a stream with m slots reshapes to shots x m."""
    groups = {}
    for k in range(n_slots):
        groups.setdefault(qpu_data.mappings[f"ro[{k}]"], []).append(k)
    cols, n = {}, None
    for name, ks in groups.items():
        v = np.asarray(qpu_data.readout_values[name].inner()).ravel()
        if len(v) % len(ks) or (n is not None and len(v) // len(ks) != n):
            raise ValueError(f"stream {name}: {len(v)} values for slots {ks}, shots so far {n}")
        n = len(v) // len(ks)
        for i, k in enumerate(ks):
            cols[k] = v.reshape(n, len(ks))[:, i]
    return np.stack([cols[k] for k in range(n_slots)], axis=1).astype(int)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--distances", default="3,5")
    ap.add_argument("--rounds", default="1,3")
    ap.add_argument("--shots", type=int, default=2000)
    ap.add_argument("--loop", action="store_true", help="in-program shot loop (for the QPU)")
    ap.add_argument("--unroll", type=int, default=1,
                    help="shots per program; the one-shot program is compiled, then copied")
    ap.add_argument("--inject", action="store_true", help="add a known error (test programs)")
    ap.add_argument("--chain", action="store_true",
                    help="branch-free: no-reset scheme, no reset between unrolled shots")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    jobs = []
    for code in CODES:
        for d in map(int, args.distances.split(",")):
            for R in map(int, args.rounds.split(",")):
                for reset in ((0,) if args.chain else (1, 0)):
                    for state in (0, 1):
                        name = f"rmid_{code}_d{d}_r{R}_{'rst' if reset else 'nor'}_s{state}"
                        name += "_inj" if args.inject else ""
                        nq = 2 * d - 1
                        with open(os.path.join(args.out, name + ".quil"), "w") as fh:
                            fh.write(to_quil(ops_mid(code, d, R, state, reset, args.inject), nq,
                                             args.loop, args.unroll > 1 and not args.chain))
                        jobs.append({"name": name, "quil": name + ".quil", "kind": "raw",
                                     "role": "repmid", "code": code, "d": d, "rounds": R,
                                     "reset": reset, "state": state, "n_qubits": nq,
                                     "n_slots": R * (d - 1) + d, "loop": args.loop,
                                     "unroll": args.unroll, "inject": args.inject,
                                     "chain": args.chain, "shots": args.shots})
    with open(os.path.join(args.out, "manifest.json"), "w") as fh:
        json.dump({"description": "repetition code, mid-circuit ancilla measurement",
                   "jobs": jobs}, fh, indent=2)
    print(f"wrote {len(jobs)} programs to {args.out}")


if __name__ == "__main__":
    main()
