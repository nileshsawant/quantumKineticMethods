"""
Choose physical qubits for a circuit from live calibration data (vendor-neutral).

The circuit is described by its interaction graph: the logical qubit pairs that share two-qubit gates
(with gate counts), the number of one-qubit gates per qubit, and the measured qubits.  Every
SWAP-free embedding of that graph into the device coupling map is a candidate layout, scored by the
estimated probability that nothing fails,

    log P = sum_edges n2 * ln f2q  +  sum_qubits n1 * ln f1q  +  sum_measured ln fro,

using the fidelities of exactly the physical qubits and couplers the layout uses.  A branch-and-bound
search returns the best layouts without enumerating all embeddings.

Calibration sources: a generic JSON file, a saved or live Rigetti QCS instruction-set architecture,
or a Qiskit ``Target`` (library use).  Interaction sources: an edge list, a Quil program, an OpenQASM
file, or a Qiskit ``QuantumCircuit``.

Examples:
    python3 qlb_port/select_qubits.py --qcs-isa isa.json --edges 0-1,1-3,2-3 --measured 0,1,2,3
    python3 qlb_port/select_qubits.py --qcs Cepheus-1-108Q --quil qcs_quil_n2/alpha_t4.quil --top 5
    python3 qlb_port/select_qubits.py --calibration cal.json --qasm circuit.qasm

Generic calibration JSON:
    {"one_qubit": {"0": 0.999, ...}, "readout": {"0": 0.98, ...}, "two_qubit": {"0-1": 0.99, ...}}
Qubits or couplers missing from the calibration are treated as unusable.
"""

import argparse
import heapq
import json
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field


@dataclass
class Calibration:
    f1q: dict                      # physical qubit -> one-qubit gate fidelity
    fro: dict                      # physical qubit -> readout fidelity
    f2q: dict                      # frozenset({a, b}) -> two-qubit gate fidelity

    def neighbors(self):
        nb = defaultdict(set)
        for e in self.f2q:
            a, b = tuple(e)
            nb[a].add(b)
            nb[b].add(a)
        return nb


@dataclass
class Interaction:
    n_qubits: int
    edges: Counter                 # frozenset({i, j}) of logical qubits -> number of two-qubit gates
    one_qubit: Counter = field(default_factory=Counter)   # logical qubit -> number of one-qubit gates
    measured: set = field(default_factory=set)
    labels: list = None            # original qubit index of each logical qubit (default: itself)

    def __post_init__(self):
        if self.labels is None:
            self.labels = list(range(self.n_qubits))


# ----------------------------------------------------------------------------- calibration loaders

def _usable(v):
    return v is not None and 0.0 < v <= 1.0


def calibration_from_json(path):
    d = json.load(open(path))
    f2q = {}
    for k, v in d.get("two_qubit", {}).items():
        a, b = (int(x) for x in re.split(r"[-,]", k))
        if _usable(v):
            f2q[frozenset((a, b))] = max(v, f2q.get(frozenset((a, b)), 0.0))
    return Calibration({int(q): v for q, v in d.get("one_qubit", {}).items() if _usable(v)},
                       {int(q): v for q, v in d.get("readout", {}).items() if _usable(v)}, f2q)


def calibration_from_qcs_isa(isa):
    """From a Rigetti QCS instruction-set architecture (dict, JSON text, or qcs_sdk object)."""
    if hasattr(isa, "json"):
        isa = isa.json()
    if isinstance(isa, str):
        isa = json.loads(isa)

    def sites(ops, op_name, char_name):
        for op in ops:
            if op["name"] == op_name:
                for s in op["sites"]:
                    for c in s["characteristics"]:
                        if c["name"] == char_name and _usable(c.get("value")):
                            yield tuple(s["node_ids"]), c["value"]

    f1q = {n[0]: v for n, v in sites(isa.get("benchmarks", []), "randomized_benchmark_1q", "fRB")}
    fro = {n[0]: v for n, v in sites(isa["instructions"], "MEASURE", "fRO")}
    f2q = {}
    for name, char in (("CZ", "fCZ"), ("ISWAP", "fISWAP"), ("XY", "fXY"), ("CPHASE", "fCPHASE")):
        for n, v in sites(isa["instructions"], name, char):
            f2q[frozenset(n)] = max(v, f2q.get(frozenset(n), 0.0))
    return Calibration(f1q, fro, f2q)


def calibration_from_qcs(device):
    """Fetch the live calibration (read-only; credentials are read by the SDK from ~/.qcs)."""
    from qcs_sdk import QCSClient
    from qcs_sdk.qpu.isa import get_instruction_set_architecture
    return calibration_from_qcs_isa(
        get_instruction_set_architecture(client=QCSClient.load(), quantum_processor_id=device))


def calibration_from_qiskit_target(target, two_qubit_gate=None, one_qubit_gate=None):
    """From a Qiskit ``Target`` (``backend.target``); fidelity = 1 - error."""
    names = set(target.operation_names)
    g2 = two_qubit_gate or next((g for g in ("cz", "ecr", "cx", "iswap", "rzz") if g in names), None)
    g1 = one_qubit_gate or next((g for g in ("sx", "x", "rx") if g in names), None)

    def fid(props):
        return None if props is None or props.error is None else 1.0 - props.error

    f2q, f1q, fro = {}, {}, {}
    for qargs, props in (target[g2].items() if g2 else []):
        v = fid(props)
        if _usable(v):
            f2q[frozenset(qargs)] = max(v, f2q.get(frozenset(qargs), 0.0))
    for qargs, props in (target[g1].items() if g1 else []):
        if _usable(fid(props)):
            f1q[qargs[0]] = fid(props)
    for qargs, props in (target["measure"].items() if "measure" in names else []):
        if _usable(fid(props)):
            fro[qargs[0]] = fid(props)
    return Calibration(f1q, fro, f2q)


# ----------------------------------------------------------------------------- interaction loaders

def interaction_from_edges(spec, measured=None, n_qubits=None):
    """``"0-1,1-3:2,2-3"``: logical pairs with optional two-qubit gate counts."""
    edges = Counter()
    for item in filter(None, spec.split(",")):
        pair, _, cnt = item.partition(":")
        a, b = (int(x) for x in pair.split("-"))
        edges[frozenset((a, b))] += int(cnt or 1)
    qubits = {q for e in edges for q in e} | set(measured or ())
    n = n_qubits or (max(qubits) + 1)
    return Interaction(n, edges, Counter({q: 1 for q in range(n)}), set(measured or ()))


def interaction_from_quil(text):
    edges, one, measured, qubits = Counter(), Counter(), set(), set()
    for line in text.splitlines():
        line = line.split("#")[0].strip()
        if not line or line.startswith(("DECLARE", "PRAGMA", "FENCE", "DELAY", "RESET", "HALT")):
            continue
        tok = re.sub(r"\(.*?\)", "", line).split()
        if tok[0] == "MEASURE":
            measured.add(int(tok[1]))
            qubits.add(int(tok[1]))
            continue
        qs = [int(t) for t in tok[1:] if t.isdigit()]
        qubits.update(qs)
        if len(qs) == 2:
            edges[frozenset(qs)] += 1
        elif len(qs) == 1:
            one[qs[0]] += 1
    labels = sorted(qubits)
    ix = {q: i for i, q in enumerate(labels)}
    return Interaction(len(labels), Counter({frozenset(ix[q] for q in e): n for e, n in edges.items()}),
                       Counter({ix[q]: n for q, n in one.items()}), {ix[q] for q in measured}, labels)


def interaction_from_qiskit(qc):
    edges, one, measured = Counter(), Counter(), set()
    for inst in qc.data:
        qs = [qc.find_bit(q).index for q in inst.qubits]
        name = inst.operation.name
        if name == "measure":
            measured.update(qs)
        elif name in ("barrier", "delay", "reset"):
            continue
        elif len(qs) == 2:
            edges[frozenset(qs)] += 1
        elif len(qs) == 1:
            one[qs[0]] += 1
        elif len(qs) > 2:
            raise ValueError(f"decompose {name} into one- and two-qubit gates before selecting qubits")
    return Interaction(qc.num_qubits, edges, one, measured)


def interaction_from_qasm(path):
    from qiskit import qasm2, qasm3
    text = open(path).read()
    qc = qasm3.loads(text) if "OPENQASM 3" in text else qasm2.loads(text)
    return interaction_from_qiskit(qc)


# ----------------------------------------------------------------------------- search

def _terms(inter):
    """Per logical qubit: the edges to earlier qubits are added when the qubit is placed."""
    adj = defaultdict(set)
    for e in inter.edges:
        a, b = tuple(e)
        adj[a].add(b)
        adj[b].add(a)
    order, seen = [], set()
    for start in sorted(range(inter.n_qubits), key=lambda q: -len(adj[q])):
        if start in seen:
            continue
        queue = [start]
        seen.add(start)
        while queue:                                   # BFS keeps each new qubit next to a placed one
            q = queue.pop(0)
            order.append(q)
            for r in sorted(adj[q], key=lambda r: -len(adj[r])):
                if r not in seen:
                    seen.add(r)
                    queue.append(r)
    return order, adj


def select_layouts(inter, cal, top=5, exclude=()):
    """Return [(log_score, layout)] best first; layout[i] = physical qubit of logical qubit i."""
    nb = cal.neighbors()
    known = set(cal.f1q) if cal.f1q else (set(cal.fro) | set(nb))
    usable = known - set(exclude)
    order, adj = _terms(inter)
    best2 = math.log(max(cal.f2q.values())) if cal.f2q else 0.0
    best1 = math.log(max(cal.f1q.values())) if cal.f1q else 0.0
    bestr = math.log(max(cal.fro.values())) if cal.fro else 0.0

    def local(q, p, placed):
        s = inter.one_qubit[q] * (math.log(cal.f1q[p]) if cal.f1q else 0.0)
        if q in inter.measured:
            if p not in cal.fro:
                return None
            s += math.log(cal.fro[p])
        for r in adj[q]:
            if r in placed:
                e = frozenset((p, placed[r]))
                if e not in cal.f2q:
                    return None
                s += inter.edges[frozenset((q, r))] * math.log(cal.f2q[e])
        return s

    # optimistic value of the terms still to be added after position k
    rest = [0.0] * (len(order) + 1)
    for k in range(len(order) - 1, -1, -1):
        q = order[k]
        later = set(order[:k])
        rest[k] = (rest[k + 1] + inter.one_qubit[q] * best1 + (bestr if q in inter.measured else 0.0)
                   + sum(inter.edges[frozenset((q, r))] for r in adj[q] if r in later) * best2)

    heap = []                                          # min-heap of the current top layouts
    placed, used = {}, set()

    def bound():
        return heap[0][0] if len(heap) == top else -math.inf

    def dfs(k, score):
        if k == len(order):
            item = (score, [placed[q] for q in range(inter.n_qubits)])
            (heapq.heappush if len(heap) < top else heapq.heappushpop)(heap, item)
            return
        q = order[k]
        anchors = [placed[r] for r in adj[q] if r in placed]
        cands = set.intersection(*(nb[a] for a in anchors)) if anchors else usable
        for p in sorted(cands & usable - used):
            s = local(q, p, placed)
            if s is None or score + s + rest[k + 1] <= bound():
                continue
            placed[q] = p
            used.add(p)
            dfs(k + 1, score + s)
            del placed[q]
            used.discard(p)

    dfs(0, 0.0)
    return sorted(heap, reverse=True)


def describe(inter, cal, layout):
    lines = []
    for e, n in sorted(inter.edges.items(), key=lambda x: sorted(x[0])):
        a, b = sorted(e)
        f = cal.f2q[frozenset((layout[a], layout[b]))]
        lines.append(f"    logical {inter.labels[a]}-{inter.labels[b]} -> physical {layout[a]}-{layout[b]}: "
                     f"f2q {f:.4f} x {n}")
    for q in range(inter.n_qubits):
        p = layout[q]
        ro = f" fro {cal.fro[p]:.4f}" if q in inter.measured else ""
        f1 = f" f1q {cal.f1q[p]:.5f}" if p in cal.f1q else ""
        lines.append(f"    logical {inter.labels[q]} -> physical {p}:{f1}{ro}")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--calibration", help="generic calibration JSON")
    src.add_argument("--qcs-isa", help="saved Rigetti QCS instruction-set architecture JSON")
    src.add_argument("--qcs", metavar="DEVICE", help="fetch the live QCS calibration, e.g. Cepheus-1-108Q")
    circ = ap.add_mutually_exclusive_group(required=True)
    circ.add_argument("--edges", help="logical interaction pairs, e.g. 0-1,1-3:2,2-3")
    circ.add_argument("--quil", help="Quil program (logical qubit indices)")
    circ.add_argument("--qasm", help="OpenQASM 2/3 file (needs qiskit)")
    ap.add_argument("--measured", help="measured logical qubits for --edges (default: all)")
    ap.add_argument("--exclude", default="", help="physical qubits to avoid")
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--save-calibration", help="write the fetched QCS calibration to this JSON file")
    a = ap.parse_args()

    if a.calibration:
        cal = calibration_from_json(a.calibration)
    elif a.qcs_isa:
        cal = calibration_from_qcs_isa(open(a.qcs_isa).read())
    else:
        from qcs_sdk import QCSClient
        from qcs_sdk.qpu.isa import get_instruction_set_architecture
        isa = get_instruction_set_architecture(client=QCSClient.load(), quantum_processor_id=a.qcs)
        if a.save_calibration:
            open(a.save_calibration, "w").write(isa.json())
        cal = calibration_from_qcs_isa(isa)

    if a.edges:
        meas = [int(x) for x in a.measured.split(",")] if a.measured else None
        inter = interaction_from_edges(a.edges, meas)
        if meas is None:
            inter.measured = set(range(inter.n_qubits))
    elif a.quil:
        inter = interaction_from_quil(open(a.quil).read())
    else:
        inter = interaction_from_qasm(a.qasm)

    exclude = [int(x) for x in a.exclude.split(",") if x]
    print(f"calibration: {len(cal.f1q | cal.fro)} qubits, {len(cal.f2q)} couplers; circuit: "
          f"{inter.n_qubits} qubits, {sum(inter.edges.values())} two-qubit gates on {len(inter.edges)} pairs, "
          f"{len(inter.measured)} measured")
    res = select_layouts(inter, cal, top=a.top, exclude=exclude)
    if not res:
        print("no SWAP-free layout exists; the circuit needs routing or a smaller interaction graph")
        return
    for rank, (score, layout) in enumerate(res, 1):
        print(f"{rank}. estimated success {math.exp(score):.4f}  --qubits {','.join(map(str, layout))}")
    if inter.labels != list(range(inter.n_qubits)):
        print(f"   (physical qubits listed for the program qubits {','.join(map(str, inter.labels))})")
    print("best layout:\n" + describe(inter, cal, res[0][1]))


if __name__ == "__main__":
    main()
