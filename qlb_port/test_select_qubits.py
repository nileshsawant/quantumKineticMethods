"""
Checks of qlb_port.select_qubits against brute force on small random devices.

Usage:
    PYTHONPATH=. python3 -m qlb_port.test_select_qubits
"""

import itertools
import math
import random
from collections import Counter

from qlb_port.select_qubits import (Calibration, Interaction, interaction_from_edges,
                                    interaction_from_quil, select_layouts)


def random_grid(rows, cols, rng, dead=()):
    q = lambda r, c: r * cols + c
    f2q = {}
    for r in range(rows):
        for c in range(cols):
            for dr, dc in ((0, 1), (1, 0)):
                if r + dr < rows and c + dc < cols:
                    f2q[frozenset((q(r, c), q(r + dr, c + dc)))] = rng.uniform(0.95, 0.997)
    n = rows * cols
    cal = Calibration({i: rng.uniform(0.995, 0.9995) for i in range(n) if i not in dead},
                      {i: rng.uniform(0.90, 0.99) for i in range(n) if i not in dead},
                      {e: f for e, f in f2q.items() if not e & set(dead)})
    return cal


def brute_force(inter, cal):
    best = []
    for perm in itertools.permutations(sorted(cal.f1q), inter.n_qubits):
        s = 0.0
        ok = True
        for e, n in inter.edges.items():
            a, b = tuple(e)
            pe = frozenset((perm[a], perm[b]))
            if pe not in cal.f2q:
                ok = False
                break
            s += n * math.log(cal.f2q[pe])
        if not ok:
            continue
        s += sum(inter.one_qubit[q] * math.log(cal.f1q[perm[q]]) for q in range(inter.n_qubits))
        s += sum(math.log(cal.fro[perm[q]]) for q in inter.measured)
        best.append((s, list(perm)))
    return sorted(best, reverse=True)


def test_matches_brute_force():
    rng = random.Random(1)
    cases = [
        interaction_from_edges("0-1,1-3,2-3", measured=[0, 1, 2, 3]),     # chain, as the 4-site QLB circuit
        interaction_from_edges("0-1:3,1-2:2", measured=[1]),              # weighted chain, one readout
        interaction_from_edges("0-1,1-2,2-3,3-0", measured=[0, 2]),       # 4-cycle
        interaction_from_edges("0-1,2-3", measured=[0, 1, 2, 3]),         # two disconnected pairs
        interaction_from_edges("0-1,0-2,0-3", measured=[1, 2, 3]),        # star
    ]
    for inter in cases:
        cal = random_grid(3, 3, rng, dead=(4,))
        got = select_layouts(inter, cal, top=3)
        ref = brute_force(inter, cal)[:3]
        assert [round(s, 12) for s, _ in got] == [round(s, 12) for s, _ in ref], (got, ref)
    print("branch-and-bound matches brute force on 5 interaction graphs")


def test_no_layout_and_quil_parse():
    rng = random.Random(2)
    cal = random_grid(2, 3, rng)
    triangle = interaction_from_edges("0-1,1-2,0-2")
    assert select_layouts(triangle, cal) == [], "a triangle cannot embed in a square grid"
    inter = interaction_from_quil("DECLARE ro BIT[2]\nRX(pi/2) 0\nCZ 0 1\nRZ(0.3) 1\nCZ 1 0\nMEASURE 1 ro[0]\n")
    assert inter.edges == Counter({frozenset((0, 1)): 2}) and inter.measured == {1}
    assert inter.one_qubit == Counter({0: 1, 1: 1}) and isinstance(inter, Interaction)
    print("triangle correctly rejected; Quil interaction graph parsed")


if __name__ == "__main__":
    test_matches_brute_force()
    test_no_layout_and_quil_parse()
