"""
Validate the ancilla (ripple-carry) streaming by emulation, and measure its
two-qubit gate reduction versus the multi-controlled-X (MCX) cascade.
================================================================================

Three checks, all by state-vector emulation (no hardware, no Aer/GPU needed):

  [1] the bare ripple-carry incrementer adds 1 (mod 2**n) and returns the carry
      ancilla to |0> for every input;
  [2] the ancilla streaming circuit reproduces the classical +/-1 shift
      permutation ``operators.streaming_reference`` on the spinor-position
      register, again leaving the ancilla clean;
  [3] the transpiled two-qubit (CX) gate count of a streaming step, MCX cascade
      vs ripple-carry ancilla, over a range of n_pos.

Run:
    module load qiskit/aer-gpu
    PYTHONPATH=. python3 -u -m qlb_port.validate_ancilla_streaming
"""

import numpy as np
from qiskit import transpile
from qiskit.quantum_info import Statevector

from . import operators as ops
from . import streaming as st

_BASIS = ["rz", "ry", "rx", "cx"]


def _perm(P):
    """perm[j] = row index of the single 1 in column j of permutation matrix P."""
    return np.argmax(np.abs(P), axis=0)


def _acts_as(qc, perm, n_reg):
    """True iff qc maps |r>|0>_anc -> |perm[r]>|0>_anc for every register basis
    state r, i.e. it realizes the permutation and leaves the ancilla in |0>."""
    dim_reg = 2 ** n_reg
    dim_tot = 2 ** qc.num_qubits
    for r in range(dim_reg):
        out = Statevector.from_int(r, dim_tot).evolve(qc).data
        exp = int(perm[r])                       # ancilla bits are 0 -> full index = exp
        if not np.isclose(abs(out[exp]), 1.0, atol=1e-9):
            return False
        if not np.allclose(np.delete(out, exp), 0.0, atol=1e-9):
            return False
    return True


def _count_2q(qc):
    t = transpile(qc, basis_gates=_BASIS, optimization_level=3)
    return sum(v for g, v in t.count_ops().items() if g == "cx")


def main():
    print("=" * 72)
    print("Ancilla ripple-carry streaming: emulation validation")
    print("=" * 72)

    ok_all = True

    # [1] incrementer: +1 mod 2^n on the register, carry ancilla restored
    print("\n[1] ripple-carry incrementer adds 1 (mod 2^n), ancilla clean:")
    for n in range(1, 6):
        qc = st.ancilla_increment_circuit(n)
        N = 2 ** n
        perm = np.array([(r + 1) % N for r in range(N)])
        ok = _acts_as(qc, perm, n)
        ok_all &= ok
        print(f"    n_pos={n}: {'PASS' if ok else 'FAIL'}  "
              f"({qc.num_qubits} qubits = {n} register + {qc.num_qubits - n} ancilla)")

    # [2] streaming vs the classical +/-1 shift permutation
    print("\n[2] ancilla streaming == classical streaming_reference permutation:")
    for axis in ("x", "y", "z"):
        for n_pos in range(1, 5):
            qc = st.ancilla_streaming_circuit(axis, n_pos)
            P = ops.streaming_reference(axis, n_pos)
            ok = _acts_as(qc, _perm(P), 2 + n_pos)
            ok_all &= ok
            print(f"    axis {axis}, n_pos={n_pos}: {'PASS' if ok else 'FAIL'}")

    # [3] two-qubit gate count: MCX cascade vs ripple-carry ancilla
    print("\n[3] two-qubit (CX) count of one streaming step (opt-level 3):")
    print("     n_pos | MCX cascade | ancilla ripple | reduction")
    print("     ------+-------------+----------------+----------")
    for n_pos in range(2, 7):
        mcx = _count_2q(st.streaming_circuit("x", n_pos))
        anc = _count_2q(st.ancilla_streaming_circuit("x", n_pos))
        print(f"     {n_pos:5d} | {mcx:11d} | {anc:14d} | {mcx / max(anc, 1):6.2f}x")

    print("\n" + ("ALL CHECKS PASSED" if ok_all else "SOME CHECKS FAILED"))
    return 0 if ok_all else 1


if __name__ == "__main__":
    raise SystemExit(main())
