"""
qlb_port — port Quantum Lattice Boltzmann (QLB) operators to emulated quantum circuits.
======================================================================================

A small, self-contained toolkit for compiling the unit operators of the QLB Dirac
solver into gate circuits and verifying them against their classical matrices on a
qiskit-aer (GPU) backend.

Layout
------
  operators.py : the QLB unit operators as NumPy matrices (single source of truth,
                 kept consistent with ``dirac_qlb_solver`` via a test).
  backend.py   : qiskit-aer helpers (GPU-first) — circuit unitary, statevector apply,
                 fidelity metrics.
  port.py      : the general "unitary -> verified circuit" harness plus convenience
                 builders for the QLB rotation and collision operators.
  streaming.py : the QLB +/-1 lattice shift as a controlled increment/decrement on a
                 position register.
  test_port.py : validation — every ported operator must match its target to ~1e-9.

Environment
-----------
Requires the HPC module ``qiskit/aer-gpu`` (qiskit-aer-gpu).  Load it with::

    module load qiskit/aer-gpu

Quick start
-----------
    from qlb_port import operators as ops
    from qlb_port.port import port_and_verify

    r = port_and_verify(ops.X_ROTATION, label="x-rotation")
    print(r["fidelity"], r["cx"], r["depth"])
"""

from . import operators, backend, port, streaming, sweep

__all__ = ["operators", "backend", "port", "streaming", "sweep"]
