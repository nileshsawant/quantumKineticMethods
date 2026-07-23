"""
Streaming as quantum circuits: the QLB per-component +/-1 lattice shift.
=======================================================================

The streaming sub-step shifts each spinor component by +/-1 site along the sweep
axis.  With the spinor encoded in qubits 0,1 the shift sign depends only on spinor
qubit 1 (components 0,1 have qubit1=0; components 2,3 have qubit1=1), so streaming is
a **controlled increment/decrement** of the position register:

    axis x, z :  qubit1 = 1  ->  +1 (increment) ,  qubit1 = 0  ->  -1 (decrement)
    axis y    :  qubit1 = 1  ->  -1 (decrement) ,  qubit1 = 0  ->  +1 (increment)

The increment adds 1 (mod 2**n_pos) with a ripple cascade of multi-controlled X
gates; the decrement is its inverse.  This realises **periodic** streaming (the
wrap-around matches the periodic y/z sweeps and the bulk of the x-sweep).  The open /
bounce-back x boundary is a separate boundary operator (not included here).

Qubit layout (matches operators.streaming_reference):
    spinor   = qubits 0,1   (qubit 0 = LSB)
    position = qubits 2..(1+n_pos)
"""

from qiskit import QuantumCircuit

from . import operators as ops

SPINOR_QUBITS = (0, 1)


def increment_gate(n_pos):
    """Gate that adds 1 (mod 2**n_pos) to an n_pos-qubit register (qubit 0 = LSB)."""
    qc = QuantumCircuit(n_pos, name="incr")
    # MSB first: lower qubits are still the pre-update value when used as controls.
    for i in range(n_pos - 1, 0, -1):
        qc.mcx(list(range(i)), i)
    qc.x(0)
    return qc.to_gate()


def decrement_gate(n_pos):
    """Gate that subtracts 1 (mod 2**n_pos); inverse of :func:`increment_gate`."""
    return increment_gate(n_pos).inverse()


def streaming_circuit(axis, n_pos):
    """
    Controlled increment/decrement circuit implementing the streaming sub-step.

    Parameters
    ----------
    axis  : 'x', 'y', or 'z'
    n_pos : number of position qubits (lattice size N = 2**n_pos)

    Returns
    -------
    QuantumCircuit on (2 + n_pos) qubits: spinor = 0,1 ; position = 2..(1+n_pos).
    """
    signs = ops.streaming_signs(axis)
    n = 2 + n_pos
    qc = QuantumCircuit(n, name=f"stream_{axis}")
    pos = list(range(2, n))

    inc = increment_gate(n_pos)
    dec = inc.inverse()

    # signs[0]==signs[1] (qubit1=0) and signs[2]==signs[3] (qubit1=1) by construction.
    g_hi = inc if signs[2] > 0 else dec      # qubit1 = 1
    g_lo = inc if signs[0] > 0 else dec      # qubit1 = 0

    qc.append(g_hi.control(1, ctrl_state="1"), [1] + pos)
    qc.append(g_lo.control(1, ctrl_state="0"), [1] + pos)
    return qc
