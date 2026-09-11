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


def reflecting_streaming_circuit(axis, n_pos):
    """
    Reflecting (bounce-back) streaming: hard walls at both ends of the axis.

    A reflecting wall is a *periodic shift on the folded 2N-site ring*
    R0 -> R1 -> ... -> R_{N-1} -> L_{N-1} -> ... -> L0 -> R0, which is a single
    2N-cycle and therefore unitary.  We realise it by unfolding the (direction,
    position) pair into an (n_pos+1)-bit coordinate p (p = x for +movers, p = 2N-1-x
    for -movers) and applying a plain +1 (mod 2^{n_pos+1}) increment: every mover
    advances by one, and a mover that reaches a wall crosses the fold and comes back
    with reversed direction (its spinor direction bit q1 flips, spin q0 preserved).

    Qubit layout: spinor = qubits 0,1 (q1 = direction); position = qubits 2..(1+n_pos).
    """
    signs = ops.streaming_signs(axis)
    plus_is_q1_one = signs[2] > 0            # is q1=1 the +mover for this axis?
    n = 2 + n_pos
    q1 = 1
    pos = list(range(2, n))
    qc = QuantumCircuit(n, name=f"reflstream_{axis}")

    if not plus_is_q1_one:                   # y-axis: make q1=1 the +mover
        qc.x(q1)
    # unfold: complement the position register for the -mover (q1=0), then set the
    # fold flag as the MSB (flag = NOT q1)
    qc.x(q1)
    for q in pos:
        qc.cx(q1, q)
    qc.x(q1)
    qc.x(q1)
    # +1 (mod 2^{n_pos+1}) on the register [pos (LSB..), q1 (MSB)]
    reg = pos + [q1]
    for i in range(len(reg) - 1, 0, -1):
        qc.mcx(reg[:i], reg[i])
    qc.x(reg[0])
    # fold back
    qc.x(q1)
    qc.x(q1)
    for q in pos:
        qc.cx(q1, q)
    qc.x(q1)
    if not plus_is_q1_one:
        qc.x(q1)
    return qc


# ================================================================================
# Ancilla-based streaming (ripple-carry incrementer)
# ================================================================================
# The increment above is a cascade of multi-controlled X gates, whose ancilla-free
# synthesis grows super-linearly with n_pos.  A ripple-carry adder instead computes
# the carry chain once into a register of clean carry ancilla, so a controlled +/-1
# shift costs O(n_pos) Toffoli gates.  The carry chain is uncomputed, so the ancilla
# start and end in |0>.  In the controlled version the carry compute/uncompute stay
# uncontrolled (they cancel when the control is off) and only the O(n_pos) bit-flips
# are gated, which is what keeps the count linear.
#
# ancilla_streaming_circuit(axis, n_pos) lives on 2 + n_pos + (n_pos-1) qubits:
#     spinor   = qubits 0,1
#     position = qubits 2..(1+n_pos)
#     carry    = qubits (2+n_pos)..(2*n_pos)   (clean before and after)

def n_stream_ancilla(n_pos):
    """Carry-ancilla count used by the ripple-carry streaming (n_pos - 1)."""
    return max(n_pos - 1, 0)


def _ripple_increment(qc, p, a, ctrl=None):
    """Append +1 (mod 2**len(p)) on register p (p[0]=LSB) with carry ancilla a
    (len a = len p - 1).  If ctrl is given only the bit-flips are controlled by it
    (state |1>); the carry chain is uncomputed either way, leaving a in |0>."""
    n = len(p)
    if n == 1:                                    # +1 mod 2 is a single flip
        if ctrl is None:
            qc.x(p[0])
        else:
            qc.cx(ctrl, p[0])
        return
    qc.cx(p[0], a[0])                             # a[0] = carry_1 = p[0]
    for k in range(2, n):                         # a[k-1] = carry_k = p[k-1] & carry_{k-1}
        qc.ccx(p[k - 1], a[k - 2], a[k - 1])
    for k in range(n - 1, 1, -1):                 # flip high bits, then uncompute carries
        if ctrl is None:
            qc.cx(a[k - 1], p[k])
        else:
            qc.ccx(ctrl, a[k - 1], p[k])
        qc.ccx(p[k - 1], a[k - 2], a[k - 1])
    if ctrl is None:
        qc.cx(a[0], p[1])
    else:
        qc.ccx(ctrl, a[0], p[1])
    qc.cx(p[0], a[0])                             # uncompute a[0]
    if ctrl is None:
        qc.x(p[0])
    else:
        qc.cx(ctrl, p[0])


def _ripple_decrement(qc, p, a, ctrl=None):
    """Append -1 (mod 2**len(p)) as X-conjugated increment (complement, +1,
    complement); the outer complements cancel when ctrl is off."""
    for q in p:
        qc.x(q)
    _ripple_increment(qc, p, a, ctrl)
    for q in p:
        qc.x(q)


def _add_ctrl_shift(qc, sign, ctrl, ctrl_state, pos, anc):
    """Append a +1 (sign>0) or -1 (sign<0) shift of pos, gated on ctrl==ctrl_state."""
    if ctrl_state == 0:
        qc.x(ctrl)
    if sign > 0:
        _ripple_increment(qc, pos, anc, ctrl)
    else:
        _ripple_decrement(qc, pos, anc, ctrl)
    if ctrl_state == 0:
        qc.x(ctrl)


def ancilla_increment_circuit(n_pos):
    """+1 (mod 2**n_pos) on qubits 0..n_pos-1 via a ripple-carry chain on the carry
    ancilla qubits n_pos..(2*n_pos-2), which are restored to |0>."""
    n_anc = n_stream_ancilla(n_pos)
    qc = QuantumCircuit(n_pos + n_anc, name="incr_anc")
    _ripple_increment(qc, list(range(n_pos)), list(range(n_pos, n_pos + n_anc)))
    return qc


def ancilla_streaming_circuit(axis, n_pos):
    """Ancilla (ripple-carry) version of :func:`streaming_circuit`: the same
    controlled +/-1 shift, but O(n_pos) Toffoli gates instead of the MCX cascade.

    Returns a QuantumCircuit on 2 + n_pos + (n_pos-1) qubits (spinor 0,1; position
    2..1+n_pos; carry ancilla 2+n_pos.., clean before and after).
    """
    signs = ops.streaming_signs(axis)
    n_anc = n_stream_ancilla(n_pos)
    n_total = 2 + n_pos + n_anc
    qc = QuantumCircuit(n_total, name=f"stream_anc_{axis}")
    q1 = 1
    pos = list(range(2, 2 + n_pos))
    anc = list(range(2 + n_pos, n_total))
    _add_ctrl_shift(qc, int(signs[2]), q1, 1, pos, anc)   # qubit1 = 1 components
    _add_ctrl_shift(qc, int(signs[0]), q1, 0, pos, anc)   # qubit1 = 0 components
    return qc
