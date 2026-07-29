The scheme and its circuits
===========================

The state is a four-component Dirac spinor :math:`\psi(\mathbf{x},t)` on a lattice.
Following Dellar (2011), the scheme uses the **Majorana form** of the Dirac equation,
in which the three matrices multiplying the spatial derivatives are real, so each
spatial gradient becomes a :math:`\pm 1` lattice shift after a fixed spinor rotation.

One time step is split into per-axis substeps. Each substep

.. math::

   \psi \;\longmapsto\; \mathbf{R}_a\,\mathrm{Stream}_a\!\left[\,\hat{Q}\,
   \mathbf{R}_a^{-1}\psi\,\right]

rotates into the characteristic frame, applies the :math:`\mathrm{SU}(2)` collision
:math:`\hat{Q}` that carries the mass and the potential, streams by one site, and
rotates back. Every factor is unitary, so the discrete :math:`\ell_2` norm is
conserved exactly -- which is precisely what a quantum circuit on :math:`n` qubits
implements.

Amplitude encoding
------------------

The four spinor components live in **2 shared spinor qubits**; each lattice axis of
:math:`N` sites is addressed by :math:`\log_2 N` **position qubits**. Writing the
spinor bits as :math:`q_0,q_1` and the position bits as :math:`p_0,\dots,p_{n-1}`,

.. math::

   c = q_0 + 2 q_1 \in \{0,1,2,3\}, \qquad
   x = \sum_{j=0}^{n-1} 2^{\,j} p_j, \qquad
   \text{basis index } i = 4x + c .

A :math:`32^3` lattice (32,768 sites) is therefore represented in
:math:`2 + 5 + 5 + 5 = 17` qubits. Implemented in :mod:`qlb_port.operators`.

The unit operations as gates
----------------------------

**Rotations and collision.** These are fixed :math:`4\times 4` unitaries on the two
spinor qubits. They are compiled with the exact Cartan (KAK) decomposition: the
rotation costs one CX gate, the collision two. See :mod:`qlb_port.port` and
:mod:`qlb_port.sweep`.

**Streaming is** :math:`+1` **on the position register.** Moving every amplitude one
site is the map :math:`|x\rangle \mapsto |x+1\rangle`, i.e. binary "add one" -- a
ripple of multi-controlled-X gates on the address bits, applied to all sites at once.
The direction qubit :math:`q_1` controls whether the register is incremented or
decremented. See :mod:`qlb_port.streaming`.

**Bounce-back is** :math:`+1` **on a folded ring.** Folding the direction qubit into
the position register as its most significant bit turns a hard-wall reflection into a
single plain increment: interior movers advance one site, and a mover that reaches a
wall crosses the fold and returns with its direction reversed. One increment does
advection and reflection at once. See
:func:`qlb_port.streaming.reflecting_streaming_circuit`.

**Position-dependent potential.** In the massless case the collision is a pure phase,
so the potential is a diagonal **phase oracle**
:math:`\mathcal{O}_V = \mathrm{diag}(\hat{a}(0),\dots,\hat{a}(N-1))` on the position
register. In the massive case it is a **position-multiplexed collision**: a uniform
vacuum collision on the spinor followed by a position-controlled correction that acts
only where the potential is nonzero. See :mod:`qlb_port.potential`.

Assembling a time step
----------------------

Single-axis sweeps (:mod:`qlb_port.sweep`) are tensored -- one position register per
axis on a shared spinor register -- into two- and three-dimensional time steps
(:mod:`qlb_port.twod`, :mod:`qlb_port.threed`). Because the position registers are
disjoint, the cost of a :math:`D`-dimensional step is the sum of the :math:`D`
single-axis sweep costs plus the constant spinor gates.

The porting harness
-------------------

All of the above is produced and checked by a single primitive
(:mod:`qlb_port.port`): given a target unitary :math:`U`, compile it to a chosen gate
set and *verify* that the compiled circuit reproduces :math:`U` on the emulator,
returning the phase-invariant fidelity
:math:`F = |\mathrm{tr}(U^\dagger U_\mathrm{circ})| / \dim`, the two-qubit gate count,
and the depth.
