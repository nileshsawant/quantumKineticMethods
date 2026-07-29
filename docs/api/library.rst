Circuit library (``qlb_port``)
==============================

.. automodule:: qlb_port
   :no-members:

Operators
---------

The QLB unit operators as NumPy matrices: the Dirac / Pauli matrices, the fixed
rotations, the collision, and the classical streaming permutation used as a reference.

.. automodule:: qlb_port.operators
   :members:

Backend
-------

Qiskit Aer helpers (GPU-first): circuit unitary, state-vector application, and fidelity
metrics.

.. automodule:: qlb_port.backend
   :members:

Porting harness
---------------

The general "unitary :math:`\rightarrow` verified circuit" primitive, plus convenience
builders for the QLB rotation and collision operators.

.. automodule:: qlb_port.port
   :members:

Streaming
---------

The QLB :math:`\pm 1` lattice shift as a controlled increment / decrement of a position
register, and the reflecting (bounce-back) streaming as a folded-ring increment.

.. automodule:: qlb_port.streaming
   :members:

Single-axis sweep
-----------------

Assemble rotation, collision, and streaming into one full single-axis substep.

.. automodule:: qlb_port.sweep
   :members:

Position-dependent potential
----------------------------

A position-dependent potential as a phase oracle (massless) or a position-multiplexed
collision (massive).

.. automodule:: qlb_port.potential
   :members:

Two dimensions
--------------

Compose x- and y-sweeps on two position registers.

.. automodule:: qlb_port.twod
   :members:

Three dimensions
----------------

Compose x-, y-, and z-sweeps on three position registers.

.. automodule:: qlb_port.threed
   :members:
