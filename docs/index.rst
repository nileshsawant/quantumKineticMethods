quantumKineticMethods
=====================

**Exact quantum circuits for the quantum lattice Boltzmann (QLB) method for the
Dirac equation.**

This project ports the three-dimensional Succi--Dellar Dirac QLB scheme, operation
by operation, to *exact* quantum circuits on qubits, and verifies on a state-vector
emulator that the circuits reproduce the classical QLB solver to machine precision.

A QLB time step is a fixed sequence of *unitary* operations -- a basis rotation, a
collision, a streaming shift, and the inverse rotation -- so it maps naturally onto a
quantum circuit. Here that mapping is made explicit and checked, layer by layer,
against a classical reference.

.. note::

   No claim of computational advantage is made. The contribution is *exact
   representability*, together with measured gate counts and a validated, reusable
   set of circuit primitives.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   installation
   scheme
   usage
   validation
   api/index

Indices and tables
-------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
