Classical solver (``dirac_qlb_solver``)
=======================================

``dirac_qlb_solver.py`` is the classical Quantum Lattice Boltzmann solver for the Dirac
equation that provides the reference every circuit is checked against. It is a runnable
script rather than an import-only library: importing it configures a simulation from
its module-level parameters and prints a summary, so it is documented here descriptively
rather than through autodoc.

.. note::

   The circuit library reproduces the *operations* of this solver exactly. The
   operators in :mod:`qlb_port.operators` are kept consistent with
   ``dirac_qlb_solver`` by a test (:mod:`qlb_port.test_port`).

What it does
------------

The solver advances a four-component Dirac spinor field on a regular lattice by the QLB
substep loop -- rotate into the characteristic frame, collide, stream by one site,
rotate back -- with the mass and any scalar potential carried by the
:math:`\mathrm{SU}(2)` collision. It implements periodic and reflecting (bounce-back)
boundaries and a position-dependent potential, and reproduces graphene Klein-tunnelling
physics in the massless limit.

Key module-level parameters
---------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Name
     - Meaning
   * - ``NX, NY, NZ``
     - lattice sizes per axis (set ``NZ = 1`` for a 2D problem)
   * - ``LX, LY, LZ``
     - physical box lengths
   * - ``DX, DY, DZ``
     - lattice spacings (``LX / NX``, ...)
   * - ``DT``
     - time step, fixed to the light-cone condition ``DX / C``
   * - ``M_PARTICLE``
     - particle rest mass (``0`` for massless graphene fermions)
   * - ``OMEGA_C``
     - Compton frequency :math:`m c^2 / \hbar`
   * - ``T_STEPS``
     - number of time steps to advance
   * - ``HBAR, C, Q_ELECTRON``
     - physical constants

Running it
----------

.. code-block:: bash

   python dirac_qlb_solver.py

Edit the parameters at the top of the file to change the lattice, mass, potential, or
number of steps. For the ported-circuit equivalents and the validation harness, see
:doc:`../usage` and :mod:`qlb_port`.

Solver tests
------------

``test_dirac_qlb_solver.py`` and ``test_qlb_validation.py`` (at the repository root)
validate the solver: matrix properties, massless advection, the dispersion relation,
convergence, and a uniform-potential global-phase check.
