Usage
=====

Port and verify a unit operator
-------------------------------

Every unit operator can be compiled to a gate set and checked against its classical
matrix in one call:

.. code-block:: python

   from qlb_port import operators as ops
   from qlb_port.port import port_and_verify

   r = port_and_verify(ops.X_ROTATION, label="x-rotation")
   print(r["fidelity"], r["cx"], r["depth"])   # -> 1.0  1  4

Build a streaming circuit
-------------------------

.. code-block:: python

   from qlb_port import streaming as st

   inc  = st.streaming_circuit("x", n_pos=5)                 # periodic +/-1 shift
   refl = st.reflecting_streaming_circuit("x", n_pos=5)      # bounce-back (folded ring)

A position-dependent potential
------------------------------

.. code-block:: python

   from qlb_port import potential as pot

   V = pot.impurity_field(n_pos=4, positions=[6, 7, 8, 9], g_value=0.5)  # a barrier
   oracle = pot.potential_oracle(V)          # massless: a diagonal phase oracle

Run the checks and regenerate the figures
-----------------------------------------

.. code-block:: bash

   # every operator and composition reproduces the classical scheme (fast, CPU)
   python -m qlb_port.test_port

   # physics validation of the classical solver against the exact Dirac equation
   python test_qlb_validation.py

   # regenerate the 1D / 2D / 3D validation figures
   python -m qlb_port.plot_validation
   python -m qlb_port.plot_validation_2d
   python -m qlb_port.plot_validation_3d

The larger circuits are emulated fastest on a GPU through ``qiskit-aer-gpu``.
