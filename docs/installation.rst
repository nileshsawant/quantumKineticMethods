Installation
============

The library is pure Python. The circuit ports are emulated with
`Qiskit Aer <https://github.com/Qiskit/qiskit-aer>`_; the classical solver needs only
NumPy and SciPy.

.. code-block:: bash

   pip install numpy scipy matplotlib qiskit qiskit-aer

For the larger circuits (for example the 17-qubit reflecting box) a GPU build of Aer
is much faster:

.. code-block:: bash

   pip install qiskit-aer-gpu     # x86_64 + CUDA

Then clone the repository and add it to your ``PYTHONPATH`` (or run from its root):

.. code-block:: bash

   git clone https://github.com/nileshsawant/quantumKineticMethods.git
   cd quantumKineticMethods

Building the documentation
--------------------------

.. code-block:: bash

   pip install -r docs/requirements.txt
   sphinx-build -b html docs docs/_build/html

The autodoc build mocks the heavy scientific stack (NumPy, SciPy, Qiskit, ...), so the
documentation builds without a GPU or a Qiskit installation.
