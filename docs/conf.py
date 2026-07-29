# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# Make the package importable for autodoc (repo root is one level up).
sys.path.insert(0, os.path.abspath(".."))

# -- Project information -----------------------------------------------------
project = "quantumKineticMethods"
author = "Nilesh Sawant"
copyright = "2026, Nilesh Sawant"
release = "0.1"
version = "0.1"

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",     # pull docstrings from the modules
    "sphinx.ext.napoleon",    # understand NumPy / Google style docstrings
    "sphinx.ext.mathjax",     # render LaTeX math
    "sphinx.ext.viewcode",    # add [source] links
    "sphinx.ext.intersphinx", # cross-link to Python / NumPy docs
]

# autodoc imports the library to read its docstrings. NumPy and SciPy are lightweight
# and are used at import time (module-level matrices in operators.py), so they are
# installed for the build; only the heavy quantum / plotting stack is mocked, so the
# docs still build anywhere (e.g. Read the Docs) without a GPU or qiskit.
autodoc_mock_imports = [
    "matplotlib",
    "qiskit",
    "qiskit_aer",
    "cupy",
]
autodoc_member_order = "bysource"
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}

napoleon_numpy_docstring = True
napoleon_google_docstring = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- HTML output -------------------------------------------------------------
# Use the Read the Docs theme if available, else fall back to the built-in one.
try:
    import sphinx_rtd_theme  # noqa: F401

    html_theme = "sphinx_rtd_theme"
except ImportError:  # pragma: no cover
    html_theme = "alabaster"

html_static_path = ["_static"]
html_title = "quantumKineticMethods"
