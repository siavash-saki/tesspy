# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import shutil
import sys
from pathlib import Path

sys.path.insert(0, os.path.abspath("../.."))

# Ensure pandoc from pypandoc_binary is on PATH (needed by nbsphinx)
try:
    import pypandoc

    _pandoc_dir = str(Path(pypandoc.get_pandoc_path()).parent)
    if _pandoc_dir not in os.environ.get("PATH", ""):
        os.environ["PATH"] = _pandoc_dir + os.pathsep + os.environ.get("PATH", "")
except ImportError:
    pass  # pandoc must be installed system-wide

from tesspy._version import __version__ as tesspy_version  # noqa: E402

# -- Sync example notebooks from Examples/ at build time ---------------------
# The single source of truth is Examples/ at the repo root.
# Sphinx requires source files inside its source tree, so we copy them here.
# The copies are .gitignored and never committed.

_examples_src = Path(__file__).resolve().parent.parent.parent / "Examples"
_examples_dst = Path(__file__).resolve().parent / "example_notebooks"
if _examples_src.exists():
    _examples_dst.mkdir(exist_ok=True)
    for _nb in _examples_src.glob("*.ipynb"):
        shutil.copy2(_nb, _examples_dst / _nb.name)

autodoc_mock_imports = [
    "numpy",
    "pandas",
    "geopandas",
    "shapely",
    "sklearn",
    "scipy",
    "hdbscan",
    "h3",
    "mercantile",
    "osmnx",
]

# -- Project information -----------------------------------------------------

project = "tesspy"
copyright = "2022, tesspy developers"
author = "tesspy developers"
release = tesspy_version

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "myst_parser",
    "nbsphinx",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

master_doc = "index"

templates_path = ["_templates"]
exclude_patterns = []

# -- nbsphinx configuration --------------------------------------------------

nbsphinx_execute = "never"

# -- MyST configuration ------------------------------------------------------

myst_enable_extensions = [
    "colon_fence",
    "deflist",
]

# -- Options for HTML output -------------------------------------------------

html_theme = "furo"
html_static_path = ["_static"]
