import os


# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------

project = "Deemian"
copyright = "2024, Muhammad Radifar"
author = "Muhammad Radifar"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ["myst_parser", "sphinxcontrib.video", "sphinx_design"]
myst_enable_extensions = ["colon_fence"]
# autoapi_dirs = ["../src"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
#
# Additional note: when using venv to build documentation venv or .venv
# should be excluded or error messages occur
# https://github.com/sphinx-doc/sphinx/issues/2066#issuecomment-474587560
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "venv", ".venv"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "furo"

html_theme_options = {
    "sidebar_hide_name": True,
    "light_logo": "Deemian_logo_web.png",
    "dark_logo": "Deemian_logo_web.png",
    "announcement": "<strong style='color:red;\
        '>Attention!</strong> Both Deemian and Deemian Viewer are currently under heavy construction.\
        Please use at your own risk.",
}

html_static_path = ["_static"]

if os.getenv("READTHEDOCS"):
    extensions.append("sphinxcontrib.googleanalytics")
    googleanalytics_id = "G-GV5XETXYHY"
