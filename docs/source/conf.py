# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usa../inputs/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usa../inputs/configuration.html#project-information

project = 'pomsimulator'
copyright = '2024, Enric Petrus'
author = 'Enric Petrus'
release = '1.0'


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usa../inputs/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc']
templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usa../inputs/configuration.html#options-for-html-output

#html_title = 'POMSimulator documentation'
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_logo = "../.img/pomsimulator_logo_reduced.png"
#html_theme_options = {
#    'logo_only': True,
#    'display_version': False,
#}