# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'FastChem'
copyright = '2019 - 2025, Daniel Kitzmann, Joachim Stock'
author = 'Daniel Kitzmann'
release = '4.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ['static']

html_theme_options = {
	'source_repository': 'https://github.com/newstrangeworlds/FastChem/',
	'source_branch': 'master',
	'source_directory': 'docs_src/',
}

html_logo = 'static/fastchem_logo_dark.png'
