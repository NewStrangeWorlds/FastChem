# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'FastChem'
copyright = '2019 - 2023, Daniel Kitzmann, Joachim Stock'
author = 'Daniel Kitzmann'
release = '3.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['static']


# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
html_theme_options = {
	'logo': 'logo_l.png',
	'description': 'FastChem - Ultrafast Equilibrium Chemistry',
	'github_button': 'true',
	'github_user': 'exoclime',
    'github_repo': 'fastchem',
	'github_type': 'watch',
	'page_width' : '1001px',
	}
