# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'xraybinaryorbit'
copyright = '2024, LAEX'
author = 'LAEX'
release = '0.2.9'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# Enable Sphinx extensions
extensions = [
    'sphinx.ext.autodoc',              # Automatically document Python modules
    'sphinx.ext.autosummary',          # Generate summary tables for modules
    'sphinx_autodoc_typehints',        # Add type hints to the documentation
    'sphinx_rtd_theme',                # Use the Read the Docs theme
]

# Paths that contain templates, relative to this directory.
templates_path = ['_templates']

# Patterns to ignore when looking for source files
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The default language for code highlighting (Python in this case)
highlight_language = 'python'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# The theme to use for HTML and HTML Help pages.
html_theme = 'sphinx_rtd_theme'  # Changed to Read the Docs theme for a modern look

# Paths that contain custom static files (such as style sheets), relative to this directory.
html_static_path = ['_static']

# Add any extra paths that contain custom JavaScript or CSS files (e.g., stylesheets).
html_css_files = [
    'custom.css',  # Add your custom CSS file if needed
]

# If true, the reST files will be treated as reStructuredText files,
# enabling support for inline markup.
# html_use_index = True  # Uncomment if you want to enable indexing
