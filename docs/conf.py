import os
import sys

# Configure the path to find your package
sys.path.insert(0, os.path.abspath('../'))

# -- Project Information ------------------------------------------------
project = 'xraybinaryorbit'
copyright = '2024, LAEX'
author = 'LAEX'
release = '0.2'

# -- General Configuration ----------------------------------------------
extensions = [
    'sphinx.ext.autodoc',     # Automatically document your code
    'sphinx.ext.napoleon',    # Support for Google style docstrings
    'recommonmark',           # Support for Markdown files
    'sphinx.ext.viewcode',    # Add links to the source code
    'sphinx.ext.linkcode',    # Link to GitHub source files
]

templates_path = ['_templates']
exclude_patterns = []

# Function to link specific code sections to GitHub
def linkcode_resolve(domain, info):
    if domain != 'py':
        return None
    if not info['module']:
        return None
    filename = info['module'].replace('.', '/')
    # Customize the URL below to match your GitHub repository structure
    return f"https://github.com/xragua/xraybinaryorbit/blob/main/{filename}.py"

# -- Options for HTML Output --------------------------------------------
html_theme = 'sphinx_rtd_theme'  # Using the Read the Docs theme

html_theme_options = {
    'canonical_url': '',
    'analytics_id': '',  # Optional: Your Google Analytics ID
    'logo_only': False,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'vcs_pageview_mode': '',
    
    # Link to your GitHub repository
    'style_nav_header_background': '#2980B9',  # Custom header color
    'navigation_depth': 4,
    'collapse_navigation': False,
    'sticky_navigation': True,
    'titles_only': False,
    'github_url': 'https://github.com/xragua/xraybinaryorbit',  # Link to your GitHub
}

html_static_path = ['_static']  # Path to your static files
