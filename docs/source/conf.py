# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('./gehong/'))

# -- Project information -----------------------------------------------------

project = 'csst-ifs-gehong'
copyright = '2025, Shuai Feng'
author = 'Shuai Feng'

# The full version, including alpha/beta/rc tags
release = '3.0.1'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

#'sphinx.ext.duration',
#              'sphinx.ext.doctest',

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.autosummary',
              'sphinx.ext.intersphinx', 
              'sphinx.ext.napoleon']
autodoc_typehints = 'description'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
#html_theme = 'alabaster'
import sphinx_rtd_theme
html_theme = "sphinx_rtd_theme"
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

math_number_all = True  # Set this option to True if you want all displayed math to be numbered. The default is False.
math_eqref_format = 'Eq.{number}'  # gets rendered as, for example, Eq.10.

# If True, displayed math equations are numbered across pages when numfig
# is enabled. The numfig_secnum_depth setting is respected. The eq, not
# numref, role must be used to reference equation numbers. Default is
# True.
math_numfig = True

# see http://www.sphinx-doc.org/en/master/usage/configuration.html#confval-numfig
# If true, figures, tables and code-blocks are automatically numbered if they have a caption.
# The numref role is enabled. Obeyed so far only by HTML and LaTeX builders. Default is False.
# The LaTeX builder always assigns numbers whether this option is enabled or not.
numfig = True
numfig_secnum_depth = 2

# A dictionary mapping 'figure', 'table', 'code-block' and 'section' to strings that are used for format of figure numbers.
# As a special character, %s will be replaced to figure number.
# Default is to use 'Fig. %s' for 'figure', 'Table %s' for 'table', 'Listing %s' for 'code-block' and 'Section' for 'section'.
numfig_format = {
    'figure': 'Fig. %s',
    'table': 'Table %s',
    'code-block': 'Listing %s',
    'section': 'Section %s',
}

numfig_format = {
    'figure': '图 %s',
    'table': '表 %s',
    'code-block': '代码 %s',
    'section': '节 %s',
}

# 添加静态资源路径
html_static_path = ['_static']

# 设置 logo（相对于 _static 路径）
html_logo = '_static/gehong_logo.png'

# 可选：为了让 logo 更好看，设置主题选项（对 alabaster、sphinx_rtd_theme 等有效）
html_theme_options = {
    'logo_only': False,        # 只显示 logo，不显示项目名（可选）
    'display_version': False  # 不显示版本号（可选）
}

# conf.py 中允许 HTML 写入
rst_prolog = """
.. role:: raw-html(raw)
   :format: html
"""

math_number_all = False