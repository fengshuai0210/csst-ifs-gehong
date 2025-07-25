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
copyright = '2025, CSST-IFS Team'
author = 'Shuai Feng'

# The full version, including alpha/beta/rc tags
release = '3.1.0'

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

# 添加静态资源路径
html_static_path = ['_static']

# 设置 logo（相对于 _static 路径）
html_logo = '_static/gehong_logo.png'

# 可选：为了让 logo 更好看，设置主题选项（对 alabaster、sphinx_rtd_theme 等有效）
html_theme_options = {
    'logo_only': True,        # 只显示 logo，不显示项目名（可选）
    'display_version': False  # 不显示版本号（可选）
}

# conf.py 中允许 HTML 写入
rst_prolog = """
.. role:: raw-html(raw)
   :format: html
"""
