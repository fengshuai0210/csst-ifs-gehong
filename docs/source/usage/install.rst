Installation
============

.. _installation:

Installing GEHONG
-----------------

To install **GEHONG**, simply use:

.. code-block:: bash

   pip install csst-ifs-gehong

This installs the core Python package. However, **spectral templates** required for modeling are not included and must be downloaded separately, as described below.

Downloading Spectral Templates
------------------------------

GEHONG relies on external datasets for stellar and emission line spectral modeling. These templates are **not bundled** in the PyPI package and must be downloaded from GitHub.

You can clone the official data repository here:

- GitHub: https://github.com/fengshuai0210/csst-ifs-gehong-data

Clone it using:

.. code-block:: bash

   git clone https://github.com/fengshuai0210/csst-ifs-gehong-data.git

This will produce a folder structure like:

.. code-block:: text

   csst-ifs-gehong-data/
   └── data/

Setting the Resource Path
--------------------------

GEHONG requires access to external resource files (e.g., SSP templates, emission line data),
located in the ``data/`` folder of the downloaded data repository.

You can configure the resource path using **any one of the following three methods**:

1. **Set an environment variable** (recommended for persistent use)

   Add the following line to your shell configuration file
   (e.g., ``~/.bash_profile``, ``~/.zshrc``, or ``~/.bashrc``):

   .. code-block:: bash

      export gehong_path='/Users/sfeng/GIT/csst-ifs-gehong-data'

   Then reload the shell environment:

   .. code-block:: bash

      source ~/.bash_profile  # or ~/.zshrc, ~/.bashrc, etc.

2. **Set the path at runtime in Python**

   For temporary or interactive use, you can configure the path manually:

   .. code-block:: python

      import gehong
      gehong.set_path("/Users/sfeng/GIT/csst-ifs-gehong-data")

3. **Use the development mode fallback**

   If you're developing GEHONG locally, you may place the ``data/`` folder directly inside the
   ``gehong/`` package directory. This fallback will be used automatically if no other path is specified.

All methods require that the target directory contains the expected subdirectory structure:

.. code-block:: text

   csst-ifs-gehong-data/
   └── data/

Verifying Installation
----------------------

To confirm that GEHONG is installed and configured correctly, you can run the following test script in Python:

.. code-block:: python

   import gehong
   import os

   # Check version
   print("GEHONG version:", gehong.__version__)

   # Get the configured resource path
   try:
       path = gehong.get_path()
   except RuntimeError as e:
       raise RuntimeError("Resource path is not set correctly:\n" + str(e))
   else:
       print("GEHONG resource path =", path)

   # Check if 'data/' subdirectory exists
   expected_path = os.path.join(path, "data")
   if not os.path.isdir(expected_path):
       raise FileNotFoundError(f"'data' folder not found at: {expected_path}")
   else:
       print("Data folder found at:", expected_path)

If all checks pass without error, your installation and configuration are complete.

Next Steps
----------

If everything is installed and configured correctly, you're ready to start using GEHONG.

👉 We recommend starting with the :doc:`Quickstart Guide <quickstart>`, which introduces example workflows using interactive Jupyter notebooks.

