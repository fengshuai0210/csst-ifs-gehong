Change Log
==========

Version 3.1.0 (2025-07-25)
--------------------------

- Enhanced ``Config`` with two predefined parameter sets for quick compatibility with **CSST-IFS-SIM** and **CSST-IFS-ETC**.
- Stellar continuum spectra can now be generated based on user-defined **star formation history (SFH)** and **chemical enrichment history (CEH)**.
- H II region spectra (ionized gas) can be generated directly from **star formation rate (SFR)** input.
- Updated **theoretical stellar spectral templates** used for single-star spectrum simulation.
- Refactored ``Map2D`` and ``Cube3D`` modules to support physical-mode modeling using **SFH**, **CEH**, and **SFR** inputs.

Version 3.0.0 (2025-01)
-----------------------

- Adjusted version number to align with other modules in the **CSST simulation system**.
- Improved ``Config`` structure with default values tuned for **CSST-IFS-ETC** compatibility.
- Codebase formatting aligned with **PEP8** standards for improved readability and maintainability.

Older Versions
--------------

For earlier versions, please refer to the `PyPI release archive <https://pypi.org/project/csst-ifs-gehong/#history>`_.
