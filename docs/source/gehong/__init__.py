"""
gehong: A Python package for generating synthetic galaxy spectra and cubes.

This package includes tools for:
- Stellar continuum modeling
- Nebular emission line modeling
- 2D and 3D spectral cube construction
- Configurable simulation of observational properties

Modules:
- config.py       : Configuration utilities
- spec1d.py       : Spectrum generation (1D)
- map2d.py        : Generation of 2D parameter maps
- cube3d.py       : Assembly of 3D datacubes

Path Management:
GEHONG uses resource files (e.g., SSP templates, emission-line libraries).
You must configure the path to these resources before use, by one of:

    1. Setting the environment variable:
       export gehong_path=/path/to/data

    2. Calling the function:
       gehong.set_path("/path/to/data")

    3. Placing the resource folder inside gehong/data/ (for development use)
"""


__version__ = '3.1.0'

# Internal global variable to store the resource path
_gehong_path = os.environ.get("gehong_path", None)


def set_path(path):
    """
    Manually set the path to GEHONG resource files.

    Parameters
    ----------
    path : str
        Absolute or relative path to the directory containing template or resource files.

    Raises
    ------
    NotADirectoryError
        If the provided path does not exist or is not a directory.
    """
    global _gehong_path
    if not os.path.isdir(path):
        raise NotADirectoryError(f"[gehong] Provided path is not a valid directory: {path}")
    _gehong_path = path


def get_path():
    """
    Get the current path to GEHONG resource files.

    Returns
    -------
    str
        Path to the resource directory.

    Raises
    ------
    RuntimeError
        If the path is not set and no fallback path is found.
    """
    global _gehong_path

    if _gehong_path is not None:
        return _gehong_path

    # Fallback: use local 'data/' directory under gehong/
    current_dir = os.path.dirname(inspect.getfile(inspect.currentframe()))
    default_path = os.path.join(current_dir, "data")
    if os.path.isdir(default_path):
        return default_path

    raise RuntimeError(
        "[gehong] Resource path is not set. Please do one of the following:\n"
        "  1. Set the environment variable 'gehong_path'\n"
        "     e.g., export gehong_path=/path/to/data\n"
        "  2. Or call gehong.set_path('/path/to/data')\n"
        "  3. Or place the resource folder inside 'gehong/data/'"
    )
