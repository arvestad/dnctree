try:
    from dnctree._version import version as __version__
except ImportError:
    # Package was not installed via setuptools-scm (e.g. running from a
    # source checkout without a build), so no _version.py was generated.
    __version__ = "0+unknown"

