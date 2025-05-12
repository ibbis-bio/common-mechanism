from importlib.metadata import version, PackageNotFoundError
try:
    __version__ = version("commec")
except ImportError, PackageNotFoundError:
    __version__ = "X.X.X"
