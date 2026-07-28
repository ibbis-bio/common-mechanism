from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("commec")
except (ImportError, PackageNotFoundError):
    __version__ = "X.X.X"
