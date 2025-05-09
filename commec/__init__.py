from importlib.metadata import version
try:
    __version__ = version(__name__)
except:
    __version__ = "X.X.X"