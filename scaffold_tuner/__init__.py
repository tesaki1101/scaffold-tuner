from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("scaffold-tuner")
except PackageNotFoundError:
    __version__ = "unknown"