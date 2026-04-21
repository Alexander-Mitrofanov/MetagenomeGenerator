"""CHIMERA (Configurable Hybrid In-silico Metagenome Emulator for Read Analysis): build simulated metagenome FASTAs from NCBI or in-house genomes."""

# Read the installed package version from metadata so the canonical source of
# truth stays in pyproject.toml. Falls back to a sentinel when the package is
# imported directly from a source tree without being installed.
from importlib.metadata import PackageNotFoundError, version as _pkg_version

try:
    __version__ = _pkg_version("chimera-metagenome-generator")
except PackageNotFoundError:  # pragma: no cover - source tree / dev checkout
    __version__ = "0.0.0+unknown"

from .download_genomes import download_genomes, load_accessions
from .chunk_genomes import build_metagenome, get_file_stats
from .genome_layout import validate_genome_dir

__all__ = [
    "__version__",
    "download_genomes",
    "load_accessions",
    "build_metagenome",
    "get_file_stats",
    "validate_genome_dir",
]
