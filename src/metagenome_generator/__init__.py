"""CHIMERA (Configurable Hybrid In-silico Metagenome Emulator for Read Analysis): build simulated metagenome FASTAs from NCBI or in-house genomes."""

__version__ = "0.1.0"

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
