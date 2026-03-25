"""
bioentro suite — Informational Analysis Suite for Biological Sequences

Tools
─────
  bioentro  : Shannon entropy, JSD, Kolmogorov complexity and IPS across four
              biological levels (genome, genes, protein, proteome-global).
  netentro  : Informational similarity network for functional prediction of
              hypothetical proteins.

Typical usage (CLI after installation):
    bioentro -i genome.fna -g genome.gff3 -m genes -o genes.tsv
    netentro -i genes.tsv  -t LOCUS_TAG   -o network.png
"""

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("bioentro-suite")
except PackageNotFoundError:          # running directly from source tree
    __version__ = "0.1.0-dev"

__all__ = ["__version__"]