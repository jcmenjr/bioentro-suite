#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bioentro v0.1.0 — Informational Zoom Suite for Biological Sequences

Computes Shannon entropy, KL divergence, Jensen-Shannon divergence (JSD),
Kolmogorov complexity proxy, and an Informational Priority Score (IPS) across
four biological levels of organization:

  genome          — whole genome vs uniform theoretical background (nucleotide)
  genes           — each CDS vs genomic background distribution (nucleotide)
  protein         — each protein vs proteome background distribution (amino acid)
  proteome-global — whole proteome vs uniform theoretical background (amino acid)

Background distributions per mode
──────────────────────────────────
  genome / proteome-global : uniform distribution (theoretical maximum)
  genes                    : k-mer distribution of the FULL GENOME
  protein                  : k-mer distribution of the FULL PROTEOME

This design is biologically meaningful:
  • Gene-vs-genome detects compositional outliers (pathogenicity islands, HGT).
  • Protein-vs-proteome detects compositional outliers within the proteome.
  • Genome/proteome-global vs uniform measures absolute informational efficiency.

Metrics
───────
  H_max      : maximum theoretical entropy (k × log2(alphabet_size))  [bits]
  H_real     : observed Shannon entropy                                [bits]
  Efficiency : H_real / H_max  (0–1 fraction)
  KL_bg      : KL(sequence || background)  — asymmetric, directional  [bits]
  JSD_bg     : JSD(sequence || background) — symmetric, bounded [0,1]
  sqrt_JSD   : √JSD_bg — proper metric distance [0,1]
  Kolmogorov : zlib compression ratio (proxy for long-range complexity)
  IPS        : Informational Priority Score = Efficiency × √JSD_bg × length_penalty
               Only for 'genes' and 'protein' modes.

Usage:
    bioentro -i genome.fna  -g genome.gff3 -m genome         -o out.tsv
    bioentro -i genome.fna  -g genome.gff3 -m genes          -o out.tsv
    bioentro -i proteins.faa               -m protein         -o out.tsv
    bioentro -i proteins.faa               -m proteome-global -o out.tsv

Requirements:
    biopython >= 1.79

Author:  bioentro contributors
License: MIT
"""

from __future__ import annotations

__version__ = "0.1.0"

# ---------------------------------------------------------------------------
# Standard library
# ---------------------------------------------------------------------------
import argparse
import csv
import logging
import math
import re
import sys
import zlib
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple, Union

# ---------------------------------------------------------------------------
# Third-party
# ---------------------------------------------------------------------------
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# ---------------------------------------------------------------------------
# Type aliases
# ---------------------------------------------------------------------------
KmerDist = Dict[str, float]   # k-mer → normalized frequency

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
PROTEIN_ALPHABET_SIZE: int = 20
DNA_ALPHABET_SIZE: int = 4

# IUPAC ambiguous characters + gap + stop codon
INVALID_AA: frozenset = frozenset("XBZ*-")
INVALID_NT: frozenset = frozenset("NRYWSKMBDHV-")

DEFAULT_K_PROTEIN: int = 1
DEFAULT_K_DNA: int = 6

ZLIB_LEVEL: int = 9

# Length scale for IPS penalty: sequences shorter than these are penalized
L0_NT: int = 300    # ~100 codons — minimum meaningful CDS
L0_AA: int = 100    # minimum meaningful protein

# Laplace smoothing for unseen k-mers when computing KL(P || bg)
KL_SMOOTHING: float = 1e-10

# ---------------------------------------------------------------------------
# Module-level logger  (never call basicConfig at import time)
# ---------------------------------------------------------------------------
logger = logging.getLogger(__name__)


def _configure_logging(verbose: bool = False) -> None:
    """Configure root logger. Called only from main()."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s — %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


# ---------------------------------------------------------------------------
# Custom exceptions  (no sys.exit inside processing functions)
# ---------------------------------------------------------------------------

class BioentroError(Exception):
    """Base class for all bioentro errors."""

class InputFileError(BioentroError):
    """Missing, unreadable, or empty input file."""

class AnnotationError(BioentroError):
    """Unrecoverable GFF3 parsing error."""

class EmptySequenceError(BioentroError):
    """No valid sequences could be processed."""


# ---------------------------------------------------------------------------
# Result dataclasses
# (field names == TSV column headers — guaranteed alignment via csv.DictWriter)
# ---------------------------------------------------------------------------

@dataclass
class GenomeResult:
    """Whole-genome metrics vs uniform theoretical background."""
    Organism: str
    Size_Mb: float
    H_max: float        # bits
    H_real: float       # bits
    Efficiency: float   # H_real / H_max  [0–1]
    JSD_uniform: float  # JSD(genome || uniform)  [0–1]
    sqrt_JSD: float     # metric distance [0–1]
    Kolmogorov: float   # zlib compression ratio


@dataclass
class GeneResult:
    """Per-CDS metrics; each gene compared against the full-genome background."""
    Gene_ID: str
    Product: str
    Len_NT: int
    H_max: float
    H_real: float
    Efficiency: float
    KL_genomebg: float   # KL(gene || genome)  — directional, bits
    JSD_genomebg: float  # JSD(gene || genome) — symmetric [0–1]
    sqrt_JSD: float
    Kolmogorov: float
    IPS: float           # Informational Priority Score [0–1]


@dataclass
class ProteinResult:
    """Per-protein metrics; each protein compared against the proteome background."""
    Prot_ID: str
    Product: str
    Len_AA: int
    H_max: float
    H_real: float
    Efficiency: float
    KL_proteomebg: float   # KL(protein || proteome)
    JSD_proteomebg: float  # JSD(protein || proteome) [0–1]
    sqrt_JSD: float
    Kolmogorov: float
    IPS: float


@dataclass
class ProteomeGlobalResult:
    """Whole-proteome metrics vs uniform theoretical background."""
    Target: str
    Size_AA: int
    H_max: float
    H_real: float
    Efficiency: float
    JSD_uniform: float
    sqrt_JSD: float
    Kolmogorov: float


# ---------------------------------------------------------------------------
# Core informational metric functions
# ---------------------------------------------------------------------------

def compute_kmer_dist(sequence: str, k: int, is_protein: bool) -> KmerDist:
    """Build a normalized k-mer frequency distribution from a sequence.

    Args:
        sequence:   Uppercase biological sequence.
        k:          K-mer size (>= 1).
        is_protein: True for amino-acid, False for nucleotide.

    Returns:
        Dict mapping k-mer string → relative frequency.
        Empty dict if no valid k-mers found.
    """
    invalid = INVALID_AA if is_protein else INVALID_NT
    kmers = [
        sequence[i : i + k]
        for i in range(len(sequence) - k + 1)
        if not any(c in invalid for c in sequence[i : i + k])
    ]
    if not kmers:
        return {}
    counts = Counter(kmers)
    total = sum(counts.values())
    return {kmer: count / total for kmer, count in counts.items()}


def calculate_h_max(k: int, is_protein: bool) -> float:
    """Maximum possible Shannon entropy for a k-mer distribution (bits).

    Uses the numerically stable form k × log2(A) instead of log2(A^k).
    """
    alphabet_size = PROTEIN_ALPHABET_SIZE if is_protein else DNA_ALPHABET_SIZE
    return k * math.log2(alphabet_size)


def calculate_shannon(p_dist: KmerDist) -> float:
    """Shannon entropy H(P) in bits.

    Args:
        p_dist: Normalized k-mer frequency distribution.

    Returns:
        H(P) = -Σ p_i log2(p_i)  [bits]
    """
    if not p_dist:
        return 0.0
    return -sum(p * math.log2(p) for p in p_dist.values() if p > 0.0)


def calculate_jsd_vs_uniform(
    p_dist: KmerDist,
    k: int,
    is_protein: bool,
) -> float:
    """JSD(P || Uniform) computed analytically — no need to enumerate all A^k k-mers.

    For k-mers NOT observed in the sequence (p_i = 0):
        M_i = (0 + 1/N) / 2 = 1/(2N)
        Contribution: 0.5 × (1/N) × log2(2)  [from Q side only]

    For k-mers observed (p_i > 0):
        M_i = (p_i + 1/N) / 2
        Contribution: 0.5 × p_i × log2(p_i/M_i) + 0.5 × (1/N) × log2((1/N)/M_i)

    Args:
        p_dist:     Observed k-mer distribution.
        k:          K-mer size.
        is_protein: Determines alphabet size N = A^k.

    Returns:
        JSD ∈ [0, 1]  (log2 base).
    """
    if not p_dist:
        return 0.0

    alphabet_size = PROTEIN_ALPHABET_SIZE if is_protein else DNA_ALPHABET_SIZE
    N = alphabet_size ** k
    p_uniform = 1.0 / N
    n_unobserved = N - len(p_dist)

    jsd = 0.0
    for p_i in p_dist.values():
        m_i = (p_i + p_uniform) / 2.0
        jsd += 0.5 * p_i * math.log2(p_i / m_i)
        jsd += 0.5 * p_uniform * math.log2(p_uniform / m_i)

    # Unobserved k-mers: p=0, so only Q side contributes
    # Q_i / M_i = p_uniform / (p_uniform / 2) = 2  →  log2(2) = 1
    jsd += n_unobserved * 0.5 * p_uniform * 1.0

    return round(max(0.0, min(1.0, jsd)), 4)


def calculate_jsd(p_dist: KmerDist, q_dist: KmerDist) -> float:
    """JSD(P || Q) between two empirical k-mer distributions.

    Symmetric: JSD(P||Q) == JSD(Q||P)
    Bounded:   JSD ∈ [0, 1]  (log2 base)
    Metric:    √JSD is a proper metric distance

    Args:
        p_dist: K-mer distribution of the target sequence.
        q_dist: K-mer distribution of the background.

    Returns:
        JSD ∈ [0, 1].
    """
    if not p_dist or not q_dist:
        return 0.0

    all_keys = set(p_dist) | set(q_dist)
    jsd = 0.0
    for key in all_keys:
        p = p_dist.get(key, 0.0)
        q = q_dist.get(key, 0.0)
        m = (p + q) / 2.0
        if p > 0.0:
            jsd += 0.5 * p * math.log2(p / m)
        if q > 0.0:
            jsd += 0.5 * q * math.log2(q / m)

    return round(max(0.0, min(1.0, jsd)), 4)


def calculate_kl(p_dist: KmerDist, q_dist: KmerDist) -> float:
    """KL(P || Q) with Laplace smoothing for k-mers absent in Q.

    Note: KL is asymmetric and unbounded — not suitable as a distance metric.
    Prefer JSD for ranking and comparison. KL is kept for its directional
    biological interpretation: how much information is lost when approximating
    P (individual sequence) with Q (background).

    Args:
        p_dist: Target sequence distribution.
        q_dist: Background distribution.

    Returns:
        KL divergence ≥ 0  [bits].
    """
    if not p_dist or not q_dist:
        return 0.0

    kl = 0.0
    for kmer, p in p_dist.items():
        q = q_dist.get(kmer, KL_SMOOTHING)
        kl += p * math.log2(p / q)

    return round(max(0.0, kl), 4)


def calculate_kolmogorov(sequence: str) -> float:
    """Kolmogorov complexity proxy via zlib compression ratio.

    Ratio → 1.0: high complexity (incompressible, random-like).
    Ratio → 0.0: high redundancy (highly compressible, repetitive).

    Complements Shannon entropy: Shannon captures local k-mer statistics;
    zlib captures long-range repetitive structure.

    Returns:
        compressed_bytes / raw_bytes, or 0.0 on error.
    """
    if not sequence:
        return 0.0
    raw = sequence.encode("ascii")
    try:
        compressed = zlib.compress(raw, level=ZLIB_LEVEL)
        return round(len(compressed) / len(raw), 4)
    except zlib.error as exc:
        logger.warning("zlib compression failed: %s", exc)
        return 0.0


def calculate_ips(
    efficiency: float,
    jsd_bg: float,
    length: int,
    l0: int,
) -> float:
    """Informational Priority Score for ranking hypothetical sequences.

    IPS = Efficiency × √JSD_bg × (1 − e^(−length / L0))

    Components:
        Efficiency    : normalized entropy [0–1] — compositional complexity.
        √JSD_bg       : metric distance from background [0–1] — how different
                        the sequence is from the genomic/proteomic context.
        length_penalty: (1 − e^(−L/L0)) — smoothly penalizes very short
                        sequences that are likely annotation artifacts.

    Higher IPS → sequence is compositionally complex AND informationally
    distinct from its background → higher priority for experimental validation.

    Args:
        efficiency: H_real / H_max  [0–1].
        jsd_bg:     JSD(sequence || background)  [0–1].
        length:     Sequence length (nt for genes, aa for proteins).
        l0:         Length scale constant (L0_NT or L0_AA).

    Returns:
        IPS ∈ [0, 1].
    """
    length_penalty = 1.0 - math.exp(-length / l0)
    ips = efficiency * math.sqrt(jsd_bg) * length_penalty
    return round(max(0.0, min(1.0, ips)), 4)


def _safe_efficiency(h_real: float, h_max: float) -> float:
    """H_real / H_max as a fraction [0–1]. Returns 0.0 if H_max == 0."""
    return round(h_real / h_max, 4) if h_max > 0.0 else 0.0


# ---------------------------------------------------------------------------
# GFF3 parsing
# ---------------------------------------------------------------------------

@dataclass
class _CdsRecord:
    """Internal GFF3 CDS feature representation."""
    chrom: str
    start: int    # 1-based, inclusive (GFF3)
    end: int      # 1-based, inclusive
    strand: str
    gene_id: str
    product: str


def _extract_attr(attributes: str, key: str, fallback: str = "") -> str:
    """Extract a GFF3 attribute value by key name."""
    match = re.search(rf"{re.escape(key)}=([^;\n]+)", attributes)
    return match.group(1).strip() if match else fallback


def _parse_gff_cds(gff_path: Path) -> Iterator[_CdsRecord]:
    """Yield CDS records from a GFF3 file.

    Args:
        gff_path: Path to the GFF3 annotation file.

    Yields:
        _CdsRecord for each CDS feature found.

    Raises:
        AnnotationError: On file I/O errors or invalid coordinate fields.
    """
    try:
        with gff_path.open("r", encoding="utf-8") as fh:
            for lineno, raw in enumerate(fh, start=1):
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 9:
                    logger.debug("GFF line %d skipped (< 9 fields).", lineno)
                    continue
                if parts[2] != "CDS":
                    continue
                try:
                    start, end = int(parts[3]), int(parts[4])
                except ValueError as exc:
                    raise AnnotationError(
                        f"GFF line {lineno}: invalid coordinates — {exc}"
                    ) from exc

                attrs = parts[8]
                yield _CdsRecord(
                    chrom=parts[0],
                    start=start,
                    end=end,
                    strand=parts[6],
                    gene_id=_extract_attr(attrs, "ID", fallback="unknown"),
                    product=_extract_attr(attrs, "product", fallback="hypothetical protein"),
                )
    except OSError as exc:
        raise AnnotationError(f"Cannot open GFF file '{gff_path}': {exc}") from exc


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def _validate_file(path: Path, suffix: Union[str, Tuple[str, ...]] = "") -> None:
    """Check that a path exists, is a regular file, and has an expected extension.

    A suffix mismatch only logs a warning — it does not abort execution.

    Raises:
        InputFileError: If path does not exist or is not a regular file.
    """
    if not path.exists():
        raise InputFileError(f"File not found: '{path}'")
    if not path.is_file():
        raise InputFileError(f"Not a regular file: '{path}'")
    if suffix:
        allowed = (suffix,) if isinstance(suffix, str) else suffix
        if path.suffix.lower() not in allowed:
            logger.warning(
                "File '%s' has unexpected extension (expected one of %s).",
                path.name, allowed,
            )


def _extract_organism_name(record: SeqRecord) -> str:
    """Extract organism name from NCBI-style FASTA description, or fall back to ID."""
    match = re.search(r"\[organism=([^\]]+)\]", record.description)
    return match.group(1) if match else record.id


def write_tsv(results: list, output_path: Path) -> None:
    """Write a list of dataclass instances to a tab-separated file.

    Column names are derived from dataclass field names — guaranteed to match data.

    Args:
        results:     Non-empty list of dataclass instances (homogeneous type).
        output_path: Destination path (parent directories created if needed).

    Raises:
        ValueError: If results is empty.
        OSError:    If the file cannot be written.
    """
    if not results:
        raise ValueError("Nothing to write — results list is empty.")

    fieldnames = list(asdict(results[0]).keys())
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for item in results:
            writer.writerow(asdict(item))


# ---------------------------------------------------------------------------
# Processing modes
# ---------------------------------------------------------------------------

def process_genome_global(
    fna_path: Path,
    gff_path: Optional[Path],   # unused here; kept for uniform signature
    k: int,
) -> List[GenomeResult]:
    """Compute whole-genome informational metrics vs uniform background.

    Background: uniform distribution (theoretical maximum entropy).
    No IPS — this mode produces a single aggregate measurement.

    Args:
        fna_path: Genome FASTA (.fna / .fa / .fasta).
        gff_path: Unused (required by uniform mode signature).
        k:        K-mer size.

    Returns:
        Single-element list of GenomeResult.

    Raises:
        InputFileError: If FASTA is empty or unreadable.
    """
    _validate_file(fna_path, suffix=(".fna", ".fa", ".fasta"))

    records = list(SeqIO.parse(fna_path, "fasta"))
    if not records:
        raise InputFileError(f"No sequences found in '{fna_path}'.")

    org_name = _extract_organism_name(records[0])
    genome_seq = "".join(str(r.seq).upper() for r in records)

    p_dist = compute_kmer_dist(genome_seq, k, is_protein=False)
    h_max = calculate_h_max(k, is_protein=False)
    h_real = calculate_shannon(p_dist)
    efficiency = _safe_efficiency(h_real, h_max)
    jsd_unif = calculate_jsd_vs_uniform(p_dist, k, is_protein=False)

    return [
        GenomeResult(
            Organism=org_name,
            Size_Mb=round(len(genome_seq) / 1e6, 4),
            H_max=round(h_max, 4),
            H_real=round(h_real, 4),
            Efficiency=efficiency,
            JSD_uniform=jsd_unif,
            sqrt_JSD=round(math.sqrt(jsd_unif), 4),
            Kolmogorov=calculate_kolmogorov(genome_seq),
        )
    ]


def process_gene_level(
    fna_path: Path,
    gff_path: Optional[Path],
    k: int,
) -> List[GeneResult]:
    """Compute per-CDS informational metrics vs the full-genome background.

    Background: k-mer distribution of the COMPLETE GENOME (all contigs
    concatenated). This makes KL_bg and JSD_bg measure how much each
    individual CDS deviates from the organism's overall nucleotide
    composition — a biologically meaningful signal for detecting
    horizontally transferred genes, pathogenicity islands, etc.

    Args:
        fna_path: Genome FASTA (.fna).
        gff_path: GFF3 annotation file.
        k:        K-mer size.

    Returns:
        List of GeneResult, one per processed CDS.

    Raises:
        InputFileError, AnnotationError, EmptySequenceError.
    """
    _validate_file(fna_path, suffix=(".fna", ".fa", ".fasta"))
    _validate_file(gff_path, suffix=(".gff", ".gff3"))

    genome = SeqIO.to_dict(SeqIO.parse(fna_path, "fasta"))
    if not genome:
        raise InputFileError(f"No sequences found in '{fna_path}'.")

    # Build genomic background distribution from the full concatenated genome
    logger.info("Computing genomic background k-mer distribution (k=%d)...", k)
    full_genome_seq = "".join(str(r.seq).upper() for r in genome.values())
    bg_dist = compute_kmer_dist(full_genome_seq, k, is_protein=False)

    h_max = calculate_h_max(k, is_protein=False)
    results: List[GeneResult] = []
    skipped = 0

    for cds in _parse_gff_cds(gff_path):
        if cds.chrom not in genome:
            logger.warning(
                "Chromosome '%s' not in FASTA — skipping gene '%s'.",
                cds.chrom, cds.gene_id,
            )
            skipped += 1
            continue

        region = genome[cds.chrom].seq[cds.start - 1 : cds.end]
        seq = str(
            region.reverse_complement() if cds.strand == "-" else region
        ).upper()

        p_dist = compute_kmer_dist(seq, k, is_protein=False)
        h_real = calculate_shannon(p_dist)
        efficiency = _safe_efficiency(h_real, h_max)
        kl_bg = calculate_kl(p_dist, bg_dist)
        jsd_bg = calculate_jsd(p_dist, bg_dist)

        results.append(
            GeneResult(
                Gene_ID=cds.gene_id,
                Product=cds.product,
                Len_NT=len(seq),
                H_max=round(h_max, 4),
                H_real=round(h_real, 4),
                Efficiency=efficiency,
                KL_genomebg=kl_bg,
                JSD_genomebg=jsd_bg,
                sqrt_JSD=round(math.sqrt(jsd_bg), 4),
                Kolmogorov=calculate_kolmogorov(seq),
                IPS=calculate_ips(efficiency, jsd_bg, len(seq), L0_NT),
            )
        )

    if skipped:
        logger.warning("%d CDS feature(s) skipped (chromosome not in FASTA).", skipped)
    if not results:
        raise EmptySequenceError("No CDS features could be processed.")

    logger.info("Processed %d CDS features.", len(results))
    return results


def process_protein_individual(
    faa_path: Path,
    k: int,
) -> List[ProteinResult]:
    """Compute per-protein informational metrics vs the full-proteome background.

    Background: k-mer distribution of ALL proteins in the .faa file
    concatenated. This makes KL_bg and JSD_bg measure how much each
    individual protein deviates from the organism's overall amino acid
    (or k-mer) composition — useful for detecting compositionally unusual
    hypothetical proteins.

    Args:
        faa_path: Protein FASTA (.faa).
        k:        K-mer size.

    Returns:
        List of ProteinResult, one per FASTA record.

    Raises:
        InputFileError, EmptySequenceError.
    """
    _validate_file(faa_path, suffix=(".faa", ".fa", ".fasta"))

    # Two-pass: first build proteome background, then process individually
    all_records = list(SeqIO.parse(faa_path, "fasta"))
    if not all_records:
        raise EmptySequenceError(f"No protein sequences found in '{faa_path}'.")

    logger.info("Computing proteome background k-mer distribution (k=%d)...", k)
    proteome_seq = "".join(str(r.seq).upper() for r in all_records)
    bg_dist = compute_kmer_dist(proteome_seq, k, is_protein=True)

    h_max = calculate_h_max(k, is_protein=True)
    results: List[ProteinResult] = []

    for record in all_records:
        seq = str(record.seq).upper()
        p_dist = compute_kmer_dist(seq, k, is_protein=True)
        h_real = calculate_shannon(p_dist)
        efficiency = _safe_efficiency(h_real, h_max)
        kl_bg = calculate_kl(p_dist, bg_dist)
        jsd_bg = calculate_jsd(p_dist, bg_dist)

        results.append(
            ProteinResult(
                Prot_ID=record.id,
                Product=record.description,
                Len_AA=len(seq),
                H_max=round(h_max, 4),
                H_real=round(h_real, 4),
                Efficiency=efficiency,
                KL_proteomebg=kl_bg,
                JSD_proteomebg=jsd_bg,
                sqrt_JSD=round(math.sqrt(jsd_bg), 4),
                Kolmogorov=calculate_kolmogorov(seq),
                IPS=calculate_ips(efficiency, jsd_bg, len(seq), L0_AA),
            )
        )

    logger.info("Processed %d protein sequences.", len(results))
    return results


def process_proteome_global(
    faa_path: Path,
    k: int,
) -> List[ProteomeGlobalResult]:
    """Compute whole-proteome informational metrics vs uniform background.

    Background: uniform distribution (theoretical maximum entropy).
    No IPS — this mode produces a single aggregate measurement.
    Analogous to genome mode but for the amino acid sequence space.

    Args:
        faa_path: Protein FASTA (.faa).
        k:        K-mer size.

    Returns:
        Single-element list of ProteomeGlobalResult.

    Raises:
        InputFileError, EmptySequenceError.
    """
    _validate_file(faa_path, suffix=(".faa", ".fa", ".fasta"))

    all_seqs = "".join(str(r.seq).upper() for r in SeqIO.parse(faa_path, "fasta"))
    if not all_seqs:
        raise EmptySequenceError(f"No sequences found in '{faa_path}'.")

    p_dist = compute_kmer_dist(all_seqs, k, is_protein=True)
    h_max = calculate_h_max(k, is_protein=True)
    h_real = calculate_shannon(p_dist)
    efficiency = _safe_efficiency(h_real, h_max)
    jsd_unif = calculate_jsd_vs_uniform(p_dist, k, is_protein=True)

    return [
        ProteomeGlobalResult(
            Target="Global_Proteome",
            Size_AA=len(all_seqs),
            H_max=round(h_max, 4),
            H_real=round(h_real, 4),
            Efficiency=efficiency,
            JSD_uniform=jsd_unif,
            sqrt_JSD=round(math.sqrt(jsd_unif), 4),
            Kolmogorov=calculate_kolmogorov(all_seqs),
        )
    ]


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="bioentro",
        description=(
            f"bioentro v{__version__} — Informational Zoom Suite\n\n"
            "Shannon entropy · KL divergence · Jensen-Shannon divergence\n"
            "Kolmogorov complexity · Informational Priority Score (IPS)\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Background distributions per mode:\n"
            "  genome          → uniform (theoretical maximum)\n"
            "  genes           → full genome k-mer distribution\n"
            "  protein         → full proteome k-mer distribution\n"
            "  proteome-global → uniform (theoretical maximum)\n\n"
            "Examples:\n"
            "  %(prog)s -i genome.fna  -g genome.gff3 -m genome         -o genome.tsv\n"
            "  %(prog)s -i genome.fna  -g genome.gff3 -m genes          -o genes.tsv\n"
            "  %(prog)s -i proteins.faa               -m protein    -k 2 -o prot.tsv\n"
            "  %(prog)s -i proteins.faa               -m proteome-global -o prot_g.tsv\n"
        ),
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Input FASTA (.fna for DNA modes, .faa for protein modes)",
    )
    parser.add_argument(
        "-g", "--gff",
        help="GFF3 annotation file (required for 'genome' and 'genes' modes)",
    )
    parser.add_argument(
        "-m", "--mode", required=True,
        choices=["genome", "genes", "protein", "proteome-global"],
        help="Analysis level",
    )
    parser.add_argument(
        "-k", "--kmer", type=int,
        help=(
            f"K-mer size "
            f"(default: {DEFAULT_K_PROTEIN} for protein modes, "
            f"{DEFAULT_K_DNA} for DNA modes)"
        ),
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output TSV file path",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Enable debug-level logging",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}",
    )
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    """CLI entry point.

    Args:
        argv: Argument list (defaults to sys.argv[1:] when None).

    Returns:
        Exit code: 0 = success, 1 = handled error, 130 = keyboard interrupt.
    """
    parser = _build_parser()
    args = parser.parse_args(argv)
    _configure_logging(verbose=args.verbose)

    input_path = Path(args.input)
    output_path = Path(args.output)
    gff_path = Path(args.gff) if args.gff else None

    # Validate k
    if args.kmer is not None:
        if args.kmer < 1:
            parser.error("--kmer must be a positive integer.")
        k = args.kmer
    else:
        k = DEFAULT_K_PROTEIN if args.mode in {"protein", "proteome-global"} else DEFAULT_K_DNA

    # GFF requirement
    if args.mode in {"genome", "genes"} and gff_path is None:
        parser.error(f"--gff is required for mode '{args.mode}'.")

    logger.info(
        "bioentro v%s | mode=%s | k=%d | input=%s",
        __version__, args.mode, k, input_path.name,
    )

    try:
        if args.mode == "genome":
            data = process_genome_global(input_path, gff_path, k)
        elif args.mode == "genes":
            data = process_gene_level(input_path, gff_path, k)
        elif args.mode == "protein":
            data = process_protein_individual(input_path, k)
        else:  # proteome-global
            data = process_proteome_global(input_path, k)

        write_tsv(data, output_path)
        logger.info(
            "Done. %d record(s) written to '%s'.", len(data), output_path
        )
        return 0

    except BioentroError as exc:
        logger.error("%s", exc)
        return 1
    except KeyboardInterrupt:
        logger.warning("Interrupted by user.")
        return 130


if __name__ == "__main__":
    sys.exit(main())