#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pangentro v0.1.1 — Pangenomic Integration Tool for bioentro-suite

Integrates bioentro informational metrics with Panaroo pangenome output
to answer the central research question:

    Do hypothetical proteins with high IPS concentrate in the shell/cloud
    genome, consistent with horizontal gene transfer or clinical niche
    adaptation in Pseudomonas aeruginosa?

Subcommands
───────────
  run-bioentro : Run bioentro protein mode on all isolates in a directory.
                 Wrapper that automates the per-isolate loop.

  integrate    : Cross bioentro TSV results with Panaroo gene_presence_absence.csv.
                 Produces a unified table with IPS + pangenomic category per cluster.
                 Two IPS modes:
                   individual  — averages per-isolate IPS values per cluster
                   pangenomic  — recalculates IPS using pan-proteome as background

  compare      : Statistical comparison of IPS distributions across pangenomic
                 categories (core, shell, cloud). Outputs summary table +
                 publication-quality figure.

Scientific rationale — IPS modes
─────────────────────────────────
  individual mode:
    Each protein's IPS is calculated against its own isolate's proteome.
    The per-cluster value is the mean IPS across all isolate representatives.
    Captures: atypicality within each organism's own compositional context.
    Limitation: backgrounds differ between isolates — means are not directly
    comparable across isolates.

  pangenomic mode (v0.1.1):
    Background = concatenated core genome sequences ONLY (clusters with freq
    >= CORE_THRESHOLD from pan_genome_reference.fa).
    Rationale: the core genome represents the stable, conserved "identity" of
    the species — genes under purifying selection that have been retained in
    all isolates. Using it as background means IPS_pan measures how different
    a protein is from the compositional baseline that defines the organism.
    Shell and cloud proteins are expected to have higher IPS_pan because they
    deviate from that stable core compositional space — consistent with HGT
    or clinical niche adaptation.
    Advantage: all IPS values are in the same space — directly comparable.
    Improvement over v0.1.0: using pan-proteome as background diluted the
    signal because shell/cloud sequences were part of their own background.

  Both modes together:
    IPS_individual high + IPS_pangenomic low  → typical in its isolate,
                                                typical in species
    IPS_individual low  + IPS_pangenomic high → typical in isolate,
                                                atypical in species (local adaptation)
    IPS_individual high + IPS_pangenomic high → atypical in both contexts
                                                (strongest HGT/novelty signal)

Pangenomic categories (Panaroo definitions)
────────────────────────────────────────────
  core      : present in >= 99% of genomes
  soft_core : present in 95–99% of genomes
  shell     : present in 15–95% of genomes
  cloud     : present in < 15% of genomes

Usage:
    pangentro run-bioentro -f bakta_outputs/ -o bioentro_results/
    pangentro integrate -p gene_presence_absence.csv -b bioentro_results/ -o metrics.tsv
    pangentro integrate -p gene_presence_absence.csv -r pan_genome_reference.fa \\
                        --ips-mode pangenomic -o metrics_pan.tsv
    pangentro compare -i metrics.tsv -o comparison/

Requirements:
    biopython >= 1.79, pandas >= 1.5, numpy >= 1.23,
    scipy >= 1.9, matplotlib >= 3.6

Author:  bioentro contributors
License: MIT
"""

from __future__ import annotations

__version__ = "0.1.1"

# ---------------------------------------------------------------------------
# Standard library
# ---------------------------------------------------------------------------
import argparse
import csv
import logging
import math
import subprocess
import sys
import zlib
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Third-party
# ---------------------------------------------------------------------------
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy import stats

# ---------------------------------------------------------------------------
# Constants — category thresholds (Panaroo defaults)
# ---------------------------------------------------------------------------
CORE_THRESHOLD: float      = 0.99   # >= 99% of genomes
SOFT_CORE_THRESHOLD: float = 0.95   # >= 95%
SHELL_THRESHOLD: float     = 0.15   # >= 15%
# cloud = < 15%

PROTEIN_ALPHABET_SIZE: int = 20
DEFAULT_K_PROTEIN: int     = 1
ZLIB_LEVEL: int            = 9
L0_AA: int                 = 100
KL_SMOOTHING: float        = 1e-10

INVALID_AA: frozenset = frozenset("XBZ*-")

# Category colors for figures — consistent with netentro palette
CATEGORY_COLORS: Dict[str, str] = {
    "core":      "#4E9AF1",   # blue
    "soft_core": "#1D9E75",   # teal
    "shell":     "#F4A460",   # orange
    "cloud":     "#E74C3C",   # red
}

CATEGORY_ORDER: List[str] = ["core", "soft_core", "shell", "cloud"]

# ---------------------------------------------------------------------------
# Logger
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
# Custom exceptions
# ---------------------------------------------------------------------------

class PangentroError(Exception):
    """Base class for pangentro errors."""

class InputError(PangentroError):
    """Missing or malformed input file."""

class InsufficientDataError(PangentroError):
    """Not enough data for the requested analysis."""


# ---------------------------------------------------------------------------
# Result dataclasses
# ---------------------------------------------------------------------------

@dataclass
class ClusterMetrics:
    """Informational metrics for one pangenomic cluster."""
    Cluster:         str
    Annotation:      str
    Category:        str     # core | soft_core | shell | cloud
    N_genomes:       int     # number of isolates where cluster is present
    Freq:            float   # N_genomes / total_genomes  [0–1]
    Is_hypothetical: bool
    # Individual mode metrics (mean across isolates)
    Mean_IPS_indiv:  float
    Std_IPS_indiv:   float
    N_isolates_indiv: int    # isolates with bioentro data for this cluster
    # Pangenomic mode metrics (recalculated against pan-proteome)
    IPS_pan:         float   # -1.0 if not calculated
    H_real_pan:      float
    Efficiency_pan:  float
    JSD_pan:         float
    sqrt_JSD_pan:    float
    Kolmogorov_pan:  float


@dataclass
class CategorySummary:
    """Statistical summary of IPS by pangenomic category."""
    Category:        str
    N_clusters:      int
    N_hypothetical:  int
    Pct_hypothetical: float
    Mean_IPS:        float
    Median_IPS:      float
    Std_IPS:         float
    Q25_IPS:         float
    Q75_IPS:         float


@dataclass
class StatTest:
    """Result of a pairwise statistical test between two categories."""
    Category_A:  str
    Category_B:  str
    Test:        str     # kruskal | mannwhitneyu
    Statistic:   float
    P_value:     float
    Significant: bool    # p < 0.05
    Effect_size: float   # rank-biserial correlation


# ---------------------------------------------------------------------------
# Informational metric functions (minimal reimplementation for pangentro)
# These mirror bioentro.py exactly — same math, no dependencies between scripts
# ---------------------------------------------------------------------------

def _compute_kmer_dist(sequence: str, k: int) -> Dict[str, float]:
    """Normalized k-mer frequency distribution (protein sequences)."""
    kmers = [
        sequence[i : i + k]
        for i in range(len(sequence) - k + 1)
        if not any(c in INVALID_AA for c in sequence[i : i + k])
    ]
    if not kmers:
        return {}
    counts = Counter(kmers)
    total = sum(counts.values())
    return {kmer: count / total for kmer, count in counts.items()}


def _calculate_shannon(p_dist: Dict[str, float]) -> float:
    """Shannon entropy H(P) in bits."""
    if not p_dist:
        return 0.0
    return -sum(p * math.log2(p) for p in p_dist.values() if p > 0.0)


def _calculate_h_max(k: int) -> float:
    """Maximum theoretical Shannon entropy for protein k-mers."""
    return k * math.log2(PROTEIN_ALPHABET_SIZE)


def _calculate_jsd(
    p_dist: Dict[str, float],
    q_dist: Dict[str, float],
) -> float:
    """JSD(P || Q) — symmetric, bounded [0,1]."""
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
    return round(max(0.0, min(1.0, jsd)), 6)


def _calculate_kolmogorov(sequence: str) -> float:
    """Kolmogorov complexity proxy via zlib compression ratio."""
    if not sequence:
        return 0.0
    raw = sequence.encode("ascii")
    try:
        compressed = zlib.compress(raw, level=ZLIB_LEVEL)
        return round(len(compressed) / len(raw), 6)
    except zlib.error:
        return 0.0


def _calculate_ips(
    efficiency: float,
    jsd_bg: float,
    length: int,
    l0: int = L0_AA,
) -> float:
    """Informational Priority Score."""
    length_penalty = 1.0 - math.exp(-length / l0)
    ips = efficiency * math.sqrt(jsd_bg) * length_penalty
    return round(max(0.0, min(1.0, ips)), 6)


def _safe_efficiency(h_real: float, h_max: float) -> float:
    """H_real / H_max. Returns 0.0 if H_max == 0."""
    return round(h_real / h_max, 6) if h_max > 0.0 else 0.0


# ---------------------------------------------------------------------------
# Pangenomic category assignment
# ---------------------------------------------------------------------------

def assign_category(freq: float) -> str:
    """Assign a pangenomic category based on presence frequency.

    Args:
        freq: Fraction of genomes where the cluster is present [0–1].

    Returns:
        One of: 'core', 'soft_core', 'shell', 'cloud'.
    """
    if freq >= CORE_THRESHOLD:
        return "core"
    if freq >= SOFT_CORE_THRESHOLD:
        return "soft_core"
    if freq >= SHELL_THRESHOLD:
        return "shell"
    return "cloud"


def is_hypothetical(annotation: str) -> bool:
    """Return True if annotation suggests an unannotated protein."""
    lower = annotation.lower()
    keywords = (
        "hypothetical", "unknown function", "uncharacterized",
        "putative", "predicted protein", "duf",
    )
    return any(kw in lower for kw in keywords)


# ---------------------------------------------------------------------------
# Panaroo matrix parser
# ---------------------------------------------------------------------------

def load_presence_absence(
    csv_path: Path,
    n_genomes: Optional[int] = None,
) -> pd.DataFrame:
    """Load and parse Panaroo gene_presence_absence.csv.

    Args:
        csv_path:  Path to gene_presence_absence.csv.
        n_genomes: Expected number of genomes. Auto-detected if None.

    Returns:
        DataFrame with columns: Cluster, Annotation, N_genomes, Freq,
        Category, Is_hypothetical, plus one column per isolate.

    Raises:
        InputError: If file is missing or malformed.
    """
    if not csv_path.exists():
        raise InputError(f"File not found: '{csv_path}'")

    try:
        df = pd.read_csv(csv_path, low_memory=False)
    except Exception as exc:
        raise InputError(f"Cannot read '{csv_path}': {exc}") from exc

    # Panaroo columns: Gene, Non-unique Gene name, Annotation, then one per isolate
    required = {"Gene", "Annotation"}
    if not required.issubset(df.columns):
        raise InputError(
            f"Expected columns {required} in gene_presence_absence.csv. "
            f"Found: {list(df.columns[:5])}..."
        )

    # Isolate columns = everything after the first 3 fixed columns
    isolate_cols = [c for c in df.columns if c not in
                    {"Gene", "Non-unique Gene name", "Annotation"}]

    if n_genomes is None:
        n_genomes = len(isolate_cols)

    logger.info(
        "Loaded pangenome matrix: %d clusters, %d isolates.",
        len(df), n_genomes,
    )

    # Count presence (non-empty cells) per cluster
    presence = df[isolate_cols].notna() & (df[isolate_cols] != "")
    df["N_genomes"] = presence.sum(axis=1).astype(int)
    df["Freq"] = df["N_genomes"] / n_genomes
    df["Category"] = df["Freq"].apply(assign_category)
    df["Is_hypothetical"] = df["Annotation"].fillna("").apply(is_hypothetical)
    df = df.rename(columns={"Gene": "Cluster"})

    return df, isolate_cols, n_genomes


# ---------------------------------------------------------------------------
# Individual mode integration
# ---------------------------------------------------------------------------

def load_bioentro_results(bioentro_dir: Path) -> Dict[str, pd.DataFrame]:
    """Load all bioentro protein TSV files from a directory.

    Args:
        bioentro_dir: Directory containing *_protein.tsv files.

    Returns:
        Dict mapping isolate_name → DataFrame.

    Raises:
        InputError: If directory is empty or no TSV files found.
    """
    if not bioentro_dir.is_dir():
        raise InputError(f"Not a directory: '{bioentro_dir}'")

    tsv_files = list(bioentro_dir.glob("*_protein.tsv"))
    if not tsv_files:
        # Also try without suffix convention
        tsv_files = list(bioentro_dir.glob("*.tsv"))

    if not tsv_files:
        raise InputError(
            f"No bioentro TSV files found in '{bioentro_dir}'. "
            "Run 'pangentro run-bioentro' first."
        )

    results: Dict[str, pd.DataFrame] = {}
    for tsv in sorted(tsv_files):
        # Derive isolate name from filename
        isolate = tsv.stem.replace("_protein", "")
        try:
            df = pd.read_csv(tsv, sep="\t")
            if "Prot_ID" not in df.columns or "IPS" not in df.columns:
                logger.warning(
                    "Skipping '%s': missing Prot_ID or IPS column.", tsv.name
                )
                continue
            results[isolate] = df
            logger.debug("Loaded bioentro results for '%s': %d proteins.", isolate, len(df))
        except Exception as exc:
            logger.warning("Cannot read '%s': %s", tsv.name, exc)

    logger.info("Loaded bioentro results for %d isolates.", len(results))
    return results


def _match_cluster_to_bioentro(
    cluster_row: pd.Series,
    isolate_cols: List[str],
    bioentro_data: Dict[str, pd.DataFrame],
) -> Tuple[float, float, int]:
    """Find IPS values for all isolate representatives of a cluster.

    For each isolate, looks up the locus tag in that isolate's bioentro TSV.

    Returns:
        (mean_IPS, std_IPS, n_found)
    """
    ips_values = []

    for isolate_col in isolate_cols:
        locus_tag = cluster_row.get(isolate_col)
        if pd.isna(locus_tag) or locus_tag == "":
            continue

        # The isolate column name should match the key in bioentro_data
        # Try exact match first, then prefix match
        df_isolate = bioentro_data.get(isolate_col)
        if df_isolate is None:
            # Try to find a matching key
            matches = [k for k in bioentro_data if isolate_col in k or k in isolate_col]
            if matches:
                df_isolate = bioentro_data[matches[0]]

        if df_isolate is None:
            continue

        # Look up the locus tag — Panaroo may list multiple loci separated by ;
        tags = str(locus_tag).split(";")
        for tag in tags:
            tag = tag.strip()
            match = df_isolate[df_isolate["Prot_ID"] == tag]
            if not match.empty:
                ips_values.append(float(match["IPS"].iloc[0]))
                break

    if not ips_values:
        return -1.0, 0.0, 0

    return (
        round(float(np.mean(ips_values)), 6),
        round(float(np.std(ips_values)), 6),
        len(ips_values),
    )


# ---------------------------------------------------------------------------
# Pangenomic mode — IPS recalculation against pan-proteome background
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Pangenomic mode — IPS recalculation against CORE GENOME background
# ---------------------------------------------------------------------------

def build_core_background(
    pan_reference_faa: Path,
    core_cluster_ids: set,
    k: int = DEFAULT_K_PROTEIN,
) -> Dict[str, float]:
    """Build k-mer background distribution from core genome sequences only.

    The core genome represents the stable, conserved compositional identity
    of the species — genes present in >= 99% of isolates under purifying
    selection. Using only core sequences as background ensures that IPS_pan
    measures atypicality relative to what is truly characteristic of the
    organism, not diluted by the compositional diversity of rare shell/cloud
    genes that are part of the signal we want to detect.

    Biological rationale:
        If we used the full pan-proteome as background (v0.1.0), shell and
        cloud sequences contributed to their own background, reducing their
        apparent JSD. By using only core sequences, we define a stable
        reference frame against which accessory genes can be meaningfully
        compared.

    Args:
        pan_reference_faa: Path to pan_genome_reference.fa from Panaroo.
                           Contains one representative sequence per cluster.
        core_cluster_ids:  Set of cluster IDs classified as core
                           (freq >= CORE_THRESHOLD). Built from the
                           gene_presence_absence.csv matrix.
        k:                 K-mer size (must match bioentro run).

    Returns:
        Normalized k-mer frequency distribution of the core proteome.

    Raises:
        InputError: If file is missing, empty, or no core sequences found.
    """
    if not pan_reference_faa.exists():
        raise InputError(f"Pan-proteome reference not found: '{pan_reference_faa}'")

    logger.info(
        "Building CORE genome background from '%s' (k=%d, %d core clusters)...",
        pan_reference_faa.name, k, len(core_cluster_ids),
    )

    core_seq_parts: List[str] = []
    n_found = 0
    n_skipped = 0

    for rec in SeqIO.parse(pan_reference_faa, "fasta"):
        if rec.id in core_cluster_ids:
            core_seq_parts.append(str(rec.seq).upper())
            n_found += 1
        else:
            n_skipped += 1

    if not core_seq_parts:
        raise InputError(
            f"No core sequences found in '{pan_reference_faa.name}'. "
            f"Checked {n_found + n_skipped} sequences against "
            f"{len(core_cluster_ids)} core cluster IDs. "
            "Verify that cluster IDs in the FASTA match those in "
            "gene_presence_absence.csv."
        )

    all_seq = "".join(core_seq_parts)
    bg = _compute_kmer_dist(all_seq, k)

    logger.info(
        "Core background: %d core sequences used, %d skipped (shell/cloud), "
        "%d residues, %d unique %d-mers.",
        n_found, n_skipped, len(all_seq), len(bg), k,
    )
    return bg


def calculate_pangenomic_ips(
    sequence: str,
    bg_dist: Dict[str, float],
    k: int = DEFAULT_K_PROTEIN,
) -> Tuple[float, float, float, float, float, float]:
    """Calculate IPS for one sequence against the pan-proteome background.

    Returns:
        (IPS, H_real, Efficiency, JSD, sqrt_JSD, Kolmogorov)
    """
    p_dist = _compute_kmer_dist(sequence, k)
    h_max = _calculate_h_max(k)
    h_real = _calculate_shannon(p_dist)
    efficiency = _safe_efficiency(h_real, h_max)
    jsd = _calculate_jsd(p_dist, bg_dist)
    sqrt_jsd = round(math.sqrt(jsd), 6)
    kolmogorov = _calculate_kolmogorov(sequence)
    ips = _calculate_ips(efficiency, jsd, len(sequence))
    return ips, h_real, efficiency, jsd, sqrt_jsd, kolmogorov


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def write_tsv(results: list, output_path: Path) -> None:
    """Write a list of dataclass instances to a TSV file."""
    if not results:
        raise ValueError("Nothing to write — results list is empty.")
    fieldnames = list(asdict(results[0]).keys())
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for item in results:
            writer.writerow(asdict(item))
    logger.info("Written %d record(s) to '%s'.", len(results), output_path)


# ---------------------------------------------------------------------------
# run-bioentro subcommand
# ---------------------------------------------------------------------------

def cmd_run_bioentro(
    faa_dir: Path,
    output_dir: Path,
    k: int = DEFAULT_K_PROTEIN,
    bioentro_cmd: str = "bioentro",
) -> None:
    """Run bioentro protein mode on all .faa files in a directory.

    Args:
        faa_dir:     Directory containing one .faa per isolate.
        output_dir:  Directory for bioentro TSV outputs.
        k:           K-mer size for bioentro.
        bioentro_cmd: Path or name of the bioentro executable.

    Raises:
        InputError: If faa_dir is empty or bioentro is not found.
    """
    if not faa_dir.is_dir():
        raise InputError(f"Not a directory: '{faa_dir}'")

    faa_files = sorted(faa_dir.glob("**/*.faa"))
    if not faa_files:
        raise InputError(f"No .faa files found in '{faa_dir}'.")

    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info("Running bioentro on %d isolates...", len(faa_files))

    ok, skipped, failed = 0, 0, 0
    for faa in faa_files:
        # Derive isolate name from directory or filename
        isolate = faa.parent.name if faa.parent != faa_dir else faa.stem
        out_tsv = output_dir / f"{isolate}_protein.tsv"

        if out_tsv.exists():
            logger.debug("Skipping '%s' — output already exists.", isolate)
            skipped += 1
            continue

        cmd = [
            bioentro_cmd,
            "-i", str(faa),
            "-m", "protein",
            "-k", str(k),
            "-o", str(out_tsv),
        ]

        logger.info("  %s...", isolate)
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            logger.error(
                "bioentro failed for '%s':\n%s", isolate, result.stderr
            )
            failed += 1
        else:
            ok += 1

    print(f"\nrun-bioentro complete:")
    print(f"  Processed : {ok}")
    print(f"  Skipped   : {skipped}  (output already existed)")
    print(f"  Failed    : {failed}")
    print(f"  Output    : '{output_dir}'")


# ---------------------------------------------------------------------------
# integrate subcommand
# ---------------------------------------------------------------------------

def cmd_integrate(
    presence_absence: Path,
    bioentro_dir: Optional[Path],
    pan_reference: Optional[Path],
    output_path: Path,
    ips_mode: str = "both",
    k: int = DEFAULT_K_PROTEIN,
) -> List[ClusterMetrics]:
    """Cross bioentro results with Panaroo pangenome matrix.

    v0.1.1 change: IPS_pan is now calculated against the CORE GENOME
    background (clusters with freq >= CORE_THRESHOLD) instead of the full
    pan-proteome. This produces higher contrast between core and shell/cloud
    proteins because the background no longer includes the signal being
    measured.

    Args:
        presence_absence: gene_presence_absence.csv from Panaroo.
        bioentro_dir:     Directory with per-isolate bioentro TSVs (individual mode).
        pan_reference:    pan_genome_reference.fa from Panaroo (pangenomic mode).
        output_path:      Output TSV path.
        ips_mode:         'individual', 'pangenomic', or 'both'.
        k:                K-mer size (must match bioentro run).

    Returns:
        List of ClusterMetrics.

    Raises:
        InputError, InsufficientDataError.
    """
    # Load pangenome matrix
    df_pan, isolate_cols, n_genomes = load_presence_absence(presence_absence)
    logger.info("Total clusters: %d | Isolates: %d", len(df_pan), n_genomes)

    # Load bioentro results if needed
    bioentro_data: Dict[str, pd.DataFrame] = {}
    if ips_mode in ("individual", "both"):
        if bioentro_dir is None:
            raise InputError(
                "--bioentro-dir is required for ips-mode 'individual' or 'both'."
            )
        bioentro_data = load_bioentro_results(bioentro_dir)

    # Build core genome background if needed (v0.1.1: core-only background)
    pan_bg: Dict[str, float] = {}
    pan_seqs: Dict[str, str] = {}
    if ips_mode in ("pangenomic", "both"):
        if pan_reference is None:
            raise InputError(
                "--pan-reference is required for ips-mode 'pangenomic' or 'both'."
            )
        # Extract core cluster IDs from the presence-absence matrix
        core_cluster_ids: set = set(
            df_pan.loc[df_pan["Category"] == "core", "Cluster"].tolist()
        )
        if not core_cluster_ids:
            logger.warning(
                "No core clusters found with threshold %.2f. "
                "Falling back to full pan-proteome background.",
                CORE_THRESHOLD,
            )
            # Fallback: use all sequences (v0.1.0 behavior)
            core_cluster_ids = set(df_pan["Cluster"].tolist())

        pan_bg = build_core_background(pan_reference, core_cluster_ids, k)
        # Index ALL sequences by cluster name for IPS calculation
        for rec in SeqIO.parse(pan_reference, "fasta"):
            pan_seqs[rec.id] = str(rec.seq).upper()

    results: List[ClusterMetrics] = []
    n_total = len(df_pan)

    for idx, row in df_pan.iterrows():
        if idx % 500 == 0:
            logger.info("  Processing cluster %d / %d...", idx, n_total)

        cluster = str(row["Cluster"])
        annotation = str(row.get("Annotation", "")) or ""
        category = str(row["Category"])
        n_gen = int(row["N_genomes"])
        freq = round(float(row["Freq"]), 4)
        is_hypo = bool(row["Is_hypothetical"])

        # ── Individual mode ───────────────────────────────────────────────
        mean_ips_i, std_ips_i, n_found = -1.0, 0.0, 0
        if ips_mode in ("individual", "both") and bioentro_data:
            mean_ips_i, std_ips_i, n_found = _match_cluster_to_bioentro(
                row, isolate_cols, bioentro_data
            )

        # ── Pangenomic mode ───────────────────────────────────────────────
        ips_pan = -1.0
        h_real_pan = efficiency_pan = jsd_pan = sqrt_jsd_pan = kolmo_pan = -1.0

        if ips_mode in ("pangenomic", "both") and pan_bg:
            # Find representative sequence for this cluster
            seq = pan_seqs.get(cluster)
            if seq:
                (ips_pan, h_real_pan, efficiency_pan,
                 jsd_pan, sqrt_jsd_pan, kolmo_pan) = calculate_pangenomic_ips(
                    seq, pan_bg, k
                )

        results.append(ClusterMetrics(
            Cluster=cluster,
            Annotation=annotation,
            Category=category,
            N_genomes=n_gen,
            Freq=freq,
            Is_hypothetical=is_hypo,
            Mean_IPS_indiv=mean_ips_i,
            Std_IPS_indiv=std_ips_i,
            N_isolates_indiv=n_found,
            IPS_pan=ips_pan,
            H_real_pan=h_real_pan,
            Efficiency_pan=efficiency_pan,
            JSD_pan=jsd_pan,
            sqrt_JSD_pan=sqrt_jsd_pan,
            Kolmogorov_pan=kolmo_pan,
        ))

    write_tsv(results, output_path)

    # Print quick summary
    df_res = pd.DataFrame([asdict(r) for r in results])
    print(f"\nIntegration complete: {len(results)} clusters")
    for cat in CATEGORY_ORDER:
        sub = df_res[df_res["Category"] == cat]
        if sub.empty:
            continue
        n_hypo = sub["Is_hypothetical"].sum()
        print(f"  {cat:<12} {len(sub):>5} clusters  |  {n_hypo:>4} hypothetical "
              f"({100*n_hypo/len(sub):.1f}%)")
    print(f"\nOutput → '{output_path}'")

    return results


# ---------------------------------------------------------------------------
# compare subcommand
# ---------------------------------------------------------------------------

def _rank_biserial(group_a: np.ndarray, group_b: np.ndarray) -> float:
    """Rank-biserial correlation as effect size for Mann-Whitney U."""
    n_a, n_b = len(group_a), len(group_b)
    if n_a == 0 or n_b == 0:
        return 0.0
    u_stat, _ = stats.mannwhitneyu(group_a, group_b, alternative="two-sided")
    return round(1.0 - (2.0 * u_stat) / (n_a * n_b), 4)


def cmd_compare(
    metrics_path: Path,
    output_dir: Path,
    ips_col: str = "IPS_pan",
    alpha: float = 0.05,
) -> Tuple[List[CategorySummary], List[StatTest]]:
    """Compare IPS distributions across pangenomic categories.

    Performs Kruskal-Wallis test (overall) and pairwise Mann-Whitney U
    tests with rank-biserial effect size. Generates summary TSVs and
    a publication-quality box/violin figure.

    Args:
        metrics_path: Output of 'integrate' subcommand.
        output_dir:   Directory for output files and figures.
        ips_col:      IPS column to use: 'IPS_pan' or 'Mean_IPS_indiv'.
        alpha:        Significance threshold (default: 0.05).

    Returns:
        (list of CategorySummary, list of StatTest)
    """
    if not metrics_path.exists():
        raise InputError(f"File not found: '{metrics_path}'")

    df = pd.read_csv(metrics_path, sep="\t")

    if ips_col not in df.columns:
        raise InputError(
            f"Column '{ips_col}' not found. "
            f"Available: {[c for c in df.columns if 'IPS' in c]}"
        )

    # Filter: only clusters with valid IPS
    df = df[df[ips_col] >= 0].copy()
    if df.empty:
        raise InsufficientDataError(
            f"No valid IPS values in column '{ips_col}'. "
            "Run integrate with the appropriate --ips-mode first."
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Per-category summary ──────────────────────────────────────────────
    summaries: List[CategorySummary] = []
    groups: Dict[str, np.ndarray] = {}

    for cat in CATEGORY_ORDER:
        sub = df[df["Category"] == cat][ips_col].dropna().values
        if len(sub) == 0:
            continue
        groups[cat] = sub
        n_hypo = int(df[(df["Category"] == cat) & df["Is_hypothetical"]][ips_col].ge(0).sum())
        summaries.append(CategorySummary(
            Category=cat,
            N_clusters=len(sub),
            N_hypothetical=n_hypo,
            Pct_hypothetical=round(100 * n_hypo / len(sub), 2) if len(sub) > 0 else 0.0,
            Mean_IPS=round(float(np.mean(sub)), 4),
            Median_IPS=round(float(np.median(sub)), 4),
            Std_IPS=round(float(np.std(sub)), 4),
            Q25_IPS=round(float(np.percentile(sub, 25)), 4),
            Q75_IPS=round(float(np.percentile(sub, 75)), 4),
        ))

    # ── Kruskal-Wallis overall test ───────────────────────────────────────
    all_groups = list(groups.values())
    if len(all_groups) >= 2:
        kw_stat, kw_p = stats.kruskal(*all_groups)
        logger.info(
            "Kruskal-Wallis H=%.4f, p=%.6f (%s)",
            kw_stat, kw_p,
            "significant" if kw_p < alpha else "not significant",
        )
    else:
        kw_stat, kw_p = 0.0, 1.0

    # ── Pairwise Mann-Whitney U tests ─────────────────────────────────────
    stat_tests: List[StatTest] = []
    cat_list = list(groups.keys())
    for i in range(len(cat_list)):
        for j in range(i + 1, len(cat_list)):
            a, b = cat_list[i], cat_list[j]
            arr_a, arr_b = groups[a], groups[b]
            if len(arr_a) < 3 or len(arr_b) < 3:
                continue
            u_stat, p_val = stats.mannwhitneyu(arr_a, arr_b, alternative="two-sided")
            effect = _rank_biserial(arr_a, arr_b)
            stat_tests.append(StatTest(
                Category_A=a,
                Category_B=b,
                Test="mannwhitneyu",
                Statistic=round(float(u_stat), 4),
                P_value=round(float(p_val), 6),
                Significant=p_val < alpha,
                Effect_size=effect,
            ))

    # ── Write tables ──────────────────────────────────────────────────────
    write_tsv(summaries, output_dir / "category_summary.tsv")
    if stat_tests:
        write_tsv(stat_tests, output_dir / "statistical_tests.tsv")

    # ── Figure ────────────────────────────────────────────────────────────
    _draw_comparison_figure(
        groups=groups,
        summaries=summaries,
        stat_tests=stat_tests,
        ips_col=ips_col,
        kw_p=kw_p,
        output_dir=output_dir,
        alpha=alpha,
    )

    # ── Print summary ─────────────────────────────────────────────────────
    print(f"\nIPS comparison across pangenomic categories ({ips_col})")
    print(f"Kruskal-Wallis: H={kw_stat:.4f}, p={kw_p:.6f} "
          f"({'*significant*' if kw_p < alpha else 'not significant'})\n")
    print(f"{'Category':<12} {'N':>6} {'N_hypo':>7} {'Mean_IPS':>10} "
          f"{'Median':>8} {'Std':>8}")
    print("─" * 55)
    for s in summaries:
        print(f"{s.Category:<12} {s.N_clusters:>6} {s.N_hypothetical:>7} "
              f"{s.Mean_IPS:>10.4f} {s.Median_IPS:>8.4f} {s.Std_IPS:>8.4f}")

    if stat_tests:
        print(f"\nPairwise tests (Mann-Whitney U, α={alpha}):")
        for t in stat_tests:
            sig = "✓" if t.Significant else "✗"
            print(f"  {sig} {t.Category_A} vs {t.Category_B}: "
                  f"p={t.P_value:.4f}, effect={t.Effect_size:.3f}")

    print(f"\nOutputs → '{output_dir}'")
    return summaries, stat_tests


def _draw_comparison_figure(
    groups: Dict[str, np.ndarray],
    summaries: List[CategorySummary],
    stat_tests: List[StatTest],
    ips_col: str,
    kw_p: float,
    output_dir: Path,
    alpha: float,
) -> None:
    """Generate boxplot + stripplot figure for IPS by category."""

    fig, axes = plt.subplots(1, 2, figsize=(14, 7))

    cats = list(groups.keys())
    colors = [CATEGORY_COLORS.get(c, "#888888") for c in cats]
    data = [groups[c] for c in cats]
    labels = [f"{c}\n(n={len(groups[c])})" for c in cats]

    # ── Left: box + strip ────────────────────────────────────────────────
    ax = axes[0]
    bp = ax.boxplot(
        data, patch_artist=True, notch=False,
        medianprops={"color": "black", "linewidth": 2},
        whiskerprops={"linewidth": 1.2},
        capprops={"linewidth": 1.2},
        flierprops={"marker": "o", "markersize": 3, "alpha": 0.3},
    )
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    # Overlay strip
    for i, (cat, arr) in enumerate(zip(cats, data)):
        x_jitter = np.random.normal(i + 1, 0.07, size=len(arr))
        ax.scatter(
            x_jitter, arr,
            alpha=0.25, s=8,
            color=CATEGORY_COLORS.get(cat, "#888"),
            zorder=3,
        )

    ax.set_xticks(range(1, len(cats) + 1))
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel("IPS", fontsize=11)
    ax.set_title(
        f"IPS distribution by pangenomic category\n"
        f"Kruskal-Wallis p = {kw_p:.4f}"
        f"{' *' if kw_p < alpha else ''}",
        fontsize=12,
    )
    ax.set_ylim(-0.02, 1.05)
    ax.grid(axis="y", alpha=0.3, linewidth=0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # ── Right: mean IPS + hypothetical % ─────────────────────────────────
    ax2 = axes[1]
    cat_labels = [s.Category for s in summaries]
    mean_ips = [s.Mean_IPS for s in summaries]
    pct_hypo = [s.Pct_hypothetical for s in summaries]
    bar_colors = [CATEGORY_COLORS.get(c, "#888") for c in cat_labels]

    x = np.arange(len(cat_labels))
    width = 0.35

    bars1 = ax2.bar(x - width / 2, mean_ips, width,
                    color=bar_colors, alpha=0.8, label="Mean IPS")
    ax2_twin = ax2.twinx()
    bars2 = ax2_twin.bar(x + width / 2, pct_hypo, width,
                         color=bar_colors, alpha=0.4,
                         hatch="///", label="% hypothetical")

    ax2.set_ylabel("Mean IPS", fontsize=11)
    ax2_twin.set_ylabel("% hypothetical proteins", fontsize=11)
    ax2.set_xticks(x)
    ax2.set_xticklabels(cat_labels, fontsize=10)
    ax2.set_ylim(0, 1.0)
    ax2_twin.set_ylim(0, 100)
    ax2.set_title("Mean IPS and % hypothetical by category", fontsize=12)
    ax2.spines["top"].set_visible(False)
    ax2_twin.spines["top"].set_visible(False)

    # Combined legend
    legend_patches = [
        mpatches.Patch(facecolor=CATEGORY_COLORS.get(c, "#888"),
                       alpha=0.8, label=c)
        for c in cat_labels
    ]
    ax2.legend(handles=legend_patches, loc="upper right", fontsize=9)

    fig.suptitle(
        f"Informational Priority Score — Pangenomic Analysis\n"
        f"pangentro v{__version__} | {ips_col}",
        fontsize=13, y=1.01,
    )
    fig.tight_layout()

    out_png = output_dir / f"ips_pangenome_comparison_{ips_col}.png"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    logger.info("Figure saved to '%s'.", out_png)
    plt.show()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pangentro",
        description=(
            f"pangentro v{__version__} — Pangenomic Integration Tool\n\n"
            "Integrates bioentro informational metrics with Panaroo pangenome\n"
            "output to analyze IPS distribution across core, shell, and cloud genomes."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Typical workflow:\n"
            "  1. pangentro run-bioentro -f bakta_outputs/ -o bioentro_results/\n"
            "  2. pangentro integrate \\\n"
            "         -p panaroo_output/gene_presence_absence.csv \\\n"
            "         -b bioentro_results/ \\\n"
            "         -r panaroo_output/pan_genome_reference.fa \\\n"
            "         --ips-mode both \\\n"
            "         -o pangenome_metrics.tsv\n"
            "  3. pangentro compare -i pangenome_metrics.tsv -o comparison/\n\n"
            "Run 'pangentro <subcommand> --help' for options."
        ),
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )

    sub = parser.add_subparsers(dest="subcommand", metavar="SUBCOMMAND")
    sub.required = True

    def _shared(p: argparse.ArgumentParser) -> None:
        p.add_argument("-v", "--verbose", action="store_true",
                       help="Enable debug-level logging")

    # ── run-bioentro ──────────────────────────────────────────────────────
    rb = sub.add_parser(
        "run-bioentro",
        help="Run bioentro protein mode on all isolates in a directory",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  pangentro run-bioentro -f bakta_outputs/ -o bioentro_results/\n"
            "  pangentro run-bioentro -f bakta_outputs/ -o results/ -k 2\n"
        ),
    )
    _shared(rb)
    rb.add_argument("-f", "--faa-dir", required=True,
                    help="Directory containing .faa files (one per isolate, "
                         "or in per-isolate subdirectories)")
    rb.add_argument("-o", "--output", required=True,
                    help="Output directory for bioentro TSV files")
    rb.add_argument("-k", "--kmer", type=int, default=DEFAULT_K_PROTEIN,
                    help=f"K-mer size (default: {DEFAULT_K_PROTEIN})")
    rb.add_argument("--bioentro-cmd", default="bioentro",
                    help="Path to bioentro executable (default: 'bioentro')")

    # ── integrate ─────────────────────────────────────────────────────────
    ig = sub.add_parser(
        "integrate",
        help="Cross bioentro results with Panaroo pangenome matrix",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "IPS modes:\n"
            "  individual  — averages per-isolate IPS; requires --bioentro-dir\n"
            "  pangenomic  — recalculates against CORE GENOME background (v0.1.1);\n"
            "                requires --pan-reference\n"
            "  both        — calculates both (default; requires both inputs)\n\n"
            "v0.1.1: IPS_pan now uses core genome as background instead of\n"
            "full pan-proteome. This increases contrast between core and\n"
            "shell/cloud proteins.\n\n"
            "Examples:\n"
            "  pangentro integrate -p gene_presence_absence.csv \\\n"
            "      -b bioentro_results/ -r pan_genome_reference.fa \\\n"
            "      --ips-mode both -o pangenome_metrics.tsv\n"
        ),
    )
    _shared(ig)
    ig.add_argument("-p", "--presence-absence", required=True,
                    help="gene_presence_absence.csv from Panaroo")
    ig.add_argument("-b", "--bioentro-dir",
                    help="Directory with per-isolate bioentro TSVs (individual mode)")
    ig.add_argument("-r", "--pan-reference",
                    help="pan_genome_reference.fa from Panaroo (pangenomic mode)")
    ig.add_argument("--ips-mode", default="both",
                    choices=["individual", "pangenomic", "both"],
                    help="IPS calculation mode (default: both)")
    ig.add_argument("-k", "--kmer", type=int, default=DEFAULT_K_PROTEIN,
                    help=f"K-mer size — must match bioentro run (default: {DEFAULT_K_PROTEIN})")
    ig.add_argument("-o", "--output", required=True,
                    help="Output TSV path")

    # ── compare ───────────────────────────────────────────────────────────
    cp = sub.add_parser(
        "compare",
        help="Statistical comparison of IPS across pangenomic categories",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Outputs:\n"
            "  category_summary.tsv     — mean, median, std IPS per category\n"
            "  statistical_tests.tsv    — pairwise Mann-Whitney U results\n"
            "  ips_pangenome_comparison_*.png — publication figure\n\n"
            "Examples:\n"
            "  pangentro compare -i pangenome_metrics.tsv -o comparison/\n"
            "  pangentro compare -i pangenome_metrics.tsv -o comparison/ "
            "--ips-col Mean_IPS_indiv\n"
        ),
    )
    _shared(cp)
    cp.add_argument("-i", "--input", required=True,
                    help="Output of 'integrate' subcommand")
    cp.add_argument("-o", "--output", required=True,
                    help="Output directory for tables and figures")
    cp.add_argument("--ips-col", default="IPS_pan",
                    choices=["IPS_pan", "Mean_IPS_indiv"],
                    help="IPS column to compare (default: IPS_pan)")
    cp.add_argument("--alpha", type=float, default=0.05,
                    help="Significance threshold (default: 0.05)")

    return parser


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> int:
    """CLI entry point.

    Returns:
        0 on success, 1 on handled error, 130 on keyboard interrupt.
    """
    parser = _build_parser()
    args = parser.parse_args(argv)
    _configure_logging(verbose=args.verbose)

    try:
        if args.subcommand == "run-bioentro":
            cmd_run_bioentro(
                faa_dir=Path(args.faa_dir),
                output_dir=Path(args.output),
                k=args.kmer,
                bioentro_cmd=args.bioentro_cmd,
            )

        elif args.subcommand == "integrate":
            cmd_integrate(
                presence_absence=Path(args.presence_absence),
                bioentro_dir=Path(args.bioentro_dir) if args.bioentro_dir else None,
                pan_reference=Path(args.pan_reference) if args.pan_reference else None,
                output_path=Path(args.output),
                ips_mode=args.ips_mode,
                k=args.kmer,
            )

        elif args.subcommand == "compare":
            cmd_compare(
                metrics_path=Path(args.input),
                output_dir=Path(args.output),
                ips_col=args.ips_col,
                alpha=args.alpha,
            )

        return 0

    except PangentroError as exc:
        logger.error("%s", exc)
        return 1
    except KeyboardInterrupt:
        logger.warning("Interrupted by user.")
        return 130


if __name__ == "__main__":
    sys.exit(main())