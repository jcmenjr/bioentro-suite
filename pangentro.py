#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pangentro v0.5.0 — anvi'o-native pangenomic prioritization for bioentro-suite

WHAT THIS DOES
──────────────
Reads the gene-clusters summary produced by `anvi-summarize` on a pan-db and
attaches, to every gene cluster (GC), a sequence-intrinsic Informational
Priority Score (IPS) plus an operational "hypothetical protein" (HP) status and
a pangenomic category (core / soft_core / shell / cloud). The output is built
to be imported straight back into anvi'o 9 ("eunice") as items-level
miscellaneous data, so that in the interactive pangenome you can sort/colour
gene clusters by IPS and HP status and decide which clusters to inspect first.

This is a *biology-blind triage* step. IPS does not predict biological
relevance — it only flags compositionally atypical proteins so that a finite
budget of downstream effort (AlphaFold, Foldseek, manual curation, pan-GWAS)
is spent where it is most likely to pay off. Hypotheses come from the full
pipeline + phenotype overlay, not from IPS alone.

    IPS  =  J'  ·  √JSD  ·  (1 − e^(−L/L₀))

  J'    = H_obs / H_max  ............ Pielou's evenness of the amino-acid
                                      composition (normalised Shannon entropy).
                                      H_max = log2(20) ≈ 4.32 bits for k = 1.
                                      Background-independent. Acts as a
                                      low-complexity filter.
  √JSD  = Jensen-Shannon *distance* against a background proteome. This is the
          primary ranking driver. √JSD (not JSD) is used because it is a true
          metric (Endres & Schindelin 2003). With log base 2, JSD ∈ [0, 1], so
          √JSD ∈ [0, 1] as well.
  L-term= length penalty; sequences much shorter than L₀ (default 100 aa) are
          down-weighted, matching the usual convention for "likely functional"
          ORFs.

NOTE ON k.  At the default k = 1 the score is purely *compositional* (residue
usage), not positional: two proteins with the same amino-acid composition but a
different order get the same IPS. This is deliberate — for single proteins k ≥ 2
makes the per-sequence k-mer distribution sparse and JSD unstable. k > 1 is
supported but flagged as experimental.

Two backgrounds are computed per gene (then averaged per GC):
  IPS_core : √JSD against the concatenated background proteome (genes in the
             categories named by --bg-categories, 'core' by default). One shared
             reference space → values comparable across all GCs.
             Caveat: a core GC's own sequences are part of this background, so
             core GCs are partially self-referential (this *deflates* their
             score, which is conservative for the "HPs concentrate in the
             accessory genome" hypothesis, but it is a real circularity — keep
             it in mind when arguing the distributional result).
  IPS_self : √JSD against each gene's OWN genome proteome. No cross-genome
             circularity; the more defensible primary evidence for the
             distributional hypothesis, but each genome is a slightly different
             reference space.

CHANGES vs v0.3.1 (the last published, Panaroo/Bakta-based release)
───────────────────────────────────────────────────────────────────
  * Pipeline pivot: input is now `anvi-summarize` output, not Panaroo's
    `pan_genome_reference.fa`. Panaroo/Bakta paths are gone; `preparo` is no
    longer upstream of pangentro. (The 0.2.0 working draft was numbered *below*
    0.3.1 — this release restores monotonic versioning.)
  * AMINO-ACID ONLY, with an explicit DNA guard (the classic v0.1.x bug fed
    nucleotides into protein math).
  * Species-agnostic: no hard-coded taxon assumptions.
  * HP status is operational (absence of annotation across declared sources):
        ORFan          — no hit in any declared annotation source
        uncharacterized— has a hit but it is COG category S / a DUF / an
                         "uncharacterized"/"unknown function" description
        annotated      — a real functional annotation
  * Pangenomic categories derived from GC occupancy across genomes; thresholds
    are explicit, validated flags for a sensitivity analysis.
  * Correctness/robustness fixes vs the 0.2.0 draft:
      - Auto-import no longer passes the (non-existent) --just-do-it to
        anvi-import-misc-data; idempotent re-import now routes through
        anvi-delete-misc-data, matching anvi'o's documented overwrite path.
      - anvi-summarize is not blindly re-run into an existing directory.
      - compare: Holm–Bonferroni correction for the pairwise tests, a fixed
        (and documented) rank-biserial sign convention, robust boolean parsing,
        and a seeded figure.
  * New, publication-oriented additions: a Monte-Carlo null model that
    calibrates IPS into a length-controlled empirical p-value and z-score
    (--null-draws); a machine-readable provenance/run report (JSON, with input
    SHA-256); an IPS percentile layer; a per-GC compression-ratio diagnostic; a
    length-bias QC; and a vectorised k = 1 scorer for large pangenomes.

THE NULL MODEL (--null-draws N)
───────────────────────────────
For each scored gene, IPS is calibrated against the null hypothesis "a
compositionally typical protein of the same length". N synthetic proteins of
that length are drawn i.i.d. from the background distribution (sampled directly
as a Multinomial over the 20 residues) and scored the same way; the observed
IPS is then reported as a one-sided empirical p-value and a z-score. Because the
length penalty is identical for a gene and its length-matched null, it cancels
out — so the empirical p isolates *compositional* atypicality and absorbs the
length-dependent sampling noise of short proteins (the very thing the length-
bias QC warns about). NOTE: a residue-*permutation* null is degenerate at k = 1
(shuffling does not change composition), which is why the null is a background
re-draw, not a shuffle of the sequence itself.

Subcommands
───────────
  integrate : summary → IPS + category + HP status per GC. Writes the full
              metrics table, the anvi'o items file, an HP-only priority list,
              an anvi'o collection of the top HP candidates, and a run report.
  compare   : statistics (Kruskal–Wallis + Holm-corrected pairwise Mann–Whitney
              U with rank-biserial effect size) of IPS across categories or HP
              status, plus an optional publication figure.

Requirements
────────────
  pandas >= 1.5, numpy >= 1.23           (always)
  scipy  >= 1.9                          (compare only)
  matplotlib >= 3.6                      (compare figure only)
  anvi'o 8/9                             (optional; only to auto-run
                                          anvi-summarize / for nicer output)

Part of bioentro-suite. License: MIT.
"""

from __future__ import annotations

__version__ = "0.5.0"

# ── standard library ───────────────────────────────────────────────────────
import argparse
import hashlib
import json
import logging
import math
import shutil
import subprocess
import sys
import zlib
from collections import Counter
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

# ── third party (hard deps) ────────────────────────────────────────────────
try:
    import numpy as np
    import pandas as pd
except ImportError:  # pragma: no cover
    sys.stderr.write(
        "pangentro needs numpy and pandas. Install them in your anvi'o env:\n"
        "    pip install numpy pandas\n"
    )
    raise

logger = logging.getLogger("pangentro")


# ===========================================================================
# Constants
# ===========================================================================
AA_ORDER: str = "ACDEFGHIKLMNPQRSTVWY"          # canonical 20, fixed order
PROTEIN_ALPHABET_SIZE: int = len(AA_ORDER)      # 20
VALID_AA: frozenset = frozenset(AA_ORDER)
AA_INDEX: Dict[str, int] = {aa: i for i, aa in enumerate(AA_ORDER)}
LOG2_20: float = math.log2(PROTEIN_ALPHABET_SIZE)

DEFAULT_K: int = 1                  # single residues — avoids sparse JSD
L0_AA: int = 100                    # length scale of the penalty (aa)
ZLIB_LEVEL: int = 9

# Default pangenomic-category thresholds (fraction of genomes a GC occurs in).
# These are NOT neutral — expose them and vary them in a sensitivity analysis.
DEFAULT_CORE: float = 0.99          # >= 99% of genomes
DEFAULT_SOFT_CORE: float = 0.95     # >= 95%
DEFAULT_SHELL: float = 0.15         # >= 15%   (cloud = below this)

CATEGORY_ORDER: List[str] = ["core", "soft_core", "shell", "cloud"]
HP_ORDER: List[str] = ["ORFan", "uncharacterized", "annotated"]

CATEGORY_COLORS: Dict[str, str] = {
    "core": "#4E9AF1", "soft_core": "#1D9E75",
    "shell": "#F4A460", "cloud": "#E74C3C",
}
HP_COLORS: Dict[str, str] = {
    "ORFan": "#8E44AD", "uncharacterized": "#E67E22", "annotated": "#95A5A6",
}

# Structural columns of an anvi-summarize gene-clusters summary. Everything
# that is NOT one of these (and not an *_ACC sibling) is a candidate functional
# annotation source.
SUMMARY_STRUCTURAL_COLS: frozenset = frozenset({
    "unique_id", "gene_cluster_id", "bin_name", "genome_name",
    "gene_callers_id", "aa_sequence", "dna_sequence",
    "num_genomes_gene_cluster_has_hits", "num_genes_in_gene_cluster",
    "max_num_paralogs", "SCG", "functional_homogeneity_index",
    "geometric_homogeneity_index", "combined_homogeneity_index",
})

# Tokens that mark a *description* as uncharacterised. NOTE: bare "putative" /
# "probable" are deliberately NOT here — "putative ABC transporter" is a real
# functional hypothesis, not an HP. We only treat a hit as unknown if the
# description itself says so.
UNCHARACTERIZED_SUBSTRINGS: Tuple[str, ...] = (
    "hypothetical", "uncharacteri", "unknown function",
    "domain of unknown function", "protein of unknown function",
    "unnamed protein",
)
# "DUF" is matched separately as a token (DUF1234) to avoid false positives on
# words that merely contain the letters d-u-f.
import re as _re  # local alias; only used for the DUF token
_DUF_RE = _re.compile(r"\bduf\s*\d", _re.IGNORECASE)

EMPTY_TOKENS: frozenset = frozenset({"", "-", "na", "n/a", "nan", "none", "null"})

# misc-data items layers pangentro writes (used for idempotent overwrite).
ITEMS_LAYER_KEYS: Tuple[str, ...] = (
    "IPS", "IPS_core_max", "efficiency", "sqrt_JSD_core", "IPS_percentile",
    "IPS_core_z", "IPS_emp_p", "pangenome_category", "HP_status",
)

# Sentinel for "not computed".
NA = float("nan")


# ===========================================================================
# anvi'o-flavoured terminal
# ===========================================================================
class Terminal:
    """Delegates to anvi'o's own Run/Progress when available, otherwise mimics
    the look so pangentro feels native inside an anvi'o session but still works
    standalone."""

    def __init__(self, quiet: bool = False):
        self.quiet = quiet
        self._run = None
        self._progress = None
        self._use_color = sys.stdout.isatty()
        self.anvio_version: Optional[str] = None
        try:
            from anvio.terminal import Run, Progress  # type: ignore
            self._run = Run(verbose=not quiet)
            self._progress = Progress(verbose=not quiet)
            self.backend = "anvio"
            try:
                import anvio  # type: ignore
                self.anvio_version = getattr(anvio, "__version__", None)
            except Exception:
                self.anvio_version = None
        except Exception:
            # anvi'o is an optional dependency; falling back is expected.
            self.backend = "builtin"

    # -- colour helpers (builtin backend only) ------------------------------
    def _c(self, text: str, code: str) -> str:
        if not self._use_color:
            return text
        return f"\033[{code}m{text}\033[0m"

    # -- key/value line -----------------------------------------------------
    def info(self, key: str, value, nl_before: int = 0, nl_after: int = 0):
        if self.quiet:
            return
        if self._run is not None:
            self._run.info(str(key), str(value), nl_before=nl_before, nl_after=nl_after)
            return
        print("\n" * nl_before, end="")
        print(f"{self._c(str(key) + ' ', '1;33')}{'.' * max(2, 34 - len(str(key)))} "
              f"{self._c(str(value), '0;36')}")
        print("\n" * nl_after, end="")

    def info_single(self, message: str, nl_before: int = 0, nl_after: int = 0,
                    level: int = 1):
        if self.quiet:
            return
        if self._run is not None:
            self._run.info_single(message, nl_before=nl_before, nl_after=nl_after,
                                  level=level)
            return
        bullet = {0: "", 1: "* ", 2: "    - ", 3: "        > "}.get(level, "* ")
        print("\n" * nl_before, end="")
        print(f"{bullet}{message}")
        print("\n" * nl_after, end="")

    def warning(self, message: str, header: str = "WARNING"):
        # warnings are shown even in quiet mode — they matter
        if self._run is not None:
            self._run.warning(message, header=header)
            return
        bar = "═" * 78
        print(self._c(f"\n{header}", "1;31"))
        print(self._c(bar, "1;31"))
        print(message)
        print(self._c(bar, "1;31") + "\n")

    # -- progress -----------------------------------------------------------
    def progress_new(self, name: str, total: Optional[int] = None):
        if self._progress is not None:
            try:
                self._progress.new(name, progress_total_items=total)
            except TypeError:
                self._progress.new(name)
            return
        self._builtin_progress_name = name
        self._builtin_progress_total = total
        if not self.quiet:
            sys.stderr.write(f"{name} ...")
            sys.stderr.flush()

    def progress_update(self, message: str):
        if self._progress is not None:
            self._progress.update(message)
            return
        if not self.quiet:
            sys.stderr.write(f"\r{getattr(self, '_builtin_progress_name', '')} "
                             f"... {message}        ")
            sys.stderr.flush()

    def progress_end(self):
        if self._progress is not None:
            self._progress.end()
            return
        if not self.quiet:
            sys.stderr.write("\r" + " " * 90 + "\r")
            sys.stderr.flush()


# ===========================================================================
# Exceptions
# ===========================================================================
class PangentroError(Exception):
    """Base class for pangentro errors."""


class InputError(PangentroError):
    """Missing or malformed input."""


class InsufficientDataError(PangentroError):
    """Not enough data for the requested analysis."""


# ===========================================================================
# Result dataclasses
# ===========================================================================
@dataclass
class GeneClusterMetrics:
    gene_cluster: str
    num_genomes: int
    total_genomes: int
    freq: float
    pangenome_category: str
    n_genes: int
    n_genes_scored: int
    mean_length_aa: float
    efficiency: float          # mean Pielou J' over scored members
    mean_compression_ratio: float   # diagnostic only (NOT part of IPS)
    sqrt_JSD_core: float       # mean √JSD vs core background
    IPS_core_mean: float
    IPS_core_max: float
    IPS_core_std: float
    IPS_self_mean: float
    IPS_core_z: float          # mean IPS z-score vs length-matched bg null
    IPS_core_emp_p: float      # median one-sided empirical p (NaN w/o --null-draws)
    IPS_percentile: float      # percentile rank (0–100) of the selected IPS col
    HP_status: str
    is_hypothetical: bool
    consensus_function: str
    consensus_cog_category: str
    sources_with_hit: str


@dataclass
class GroupSummary:
    group: str
    n_clusters: int
    n_hypothetical: int
    pct_hypothetical: float
    mean_IPS: float
    median_IPS: float
    std_IPS: float
    q25_IPS: float
    q75_IPS: float


@dataclass
class StatTest:
    group_a: str
    group_b: str
    test: str
    statistic: float
    p_value: float
    p_adjusted: float          # Holm–Bonferroni across all pairwise tests
    significant: bool          # based on p_adjusted < alpha
    effect_size: float         # rank-biserial; > 0 ⇒ group_a tends to higher IPS


# ===========================================================================
# Informational metric functions  (AMINO ACIDS ONLY)
#
# These scalar/dict implementations are the readable reference. The default
# k = 1 path is computed by the vectorised scorer below; a unit test checks the
# two agree. k > 1 uses these directly.
# ===========================================================================
def kmer_counter(sequence: str, k: int) -> Counter:
    """Count valid amino-acid k-mers. A k-mer is valid only if every residue is
    one of the canonical 20 (whitelist). Silently drops X/B/Z/J/U/O/*/-, gaps
    and anything non-standard."""
    n = len(sequence)
    if n < k:
        return Counter()
    if k == 1:
        return Counter(c for c in sequence if c in VALID_AA)
    counts: Counter = Counter()
    for i in range(n - k + 1):
        window = sequence[i:i + k]
        if all(c in VALID_AA for c in window):
            counts[window] += 1
    return counts


def normalise(counts: Counter) -> Dict[str, float]:
    total = sum(counts.values())
    if total == 0:
        return {}
    return {kmer: c / total for kmer, c in counts.items()}


def shannon(dist: Dict[str, float]) -> float:
    """H(P) in bits."""
    if not dist:
        return 0.0
    return -sum(p * math.log2(p) for p in dist.values() if p > 0.0)


def h_max(k: int) -> float:
    """Maximum Shannon entropy for amino-acid k-mers: k · log2(20)."""
    return k * LOG2_20


def pielou_evenness(dist: Dict[str, float], k: int) -> float:
    """J' = H_obs / H_max, clamped to [0, 1]."""
    hm = h_max(k)
    if hm <= 0.0:
        return 0.0
    return min(1.0, max(0.0, shannon(dist) / hm))


def jensen_shannon_divergence(p: Dict[str, float], q: Dict[str, float]) -> float:
    """JSD(P || Q) in bits, symmetric, bounded [0, 1]."""
    if not p or not q:
        return 0.0
    keys = set(p) | set(q)
    jsd = 0.0
    for key in keys:
        pi = p.get(key, 0.0)
        qi = q.get(key, 0.0)
        m = (pi + qi) / 2.0
        if pi > 0.0:
            jsd += 0.5 * pi * math.log2(pi / m)
        if qi > 0.0:
            jsd += 0.5 * qi * math.log2(qi / m)
    return min(1.0, max(0.0, jsd))


def length_penalty(length: int, l0: int) -> float:
    return 1.0 - math.exp(-length / l0)


def ips(evenness: float, jsd: float, length: int, l0: int) -> float:
    """IPS = J' · √JSD · (1 − e^(−L/L₀))."""
    val = evenness * math.sqrt(max(0.0, jsd)) * length_penalty(length, l0)
    return min(1.0, max(0.0, val))


def compression_ratio(sequence: str) -> float:
    """Compression-ratio proxy for Kolmogorov complexity (diagnostic only).
    Low values flag highly repetitive / low-information sequences."""
    if not sequence:
        return NA
    raw = sequence.encode("ascii", errors="ignore")
    if not raw:
        return NA
    try:
        return len(zlib.compress(raw, level=ZLIB_LEVEL)) / len(raw)
    except zlib.error:  # pragma: no cover
        return NA


# ===========================================================================
# Vectorised k = 1 scorer (numpy over the fixed 20-letter alphabet)
# ===========================================================================
def _aa_count_matrix(seqs: Sequence[str]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (counts[n,20] int64, lengths[n] int64, comp[n] float).
    `lengths` is the full sequence length; `comp` is the compression ratio."""
    n = len(seqs)
    counts = np.zeros((n, PROTEIN_ALPHABET_SIZE), dtype=np.int64)
    lengths = np.zeros(n, dtype=np.int64)
    comp = np.full(n, NA, dtype=np.float64)
    for i, seq in enumerate(seqs):
        lengths[i] = len(seq)
        comp[i] = compression_ratio(seq)
        if not seq:
            continue
        cnt = Counter(seq)  # C-level; far faster than a Python char loop
        row = counts[i]
        for aa, idx in AA_INDEX.items():
            c = cnt.get(aa)
            if c:
                row[idx] = c
    return counts, lengths, comp


def _jsd_rows(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Row-wise JSD (bits) between probability rows P[n,20] and Q, where Q is
    either a single background vector [20] (broadcast) or a per-row matrix
    [n,20]. Result is clipped to [0, 1]. Invalid (all-zero) rows of P give
    meaningless values here and must be masked by the caller."""
    M = (P + Q) / 2.0
    with np.errstate(divide="ignore", invalid="ignore"):
        term_p = np.where(P > 0.0, P * np.log2(np.where(P > 0.0, P / M, 1.0)), 0.0)
        term_q = np.where(Q > 0.0, Q * np.log2(np.where(Q > 0.0, Q / M, 1.0)), 0.0)
    jsd = 0.5 * term_p.sum(axis=1) + 0.5 * term_q.sum(axis=1)
    return np.clip(jsd, 0.0, 1.0)


def _null_ips_by_length(
    uniq_lengths: np.ndarray, q: np.ndarray, l0: int, n_draws: int,
    rng: "np.random.Generator",
) -> Dict[int, Tuple[np.ndarray, float, float]]:
    """For each length L, draw `n_draws` synthetic proteins of length L whose
    residues are i.i.d. from the background distribution `q`, and return the
    sorted null IPS_core distribution plus its mean and std. The composition is
    sampled directly as Multinomial(L, q), which avoids materialising sequences
    and is independent of L in cost.

    A residue *permutation* null is degenerate at k = 1 (shuffling does not
    change composition, hence not IPS). The meaningful null is therefore
    'a compositionally typical protein of the same length', which also absorbs
    the length-dependent sampling noise of short proteins."""
    q = np.asarray(q, dtype=np.float64)
    out: Dict[int, Tuple[np.ndarray, float, float]] = {}
    for L in uniq_lengths:
        Li = int(L)
        draws = rng.multinomial(Li, q, size=n_draws).astype(np.float64)  # (n_draws,20)
        Pn = draws / float(Li)
        with np.errstate(divide="ignore", invalid="ignore"):
            logPn = np.where(Pn > 0.0, np.log2(np.where(Pn > 0.0, Pn, 1.0)), 0.0)
        effn = np.clip(-(Pn * logPn).sum(axis=1) / LOG2_20, 0.0, 1.0)
        jsdn = _jsd_rows(Pn, q)
        lp = 1.0 - math.exp(-Li / l0)
        ipsn = np.clip(effn * np.sqrt(jsdn) * lp, 0.0, 1.0)
        ipsn.sort()
        out[Li] = (ipsn, float(ipsn.mean()), float(ipsn.std()))
    return out


def _score_k1_vectorised(
    seqs: List[str], genomes: List[str], gcs: List[str],
    cat_by_gc: Dict[str, str], bg_set: set, l0: int, term: Terminal,
    null_draws: int = 0, seed: int = 0,
) -> Dict[str, np.ndarray]:
    n = len(seqs)
    term.progress_new("Scoring genes (k=1, vectorised)", total=n)
    term.progress_update("counting residues")
    counts, lengths, comp = _aa_count_matrix(seqs)
    totals = counts.sum(axis=1)
    valid = totals > 0

    P = np.zeros((n, PROTEIN_ALPHABET_SIZE), dtype=np.float64)
    P[valid] = counts[valid] / totals[valid][:, None]

    # Pielou evenness
    with np.errstate(divide="ignore", invalid="ignore"):
        logP = np.where(P > 0.0, np.log2(np.where(P > 0.0, P, 1.0)), 0.0)
    H = -(P * logP).sum(axis=1)
    eff = np.clip(H / LOG2_20, 0.0, 1.0)
    eff[~valid] = NA

    # core / chosen-categories background
    term.progress_update("building backgrounds")
    in_bg = np.array([cat_by_gc.get(g, "cloud") in bg_set for g in gcs]) & valid
    core_counts = counts[in_bg].sum(axis=0).astype(np.float64)
    bg_fellback = False
    if core_counts.sum() == 0:
        bg_fellback = True
        core_counts = counts[valid].sum(axis=0).astype(np.float64)
    q_core = core_counts / core_counts.sum() if core_counts.sum() else core_counts

    # per-genome background
    Q_self = np.zeros((n, PROTEIN_ALPHABET_SIZE), dtype=np.float64)
    genome_series = pd.Series(genomes)
    self_ok = np.zeros(n, dtype=bool)
    for _, idx in genome_series.groupby(genome_series).indices.items():
        sub = counts[idx].sum(axis=0).astype(np.float64)
        tot = sub.sum()
        if tot > 0:
            Q_self[idx] = sub / tot
            self_ok[idx] = True

    term.progress_update("computing IPS")
    lenpen = 1.0 - np.exp(-lengths / float(l0))
    jsd_core = _jsd_rows(P, q_core)
    jsd_self = _jsd_rows(P, Q_self)

    sqrtjsd_core = np.full(n, NA)
    ips_core = np.full(n, NA)
    ips_self = np.full(n, NA)
    sqrtjsd_core[valid] = np.sqrt(jsd_core[valid])
    ips_core[valid] = np.clip(eff[valid] * sqrtjsd_core[valid] * lenpen[valid], 0.0, 1.0)
    sok = valid & self_ok
    ips_self[sok] = np.clip(eff[sok] * np.sqrt(jsd_self[sok]) * lenpen[sok], 0.0, 1.0)
    term.progress_end()

    if bg_fellback:
        term.warning(
            "No sequences fell into the background categories — fell back to the "
            "whole-pangenome background. Consider relaxing --core or widening "
            "--bg-categories.",
            header="EMPTY BACKGROUND",
        )
    term.info("Residues in background", f"{int(core_counts.sum()):,}")

    # Monte-Carlo empirical calibration of IPS_core against the background null.
    ips_core_z = np.full(n, NA)
    ips_core_emp_p = np.full(n, NA)
    if null_draws > 0 and core_counts.sum() > 0:
        term.progress_new("Calibrating IPS against a length-matched null")
        rng = np.random.default_rng(seed)
        uniq = np.unique(lengths[valid])
        term.progress_update(f"{len(uniq):,} unique lengths × {null_draws} draws")
        null_by_len = _null_ips_by_length(uniq, q_core, l0, null_draws, rng)
        for i in np.where(valid)[0]:
            arr, mean, std = null_by_len[int(lengths[i])]
            obs = ips_core[i]
            n_ge = arr.size - int(np.searchsorted(arr, obs, side="left"))
            ips_core_emp_p[i] = (1.0 + n_ge) / (null_draws + 1.0)
            if std > 0.0:
                ips_core_z[i] = (obs - mean) / std
        term.progress_end()
        n_sig = int(np.nansum(ips_core_emp_p < 0.05))
        term.info("Genes with empirical p < 0.05",
                  f"{n_sig:,} / {int(valid.sum()):,}")

    return {
        "_eff": eff, "_len": lengths, "_comp": comp,
        "_sqrtjsd_core": sqrtjsd_core, "_ips_core": ips_core, "_ips_self": ips_self,
        "_ips_core_z": ips_core_z, "_ips_core_emp_p": ips_core_emp_p,
    }


def _score_general(
    seqs: List[str], genomes: List[str], gcs: List[str],
    cat_by_gc: Dict[str, str], bg_set: set, k: int, l0: int, term: Terminal,
    null_draws: int = 0,
) -> Dict[str, np.ndarray]:
    """Dict-based scorer for k > 1 (slower; kept simple and obviously correct).
    The Monte-Carlo null is k = 1 only, so it is skipped here."""
    if null_draws > 0:
        term.warning("The Monte-Carlo null model (--null-draws) is implemented "
                     "for k = 1 only; skipping it for k > 1.",
                     header="NULL MODEL SKIPPED")
    n = len(seqs)
    row_dist: List[Dict[str, float]] = [dict() for _ in range(n)]
    eff = np.full(n, NA)
    lengths = np.zeros(n, dtype=np.int64)
    comp = np.full(n, NA)
    genome_counts: Dict[str, Counter] = {}
    bg_counter: Counter = Counter()

    term.progress_new(f"Scoring genes (k={k}, pass 1/2)", total=n)
    for i in range(n):
        seq = seqs[i]
        lengths[i] = len(seq)
        comp[i] = compression_ratio(seq)
        counts = kmer_counter(seq, k)
        if counts:
            dist = normalise(counts)
            row_dist[i] = dist
            eff[i] = pielou_evenness(dist, k)
            genome_counts.setdefault(genomes[i], Counter()).update(counts)
            if cat_by_gc.get(gcs[i], "cloud") in bg_set:
                bg_counter.update(counts)
        if i % 2000 == 0:
            term.progress_update(f"{i:,} / {n:,}")
    term.progress_end()

    if sum(bg_counter.values()) == 0:
        term.warning(
            "No sequences fell into the background categories — fell back to the "
            "whole-pangenome background.", header="EMPTY BACKGROUND")
        for c in genome_counts.values():
            bg_counter.update(c)
    core_bg = normalise(bg_counter)
    genome_bg = {g: normalise(c) for g, c in genome_counts.items()}
    term.info("k-mers in background", f"{len(core_bg):,}")

    sqrtjsd_core = np.full(n, NA)
    ips_core = np.full(n, NA)
    ips_self = np.full(n, NA)
    term.progress_new(f"Scoring genes (k={k}, pass 2/2)", total=n)
    for i in range(n):
        dist = row_dist[i]
        if not dist:
            continue
        e = eff[i]
        L = int(lengths[i])
        jsd_c = jensen_shannon_divergence(dist, core_bg)
        sqrtjsd_core[i] = math.sqrt(max(0.0, jsd_c))
        ips_core[i] = ips(e, jsd_c, L, l0)
        gbg = genome_bg.get(genomes[i])
        if gbg:
            ips_self[i] = ips(e, jensen_shannon_divergence(dist, gbg), L, l0)
        if i % 2000 == 0:
            term.progress_update(f"{i:,} / {n:,}")
    term.progress_end()

    return {
        "_eff": eff, "_len": lengths, "_comp": comp,
        "_sqrtjsd_core": sqrtjsd_core, "_ips_core": ips_core, "_ips_self": ips_self,
        "_ips_core_z": np.full(n, NA), "_ips_core_emp_p": np.full(n, NA),
    }


def score_genes(
    df: pd.DataFrame, occ: pd.DataFrame, k: int, l0: int,
    bg_categories: Sequence[str], term: Terminal,
    null_draws: int = 0, seed: int = 0,
) -> pd.DataFrame:
    """Add per-gene columns: _eff, _len, _comp, _sqrtjsd_core, _ips_core,
    _ips_self, _ips_core_z, _ips_core_emp_p. Uses the vectorised path for k = 1,
    dict path otherwise."""
    cat_by_gc = dict(zip(occ["gene_cluster_id"], occ["pangenome_category"]))
    seqs = df["aa_sequence"].tolist()
    genomes = df["genome_name"].tolist()
    gcs = df["gene_cluster_id"].tolist()
    bg_set = set(bg_categories)

    if k == 1:
        cols = _score_k1_vectorised(seqs, genomes, gcs, cat_by_gc, bg_set, l0,
                                    term, null_draws, seed)
    else:
        cols = _score_general(seqs, genomes, gcs, cat_by_gc, bg_set, k, l0, term,
                              null_draws)

    out = df.copy()
    for name, arr in cols.items():
        out[name] = arr
    return out


# ===========================================================================
# Input parsing
# ===========================================================================
def load_gene_clusters_summary(path: Path, term: Terminal) -> pd.DataFrame:
    """Load an anvi-summarize gene-clusters summary (.txt or .txt.gz)."""
    if not path.exists():
        raise InputError(f"Gene-clusters summary not found: '{path}'")

    term.progress_new("Reading gene-clusters summary")
    term.progress_update("parsing TSV (this can take a moment for big pangenomes)")
    try:
        df = pd.read_csv(path, sep="\t", compression="infer", low_memory=False,
                         dtype=str)
    except Exception as exc:
        term.progress_end()
        raise InputError(f"Could not read '{path}' as TSV: {exc}") from exc
    term.progress_end()

    required = {"gene_cluster_id", "genome_name", "aa_sequence"}
    missing = required - set(df.columns)
    if missing:
        raise InputError(
            f"'{path}' does not look like an anvi-summarize gene-clusters "
            f"summary. Missing required column(s): {sorted(missing)}.\n"
            f"Generate it with:\n"
            f"    anvi-script-add-default-collection -p PAN.db -C DEFAULT\n"
            f"    anvi-summarize -p PAN.db -g GENOMES.db -C DEFAULT -o SUMMARY/\n"
            f"and point --gene-clusters-summary at "
            f"SUMMARY/*_gene_clusters_summary.txt.gz"
        )

    df["aa_sequence"] = (df["aa_sequence"].fillna("").astype(str)
                         .str.strip().str.rstrip("*").str.upper())
    df["genome_name"] = df["genome_name"].fillna("").astype(str)
    df["gene_cluster_id"] = df["gene_cluster_id"].astype(str)

    df = df[df["aa_sequence"] != ""].reset_index(drop=True)
    if df.empty:
        raise InsufficientDataError("No amino-acid sequences in the summary.")

    _guard_against_nucleotides(df, term)
    return df


def _guard_against_nucleotides(df: pd.DataFrame, term: Terminal) -> None:
    """Catch the classic mistake: DNA sequences exported instead of protein. If
    almost every residue across a sample is in {A,C,G,T,U,N}, these are almost
    certainly nucleotides (real proteins average ~25% A/C/G/T/N residues)."""
    sample = df["aa_sequence"].head(200)
    nt_like = set("ACGTUN")
    total = nt = 0
    for seq in sample:
        for c in seq:
            total += 1
            if c in nt_like:
                nt += 1
        if total > 50_000:
            break
    if total > 0 and (nt / total) > 0.95:
        term.warning(
            "The sequences in this summary look like NUCLEOTIDES, not amino "
            "acids (>95% of residues are A/C/G/T/U/N). pangentro computes IPS "
            "in amino-acid space (H_max = log2(20)); feeding it DNA produces "
            "meaningless scores.\nRe-export protein sequences: make sure you did "
            "NOT pass --report-DNA-sequences to anvi-summarize.",
            header="LOOKS LIKE DNA",
        )


def detect_annotation_sources(
    df: pd.DataFrame,
    user_sources: Optional[List[str]],
    user_category_col: Optional[str],
    term: Terminal,
) -> Tuple[List[str], Optional[str]]:
    """Decide which columns hold functional annotations and which holds COG
    categories. eggNOG-mapper output imported into anvi'o lands under a
    user-chosen source name (plus a `<source>_ACC` sibling), so detection is
    heuristic and overridable."""
    cols = list(df.columns)

    # category column
    category_col = user_category_col
    if category_col is None:
        cands = [c for c in cols if c.upper().endswith("CATEGORY")
                 and not c.upper().endswith("CATEGORY_ACC")]
        cog_cands = [c for c in cands if "COG" in c.upper()]
        category_col = (cog_cands or cands or [None])[0]
    elif category_col not in cols:
        term.warning(f"--category-source '{category_col}' not found; "
                     f"COG-category resolution disabled.")
        category_col = None

    # function-description sources
    if user_sources:
        source_cols = [c for c in user_sources if c in cols]
        unknown = [c for c in user_sources if c not in cols]
        if unknown:
            term.warning(f"These --annotation-sources are not in the summary "
                         f"and will be ignored: {unknown}")
    else:
        acc_bases = {c[:-4] for c in cols if c.endswith("_ACC")}
        source_cols = []
        for c in cols:
            if c in SUMMARY_STRUCTURAL_COLS or c.endswith("_ACC"):
                continue
            if category_col is not None and c == category_col:
                continue
            cu = c.upper()
            if cu.endswith("CATEGORY") or cu.endswith("CLASS"):
                continue
            known = any(tok in cu for tok in
                        ("COG", "KOFAM", "KEGG", "PFAM", "EGGNOG", "EGG_NOG",
                         "EC_NUMBER", "GO_TERM", "CAZY", "AMRFINDER",
                         "PRODUCT", "FUNCTION", "ANNOTATION", "DESCRIPTION"))
            if c in acc_bases or known:
                source_cols.append(c)

    if not source_cols:
        term.warning(
            "No functional annotation sources were detected in the summary. "
            "Every gene cluster will be called an ORFan, which makes the HP "
            "split meaningless. Did you import functions (e.g. eggNOG-mapper, "
            "COG20, Pfam, KOfam) into the contigs-dbs BEFORE building the "
            "genomes storage? You can also name the columns explicitly with "
            "--annotation-sources.",
            header="NO ANNOTATION SOURCES",
        )

    term.info("Functional annotation sources",
              ", ".join(source_cols) if source_cols else "none")
    term.info("COG-category column", category_col or "none")
    return source_cols, category_col


# ===========================================================================
# Per-GC occupancy and category
# ===========================================================================
def assign_category(freq: float, core: float, soft_core: float,
                    shell: float) -> str:
    if freq >= core:
        return "core"
    if freq >= soft_core:
        return "soft_core"
    if freq >= shell:
        return "shell"
    return "cloud"


def compute_occupancy(
    df: pd.DataFrame, core: float, soft_core: float, shell: float,
) -> Tuple[pd.DataFrame, int]:
    """Return per-GC occupancy table and the total number of genomes."""
    total_genomes = df["genome_name"].nunique()
    if total_genomes == 0:
        raise InsufficientDataError("No genomes found in the summary.")

    occ = (df.groupby("gene_cluster_id")["genome_name"]
           .nunique().rename("num_genomes").reset_index())
    occ["total_genomes"] = total_genomes
    occ["freq"] = occ["num_genomes"] / total_genomes
    occ["pangenome_category"] = occ["freq"].apply(
        lambda f: assign_category(f, core, soft_core, shell))
    return occ, total_genomes


# ===========================================================================
# HP classification
# ===========================================================================
def _is_empty(value) -> bool:
    if value is None:
        return True
    if isinstance(value, float) and math.isnan(value):
        return True
    return str(value).strip().lower() in EMPTY_TOKENS


def _looks_uncharacterized(description: str) -> bool:
    d = description.lower()
    if any(tok in d for tok in UNCHARACTERIZED_SUBSTRINGS):
        return True
    return bool(_DUF_RE.search(d))


def _most_common_nonempty(series: pd.Series) -> str:
    vals = [str(v).strip() for v in series if not _is_empty(v)]
    if not vals:
        return ""
    return Counter(vals).most_common(1)[0][0]


def classify_hp_per_gc(
    df: pd.DataFrame,
    source_cols: List[str],
    category_col: Optional[str],
    weak_categories: frozenset,
) -> pd.DataFrame:
    """Per GC, build a consensus annotation per source (majority vote of
    non-empty member values) and derive HP status (ORFan / uncharacterized /
    annotated)."""
    grp = df.groupby("gene_cluster_id")

    consensus: Dict[str, pd.Series] = {
        col: grp[col].apply(_most_common_nonempty) for col in source_cols
    }
    cat_consensus = (grp[category_col].apply(_most_common_nonempty)
                     if category_col else None)

    rows = []
    for gc in grp.groups.keys():
        per_source = {col: consensus[col].loc[gc] for col in source_cols}
        hits = {col: v for col, v in per_source.items() if not _is_empty(v)}
        has_any_hit = len(hits) > 0

        strong = any(not _looks_uncharacterized(v) for v in hits.values())

        cog_cat = cat_consensus.loc[gc] if cat_consensus is not None else ""
        cog_cat = "" if _is_empty(cog_cat) else str(cog_cat).strip()
        cat_non_weak = bool(cog_cat) and any(
            ch not in weak_categories for ch in cog_cat)

        if not has_any_hit:
            status = "ORFan"
        elif strong or cat_non_weak:
            status = "annotated"
        else:
            status = "uncharacterized"

        # human-readable representative function for the priority table
        rep = next((v for v in hits.values() if not _looks_uncharacterized(v)), "")
        if not rep and hits:
            rep = next(iter(hits.values()))

        rows.append({
            "gene_cluster_id": gc,
            "HP_status": status,
            "is_hypothetical": status in ("ORFan", "uncharacterized"),
            "consensus_function": rep,
            "consensus_cog_category": cog_cat,
            "sources_with_hit": ";".join(sorted(hits.keys())),
        })
    return pd.DataFrame(rows)


# ===========================================================================
# Aggregate per-GC metrics
# ===========================================================================
def aggregate_metrics(
    scored: pd.DataFrame, occ: pd.DataFrame, hp: pd.DataFrame, ips_column: str,
) -> pd.DataFrame:
    grp = scored.groupby("gene_cluster_id")
    agg = grp.agg(
        n_genes=("aa_sequence", "size"),
        n_genes_scored=("_ips_core", lambda s: int(s.notna().sum())),
        mean_length_aa=("_len", "mean"),
        efficiency=("_eff", "mean"),
        mean_compression_ratio=("_comp", "mean"),
        sqrt_JSD_core=("_sqrtjsd_core", "mean"),
        IPS_core_mean=("_ips_core", "mean"),
        IPS_core_max=("_ips_core", "max"),
        IPS_core_std=("_ips_core", "std"),
        IPS_self_mean=("_ips_self", "mean"),
        IPS_core_z=("_ips_core_z", "mean"),
        IPS_core_emp_p=("_ips_core_emp_p", "median"),
    ).reset_index()

    metrics = (occ.merge(agg, on="gene_cluster_id", how="left")
               .merge(hp, on="gene_cluster_id", how="left"))

    metrics["IPS_core_std"] = metrics["IPS_core_std"].fillna(0.0)
    metrics["is_hypothetical"] = metrics["is_hypothetical"].fillna(False).astype(bool)
    metrics["HP_status"] = metrics["HP_status"].fillna("ORFan")
    for c in ("consensus_function", "consensus_cog_category", "sources_with_hit"):
        metrics[c] = metrics[c].fillna("")

    # percentile rank (0–100) of the IPS column actually used for ranking
    metrics["IPS_percentile"] = (
        pd.to_numeric(metrics[ips_column], errors="coerce").rank(pct=True) * 100.0)

    for c in ("freq", "mean_length_aa", "efficiency", "mean_compression_ratio",
              "sqrt_JSD_core", "IPS_core_mean", "IPS_core_max", "IPS_core_std",
              "IPS_self_mean", "IPS_core_z", "IPS_core_emp_p", "IPS_percentile"):
        metrics[c] = pd.to_numeric(metrics[c], errors="coerce").round(6)

    metrics = metrics.rename(columns={"gene_cluster_id": "gene_cluster"})
    ordered = [f.name for f in GeneClusterMetrics.__dataclass_fields__.values()]
    return metrics[ordered]


# ===========================================================================
# Writers
# ===========================================================================
def write_metrics(metrics: pd.DataFrame, path: Path, term: Terminal) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(path, sep="\t", index=False)
    term.info("Full per-GC metrics", str(path))


def write_anvio_items(
    metrics: pd.DataFrame, path: Path, ips_column: str, term: Terminal,
) -> None:
    """Write an anvi'o misc-data-items file (one row per gene cluster). The
    first column holds gene-cluster IDs (used as keys); the rest become layers.
    Import with:
        anvi-import-misc-data FILE -p PAN.db --target-data-table items
    """
    items = pd.DataFrame({
        "gene_cluster": metrics["gene_cluster"],
        "IPS": metrics[ips_column],
        "IPS_core_max": metrics["IPS_core_max"],
        "efficiency": metrics["efficiency"],
        "sqrt_JSD_core": metrics["sqrt_JSD_core"],
        "IPS_percentile": metrics["IPS_percentile"],
        "IPS_core_z": metrics["IPS_core_z"],
        "IPS_emp_p": metrics["IPS_core_emp_p"],
        "pangenome_category": metrics["pangenome_category"],
        "HP_status": metrics["HP_status"],
    })
    # numerical layers: empty cells for NaN so anvi'o reads them as missing
    for col in ("IPS", "IPS_core_max", "efficiency", "sqrt_JSD_core",
                "IPS_percentile", "IPS_core_z", "IPS_emp_p"):
        items[col] = items[col].map(lambda x: "" if pd.isna(x) else f"{float(x):.6f}")
    path.parent.mkdir(parents=True, exist_ok=True)
    items.to_csv(path, sep="\t", index=False)
    term.info("anvi'o items file", str(path))


def write_hp_priority(
    metrics: pd.DataFrame, path: Path, ips_column: str, term: Terminal,
) -> int:
    hp = metrics[metrics["is_hypothetical"]].copy()
    hp = hp.sort_values(by=ips_column, ascending=False, kind="mergesort")
    path.parent.mkdir(parents=True, exist_ok=True)
    hp.to_csv(path, sep="\t", index=False)
    term.info("HP priority list", f"{path}  ({len(hp):,} hypothetical GCs)")
    return len(hp)


def write_collection(
    metrics: pd.DataFrame, path: Path, ips_column: str, top: int, term: Terminal,
) -> int:
    """Write an anvi'o collection of the top-`top` HP gene clusters (by IPS),
    binned by HP_status. Import with:
        anvi-import-collection FILE -p PAN.db -C pangentro_HP
    """
    hp = metrics[metrics["is_hypothetical"]].copy()
    hp = hp.sort_values(by=ips_column, ascending=False, kind="mergesort")
    if top and top > 0:
        hp = hp.head(top)
    if hp.empty:
        term.warning("No hypothetical gene clusters to put in a collection.")
        return 0
    coll = hp[["gene_cluster", "HP_status"]].copy()
    path.parent.mkdir(parents=True, exist_ok=True)
    coll.to_csv(path, sep="\t", index=False, header=False)
    term.info("anvi'o collection (top HP)", f"{path}  ({len(coll):,} GCs)")
    return len(coll)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_run_report(
    path: Path, args: argparse.Namespace, summary_path: Path, df: pd.DataFrame,
    occ: pd.DataFrame, metrics: pd.DataFrame, source_cols: List[str],
    category_col: Optional[str], term: Terminal, length_bias: float,
) -> None:
    """Machine-readable provenance + run summary (publication reproducibility)."""
    report = {
        "tool": "pangentro",
        "version": __version__,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "command_line": " ".join(sys.argv),
        "input_summary": str(summary_path),
        "input_sha256": _sha256(summary_path),
        "environment": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "anvio_backend": term.backend,
            "anvio_version": term.anvio_version,
        },
        "parameters": {
            "k": args.k, "l0": args.l0, "ips_column": args.ips_column,
            "core": args.core, "soft_core": args.soft_core, "shell": args.shell,
            "bg_categories": list(args.bg_categories),
            "weak_cog_categories": list(args.weak_cog_categories),
            "null_draws": args.null_draws, "seed": args.seed,
            "detected_annotation_sources": source_cols,
            "detected_category_column": category_col,
        },
        "counts": {
            "gene_calls": int(len(df)),
            "genomes": int(df["genome_name"].nunique()),
            "gene_clusters": int(df["gene_cluster_id"].nunique()),
            "pangenome_category": occ["pangenome_category"].value_counts().to_dict(),
            "HP_status": metrics["HP_status"].value_counts().to_dict(),
        },
        "qc": {
            "length_bias_spearman_rho": (None if length_bias is None
                                         else round(float(length_bias), 4)),
        },
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, ensure_ascii=False)
    term.info("Run report (provenance)", str(path))


# ===========================================================================
# Optional: drive anvi'o directly
# ===========================================================================
def run_anvi_summarize(
    pan_db: Path, genomes_db: Path, out_dir: Path, term: Terminal,
) -> Path:
    """Create a default collection if needed and run anvi-summarize, returning
    the path to the gene-clusters summary. If `out_dir` already contains a
    summary, reuse it rather than letting anvi-summarize abort on an existing
    directory."""
    existing = (sorted(out_dir.glob("*gene_clusters_summary.txt.gz")) +
                sorted(out_dir.glob("*gene_clusters_summary.txt")))
    if existing:
        term.info_single(f"Reusing existing summary in '{out_dir}'.", level=1)
        return existing[0]

    for tool in ("anvi-script-add-default-collection", "anvi-summarize"):
        if shutil.which(tool) is None:
            raise InputError(
                f"'{tool}' is not on PATH. Activate your anvi'o environment, "
                f"or run anvi-summarize yourself and pass "
                f"--gene-clusters-summary.")
    if out_dir.exists():
        raise InputError(
            f"Summary dir '{out_dir}' exists but has no gene-clusters summary. "
            f"anvi-summarize refuses to write into an existing directory — "
            f"remove it or pass a fresh --summary-dir.")

    term.info_single("Adding a default collection (idempotent) ...", level=1)
    subprocess.run(
        ["anvi-script-add-default-collection", "-p", str(pan_db), "-C", "DEFAULT"],
        check=False)
    term.info_single("Running anvi-summarize ...", level=1)
    subprocess.run(
        ["anvi-summarize", "-p", str(pan_db), "-g", str(genomes_db),
         "-C", "DEFAULT", "-o", str(out_dir)],
        check=True)
    hits = (sorted(out_dir.glob("*gene_clusters_summary.txt.gz")) +
            sorted(out_dir.glob("*gene_clusters_summary.txt")))
    if not hits:
        raise InputError(f"anvi-summarize ran but no gene-clusters summary was "
                         f"found in '{out_dir}'.")
    return hits[0]


def maybe_import_into_anvio(
    items_file: Path, collection_file: Path, pan_db: Path,
    collection_name: str, overwrite: bool, term: Terminal,
) -> None:
    if shutil.which("anvi-import-misc-data") is None:
        term.warning("anvi-import-misc-data not on PATH; skipping --import.")
        return

    if overwrite and shutil.which("anvi-delete-misc-data"):
        # anvi-import-misc-data has no --just-do-it; the documented overwrite
        # path is to delete the colliding keys first.
        term.info_single("Clearing any previous pangentro items layers ...", level=1)
        subprocess.run(
            ["anvi-delete-misc-data", "-p", str(pan_db),
             "--target-data-table", "items",
             "--keys-to-remove", ",".join(ITEMS_LAYER_KEYS)],
            check=False)

    term.info_single("Importing items data into the pan-db ...", level=1)
    try:
        subprocess.run(
            ["anvi-import-misc-data", str(items_file), "-p", str(pan_db),
             "--target-data-table", "items"],
            check=True)
    except subprocess.CalledProcessError:
        term.warning(
            "anvi-import-misc-data failed — most likely some of these layer "
            "keys already exist. Re-run with --overwrite, or clear them first:\n"
            f"    anvi-delete-misc-data -p {pan_db} --target-data-table items "
            f"--keys-to-remove {','.join(ITEMS_LAYER_KEYS)}",
            header="IMPORT FAILED")
        return

    if collection_file.exists() and shutil.which("anvi-import-collection"):
        term.info_single("Importing the top-HP collection ...", level=1)
        subprocess.run(
            ["anvi-import-collection", str(collection_file), "-p", str(pan_db),
             "-C", collection_name],
            check=False)


# ===========================================================================
# Welcome banner
# ===========================================================================
def welcome(term: Terminal) -> None:
    if term.quiet:
        return
    line = "═" * 78
    term.info_single(line, level=0)
    term.info_single(f"pangentro v{__version__}  ·  IPS × anvi'o pangenomics", level=0)
    term.info_single("Biology-blind triage of hypothetical proteins by "
                     "compositional atypicality.", level=0)
    term.info_single("IPS = J' · √JSD · (1 − e^(−L/L₀))", level=0)
    term.info_single(line, level=0, nl_after=1)


# ===========================================================================
# Argument validation
# ===========================================================================
def _validate_integrate_args(args: argparse.Namespace) -> None:
    if args.k < 1:
        raise InputError("--k must be >= 1.")
    if args.k > 1:
        logger.warning("k > 1 is experimental: per-sequence k-mer distributions "
                       "are sparse and JSD becomes unstable for short proteins.")
    if args.l0 < 1:
        raise InputError("--l0 must be >= 1.")
    if args.null_draws < 0:
        raise InputError("--null-draws must be >= 0 (0 disables the null model).")
    for name, val in (("core", args.core), ("soft-core", args.soft_core),
                      ("shell", args.shell)):
        if not (0.0 <= val <= 1.0):
            raise InputError(f"--{name} must be a fraction in [0, 1] (got {val}).")
    if not (args.shell <= args.soft_core <= args.core):
        raise InputError(
            f"Category thresholds must satisfy shell <= soft_core <= core "
            f"(got shell={args.shell}, soft_core={args.soft_core}, "
            f"core={args.core}).")
    unknown_bg = [c for c in args.bg_categories if c not in CATEGORY_ORDER]
    if unknown_bg:
        raise InputError(
            f"--bg-categories contains unknown categories {unknown_bg}; "
            f"choose from {CATEGORY_ORDER}.")
    if args.top_hp < 0:
        raise InputError("--top-hp must be >= 0 (0 = all HP).")


# ===========================================================================
# Subcommand: integrate
# ===========================================================================
def cmd_integrate(args: argparse.Namespace, term: Terminal) -> int:
    welcome(term)
    _validate_integrate_args(args)

    if args.gene_clusters_summary:
        summary_path = Path(args.gene_clusters_summary)
    elif args.pan_db and args.genomes_storage:
        summary_path = run_anvi_summarize(
            Path(args.pan_db), Path(args.genomes_storage),
            Path(args.summary_dir), term)
    else:
        raise InputError(
            "Provide either --gene-clusters-summary FILE, or both --pan-db and "
            "--genomes-storage (to let pangentro run anvi-summarize).")

    term.info("Input summary", str(summary_path))
    term.info("k (k-mer size)", args.k)
    term.info("L0 (length scale, aa)", args.l0)
    term.info("Category thresholds (core/soft/shell)",
              f"{args.core} / {args.soft_core} / {args.shell}")
    term.info("Background categories", ", ".join(args.bg_categories))
    term.info("Null-model draws", args.null_draws or "off")
    term.info("Layout / ranking IPS column", args.ips_column)

    df = load_gene_clusters_summary(summary_path, term)
    term.info("Genes (gene calls) in summary", f"{len(df):,}", nl_before=1)
    term.info("Genomes", f"{df['genome_name'].nunique():,}")
    term.info("Gene clusters", f"{df['gene_cluster_id'].nunique():,}")

    source_cols, category_col = detect_annotation_sources(
        df, args.annotation_sources, args.category_source, term)

    occ, total_genomes = compute_occupancy(
        df, args.core, args.soft_core, args.shell)

    # report the partition
    counts = occ["pangenome_category"].value_counts().to_dict()
    singletons = int((occ["num_genomes"] == 1).sum())
    singleton_cat = (assign_category(1.0 / total_genomes, args.core,
                                     args.soft_core, args.shell)
                     if total_genomes else "cloud")
    term.info_single("Pangenomic partition:", level=1, nl_before=1)
    for cat in CATEGORY_ORDER:
        term.info_single(f"{cat:<10} {counts.get(cat, 0):>6} gene clusters", level=2)
    if singletons:
        term.info_single(f"(of which {singletons:,} are singletons → "
                         f"category '{singleton_cat}' at these thresholds)", level=2)
    if total_genomes < 10 and counts.get("cloud", 0) == 0 and singletons:
        term.info_single(
            f"note: with only {total_genomes} genomes a singleton sits at "
            f"freq={1.0/total_genomes:.2f}, ≥ --shell ({args.shell}), so the "
            f"'cloud' bin is empty. Lower --shell if you want singletons "
            f"separated out (one for the sensitivity analysis).", level=2)

    scored = score_genes(df, occ, args.k, args.l0, args.bg_categories, term,
                         null_draws=args.null_draws, seed=args.seed)
    hp = classify_hp_per_gc(scored, source_cols, category_col,
                            frozenset(args.weak_cog_categories))
    metrics = aggregate_metrics(scored, occ, hp, args.ips_column)

    # ---- HP / IPS summary -------------------------------------------------
    hp_counts = metrics["HP_status"].value_counts().to_dict()
    term.info_single("Hypothetical-protein status (per GC):", level=1, nl_before=1)
    for st in HP_ORDER:
        term.info_single(f"{st:<16} {hp_counts.get(st, 0):>6} gene clusters", level=2)

    hyp = metrics[metrics["is_hypothetical"]]
    if not hyp.empty:
        term.info("Median IPS (hypothetical GCs)",
                  round(float(hyp[args.ips_column].median()), 4), nl_before=1)
        ann = metrics[~metrics["is_hypothetical"]]
        if not ann.empty:
            term.info("Median IPS (annotated GCs)",
                      round(float(ann[args.ips_column].median()), 4))

    # ---- length-bias QC ---------------------------------------------------
    qc = metrics[[args.ips_column, "mean_length_aa"]].apply(pd.to_numeric,
                                                            errors="coerce")
    length_bias = (qc[args.ips_column].corr(qc["mean_length_aa"], method="spearman")
                   if qc.dropna().shape[0] > 2 else None)
    if length_bias is not None:
        term.info("Length-bias (Spearman IPS vs length)", round(float(length_bias), 3))
        if abs(length_bias) > 0.6:
            term.warning(
                f"IPS correlates strongly with length (Spearman ρ={length_bias:.2f}). "
                f"The (1 − e^(−L/L₀)) term may be dominating the ranking; consider "
                f"reporting IPS within length strata, or lowering --l0.",
                header="LENGTH BIAS")

    # ---- write outputs ----------------------------------------------------
    out_dir = Path(args.output)
    proj = args.project_name
    metrics_path = out_dir / f"{proj}_pangentro_metrics.txt"
    items_path = out_dir / f"{proj}_items_for_anvio.txt"
    hp_path = out_dir / f"{proj}_HP_prioritized.txt"
    coll_path = out_dir / f"{proj}_HP_collection.txt"
    report_path = out_dir / f"{proj}_run_report.json"

    term.info_single("Writing outputs:", level=1, nl_before=1)
    write_metrics(metrics, metrics_path, term)
    write_anvio_items(metrics, items_path, args.ips_column, term)
    write_hp_priority(metrics, hp_path, args.ips_column, term)
    write_collection(metrics, coll_path, args.ips_column, args.top_hp, term)
    write_run_report(report_path, args, summary_path, df, occ, metrics,
                     source_cols, category_col, term, length_bias)

    # ---- next steps -------------------------------------------------------
    pan_for_msg = args.pan_db or "PAN.db"
    coll_name = args.collection_name
    term.info_single("Next, decorate your pangenome and explore it in anvi'o:",
                     level=1, nl_before=1)
    term.info_single(f"anvi-import-misc-data {items_path} -p {pan_for_msg} "
                     f"--target-data-table items", level=2)
    term.info_single(f"anvi-import-collection {coll_path} -p {pan_for_msg} "
                     f"-C {coll_name}", level=2)
    term.info_single(f"anvi-display-pan -p {pan_for_msg} -g GENOMES.db", level=2)
    term.info_single("In the interface, order items by IPS and colour by "
                     "HP_status to see the highest-IPS hypothetical clusters; "
                     "open a GC to inspect its aligned homologues.", level=2)

    if args.import_into_anvio:
        if not args.pan_db:
            term.warning("--import needs --pan-db; skipping import.")
        else:
            maybe_import_into_anvio(items_path, coll_path, Path(args.pan_db),
                                    coll_name, args.overwrite, term)

    term.info_single("Done.", level=0, nl_before=1)
    return 0


# ===========================================================================
# Subcommand: compare
# ===========================================================================
def _coerce_bool(series: pd.Series) -> pd.Series:
    """Robustly parse an is_hypothetical column that may come back from a TSV as
    bools, ints, or the strings 'True'/'False'/'1'/'0'."""
    if series.dtype == bool:
        return series.fillna(False)
    return (series.astype(str).str.strip().str.lower()
            .isin({"true", "1", "yes", "t"}))


def _rank_biserial(u1: float, na: int, nb: int) -> float:
    """Kerby rank-biserial correlation from the Mann–Whitney U of group a.
    Convention: r > 0 ⇒ group a is stochastically larger. r ∈ [−1, 1]."""
    if na == 0 or nb == 0:
        return 0.0
    return round(2.0 * u1 / (na * nb) - 1.0, 4)


def _holm_bonferroni(pvals: List[float]) -> List[float]:
    """Holm–Bonferroni step-down adjusted p-values."""
    m = len(pvals)
    if m == 0:
        return []
    order = sorted(range(m), key=lambda i: pvals[i])
    adj = [0.0] * m
    running = 0.0
    for rank, idx in enumerate(order):
        running = max(running, (m - rank) * pvals[idx])
        adj[idx] = min(1.0, running)
    return adj


def cmd_compare(args: argparse.Namespace, term: Terminal) -> int:
    welcome(term)
    try:
        from scipy import stats as stats_mod
    except ImportError:
        raise InputError("compare needs scipy:  pip install scipy")

    path = Path(args.input)
    if not path.exists():
        raise InputError(f"Metrics file not found: '{path}'")
    df = pd.read_csv(path, sep="\t")

    if args.ips_column not in df.columns:
        raise InputError(f"Column '{args.ips_column}' not found. IPS columns "
                         f"available: {[c for c in df.columns if 'IPS' in c]}")
    group_col = args.group_by
    if group_col not in df.columns:
        raise InputError(f"Grouping column '{group_col}' not found.")

    if "is_hypothetical" in df.columns:
        df["is_hypothetical"] = _coerce_bool(df["is_hypothetical"])
    if args.hypothetical_only and "is_hypothetical" in df.columns:
        df = df[df["is_hypothetical"]].copy()
        term.info("Restricting to hypothetical GCs", f"{len(df):,} rows")

    df = df[pd.to_numeric(df[args.ips_column], errors="coerce").notna()].copy()
    df[args.ips_column] = df[args.ips_column].astype(float)
    if df.empty:
        raise InsufficientDataError("No usable IPS values to compare.")

    order = (CATEGORY_ORDER if group_col == "pangenome_category"
             else HP_ORDER if group_col == "HP_status"
             else sorted(df[group_col].dropna().astype(str).unique()))
    groups: Dict[str, np.ndarray] = {}
    summaries: List[GroupSummary] = []
    for g in order:
        sub = df[df[group_col].astype(str) == str(g)]
        vals = sub[args.ips_column].dropna().values
        if len(vals) == 0:
            continue
        groups[g] = vals
        n_hyp = int(sub["is_hypothetical"].sum()) if "is_hypothetical" in df.columns else 0
        summaries.append(GroupSummary(
            group=g, n_clusters=len(vals), n_hypothetical=n_hyp,
            pct_hypothetical=round(100 * n_hyp / len(vals), 2),
            mean_IPS=round(float(np.mean(vals)), 4),
            median_IPS=round(float(np.median(vals)), 4),
            std_IPS=round(float(np.std(vals)), 4),
            q25_IPS=round(float(np.percentile(vals, 25)), 4),
            q75_IPS=round(float(np.percentile(vals, 75)), 4),
        ))

    if len(groups) < 2:
        raise InsufficientDataError(
            f"Need at least 2 non-empty groups in '{group_col}' to compare.")

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Kruskal–Wallis
    kw_stat, kw_p = stats_mod.kruskal(*groups.values())
    term.info("Kruskal-Wallis H", round(float(kw_stat), 4), nl_before=1)
    term.info("Kruskal-Wallis p", f"{kw_p:.3e} "
              f"({'significant' if kw_p < args.alpha else 'n.s.'})")

    # pairwise Mann–Whitney U with Holm correction
    keys = list(groups.keys())
    raw = []
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            a, b = groups[keys[i]], groups[keys[j]]
            u, p = stats_mod.mannwhitneyu(a, b, alternative="two-sided")
            raw.append((keys[i], keys[j], float(u), float(p), len(a), len(b)))
    adj = _holm_bonferroni([r[3] for r in raw])
    tests: List[StatTest] = []
    for (ga, gb, u, p, na, nb), padj in zip(raw, adj):
        tests.append(StatTest(
            group_a=ga, group_b=gb, test="mannwhitneyu",
            statistic=round(u, 2), p_value=p, p_adjusted=padj,
            significant=bool(padj < args.alpha),
            effect_size=_rank_biserial(u, na, nb)))

    summ_path = out_dir / "group_summary.tsv"
    tests_path = out_dir / "statistical_tests.tsv"
    pd.DataFrame([asdict(s) for s in summaries]).to_csv(summ_path, sep="\t", index=False)
    pd.DataFrame([asdict(t) for t in tests]).to_csv(tests_path, sep="\t", index=False)
    term.info("Per-group summary", str(summ_path), nl_before=1)
    term.info("Pairwise tests (Holm-corrected)", str(tests_path))

    if not args.no_figure:
        _compare_figure(groups, kw_p, args, out_dir, group_col, term)

    term.info_single("Done.", level=0, nl_before=1)
    return 0


def _compare_figure(groups, kw_p, args, out_dir, group_col, term):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        term.warning("matplotlib not available; skipping the figure "
                     "(tables were still written).")
        return

    rng = np.random.default_rng(args.seed)
    palette = (CATEGORY_COLORS if group_col == "pangenome_category"
               else HP_COLORS if group_col == "HP_status" else {})
    cats = list(groups.keys())
    data = [groups[c] for c in cats]
    colors = [palette.get(c, "#888888") for c in cats]
    labels = [f"{c}\n(n={len(groups[c])})" for c in cats]

    fig, ax = plt.subplots(figsize=(max(6, 1.8 * len(cats)), 6))
    bp = ax.boxplot(data, patch_artist=True, showfliers=False,
                    medianprops={"color": "black", "linewidth": 2})
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.65)
    for i, arr in enumerate(data):
        x = rng.normal(i + 1, 0.06, size=len(arr))
        ax.scatter(x, arr, s=8, alpha=0.2, color=colors[i], zorder=3)
    ax.set_xticks(range(1, len(cats) + 1))
    ax.set_xticklabels(labels)
    ax.set_ylabel(args.ips_column)
    ax.set_ylim(-0.02, 1.02)
    ax.set_title(f"IPS by {group_col}  ·  Kruskal-Wallis p = {kw_p:.2e}"
                 f"{' *' if kw_p < args.alpha else ''}")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    out_png = out_dir / f"ips_by_{group_col}.png"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    term.info("Figure", str(out_png))


# ===========================================================================
# CLI
# ===========================================================================
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="pangentro",
        description=f"pangentro v{__version__} — anvi'o-native IPS "
                    f"prioritisation of hypothetical proteins.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Typical workflow (upstream is plain anvi'o 9):\n"
            "  # 1. gene calls + functions (BEFORE the genomes storage!)\n"
            "  anvi-gen-contigs-database -f assembly.fa -o sample.db\n"
            "  anvi-run-ncbi-cogs -c sample.db            # COG20 + categories\n"
            "  #   eggNOG path: emapper.py -i sample.faa --itype proteins -o s; \n"
            "  #   parse s.emapper.annotations to gene_callers_id/source/accession/\n"
            "  #   function/e_value and  anvi-import-functions -c sample.db -i fns.txt\n"
            "  anvi-gen-genomes-storage -e external.txt -o MY-GENOMES.db\n"
            "  anvi-pan-genome -g MY-GENOMES.db -n MYPROJECT\n"
            "  anvi-script-add-default-collection -p MYPROJECT/MYPROJECT-PAN.db "
            "-C DEFAULT\n"
            "  anvi-summarize -p MYPROJECT/MYPROJECT-PAN.db -g MY-GENOMES.db "
            "-C DEFAULT -o SUMMARY/\n\n"
            "  # 2. pangentro\n"
            "  pangentro integrate \\\n"
            "      --gene-clusters-summary SUMMARY/*_gene_clusters_summary.txt.gz \\\n"
            "      --pan-db MYPROJECT/MYPROJECT-PAN.db -o pangentro_out/\n"
            "  pangentro compare -i pangentro_out/MYPROJECT_pangentro_metrics.txt "
            "-o cmp/\n"
        ),
    )
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = p.add_subparsers(dest="subcommand", metavar="SUBCOMMAND")
    sub.required = True

    def common(sp):
        sp.add_argument("-v", "--verbose", action="store_true",
                        help="debug-level logging")
        sp.add_argument("--quiet", action="store_true",
                        help="suppress the banner and progress chatter "
                             "(warnings still print)")
        sp.add_argument("--seed", type=int, default=0,
                        help="RNG seed for any stochastic step (default: 0)")

    # ---- integrate --------------------------------------------------------
    ig = sub.add_parser(
        "integrate",
        help="summary → IPS + pangenomic category + HP status per gene cluster",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    common(ig)
    src = ig.add_argument_group("input (give a summary, or a pan-db+storage)")
    src.add_argument("--gene-clusters-summary",
                     help="*_gene_clusters_summary.txt(.gz) from anvi-summarize")
    src.add_argument("-p", "--pan-db",
                     help="pan-db. Used in printed commands; with "
                          "--genomes-storage lets pangentro run anvi-summarize")
    src.add_argument("-g", "--genomes-storage",
                     help="genomes storage db (for auto anvi-summarize)")
    src.add_argument("--summary-dir", default="SUMMARY",
                     help="output dir for auto anvi-summarize (default: SUMMARY)")

    out = ig.add_argument_group("output")
    out.add_argument("-o", "--output", required=True, help="output directory")
    out.add_argument("--project-name", default="pangentro",
                     help="prefix for output files (default: pangentro)")
    out.add_argument("--ips-column", default="IPS_core_mean",
                     choices=["IPS_core_mean", "IPS_self_mean", "IPS_core_max"],
                     help="IPS variant written to the anvi'o layout / used for "
                          "ranking (default: IPS_core_mean)")
    out.add_argument("--top-hp", type=int, default=50,
                     help="how many top HP clusters go into the collection "
                          "(0 = all HP; default: 50)")
    out.add_argument("--collection-name", default="pangentro_HP",
                     help="name of the exported collection (default: pangentro_HP)")
    out.add_argument("--import", dest="import_into_anvio", action="store_true",
                     help="after writing, run anvi-import-misc-data / "
                          "anvi-import-collection (requires --pan-db)")
    out.add_argument("--overwrite", action="store_true",
                     help="with --import, first delete any existing pangentro "
                          "items layers so re-import is idempotent")

    sc = ig.add_argument_group("scoring")
    sc.add_argument("-k", "--k", type=int, default=DEFAULT_K,
                    help=f"k-mer size (default: {DEFAULT_K}; k>1 is experimental)")
    sc.add_argument("--l0", type=int, default=L0_AA,
                    help=f"length-penalty scale in aa (default: {L0_AA})")
    sc.add_argument("--bg-categories", default="core",
                    help="comma-separated categories used to build the "
                         "background proteome (default: core). e.g. "
                         "'core,soft_core'")
    sc.add_argument("--null-draws", type=int, default=0,
                    help="Monte-Carlo draws per sequence length to calibrate IPS "
                         "into a length-controlled empirical p-value and z-score "
                         "against the background (0 = off; try 100–1000). "
                         "k=1 only; reproducible via --seed.")

    cat = ig.add_argument_group("pangenomic categories")
    cat.add_argument("--core", type=float, default=DEFAULT_CORE,
                     help=f"core frequency threshold (default: {DEFAULT_CORE})")
    cat.add_argument("--soft-core", type=float, default=DEFAULT_SOFT_CORE,
                     help=f"soft-core threshold (default: {DEFAULT_SOFT_CORE})")
    cat.add_argument("--shell", type=float, default=DEFAULT_SHELL,
                     help=f"shell threshold; below this is cloud "
                          f"(default: {DEFAULT_SHELL})")

    ann = ig.add_argument_group("annotation / HP definition")
    ann.add_argument("--annotation-sources", type=lambda s: s.split(","),
                     default=None,
                     help="comma-separated summary columns to treat as "
                          "functional annotation sources (default: auto-detect)")
    ann.add_argument("--category-source", default=None,
                     help="summary column holding COG functional categories "
                          "(default: auto-detect a *_CATEGORY column)")
    ann.add_argument("--weak-cog-categories", type=lambda s: list(s),
                     default=["S"],
                     help="COG category letters treated as 'unknown' "
                          "(default: S). e.g. 'SR' to also include R")
    ig.set_defaults(func=cmd_integrate)

    # ---- compare ----------------------------------------------------------
    cp = sub.add_parser(
        "compare",
        help="statistics of IPS across categories or HP status",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    common(cp)
    cp.add_argument("-i", "--input", required=True,
                    help="*_pangentro_metrics.txt from integrate")
    cp.add_argument("-o", "--output", required=True, help="output directory")
    cp.add_argument("--ips-column", default="IPS_core_mean",
                    help="IPS column to compare (default: IPS_core_mean)")
    cp.add_argument("--group-by", default="pangenome_category",
                    choices=["pangenome_category", "HP_status"],
                    help="grouping variable (default: pangenome_category)")
    cp.add_argument("--hypothetical-only", action="store_true",
                    help="restrict to hypothetical gene clusters")
    cp.add_argument("--alpha", type=float, default=0.05,
                    help="significance threshold (default: 0.05)")
    cp.add_argument("--no-figure", action="store_true", help="skip the figure")
    cp.set_defaults(func=cmd_compare)

    return p


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    level = logging.DEBUG if getattr(args, "verbose", False) else logging.WARNING
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s — %(message)s",
        datefmt="%H:%M:%S")

    # normalise comma-separated category list
    if getattr(args, "bg_categories", None) is not None:
        args.bg_categories = [c.strip() for c in args.bg_categories.split(",")
                              if c.strip()]

    term = Terminal(quiet=getattr(args, "quiet", False))
    try:
        return args.func(args, term)
    except KeyboardInterrupt:
        sys.stderr.write("\nInterrupted.\n")
        return 130
    except PangentroError as exc:
        term.warning(str(exc), header="ERROR")
        return 1


if __name__ == "__main__":
    sys.exit(main())
