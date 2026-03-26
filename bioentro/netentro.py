#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
netentro v0.2.0 — Informational Similarity Network for Hypothetical Proteins

Subcommands
───────────
  network  : Build and visualize an informational similarity network (v0.1.0)
  predict  : Predict functional class of a target protein by centroid distance
  validate : Leave-one-out cross-validation of the predict method on annotated
             proteins — quantifies method accuracy per functional class.

Scientific rationale — predict mode
────────────────────────────────────
For each known functional class C_k (e.g. Oxidoreductase, Transport), we
compute its CENTROID as the mean feature vector over all annotated members:

    centroid_k = mean( [sqrt_JSD_i, Efficiency_i, Kolmogorov_i] )
                 for all i in C_k

Distance from the target protein t to each centroid:

    d(t, C_k) = weighted_euclidean(feature_t, centroid_k)

Weights: sqrt_JSD=0.6, Efficiency=0.2, Kolmogorov=0.2
(same as build_distance_matrix — methodologically consistent)

Confidence is computed relative to the intra-class dispersion (σ_Ck),
defined as the mean distance of class members to their own centroid:

    Confidence(C_k) = max(0,  1 - d(t, C_k) / (σ_Ck + ε))

Interpretation:
  Confidence = 1.0 → target is AT the class centroid (identical composition)
  Confidence = 0.5 → target is as far as the average class member
  Confidence = 0.0 → target is farther than any typical class member

This formula can be LOW even for the top-ranked class, which correctly
signals that the target does not compositionally resemble any known class.
The previous ranking-based normalization always gave 1.0 to rank 1,
masking exactly this situation.

Predicted class = argmin_k d(t, C_k)

Why centroids instead of individual proteins
────────────────────────────────────────────
Comparing a target to individual proteins has high variance — a single
atypical protein in a class can distort the signal. Centroids aggregate
the compositional signal of the whole class, cancelling individual noise.
This is the same rationale behind k-means and prototype-based classifiers.

Validation — leave-one-out (LOO)
─────────────────────────────────
For each annotated protein p with a known class C_k:
  1. Remove p from the dataset.
  2. Recompute centroids without p.
  3. Predict the class of p.
  4. Compare prediction to true class.
Aggregate: precision, recall, F1 per class. Micro-averaged accuracy overall.
This gives a rigorous, unbiased estimate of method performance.

Usage:
    netentro network  -i protein.tsv -t NFOBNJ_05247 -o network.png
    netentro predict  -i protein.tsv -t NFOBNJ_05247 -o prediction.tsv
    netentro validate -i protein.tsv -o validation.tsv
    netentro network  -i protein.tsv --list-ids

Requirements:
    pandas, numpy, networkx, matplotlib, scipy

Author:  bioentro contributors
License: MIT
"""

from __future__ import annotations

__version__ = "0.2.0"

# ---------------------------------------------------------------------------
# Standard library
# ---------------------------------------------------------------------------
import argparse
import csv
import logging
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Third-party
# ---------------------------------------------------------------------------
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
REQUIRED_COLS: List[str] = [
    "Prot_ID", "Product", "IPS", "sqrt_JSD",
    "Efficiency", "Kolmogorov", "JSD_proteomebg",
]

# Feature columns used for distance computation (order matters)
FEATURE_COLS: List[str] = ["sqrt_JSD", "Efficiency", "Kolmogorov"]

# Weights for weighted Euclidean distance
# sqrt_JSD is the primary informational metric; the others are complementary
FEATURE_WEIGHTS: Dict[str, float] = {
    "sqrt_JSD":   0.60,
    "Efficiency": 0.20,
    "Kolmogorov": 0.20,
}

# Minimum number of annotated members for a class to be used in prediction
MIN_CLASS_SIZE: int = 5

# Keywords that flag a protein as hypothetical (case-insensitive)
HYPOTHETICAL_KEYWORDS: Tuple[str, ...] = (
    "hypothetical", "unknown function", "uncharacterized",
    "putative", "predicted protein",
)

# Broad functional category keywords — order matters: first match wins
FUNCTIONAL_CATEGORIES: List[Tuple[str, str, str]] = [
    # (keyword, display_label, hex_color)
    ("transport",       "Transport",        "#4E9AF1"),
    ("regulator",       "Regulator",        "#F4A460"),
    ("protease",        "Protease",         "#9B59B6"),
    ("oxidoreductase",  "Oxidoreductase",   "#2ECC71"),
    ("transferase",     "Transferase",      "#1ABC9C"),
    ("kinase",          "Kinase",           "#F39C12"),
    ("synthetase",      "Synthetase",       "#E74C3C"),
    ("reductase",       "Reductase",        "#27AE60"),
    ("dehydrogenase",   "Dehydrogenase",    "#16A085"),
    ("binding",         "Binding protein",  "#8E44AD"),
    ("membrane",        "Membrane protein", "#D35400"),
    ("secreted",        "Secreted",         "#C0392B"),
]

HYPOTHETICAL_COLOR: str     = "#F0A500"
TARGET_COLOR: str           = "red"
ANNOTATED_DEFAULT_COLOR: str = "#95A5A6"

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

class NetentroError(Exception):
    """Base class for netentro errors."""

class InputError(NetentroError):
    """Missing or malformed input file or column."""

class TargetNotFoundError(NetentroError):
    """Target protein ID not found in the dataset."""

class InsufficientDataError(NetentroError):
    """Not enough annotated proteins to build class profiles."""


# ---------------------------------------------------------------------------
# Result dataclasses
# ---------------------------------------------------------------------------

@dataclass
class PredictionResult:
    """One row per functional class, ranked by distance to target."""
    Target_ID: str
    Functional_Class: str
    N_members: int
    Centroid_dist: float    # weighted Euclidean distance to class centroid
    Confidence: float       # (d_max - d_k) / (d_max - d_min)  ∈ [0, 1]
    Rank: int               # 1 = closest (predicted class)
    Min_class_size_flag: bool  # True if N_members < MIN_CLASS_SIZE


@dataclass
class ValidationResult:
    """One row per annotated protein used in LOO validation."""
    Prot_ID: str
    True_Class: str
    Predicted_Class: str
    Correct: bool
    Confidence: float
    Centroid_dist: float


@dataclass
class ValidationSummary:
    """Per-class accuracy summary from LOO validation."""
    Functional_Class: str
    N_total: int
    N_correct: int
    Precision: float
    Recall: float
    F1: float


# ---------------------------------------------------------------------------
# Data loading and validation
# ---------------------------------------------------------------------------

def load_and_validate(tsv_path: Path) -> pd.DataFrame:
    """Load bioentro protein TSV and validate required columns.

    Args:
        tsv_path: Path to the bioentro output TSV.

    Returns:
        Validated DataFrame with reset index.

    Raises:
        InputError: If file is missing, unreadable, or missing required columns.
    """
    if not tsv_path.exists():
        raise InputError(f"File not found: '{tsv_path}'")

    try:
        df = pd.read_csv(tsv_path, sep="\t")
    except Exception as exc:
        raise InputError(f"Cannot read '{tsv_path}': {exc}") from exc

    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise InputError(
            f"Missing required columns: {missing}\n"
            f"Found: {list(df.columns)}\n"
            "Make sure this is a bioentro 'protein' mode output."
        )

    before = len(df)
    df = df.dropna(subset=REQUIRED_COLS).copy()
    dropped = before - len(df)
    if dropped:
        logger.warning("Dropped %d rows with NaN values.", dropped)
    if df.empty:
        raise InputError("DataFrame is empty after removing NaN rows.")

    logger.info("Loaded %d protein records from '%s'.", len(df), tsv_path.name)
    return df.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Protein classification helpers
# ---------------------------------------------------------------------------

def is_hypothetical(product: str) -> bool:
    """Return True if the product description suggests an unannotated protein."""
    lower = product.lower()
    return any(kw in lower for kw in HYPOTHETICAL_KEYWORDS)


def classify_product(product: str) -> Tuple[str, str]:
    """Map a product description to a broad functional category and color.

    Returns:
        (display_label, hex_color). Falls back to
        ("Annotated (other)", ANNOTATED_DEFAULT_COLOR) if no keyword matches.
    """
    lower = product.lower()
    for keyword, label, color in FUNCTIONAL_CATEGORIES:
        if keyword in lower:
            return label, color
    return "Annotated (other)", ANNOTATED_DEFAULT_COLOR


# ---------------------------------------------------------------------------
# Feature normalization (shared by network and predict modes)
# ---------------------------------------------------------------------------

def _compute_normalization_params(df: pd.DataFrame) -> Dict[str, Tuple[float, float]]:
    """Compute min/max for each feature column over the full dataset.

    Using dataset-wide statistics ensures that the normalization is the same
    whether we process all proteins or a single target — critical for predict
    mode to be consistent with distance_matrix mode.

    Returns:
        Dict mapping column name → (min, max).
    """
    params = {}
    for col in FEATURE_COLS:
        lo, hi = float(df[col].min()), float(df[col].max())
        params[col] = (lo, hi)
    return params


def _normalize_row(
    row: pd.Series,
    norm_params: Dict[str, Tuple[float, float]],
) -> np.ndarray:
    """Normalize a single protein's features to [0,1] and apply weights.

    Args:
        row:         DataFrame row with FEATURE_COLS columns.
        norm_params: Output of _compute_normalization_params.

    Returns:
        1-D weighted feature vector, shape (len(FEATURE_COLS),).
    """
    vec = []
    for col in FEATURE_COLS:
        lo, hi = norm_params[col]
        val = float(row[col])
        normalized = (val - lo) / (hi - lo) if hi > lo else 0.0
        vec.append(FEATURE_WEIGHTS[col] * normalized)
    return np.array(vec)


def _build_feature_matrix(
    df: pd.DataFrame,
    norm_params: Dict[str, Tuple[float, float]],
) -> np.ndarray:
    """Build weighted, normalized feature matrix for all proteins in df.

    Args:
        df:          DataFrame with FEATURE_COLS columns.
        norm_params: Output of _compute_normalization_params (from FULL dataset).

    Returns:
        Array of shape (n_proteins, n_features).
    """
    return np.vstack([_normalize_row(row, norm_params) for _, row in df.iterrows()])


# ---------------------------------------------------------------------------
# Distance computation (network mode — unchanged from v0.1.0)
# ---------------------------------------------------------------------------

def build_distance_matrix(df: pd.DataFrame) -> np.ndarray:
    """Build pairwise weighted Euclidean distance matrix from bioentro features.

    Weights: sqrt_JSD=0.6 (primary metric distance), Efficiency=0.2,
    Kolmogorov=0.2. Features are min-max normalized using dataset range
    before weighting, so weight differences reflect importance, not scale.

    Args:
        df: Validated DataFrame with bioentro protein columns.

    Returns:
        Square distance matrix (n × n), dtype float64.
    """
    norm_params = _compute_normalization_params(df)
    X = _build_feature_matrix(df, norm_params)
    dist_matrix = squareform(pdist(X, metric="euclidean"))
    logger.debug(
        "Distance matrix shape: %s, range [%.4f, %.4f]",
        dist_matrix.shape, dist_matrix.min(), dist_matrix.max(),
    )
    return dist_matrix


def calibrate_threshold(
    df: pd.DataFrame,
    dist_matrix: np.ndarray,
    k: int,
) -> float:
    """Estimate a biologically grounded connection threshold.

    For each annotated protein, find its k nearest annotated neighbors
    and record the maximum distance. Threshold = mean of these values.

    This defines "similar" as what annotated proteins of known function
    consider similar to each other — a biologically grounded cutoff.

    Args:
        df:          Full DataFrame.
        dist_matrix: Pairwise distance matrix aligned with df rows.
        k:           Number of neighbors for calibration.

    Returns:
        Calibrated threshold, or median of pairwise distances as fallback.
    """
    annotated_idx = [
        i for i, prod in enumerate(df["Product"])
        if not is_hypothetical(str(prod))
    ]

    if len(annotated_idx) < k + 1:
        fallback = float(np.median(dist_matrix[dist_matrix > 0]))
        logger.warning(
            "Too few annotated proteins (%d) for calibration (need >%d). "
            "Using median pairwise distance %.4f as threshold.",
            len(annotated_idx), k, fallback,
        )
        return fallback

    max_knn_dists = []
    ann_sub = dist_matrix[np.ix_(annotated_idx, annotated_idx)]
    for i in range(len(annotated_idx)):
        row = ann_sub[i].copy()
        row[i] = np.inf
        knn_dists = np.sort(row)[:k]
        max_knn_dists.append(knn_dists[-1])

    threshold = float(np.mean(max_knn_dists))
    logger.info(
        "Calibrated threshold: %.4f (mean max-%d-NN among %d annotated proteins).",
        threshold, k, len(annotated_idx),
    )
    return threshold


# ---------------------------------------------------------------------------
# Predict mode — class profile construction
# ---------------------------------------------------------------------------

def build_class_profiles(
    df: pd.DataFrame,
    norm_params: Dict[str, Tuple[float, float]],
    exclude_id: Optional[str] = None,
) -> Dict[str, Tuple[np.ndarray, float]]:
    """Compute centroid and intra-class dispersion for each functional class.

    Only annotated, non-hypothetical proteins are used. Classes with fewer
    than MIN_CLASS_SIZE members are retained but flagged in predict output.

    Intra-class dispersion (σ_Ck) is the mean distance of each class member
    to the class centroid. It defines the "natural radius" of the class in
    informational space — used by predict_function to compute a confidence
    score that can be genuinely low when the target does not resemble any class.

    Args:
        df:          Full DataFrame.
        norm_params: Normalization parameters from the FULL dataset.
                     Must NOT be recomputed per-subset (LOO or otherwise)
                     to keep all distances in the same space.
        exclude_id:  Prot_ID to exclude before building profiles (for LOO).

    Returns:
        Dict mapping class label → (centroid_vector, intra_class_std).
        centroid_vector : 1-D numpy array, shape (n_features,).
        intra_class_std : mean distance of members to centroid (float ≥ 0).
    """
    df_ann = df[
        df["Product"].apply(lambda p: not is_hypothetical(str(p)))
    ].copy()

    if exclude_id is not None:
        df_ann = df_ann[df_ann["Prot_ID"] != exclude_id]

    raw: Dict[str, List[np.ndarray]] = {}
    for _, row in df_ann.iterrows():
        label, _ = classify_product(str(row["Product"]))
        if label == "Annotated (other)":
            continue
        vec = _normalize_row(row, norm_params)
        raw.setdefault(label, []).append(vec)

    profiles: Dict[str, Tuple[np.ndarray, float]] = {}
    for label, vecs in raw.items():
        matrix = np.vstack(vecs)                          # (n_members, n_features)
        centroid = matrix.mean(axis=0)
        # Mean distance of each member to the centroid
        member_dists = np.linalg.norm(matrix - centroid, axis=1)
        sigma = float(member_dists.mean()) if len(vecs) > 1 else 0.0
        profiles[label] = (centroid, sigma)

    return profiles


def _class_sizes(
    df: pd.DataFrame,
    exclude_id: Optional[str] = None,
) -> Dict[str, int]:
    """Count annotated members per functional class."""
    sizes: Dict[str, int] = {}
    for _, row in df.iterrows():
        if is_hypothetical(str(row["Product"])):
            continue
        if exclude_id and row["Prot_ID"] == exclude_id:
            continue
        label, _ = classify_product(str(row["Product"]))
        if label == "Annotated (other)":
            continue
        sizes[label] = sizes.get(label, 0) + 1
    return sizes


def predict_function(
    target_id: str,
    df: pd.DataFrame,
    norm_params: Dict[str, Tuple[float, float]],
    profiles: Dict[str, Tuple[np.ndarray, float]],
    sizes: Dict[str, int],
) -> List[PredictionResult]:
    """Rank functional classes by centroid distance from a target protein.

    Distance metric: weighted Euclidean (same weights as distance matrix).

    Confidence formula:

        Confidence(C_k) = max(0,  1 - d(t, C_k) / (σ_Ck + ε))

    where σ_Ck is the mean distance of class members to their centroid
    (intra-class dispersion). A confidence of 0.5 means the target is
    exactly as far from the centroid as the average class member.
    A confidence below 0 (clipped to 0) means the target is farther than
    any typical class member — the prediction is unreliable.

    This formula produces genuinely low confidence when the target does
    not resemble any class — unlike ranking-based normalization which
    always gives 1.0 to rank 1 regardless of absolute distance.

    Args:
        target_id:   Prot_ID of the target protein.
        df:          Full DataFrame.
        norm_params: Normalization parameters from the full dataset.
        profiles:    Output of build_class_profiles() — (centroid, sigma) per class.
        sizes:       Class member counts from _class_sizes().

    Returns:
        List of PredictionResult sorted by distance (closest first).

    Raises:
        TargetNotFoundError:    If target_id is not in df.
        InsufficientDataError:  If no class profiles could be built.
    """
    id_list = df["Prot_ID"].tolist()
    if target_id not in id_list:
        raise TargetNotFoundError(
            f"Target ID '{target_id}' not found.\n"
            f"First 10 IDs: {id_list[:10]}"
        )
    if not profiles:
        raise InsufficientDataError(
            "No functional class profiles could be built. "
            "Check that the dataset contains annotated proteins with "
            "recognizable functional keywords."
        )

    _EPS = 1e-9   # prevents division by zero for singleton classes

    target_row = df[df["Prot_ID"] == target_id].iloc[0]
    target_vec = _normalize_row(target_row, norm_params)

    # Distance from target to each class centroid
    dists: Dict[str, float] = {
        label: float(np.linalg.norm(target_vec - centroid))
        for label, (centroid, _sigma) in profiles.items()
    }

    results = []
    for rank, (label, dist) in enumerate(
        sorted(dists.items(), key=lambda x: x[1]), start=1
    ):
        _centroid, sigma = profiles[label]
        # Confidence relative to intra-class dispersion
        confidence = max(0.0, 1.0 - dist / (sigma + _EPS))
        n = sizes.get(label, 0)
        results.append(
            PredictionResult(
                Target_ID=target_id,
                Functional_Class=label,
                N_members=n,
                Centroid_dist=round(dist, 6),
                Confidence=round(confidence, 4),
                Rank=rank,
                Min_class_size_flag=n < MIN_CLASS_SIZE,
            )
        )

    top = results[0]
    logger.info(
        "Prediction for '%s': %s (confidence=%.4f, centroid_dist=%.6f, σ=%.6f, n=%d)",
        target_id, top.Functional_Class, top.Confidence, top.Centroid_dist,
        profiles[top.Functional_Class][1], top.N_members,
    )
    return results


# ---------------------------------------------------------------------------
# Validate mode — leave-one-out cross-validation
# ---------------------------------------------------------------------------

def run_loo_validation(
    df: pd.DataFrame,
    norm_params: Dict[str, Tuple[float, float]],
) -> Tuple[List[ValidationResult], List[ValidationSummary]]:
    """Leave-one-out cross-validation of the predict method.

    For each annotated protein with a recognized functional class:
      1. Remove it from the dataset.
      2. Rebuild class profiles without it (norm_params unchanged).
      3. Predict its class.
      4. Compare prediction vs true class.

    This gives an unbiased estimate of method accuracy because the
    protein being predicted is never part of its own reference set.

    Args:
        df:          Full validated DataFrame.
        norm_params: Normalization parameters (from full dataset, never recomputed).

    Returns:
        Tuple (per-protein results, per-class summary with precision/recall/F1).
    """
    # Build the pool: annotated proteins with a recognized class
    pool = []
    for _, row in df.iterrows():
        if is_hypothetical(str(row["Product"])):
            continue
        label, _ = classify_product(str(row["Product"]))
        if label == "Annotated (other)":
            continue
        pool.append((row["Prot_ID"], label))

    if len(pool) < MIN_CLASS_SIZE * 2:
        raise InsufficientDataError(
            f"Only {len(pool)} annotated proteins with recognized classes found. "
            f"Need at least {MIN_CLASS_SIZE * 2} for meaningful LOO validation."
        )

    logger.info("Running LOO validation on %d annotated proteins...", len(pool))

    per_protein: List[ValidationResult] = []

    for prot_id, true_class in pool:
        # Rebuild profiles excluding this protein
        profiles_loo = build_class_profiles(df, norm_params, exclude_id=prot_id)
        sizes_loo = _class_sizes(df, exclude_id=prot_id)

        # Skip if true class disappears entirely after exclusion (singleton class)
        if true_class not in profiles_loo:
            logger.debug(
                "Skipping '%s': class '%s' has only 1 member.", prot_id, true_class
            )
            continue

        try:
            preds = predict_function(prot_id, df, norm_params, profiles_loo, sizes_loo)
        except NetentroError as exc:
            logger.warning("LOO prediction failed for '%s': %s", prot_id, exc)
            continue

        top = preds[0]
        per_protein.append(
            ValidationResult(
                Prot_ID=prot_id,
                True_Class=true_class,
                Predicted_Class=top.Functional_Class,
                Correct=top.Functional_Class == true_class,
                Confidence=top.Confidence,
                Centroid_dist=top.Centroid_dist,
            )
        )

    # ── Per-class summary ────────────────────────────────────────────────────
    # True positives, false positives, false negatives per class
    all_classes = sorted({r.True_Class for r in per_protein})
    summaries: List[ValidationSummary] = []

    for cls in all_classes:
        tp = sum(1 for r in per_protein if r.True_Class == cls and r.Correct)
        fp = sum(1 for r in per_protein if r.Predicted_Class == cls and not r.Correct)
        fn = sum(1 for r in per_protein if r.True_Class == cls and not r.Correct)
        n_total = tp + fn

        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall    = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1        = (
            2 * precision * recall / (precision + recall)
            if (precision + recall) > 0 else 0.0
        )

        summaries.append(ValidationSummary(
            Functional_Class=cls,
            N_total=n_total,
            N_correct=tp,
            Precision=round(precision, 4),
            Recall=round(recall, 4),
            F1=round(f1, 4),
        ))

    # Overall accuracy
    overall = sum(r.Correct for r in per_protein) / len(per_protein) if per_protein else 0.0
    logger.info(
        "LOO validation complete: %d proteins, overall accuracy=%.4f",
        len(per_protein), overall,
    )

    return per_protein, summaries


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def write_tsv(results: list, output_path: Path) -> None:
    """Write a list of dataclass instances to a tab-separated file.

    Args:
        results:     Non-empty list of homogeneous dataclass instances.
        output_path: Destination path (parent directories created as needed).

    Raises:
        ValueError: If results is empty.
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

    logger.info("Written %d record(s) to '%s'.", len(results), output_path)


# ---------------------------------------------------------------------------
# Node / edge selection (network mode — unchanged)
# ---------------------------------------------------------------------------

def select_subset(
    df: pd.DataFrame,
    dist_matrix: np.ndarray,
    target_id: str,
    k_neighbors: int,
    top_ips: int,
) -> pd.DataFrame:
    """Select proteins to include in the network visualization.

    Includes: target + k nearest neighbors (full dataset) + top-N by IPS.
    Full-dataset distances are used so neighbor relationships are not
    distorted by re-scaling a subset.

    Raises:
        TargetNotFoundError: If target_id is not in df.
    """
    id_list = df["Prot_ID"].tolist()
    if target_id not in id_list:
        raise TargetNotFoundError(
            f"Target ID '{target_id}' not found.\n"
            f"First 10 IDs: {id_list[:10]}"
        )

    target_pos = id_list.index(target_id)
    dists_to_target = dist_matrix[target_pos]

    sorted_idx = np.argsort(dists_to_target)
    neighbors = [i for i in sorted_idx if i != target_pos][:k_neighbors]
    top_idx = df["IPS"].nlargest(top_ips).index.tolist()
    selected = list({target_pos} | set(neighbors) | set(top_idx))

    df_sub = df.iloc[selected].copy().reset_index(drop=True)
    logger.info(
        "Subset: %d proteins (%d neighbors + %d top-IPS + target).",
        len(df_sub), len(neighbors), top_ips,
    )
    return df_sub


# ---------------------------------------------------------------------------
# Graph construction (network mode — unchanged)
# ---------------------------------------------------------------------------

def build_graph(
    df_sub: pd.DataFrame,
    dist_matrix_full: np.ndarray,
    df_full: pd.DataFrame,
    threshold: float,
    target_id: str,
) -> nx.Graph:
    """Build an annotated NetworkX graph from the protein subset.

    Edges use full-dataset distances (not subset-recalculated) to avoid
    distortion from re-scaling.

    Node attributes: product, is_hypo, func_label, color, ips, efficiency, kolmogorov
    Edge attributes: distance, weight (= 1 - distance)
    """
    G = nx.Graph()
    full_ids = df_full["Prot_ID"].tolist()

    for _, row in df_sub.iterrows():
        hypo = is_hypothetical(str(row["Product"]))
        func_label, color = classify_product(str(row["Product"]))

        if row["Prot_ID"] == target_id:
            color = TARGET_COLOR
            func_label = "Target (hypothetical)"
        elif hypo:
            color = HYPOTHETICAL_COLOR
            func_label = "Hypothetical"

        G.add_node(
            row["Prot_ID"],
            product=row["Product"],
            is_hypo=hypo,
            func_label=func_label,
            color=color,
            ips=float(row["IPS"]),
            efficiency=float(row["Efficiency"]),
            kolmogorov=float(row["Kolmogorov"]),
        )

    sub_ids = df_sub["Prot_ID"].tolist()
    edge_count = 0
    for i, id_i in enumerate(sub_ids):
        for j, id_j in enumerate(sub_ids):
            if j <= i:
                continue
            pos_i = full_ids.index(id_i)
            pos_j = full_ids.index(id_j)
            dist = dist_matrix_full[pos_i, pos_j]
            if dist < threshold:
                G.add_edge(id_i, id_j, distance=dist, weight=1.0 - dist)
                edge_count += 1

    logger.info(
        "Graph: %d nodes, %d edges (threshold=%.4f).",
        G.number_of_nodes(), edge_count, threshold,
    )
    if edge_count == 0:
        logger.warning(
            "No edges created. Consider increasing --threshold-percentile "
            "or --calibration-k."
        )
    return G


# ---------------------------------------------------------------------------
# Visualization (network mode — unchanged)
# ---------------------------------------------------------------------------

def _make_label(prot_id: str, product: str, max_chars: int = 12) -> str:
    """Generate a readable two-line node label: ID suffix + first keyword."""
    short_id = prot_id.split("_")[-1]
    skip = {"protein", "the", "of", "a", "an", "and", "or", "with"}
    words = [w for w in product.lower().split() if w not in skip and len(w) > 2]
    keyword = words[0][:max_chars] if words else "?"
    return f"{short_id}\n({keyword})"


def draw_network(
    G: nx.Graph,
    target_id: str,
    config: "NetworkConfig",
    output_path: Path,
    threshold: float,
) -> None:
    """Render and save the informational similarity network.

    Node color  : functional category
    Node size   : IPS score
    Edge width  : inversely proportional to distance (thicker = closer)
    """
    if G.number_of_nodes() == 0:
        logger.error("Graph is empty — nothing to draw.")
        return

    pos = nx.spring_layout(
        G,
        k=config.layout_k,
        iterations=config.layout_iterations,
        seed=config.layout_seed,
    )

    fig, ax = plt.subplots(figsize=config.figsize)

    # ── Edges ────────────────────────────────────────────────────────────────
    edges = list(G.edges(data=True))
    if edges:
        max_dist = max(d["distance"] for _, _, d in edges) + 1e-9
        for u, v, data in edges:
            proximity = 1.0 - data["distance"] / max_dist
            alpha = config.edge_alpha_base + config.edge_alpha_scale * proximity
            width = 0.5 + 2.5 * proximity
            nx.draw_networkx_edges(
                G, pos, edgelist=[(u, v)], ax=ax,
                alpha=alpha, edge_color="gray", width=width,
            )

    # ── Non-target nodes ─────────────────────────────────────────────────────
    non_targets = [n for n in G.nodes if n != target_id]
    if non_targets:
        colors = [G.nodes[n]["color"] for n in non_targets]
        sizes = [
            config.node_size_base + config.node_size_ips_scale * G.nodes[n]["ips"]
            for n in non_targets
        ]
        nx.draw_networkx_nodes(
            G, pos, nodelist=non_targets, ax=ax,
            node_color=colors, node_size=sizes,
            alpha=0.85, edgecolors="black", linewidths=0.8,
        )

    # ── Target node ──────────────────────────────────────────────────────────
    if target_id in G.nodes:
        target_ips = G.nodes[target_id]["ips"]
        nx.draw_networkx_nodes(
            G, pos, nodelist=[target_id], ax=ax,
            node_color=TARGET_COLOR,
            node_size=config.node_size_base + config.node_size_ips_scale * target_ips * 1.5,
            node_shape="*", edgecolors="black", linewidths=1.2,
        )

    # ── Labels ───────────────────────────────────────────────────────────────
    labels = {n: _make_label(n, G.nodes[n]["product"]) for n in G.nodes}
    nx.draw_networkx_labels(G, pos, labels, ax=ax, font_size=7.5, font_weight="bold")

    # ── Legend ───────────────────────────────────────────────────────────────
    legend_items: List[mpatches.Patch] = [
        mpatches.Patch(color=TARGET_COLOR, label=f"Target: {target_id}"),
        mpatches.Patch(color=HYPOTHETICAL_COLOR, label="Hypothetical protein"),
        mpatches.Patch(color=ANNOTATED_DEFAULT_COLOR, label="Annotated (other)"),
    ]
    seen_labels: set = set()
    for n in G.nodes:
        lbl = G.nodes[n]["func_label"]
        col = G.nodes[n]["color"]
        if lbl not in seen_labels and col not in (
            TARGET_COLOR, HYPOTHETICAL_COLOR, ANNOTATED_DEFAULT_COLOR
        ):
            legend_items.append(mpatches.Patch(color=col, label=lbl))
            seen_labels.add(lbl)

    ax.legend(
        handles=legend_items + [
            plt.scatter(
                [], [],
                s=config.node_size_base + config.node_size_ips_scale * v,
                c="gray", alpha=0.5, edgecolors="black", linewidths=0.8,
                label=f"IPS={v:.2f}",
            )
            for v in [0.25, 0.50, 0.75]
        ],
        loc="upper left", fontsize=8, framealpha=0.85,
        title="Node color = functional class\nNode size ∝ IPS",
    )

    # ── Title ────────────────────────────────────────────────────────────────
    n_hypo = sum(1 for n in G.nodes if G.nodes[n]["is_hypo"] and n != target_id)
    n_ann  = G.number_of_nodes() - n_hypo - (1 if target_id in G.nodes else 0)
    ax.set_title(
        f"Informational Similarity Network — Focus: {target_id}\n"
        f"netentro v{__version__} | "
        f"Edges: √JSD distance < {threshold:.4f} | "
        f"Nodes: {G.number_of_nodes()} ({n_ann} annotated, {n_hypo} hypothetical)",
        fontsize=13, pad=14,
    )
    ax.axis("off")
    fig.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
    logger.info("Network saved to '%s'.", output_path)
    plt.show()


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class NetworkConfig:
    """Tunable parameters for network mode — nothing hardcoded in logic."""
    target_id: str
    k_neighbors: int = 25
    top_ips: int = 10
    threshold_percentile: float = 20.0
    use_calibrated_threshold: bool = True
    calibration_k: int = 5
    layout_k: float = 2.5
    layout_iterations: int = 150
    layout_seed: int = 42
    figsize: Tuple[int, int] = (15, 12)
    node_size_base: float = 400.0
    node_size_ips_scale: float = 1500.0
    edge_alpha_base: float = 0.15
    edge_alpha_scale: float = 0.6
    dpi: int = 300


# ---------------------------------------------------------------------------
# CLI — subcommand architecture
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="netentro",
        description=(
            f"netentro v{__version__} — Informational Similarity Network\n\n"
            "Subcommands:\n"
            "  network   Build and visualize an informational similarity network\n"
            "  predict   Predict functional class of a target protein\n"
            "  validate  LOO cross-validation of predict on annotated proteins"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  netentro network  -i protein.tsv -t NFOBNJ_05247 -o network.png\n"
            "  netentro network  -i protein.tsv -t NFOBNJ_05247 -k 30 --top-ips 15 -o net.png\n"
            "  netentro network  -i protein.tsv --list-ids\n"
            "  netentro predict  -i protein.tsv -t NFOBNJ_05247 -o prediction.tsv\n"
            "  netentro validate -i protein.tsv -o validation.tsv\n\n"
            "Run 'netentro <subcommand> --help' for subcommand-specific options."
        ),
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )

    sub = parser.add_subparsers(dest="subcommand", metavar="SUBCOMMAND")
    sub.required = True

    # ── Shared arguments factory ──────────────────────────────────────────────
    def _add_shared(p: argparse.ArgumentParser) -> None:
        p.add_argument("-i", "--input", required=True,
                       help="bioentro 'protein' mode TSV output")
        p.add_argument("-v", "--verbose", action="store_true",
                       help="Enable debug-level logging")

    # ── network subcommand ────────────────────────────────────────────────────
    net = sub.add_parser(
        "network",
        help="Build and visualize an informational similarity network",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  netentro network -i protein.tsv -t NFOBNJ_05247 -o network.png\n"
            "  netentro network -i protein.tsv -t NFOBNJ_05247 -k 30 --top-ips 15\n"
            "  netentro network -i protein.tsv --list-ids\n"
        ),
    )
    _add_shared(net)
    net.add_argument("-t", "--target",
                     help="Prot_ID of the target hypothetical protein")
    net.add_argument("-o", "--output", default="network.png",
                     help="Output image path (default: network.png)")
    net.add_argument("-k", "--k-neighbors", type=int, default=25,
                     help="Nearest informational neighbors (default: 25)")
    net.add_argument("--top-ips", type=int, default=10,
                     help="Additional top-IPS proteins to include (default: 10)")
    net.add_argument("--threshold-percentile", type=float, default=20.0,
                     help="Fallback connection threshold percentile (default: 20.0)")
    net.add_argument("--calibration-k", type=int, default=5,
                     help="k for threshold calibration (default: 5)")
    net.add_argument("--no-calibration", action="store_true",
                     help="Skip calibration, use --threshold-percentile directly")
    net.add_argument("--layout-k", type=float, default=2.5,
                     help="Spring layout repulsion constant (default: 2.5)")
    net.add_argument("--layout-seed", type=int, default=42,
                     help="Random seed for reproducible layout (default: 42)")
    net.add_argument("--dpi", type=int, default=300,
                     help="Output image DPI (default: 300)")
    net.add_argument("--list-ids", action="store_true",
                     help="Print first 50 Prot_IDs and exit")

    # ── predict subcommand ────────────────────────────────────────────────────
    pred = sub.add_parser(
        "predict",
        help="Predict functional class of a target protein by centroid distance",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Output columns:\n"
            "  Functional_Class    — class name\n"
            "  N_members           — annotated proteins in that class\n"
            "  Centroid_dist       — weighted Euclidean distance to class centroid\n"
            "  Confidence          — max(0, 1 - dist/σ_class)  ∈ [0, 1]\n"
            "                        0.5 = target as far as avg class member\n"
            "                        0.0 = target outside typical class range\n"
            "  Rank                — 1 = predicted class (closest centroid)\n"
            "  Min_class_size_flag — True if N_members < 5 (interpret with caution)\n\n"
            "Examples:\n"
            "  netentro predict -i protein.tsv -t NFOBNJ_05247 -o prediction.tsv\n"
            "  netentro predict -i protein.tsv -t NFOBNJ_05247  # stdout only\n"
        ),
    )
    _add_shared(pred)
    pred.add_argument("-t", "--target", required=True,
                      help="Prot_ID of the target hypothetical protein")
    pred.add_argument("-o", "--output", default="prediction.tsv",
                      help="Output TSV path (default: prediction.tsv)")

    # ── validate subcommand ───────────────────────────────────────────────────
    val = sub.add_parser(
        "validate",
        help="LOO cross-validation of predict accuracy on annotated proteins",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Outputs two files:\n"
            "  <output>            — per-protein LOO results\n"
            "  <output>.summary    — per-class precision / recall / F1\n\n"
            "Example:\n"
            "  netentro validate -i protein.tsv -o validation.tsv\n"
        ),
    )
    _add_shared(val)
    val.add_argument("-o", "--output", default="validation.tsv",
                     help="Output TSV path (default: validation.tsv)")

    return parser


# ---------------------------------------------------------------------------
# Subcommand handlers
# ---------------------------------------------------------------------------

def _run_network(args: argparse.Namespace, df: pd.DataFrame) -> int:
    """Execute the 'network' subcommand."""
    if args.list_ids:
        print("\n".join(df["Prot_ID"].head(50).tolist()))
        return 0

    if not args.target:
        print("error: --target is required unless --list-ids is used.", file=sys.stderr)
        return 1

    config = NetworkConfig(
        target_id=args.target,
        k_neighbors=args.k_neighbors,
        top_ips=args.top_ips,
        threshold_percentile=args.threshold_percentile,
        use_calibrated_threshold=not args.no_calibration,
        calibration_k=args.calibration_k,
        layout_k=args.layout_k,
        layout_seed=args.layout_seed,
        dpi=args.dpi,
    )

    dist_matrix = build_distance_matrix(df)

    if config.use_calibrated_threshold:
        threshold = calibrate_threshold(df, dist_matrix, config.calibration_k)
    else:
        threshold = float(np.percentile(
            dist_matrix[dist_matrix > 0], config.threshold_percentile
        ))
        logger.info("Using percentile-%.1f threshold: %.4f",
                    config.threshold_percentile, threshold)

    df_sub = select_subset(
        df, dist_matrix, config.target_id, config.k_neighbors, config.top_ips
    )
    G = build_graph(df_sub, dist_matrix, df, threshold, config.target_id)
    draw_network(G, config.target_id, config, Path(args.output), threshold)
    return 0


def _run_predict(args: argparse.Namespace, df: pd.DataFrame) -> int:
    """Execute the 'predict' subcommand."""
    norm_params = _compute_normalization_params(df)
    profiles = build_class_profiles(df, norm_params)
    sizes = _class_sizes(df)

    results = predict_function(args.target, df, norm_params, profiles, sizes)
    write_tsv(results, Path(args.output))

    # Print a human-readable summary to stdout
    print(f"\nPrediction for {args.target}")
    print(f"{'Rank':<6}{'Class':<22}{'N':<8}{'Dist':<12}{'Confidence'}")
    print("─" * 54)
    for r in results:
        flag = " ⚠ low N" if r.Min_class_size_flag else ""
        print(
            f"{r.Rank:<6}{r.Functional_Class:<22}{r.N_members:<8}"
            f"{r.Centroid_dist:<12.6f}{r.Confidence:.4f}{flag}"
        )
    print(f"\n→ Predicted class: {results[0].Functional_Class} "
          f"(confidence={results[0].Confidence:.4f})")
    print(f"  Full results written to '{args.output}'")
    return 0


def _run_validate(args: argparse.Namespace, df: pd.DataFrame) -> int:
    """Execute the 'validate' subcommand."""
    norm_params = _compute_normalization_params(df)
    per_protein, summaries = run_loo_validation(df, norm_params)

    output_path = Path(args.output)
    summary_path = output_path.with_suffix(".summary.tsv")

    write_tsv(per_protein, output_path)
    write_tsv(summaries, summary_path)

    # Print summary table to stdout
    overall = sum(r.Correct for r in per_protein) / len(per_protein) if per_protein else 0.0
    print(f"\nLOO Validation Summary ({len(per_protein)} proteins)")
    print(f"Overall accuracy: {overall:.4f}\n")
    print(f"{'Class':<22}{'N':<8}{'Correct':<10}{'Precision':<12}{'Recall':<10}{'F1'}")
    print("─" * 66)
    for s in sorted(summaries, key=lambda x: x.F1, reverse=True):
        print(
            f"{s.Functional_Class:<22}{s.N_total:<8}{s.N_correct:<10}"
            f"{s.Precision:<12.4f}{s.Recall:<10.4f}{s.F1:.4f}"
        )
    print(f"\n  Per-protein results → '{output_path}'")
    print(f"  Per-class summary   → '{summary_path}'")
    return 0


# ---------------------------------------------------------------------------
# Main entry point
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
        df = load_and_validate(Path(args.input))

        if args.subcommand == "network":
            return _run_network(args, df)
        elif args.subcommand == "predict":
            return _run_predict(args, df)
        else:  # validate
            return _run_validate(args, df)

    except NetentroError as exc:
        logger.error("%s", exc)
        return 1
    except KeyboardInterrupt:
        logger.warning("Interrupted by user.")
        return 130


if __name__ == "__main__":
    sys.exit(main())