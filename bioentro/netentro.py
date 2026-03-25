#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
netentro v0.1.0 — Informational Similarity Network for Hypothetical Proteins

Builds a protein similarity network where edges represent informational
closeness (√JSD distance) between proteins. Designed to complement
bioentro output and support functional prediction of hypothetical proteins.

Scientific rationale
────────────────────
Distance metric: √JSD(P_i || P_j) computed directly from the k-mer
distributions output by bioentro (sqrt_JSD column). This is a proper
metric distance (symmetric, bounded [0,1], satisfies triangle inequality),
mathematically justified for use in network construction and ranking.

Node selection strategy
───────────────────────
For a target hypothetical protein, we select:
  1. Its k nearest informational neighbors from the FULL dataset
     (using sqrt_JSD from bioentro output — no re-scaling).
  2. Optionally, top-N proteins by IPS for global context.

Threshold calibration
─────────────────────
The connection threshold is derived from the within-class mean distance
of annotated proteins to their k nearest annotated neighbors (i.e., what
"similar" means for proteins of known function). Connections below this
threshold are drawn. This gives a biologically grounded cutoff rather than
an arbitrary percentile of a re-scaled subset.

Node visual encoding
────────────────────
  Color : functional annotation class
    - Red star    : target hypothetical protein
    - Orange      : other hypothetical proteins
    - Blue shades : annotated proteins (color by broad functional category)
  Size  : IPS score (larger = higher priority)
  Edge  : √JSD distance (thicker/darker = more similar)

Usage:
    python netentro.py -i protein.tsv -t NFOBNJ_05247 -o network.png
    python netentro.py -i protein.tsv -t NFOBNJ_05247 -k 30 --top-ips 15 -o net.png
    python netentro.py -i protein.tsv -t NFOBNJ_05247 --list-ids  # inspect IDs

Requirements:
    pandas, numpy, networkx, matplotlib, scipy, scikit-learn

Author:  bioentro contributors
License: MIT
"""

from __future__ import annotations

__version__ = "0.1.0"

# ---------------------------------------------------------------------------
# Standard library
# ---------------------------------------------------------------------------
import argparse
import logging
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple

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
REQUIRED_COLS: List[str] = ["Prot_ID", "Product", "IPS", "sqrt_JSD",
                             "Efficiency", "Kolmogorov", "JSD_proteomebg"]

# Keywords that flag a protein as hypothetical (case-insensitive)
HYPOTHETICAL_KEYWORDS: Tuple[str, ...] = (
    "hypothetical", "unknown function", "uncharacterized",
    "putative", "predicted protein",
)

# Broad functional category keywords for color encoding
# Order matters: first match wins
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

HYPOTHETICAL_COLOR: str = "#F0A500"   # orange
TARGET_COLOR: str = "red"
ANNOTATED_DEFAULT_COLOR: str = "#95A5A6"   # gray for uncategorized annotated

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


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class NetworkConfig:
    """All tunable parameters in one place — nothing hardcoded in logic."""
    target_id: str
    k_neighbors: int = 25          # k nearest neighbors of target
    top_ips: int = 10              # additional top-IPS proteins for global context
    threshold_percentile: float = 20.0  # fallback if calibration fails
    use_calibrated_threshold: bool = True
    calibration_k: int = 5        # k used when estimating within-class distance
    layout_k: float = 2.5         # spring layout repulsion constant
    layout_iterations: int = 150
    layout_seed: int = 42
    figsize: Tuple[int, int] = (15, 12)
    node_size_base: float = 400.0
    node_size_ips_scale: float = 1500.0
    edge_alpha_base: float = 0.15
    edge_alpha_scale: float = 0.6
    dpi: int = 300


# ---------------------------------------------------------------------------
# Data loading and validation
# ---------------------------------------------------------------------------

def load_and_validate(tsv_path: Path) -> pd.DataFrame:
    """Load bioentro protein TSV and validate required columns.

    Args:
        tsv_path: Path to the bioentro output TSV.

    Returns:
        Validated DataFrame.

    Raises:
        InputError: If file is missing or required columns are absent.
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

    Args:
        product: Full product description string.

    Returns:
        (display_label, hex_color) tuple.
        Returns ("Annotated (other)", ANNOTATED_DEFAULT_COLOR) if no keyword matches.
    """
    lower = product.lower()
    for keyword, label, color in FUNCTIONAL_CATEGORIES:
        if keyword in lower:
            return label, color
    return "Annotated (other)", ANNOTATED_DEFAULT_COLOR


# ---------------------------------------------------------------------------
# Distance computation
# ---------------------------------------------------------------------------

def build_distance_matrix(df: pd.DataFrame) -> np.ndarray:
    """Build pairwise distance matrix using sqrt_JSD from bioentro output.

    sqrt_JSD is a proper metric distance (symmetric, bounded [0,1],
    triangle inequality). We use it directly — no re-scaling needed,
    as it already lives in a mathematically defined space.

    For multi-feature distance (Efficiency, Kolmogorov, sqrt_JSD):
    We use a weighted Euclidean distance in the normalized original space.
    Weights reflect theoretical importance:
      sqrt_JSD   : 0.6 — primary informational distance (metric distance)
      Efficiency  : 0.2 — complexity contribution
      Kolmogorov  : 0.2 — long-range structure contribution

    This is documented and justified — not an arbitrary Euclidean in scaled space.

    Args:
        df: Validated DataFrame with bioentro columns.

    Returns:
        Square distance matrix (n × n), dtype float64.
    """
    # Use sqrt_JSD as the primary axis — it's already a proper distance
    # Normalize Efficiency and Kolmogorov to [0,1] using dataset range
    # so they contribute proportionally, not because of scale differences
    def minmax(col: pd.Series) -> pd.Series:
        lo, hi = col.min(), col.max()
        return (col - lo) / (hi - lo) if hi > lo else col * 0.0

    w_jsd = 0.60
    w_eff = 0.20
    w_kol = 0.20

    jsd_n = minmax(df["sqrt_JSD"]).values
    eff_n = minmax(df["Efficiency"]).values
    kol_n = minmax(df["Kolmogorov"]).values

    # Weighted feature matrix
    X = np.column_stack([
        w_jsd * jsd_n,
        w_eff * eff_n,
        w_kol * kol_n,
    ])

    dist_matrix = squareform(pdist(X, metric="euclidean"))
    logger.debug("Distance matrix shape: %s, range [%.4f, %.4f]",
                 dist_matrix.shape, dist_matrix.min(), dist_matrix.max())
    return dist_matrix


def calibrate_threshold(
    df: pd.DataFrame,
    dist_matrix: np.ndarray,
    k: int,
) -> float:
    """Estimate a biologically grounded connection threshold.

    Strategy: for each annotated (non-hypothetical) protein, find its
    k nearest annotated neighbors and record the maximum distance.
    The threshold is the mean of these values across all annotated proteins.

    Rationale: this defines "similar" as what annotated proteins of known
    function consider similar to each other. Connections below this threshold
    are therefore as tight as connections within known functional groups.

    Args:
        df:          Full DataFrame.
        dist_matrix: Pairwise distance matrix aligned with df rows.
        k:           Number of neighbors for calibration.

    Returns:
        Calibrated distance threshold.
        Falls back to median of all pairwise distances if too few annotated.
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
        row[i] = np.inf   # exclude self
        knn_dists = np.sort(row)[:k]
        max_knn_dists.append(knn_dists[-1])

    threshold = float(np.mean(max_knn_dists))
    logger.info(
        "Calibrated threshold: %.4f (mean max-%d-NN distance among %d annotated proteins).",
        threshold, k, len(annotated_idx),
    )
    return threshold


# ---------------------------------------------------------------------------
# Node / edge selection
# ---------------------------------------------------------------------------

def select_subset(
    df: pd.DataFrame,
    dist_matrix: np.ndarray,
    target_id: str,
    k_neighbors: int,
    top_ips: int,
) -> pd.DataFrame:
    """Select the proteins to include in the network.

    Selection:
        - The target protein itself.
        - Its k nearest informational neighbors (full dataset distances).
        - Top-N proteins by IPS (global priority context).

    Using distances from the FULL dataset ensures that neighbor
    relationships are not distorted by re-scaling a subset.

    Args:
        df:          Full validated DataFrame.
        dist_matrix: Full pairwise distance matrix.
        target_id:   Prot_ID of the target hypothetical protein.
        k_neighbors: Number of nearest neighbors to include.
        top_ips:     Number of top-IPS proteins to include.

    Returns:
        Subset DataFrame (reset index), aligned with a new distance matrix.

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

    # Nearest neighbors (excluding self: distance = 0)
    sorted_idx = np.argsort(dists_to_target)
    neighbors = [i for i in sorted_idx if i != target_pos][:k_neighbors]

    # Top IPS
    top_idx = df["IPS"].nlargest(top_ips).index.tolist()

    # Union, always include target
    selected = list({target_pos} | set(neighbors) | set(top_idx))
    df_sub = df.iloc[selected].copy().reset_index(drop=True)

    logger.info(
        "Subset: %d proteins (%d neighbors + %d top-IPS + target).",
        len(df_sub), len(neighbors), top_ips,
    )
    return df_sub


# ---------------------------------------------------------------------------
# Graph construction
# ---------------------------------------------------------------------------

def build_graph(
    df_sub: pd.DataFrame,
    dist_matrix_full: np.ndarray,
    df_full: pd.DataFrame,
    threshold: float,
    target_id: str,
) -> nx.Graph:
    """Build a NetworkX graph from the subset with biologically annotated nodes.

    Edges are added only when the full-dataset distance between two proteins
    is below the threshold — not the subset-recalculated distance.

    Node attributes:
        product, is_hypo, func_label, color, ips, efficiency, kolmogorov

    Edge attributes:
        distance, weight (= 1 - distance, for layout algorithms)

    Args:
        df_sub:          Subset DataFrame.
        dist_matrix_full: Full pairwise distance matrix.
        df_full:          Full DataFrame (to look up row positions).
        threshold:        Maximum distance for edge creation.
        target_id:        Prot_ID of the target.

    Returns:
        Annotated NetworkX Graph.
    """
    G = nx.Graph()
    full_ids = df_full["Prot_ID"].tolist()

    # Add nodes with attributes
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

    # Add edges using FULL dataset distances (no re-scaling)
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
# Visualization
# ---------------------------------------------------------------------------

def _make_label(prot_id: str, product: str, max_chars: int = 12) -> str:
    """Generate a readable two-line node label.

    Format: short numeric ID suffix + truncated product keyword.
    Avoids splitting on spaces which can produce biologically meaningless tokens.

    Args:
        prot_id:   Full protein ID string.
        product:   Full product description.
        max_chars: Maximum characters for the product fragment.

    Returns:
        Two-line label string.
    """
    short_id = prot_id.split("_")[-1]   # last numeric suffix
    # Take the first meaningful word after dropping generic terms
    skip = {"protein", "the", "of", "a", "an", "and", "or", "with"}
    words = [w for w in product.lower().split() if w not in skip and len(w) > 2]
    keyword = words[0][:max_chars] if words else "?"
    return f"{short_id}\n({keyword})"


def draw_network(
    G: nx.Graph,
    target_id: str,
    config: NetworkConfig,
    output_path: Path,
    threshold: float,
) -> None:
    """Render and save the informational similarity network.

    Visual encoding:
        Node color   : functional category (red=target, orange=hypothetical,
                       blue shades=annotated by class, gray=uncategorized)
        Node size    : IPS score (larger = higher informational priority)
        Edge width   : inversely proportional to √JSD distance (thicker = closer)
        Edge opacity : same as width

    Args:
        G:           Annotated NetworkX graph.
        target_id:   Prot_ID to mark as the focal node.
        config:      NetworkConfig instance.
        output_path: Destination PNG/SVG path.
        threshold:   Used in the plot subtitle for interpretability.
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
        # Width and alpha inversely proportional to distance (closer = more visible)
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
    labels = {
        n: _make_label(n, G.nodes[n]["product"])
        for n in G.nodes
    }
    nx.draw_networkx_labels(G, pos, labels, ax=ax, font_size=7.5, font_weight="bold")

    # ── Legend ───────────────────────────────────────────────────────────────
    legend_items: List[mpatches.Patch] = [
        mpatches.Patch(color=TARGET_COLOR, label=f"Target: {target_id}"),
        mpatches.Patch(color=HYPOTHETICAL_COLOR, label="Hypothetical protein"),
        mpatches.Patch(color=ANNOTATED_DEFAULT_COLOR, label="Annotated (other)"),
    ]
    seen_labels = set()
    for n in G.nodes:
        lbl = G.nodes[n]["func_label"]
        col = G.nodes[n]["color"]
        if lbl not in seen_labels and col not in (TARGET_COLOR, HYPOTHETICAL_COLOR, ANNOTATED_DEFAULT_COLOR):
            legend_items.append(mpatches.Patch(color=col, label=lbl))
            seen_labels.add(lbl)

    ax.legend(handles=legend_items, loc="upper left", fontsize=9,
              framealpha=0.85, title="Functional category")

    # ── IPS size guide ───────────────────────────────────────────────────────
    for ips_val, label in [(0.25, "IPS=0.25"), (0.5, "IPS=0.50"), (0.75, "IPS=0.75")]:
        size = config.node_size_base + config.node_size_ips_scale * ips_val
        ax.scatter([], [], s=size, c="gray", alpha=0.5,
                   edgecolors="black", linewidths=0.8, label=label)
    ax.legend(
        handles=legend_items + [
            plt.scatter([], [], s=config.node_size_base + config.node_size_ips_scale * v,
                        c="gray", alpha=0.5, edgecolors="black", linewidths=0.8,
                        label=f"IPS={v:.2f}")
            for v in [0.25, 0.50, 0.75]
        ],
        loc="upper left", fontsize=8, framealpha=0.85,
        title="Node color = functional class\nNode size ∝ IPS",
    )

    # ── Title and metadata ────────────────────────────────────────────────────
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
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="netentro",
        description=(
            f"netentro v{__version__} — Informational Similarity Network\n\n"
            "Builds a protein network where edges represent informational\n"
            "proximity (weighted √JSD distance) from bioentro output."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  %(prog)s -i protein.tsv -t NFOBNJ_05247 -o network.png\n"
            "  %(prog)s -i protein.tsv -t NFOBNJ_05247 -k 30 --top-ips 15 -o net.png\n"
            "  %(prog)s -i protein.tsv --list-ids\n"
        ),
    )
    parser.add_argument("-i", "--input", required=True,
                        help="bioentro 'protein' mode TSV output")
    parser.add_argument("-t", "--target",
                        help="Prot_ID of the target hypothetical protein")
    parser.add_argument("-o", "--output", default="network.png",
                        help="Output image path (default: network.png)")
    parser.add_argument("-k", "--k-neighbors", type=int, default=25,
                        help="Number of nearest informational neighbors (default: 25)")
    parser.add_argument("--top-ips", type=int, default=10,
                        help="Additional top-IPS proteins to include (default: 10)")
    parser.add_argument("--threshold-percentile", type=float, default=20.0,
                        help="Fallback connection threshold percentile (default: 20.0)")
    parser.add_argument("--calibration-k", type=int, default=5,
                        help="k for threshold calibration from annotated proteins (default: 5)")
    parser.add_argument("--no-calibration", action="store_true",
                        help="Skip calibration, use --threshold-percentile directly")
    parser.add_argument("--layout-k", type=float, default=2.5,
                        help="Spring layout repulsion constant (default: 2.5)")
    parser.add_argument("--layout-seed", type=int, default=42,
                        help="Random seed for reproducible layout (default: 42)")
    parser.add_argument("--dpi", type=int, default=300,
                        help="Output image DPI (default: 300)")
    parser.add_argument("--list-ids", action="store_true",
                        help="Print first 50 Prot_IDs and exit")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable debug-level logging")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    """CLI entry point.

    Returns:
        0 on success, 1 on handled error, 130 on keyboard interrupt.
    """
    parser = _build_parser()
    args = parser.parse_args(argv)
    _configure_logging(verbose=args.verbose)

    input_path = Path(args.input)

    try:
        df = load_and_validate(input_path)
    except NetentroError as exc:
        logger.error("%s", exc)
        return 1

    # Utility mode: list IDs and exit
    if args.list_ids:
        print("\n".join(df["Prot_ID"].head(50).tolist()))
        return 0

    if not args.target:
        parser.error("--target is required unless --list-ids is used.")

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

    try:
        dist_matrix = build_distance_matrix(df)

        if config.use_calibrated_threshold:
            threshold = calibrate_threshold(df, dist_matrix, config.calibration_k)
        else:
            threshold = float(np.percentile(dist_matrix[dist_matrix > 0],
                                             config.threshold_percentile))
            logger.info("Using percentile-%.1f threshold: %.4f",
                        config.threshold_percentile, threshold)

        df_sub = select_subset(
            df, dist_matrix, config.target_id,
            config.k_neighbors, config.top_ips,
        )
        G = build_graph(df_sub, dist_matrix, df, threshold, config.target_id)
        draw_network(G, config.target_id, config, Path(args.output), threshold)
        return 0

    except NetentroError as exc:
        logger.error("%s", exc)
        return 1
    except KeyboardInterrupt:
        logger.warning("Interrupted by user.")
        return 130


if __name__ == "__main__":
    sys.exit(main())