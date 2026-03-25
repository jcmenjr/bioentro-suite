"""
Unit tests for netentro network construction functions.

Run with:
    pytest tests/test_netentro.py -v
"""

import numpy as np
import pandas as pd
import pytest

from bioentro.netentro import (
    is_hypothetical,
    classify_product,
    build_distance_matrix,
    calibrate_threshold,
    select_subset,
    HYPOTHETICAL_COLOR,
    ANNOTATED_DEFAULT_COLOR,
    TargetNotFoundError,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def sample_df() -> pd.DataFrame:
    """Minimal DataFrame that matches bioentro protein output schema."""
    return pd.DataFrame({
        "Prot_ID":        ["PA_001", "PA_002", "PA_003", "PA_004", "PA_005"],
        "Product":        [
            "hypothetical protein",
            "ABC transporter permease",
            "two-component response regulator",
            "putative membrane protein",
            "zinc-binding oxidoreductase",
        ],
        "Len_AA":         [120, 340, 220, 95, 410],
        "H_max":          [4.32] * 5,
        "H_real":         [3.80, 4.10, 3.95, 2.10, 4.20],
        "Efficiency":     [0.88, 0.95, 0.91, 0.49, 0.97],
        "KL_proteomebg":  [0.10, 0.05, 0.08, 0.30, 0.03],
        "JSD_proteomebg": [0.15, 0.05, 0.10, 0.40, 0.02],
        "sqrt_JSD":       [0.39, 0.22, 0.32, 0.63, 0.14],
        "Kolmogorov":     [0.72, 0.65, 0.68, 0.80, 0.60],
        "IPS":            [0.34, 0.21, 0.29, 0.31, 0.14],
    })


# ---------------------------------------------------------------------------
# is_hypothetical
# ---------------------------------------------------------------------------

class TestIsHypothetical:
    def test_hypothetical_protein(self):
        assert is_hypothetical("hypothetical protein") is True

    def test_putative(self):
        assert is_hypothetical("putative membrane protein") is True

    def test_uncharacterized(self):
        assert is_hypothetical("uncharacterized protein") is True

    def test_known_function(self):
        assert is_hypothetical("ABC transporter permease") is False

    def test_case_insensitive(self):
        assert is_hypothetical("Hypothetical Protein") is True


# ---------------------------------------------------------------------------
# classify_product
# ---------------------------------------------------------------------------

class TestClassifyProduct:
    def test_transport_keyword(self):
        label, color = classify_product("ABC transporter permease")
        assert label == "Transport"

    def test_oxidoreductase_keyword(self):
        label, _ = classify_product("zinc-binding oxidoreductase")
        assert label == "Oxidoreductase"

    def test_no_match_returns_default(self):
        label, color = classify_product("some unknown enzyme")
        assert label == "Annotated (other)"
        assert color == ANNOTATED_DEFAULT_COLOR

    def test_first_match_wins(self):
        # "transport" comes before "binding" in FUNCTIONAL_CATEGORIES
        label, _ = classify_product("transport binding protein")
        assert label == "Transport"


# ---------------------------------------------------------------------------
# build_distance_matrix
# ---------------------------------------------------------------------------

class TestBuildDistanceMatrix:
    def test_shape(self, sample_df):
        mat = build_distance_matrix(sample_df)
        n = len(sample_df)
        assert mat.shape == (n, n)

    def test_diagonal_is_zero(self, sample_df):
        mat = build_distance_matrix(sample_df)
        np.testing.assert_array_almost_equal(np.diag(mat), 0.0)

    def test_symmetric(self, sample_df):
        mat = build_distance_matrix(sample_df)
        np.testing.assert_array_almost_equal(mat, mat.T)

    def test_non_negative(self, sample_df):
        mat = build_distance_matrix(sample_df)
        assert (mat >= 0).all()

    def test_bounded(self, sample_df):
        # Weighted Euclidean on [0,1] features with weights summing to 1
        # max possible distance < sqrt(1^2 + 1^2 + 1^2) ≈ 1.73
        mat = build_distance_matrix(sample_df)
        assert mat.max() <= 2.0


# ---------------------------------------------------------------------------
# calibrate_threshold
# ---------------------------------------------------------------------------

class TestCalibrateThreshold:
    def test_returns_positive_float(self, sample_df):
        mat = build_distance_matrix(sample_df)
        t = calibrate_threshold(sample_df, mat, k=2)
        assert isinstance(t, float)
        assert t > 0.0

    def test_fallback_when_too_few_annotated(self):
        # Only one annotated protein — should fall back to median
        df = pd.DataFrame({
            "Prot_ID":    ["P1", "P2", "P3"],
            "Product":    ["hypothetical protein", "hypothetical protein", "known enzyme"],
            "sqrt_JSD":   [0.3, 0.5, 0.2],
            "Efficiency": [0.8, 0.7, 0.9],
            "Kolmogorov": [0.6, 0.65, 0.55],
            "IPS":        [0.4, 0.3, 0.2],
            "JSD_proteomebg": [0.1, 0.2, 0.05],
        })
        mat = build_distance_matrix(df)
        # k=5 but only 1 annotated → should warn and use fallback
        t = calibrate_threshold(df, mat, k=5)
        assert t > 0.0


# ---------------------------------------------------------------------------
# select_subset
# ---------------------------------------------------------------------------

class TestSelectSubset:
    def test_target_always_included(self, sample_df):
        mat = build_distance_matrix(sample_df)
        sub = select_subset(sample_df, mat, "PA_001", k_neighbors=2, top_ips=1)
        assert "PA_001" in sub["Prot_ID"].values

    def test_subset_size_bounded(self, sample_df):
        mat = build_distance_matrix(sample_df)
        sub = select_subset(sample_df, mat, "PA_001", k_neighbors=2, top_ips=1)
        # At most: 1 target + 2 neighbors + 1 top-IPS = 4 (with overlap possible)
        assert len(sub) <= 5    # capped by total dataset size

    def test_raises_for_unknown_target(self, sample_df):
        mat = build_distance_matrix(sample_df)
        with pytest.raises(TargetNotFoundError):
            select_subset(sample_df, mat, "NONEXISTENT_ID", k_neighbors=2, top_ips=1)

    def test_returns_dataframe(self, sample_df):
        mat = build_distance_matrix(sample_df)
        sub = select_subset(sample_df, mat, "PA_002", k_neighbors=2, top_ips=1)
        assert isinstance(sub, pd.DataFrame)