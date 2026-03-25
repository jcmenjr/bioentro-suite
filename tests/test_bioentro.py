"""
Unit tests for bioentro core informational metric functions.

These tests verify mathematical correctness of each metric independently.
They do NOT require any biological input files — all inputs are synthetic.

Run with:
    pytest tests/test_bioentro.py -v
"""

import math
import pytest

from bioentro.bioentro import (
    compute_kmer_dist,
    calculate_shannon,
    calculate_h_max,
    calculate_jsd,
    calculate_jsd_vs_uniform,
    calculate_kl,
    calculate_kolmogorov,
    calculate_ips,
    _safe_efficiency,
    INVALID_AA,
    INVALID_NT,
)


# ---------------------------------------------------------------------------
# compute_kmer_dist
# ---------------------------------------------------------------------------

class TestComputeKmerDist:
    def test_basic_protein_k1(self):
        dist = compute_kmer_dist("ACDEFGH", k=1, is_protein=True)
        assert set(dist.keys()) == set("ACDEFGH")
        assert abs(sum(dist.values()) - 1.0) < 1e-9

    def test_basic_dna_k1(self):
        dist = compute_kmer_dist("AACCGGTT", k=1, is_protein=False)
        assert set(dist.keys()) == {"A", "C", "G", "T"}
        assert abs(sum(dist.values()) - 1.0) < 1e-9

    def test_invalid_aa_excluded(self):
        # 'X' and '*' should be excluded from k-mer construction
        dist = compute_kmer_dist("AXCX*D", k=1, is_protein=True)
        assert "X" not in dist
        assert "*" not in dist
        assert "A" in dist and "C" in dist and "D" in dist

    def test_invalid_nt_excluded(self):
        dist = compute_kmer_dist("AANTCG", k=1, is_protein=False)
        assert "N" not in dist

    def test_empty_sequence_returns_empty(self):
        assert compute_kmer_dist("", k=1, is_protein=True) == {}

    def test_sequence_shorter_than_k_returns_empty(self):
        assert compute_kmer_dist("AC", k=6, is_protein=False) == {}

    def test_uniform_sequence_single_kmer(self):
        dist = compute_kmer_dist("AAAA", k=1, is_protein=False)
        assert dist == {"A": 1.0}

    def test_normalized_sum_k2(self):
        dist = compute_kmer_dist("ACGTACGT", k=2, is_protein=False)
        assert abs(sum(dist.values()) - 1.0) < 1e-9


# ---------------------------------------------------------------------------
# calculate_h_max
# ---------------------------------------------------------------------------

class TestCalculateHMax:
    def test_protein_k1(self):
        # k=1, A=20 → log2(20) ≈ 4.3219
        h = calculate_h_max(k=1, is_protein=True)
        assert abs(h - math.log2(20)) < 1e-9

    def test_dna_k6(self):
        # k=6, A=4 → 6 * log2(4) = 6 * 2 = 12.0
        h = calculate_h_max(k=6, is_protein=False)
        assert abs(h - 12.0) < 1e-9

    def test_numerically_stable(self):
        # k * log2(A) should equal log2(A^k) for reasonable k
        h_stable = calculate_h_max(k=10, is_protein=False)
        h_naive = math.log2(4 ** 10)
        assert abs(h_stable - h_naive) < 1e-9


# ---------------------------------------------------------------------------
# calculate_shannon
# ---------------------------------------------------------------------------

class TestCalculateShannon:
    def test_uniform_distribution_is_h_max(self):
        # 4 equiprobable k-mers → H = log2(4) = 2.0
        dist = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}
        h = calculate_shannon(dist)
        assert abs(h - 2.0) < 1e-9

    def test_deterministic_distribution_is_zero(self):
        dist = {"A": 1.0}
        assert calculate_shannon(dist) == 0.0

    def test_empty_distribution_is_zero(self):
        assert calculate_shannon({}) == 0.0

    def test_two_equiprobable(self):
        dist = {"A": 0.5, "C": 0.5}
        assert abs(calculate_shannon(dist) - 1.0) < 1e-9


# ---------------------------------------------------------------------------
# calculate_jsd
# ---------------------------------------------------------------------------

class TestCalculateJSD:
    def test_identical_distributions_zero(self):
        dist = {"A": 0.5, "C": 0.5}
        assert calculate_jsd(dist, dist) == 0.0

    def test_completely_disjoint_distributions_max(self):
        # No overlap → JSD should be 1.0
        p = {"A": 1.0}
        q = {"C": 1.0}
        jsd = calculate_jsd(p, q)
        assert abs(jsd - 1.0) < 1e-4

    def test_symmetric(self):
        p = {"A": 0.7, "C": 0.3}
        q = {"A": 0.2, "G": 0.8}
        assert abs(calculate_jsd(p, q) - calculate_jsd(q, p)) < 1e-9

    def test_bounded_between_0_and_1(self):
        p = {"A": 0.6, "C": 0.4}
        q = {"G": 0.5, "T": 0.5}
        jsd = calculate_jsd(p, q)
        assert 0.0 <= jsd <= 1.0

    def test_sqrt_jsd_triangle_inequality(self):
        # Verify √JSD satisfies the triangle inequality for three distributions
        p = {"A": 0.8, "C": 0.2}
        q = {"A": 0.5, "C": 0.5}
        r = {"C": 0.9, "G": 0.1}
        d_pq = math.sqrt(calculate_jsd(p, q))
        d_qr = math.sqrt(calculate_jsd(q, r))
        d_pr = math.sqrt(calculate_jsd(p, r))
        assert d_pr <= d_pq + d_qr + 1e-9

    def test_empty_inputs_return_zero(self):
        assert calculate_jsd({}, {"A": 1.0}) == 0.0
        assert calculate_jsd({"A": 1.0}, {}) == 0.0


# ---------------------------------------------------------------------------
# calculate_jsd_vs_uniform
# ---------------------------------------------------------------------------

class TestCalculateJSDVsUniform:
    def test_uniform_input_returns_zero(self):
        # A perfectly uniform distribution should have JSD=0 against uniform
        # With k=1, DNA: 4 k-mers each at p=0.25
        dist = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}
        jsd = calculate_jsd_vs_uniform(dist, k=1, is_protein=False)
        assert jsd < 1e-4

    def test_empty_input_returns_zero(self):
        assert calculate_jsd_vs_uniform({}, k=1, is_protein=True) == 0.0

    def test_bounded(self):
        dist = {"A": 0.9, "C": 0.1}
        jsd = calculate_jsd_vs_uniform(dist, k=1, is_protein=False)
        assert 0.0 <= jsd <= 1.0


# ---------------------------------------------------------------------------
# calculate_kl
# ---------------------------------------------------------------------------

class TestCalculateKL:
    def test_identical_distributions_zero(self):
        dist = {"A": 0.5, "C": 0.5}
        assert calculate_kl(dist, dist) == 0.0

    def test_non_negative(self):
        p = {"A": 0.7, "C": 0.3}
        q = {"A": 0.4, "C": 0.6}
        assert calculate_kl(p, q) >= 0.0

    def test_asymmetric(self):
        p = {"A": 0.8, "C": 0.2}
        q = {"A": 0.3, "C": 0.7}
        assert calculate_kl(p, q) != calculate_kl(q, p)

    def test_smoothing_for_unseen_kmers(self):
        # q does not contain "G" — should not raise, should use smoothing
        p = {"A": 0.5, "G": 0.5}
        q = {"A": 1.0}
        result = calculate_kl(p, q)
        assert result > 0.0 and math.isfinite(result)

    def test_empty_inputs_return_zero(self):
        assert calculate_kl({}, {"A": 1.0}) == 0.0
        assert calculate_kl({"A": 1.0}, {}) == 0.0


# ---------------------------------------------------------------------------
# calculate_kolmogorov
# ---------------------------------------------------------------------------

class TestCalculateKolmogorov:
    def test_empty_sequence_returns_zero(self):
        assert calculate_kolmogorov("") == 0.0

    def test_highly_repetitive_is_low(self):
        # A very repetitive sequence should compress well → low ratio
        ratio = calculate_kolmogorov("ATCG" * 500)
        assert ratio < 0.5

    def test_ratio_between_0_and_1(self):
        ratio = calculate_kolmogorov("ATCGATCGATCG")
        assert 0.0 < ratio <= 1.0

    def test_returns_float(self):
        assert isinstance(calculate_kolmogorov("ACGT"), float)


# ---------------------------------------------------------------------------
# calculate_ips
# ---------------------------------------------------------------------------

class TestCalculateIPS:
    def test_zero_efficiency_gives_zero(self):
        assert calculate_ips(efficiency=0.0, jsd_bg=0.5, length=500, l0=100) == 0.0

    def test_zero_jsd_gives_zero(self):
        assert calculate_ips(efficiency=0.8, jsd_bg=0.0, length=500, l0=100) == 0.0

    def test_very_short_sequence_is_penalized(self):
        ips_short = calculate_ips(efficiency=0.9, jsd_bg=0.9, length=10,  l0=100)
        ips_long  = calculate_ips(efficiency=0.9, jsd_bg=0.9, length=500, l0=100)
        assert ips_short < ips_long

    def test_bounded_0_to_1(self):
        ips = calculate_ips(efficiency=1.0, jsd_bg=1.0, length=1000, l0=100)
        assert 0.0 <= ips <= 1.0

    def test_perfect_values_approach_1(self):
        ips = calculate_ips(efficiency=1.0, jsd_bg=1.0, length=10_000, l0=100)
        assert ips > 0.99


# ---------------------------------------------------------------------------
# _safe_efficiency
# ---------------------------------------------------------------------------

class TestSafeEfficiency:
    def test_normal_case(self):
        assert abs(_safe_efficiency(3.0, 4.0) - 0.75) < 1e-9

    def test_zero_h_max_returns_zero(self):
        assert _safe_efficiency(1.0, 0.0) == 0.0

    def test_perfect_efficiency(self):
        assert abs(_safe_efficiency(4.0, 4.0) - 1.0) < 1e-9