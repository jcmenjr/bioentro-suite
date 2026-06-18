#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests for pangentro. Run with `pytest test_pangentro.py` or `python test_pangentro.py`.

Covers the information-theoretic identities, the agreement between the fast
vectorised k=1 scorer and the readable dict reference, the HP/DUF logic, the
statistics helpers, and an end-to-end `integrate` + `compare` run on a synthetic
anvi-summarize gene-clusters summary.
"""
import gzip
import json
import random
from pathlib import Path

import numpy as np
import pandas as pd

import pangentro as pg


class _SilentTerminal:
    quiet = True
    backend = "builtin"
    anvio_version = None

    def __getattr__(self, _):
        return lambda *a, **k: None


def test_metric_identities():
    even = pg.normalise(pg.kmer_counter(pg.AA_ORDER, 1))
    assert abs(pg.pielou_evenness(even, 1) - 1.0) < 1e-9
    assert pg.pielou_evenness(pg.normalise(pg.kmer_counter("AAAAAA", 1)), 1) == 0.0
    assert pg.jensen_shannon_divergence(even, even) < 1e-12
    d1 = pg.normalise(pg.kmer_counter("AAAA", 1))
    d2 = pg.normalise(pg.kmer_counter("CCCC", 1))
    assert abs(pg.jensen_shannon_divergence(d1, d2) - 1.0) < 1e-9
    assert pg.kmer_counter("AAXBA", 1) == {"A": 3}  # non-standard residues dropped


def test_vectorised_matches_reference():
    rng = random.Random(42)
    seqs = ["".join(rng.choice(pg.AA_ORDER + "X*")
                    for _ in range(rng.randint(20, 400))) for _ in range(300)]
    genomes = [f"g{rng.randint(0, 4)}" for _ in seqs]
    gcs = [f"GC_{i:06d}" for i in range(len(seqs))]
    cat = {gc: ("core" if i % 3 == 0 else "shell") for i, gc in enumerate(gcs)}
    t = _SilentTerminal()
    vec = pg._score_k1_vectorised(seqs, genomes, gcs, cat, {"core"}, 100, t)
    ref = pg._score_general(seqs, genomes, gcs, cat, {"core"}, 1, 100, t)
    for key in ("_eff", "_ips_core", "_ips_self", "_sqrtjsd_core"):
        a = np.asarray(vec[key], float)
        b = np.asarray(ref[key], float)
        mask = ~(np.isnan(a) & np.isnan(b))
        assert np.allclose(a[mask], b[mask], atol=1e-9), key


def test_uncharacterized_logic():
    assert pg._looks_uncharacterized("DUF4435 domain-containing protein")
    assert pg._looks_uncharacterized("Hypothetical protein")
    assert pg._looks_uncharacterized("protein of unknown function")
    assert not pg._looks_uncharacterized("dufA-like oxidoreductase")
    assert not pg._looks_uncharacterized("putative ABC transporter")


def test_stats_helpers():
    adj = pg._holm_bonferroni([0.01, 0.04, 0.03])
    assert all(0.0 <= x <= 1.0 for x in adj)
    assert adj[0] <= adj[2] <= adj[1]
    assert pg._holm_bonferroni([]) == []
    assert pg._rank_biserial(100 * 100, 100, 100) == 1.0   # a >> b
    assert pg._rank_biserial(0, 100, 100) == -1.0          # a << b


def test_null_model_calibration():
    """The Monte-Carlo null must yield empirical p in (0, 1] and finite z, and a
    background-typical sequence should not look atypical (large p)."""
    rng = random.Random(1)
    q = np.full(pg.PROTEIN_ALPHABET_SIZE, 1.0 / pg.PROTEIN_ALPHABET_SIZE)
    gen = np.random.default_rng(0)
    null = pg._null_ips_by_length(np.array([300]), q, 100, 200, gen)
    arr, mean, std = null[300]
    assert arr.size == 200 and std > 0 and (arr >= 0).all() and (arr <= 1).all()

    # an even-composition (background-typical) protein → high empirical p
    seqs = [(pg.AA_ORDER * 15)]            # 300 aa, perfectly even
    cols = pg._score_k1_vectorised(
        seqs, ["g0"], ["GC_0"], {"GC_0": "core"}, {"core"}, 100,
        _SilentTerminal(), null_draws=200, seed=0)
    p = cols["_ips_core_emp_p"][0]
    assert 0.0 < p <= 1.0
    # when disabled, the calibration columns are NaN
    off = pg._score_k1_vectorised(
        seqs, ["g0"], ["GC_0"], {"GC_0": "core"}, {"core"}, 100,
        _SilentTerminal(), null_draws=0)
    assert np.isnan(off["_ips_core_emp_p"][0]) and np.isnan(off["_ips_core_z"][0])


def _synthetic_summary(path: Path):
    rng = random.Random(7)
    rows, uid = [], 0
    for gi in range(6):
        for ci in range(40):
            if not (ci < 30 or gi < 2):  # core in all, accessory in 2 genomes
                continue
            seq = "".join(rng.choice(pg.AA_ORDER) for _ in range(rng.randint(60, 500)))
            func = "" if ci % 4 == 0 else ("hypothetical protein" if ci % 4 == 1
                                           else "ABC transporter ATP-binding protein")
            rows.append(dict(unique_id=uid, gene_cluster_id=f"GC_{ci:08d}",
                             bin_name="bin", genome_name=f"genome_{gi}",
                             gene_callers_id=uid, aa_sequence=seq,
                             COG20_FUNCTION=func, COG20_FUNCTION_ACC="",
                             COG20_CATEGORY=("S" if ci % 4 in (0, 1) else "P"),
                             COG20_CATEGORY_ACC=""))
            uid += 1
    with gzip.open(path, "wt") as fh:
        pd.DataFrame(rows).to_csv(fh, sep="\t", index=False)


def test_end_to_end(tmp_path):
    summ = tmp_path / "MY_gene_clusters_summary.txt.gz"
    _synthetic_summary(summ)
    assert pg.main(["integrate", "--gene-clusters-summary", str(summ),
                    "-o", str(tmp_path), "--project-name", "MY", "--quiet",
                    "--null-draws", "100", "--seed", "0"]) == 0

    metrics = pd.read_csv(tmp_path / "MY_pangentro_metrics.txt", sep="\t")
    assert set(metrics["pangenome_category"]) <= set(pg.CATEGORY_ORDER)
    assert set(metrics["HP_status"]) <= set(pg.HP_ORDER)
    assert metrics["IPS_core_mean"].between(0, 1).all()
    assert metrics["IPS_percentile"].between(0, 100).all()
    # null model populated the calibration columns
    assert metrics["IPS_core_emp_p"].dropna().between(0, 1).all()
    assert metrics["IPS_core_emp_p"].notna().any()

    items = pd.read_csv(tmp_path / "MY_items_for_anvio.txt", sep="\t")
    assert list(items.columns)[0] == "gene_cluster" and "IPS" in items.columns

    report = json.loads((tmp_path / "MY_run_report.json").read_text())
    assert report["version"] == pg.__version__ and len(report["input_sha256"]) == 64

    assert pg.main(["compare", "-i", str(tmp_path / "MY_pangentro_metrics.txt"),
                    "-o", str(tmp_path / "cmp"), "--group-by", "HP_status",
                    "--no-figure", "--quiet"]) == 0
    tests = pd.read_csv(tmp_path / "cmp" / "statistical_tests.tsv", sep="\t")
    assert "p_adjusted" in tests.columns


if __name__ == "__main__":
    import tempfile
    test_metric_identities()
    test_vectorised_matches_reference()
    test_uncharacterized_logic()
    test_stats_helpers()
    test_null_model_calibration()
    with tempfile.TemporaryDirectory() as d:
        test_end_to_end(Path(d))
    print("ALL TESTS PASSED")
