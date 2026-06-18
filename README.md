# bioentro-suite · pangentro

[![CI](https://github.com/jcmenjr/bioentro-suite/actions/workflows/ci.yml/badge.svg)](https://github.com/jcmenjr/bioentro-suite/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/)

**Informational triage of hypothetical proteins, native to [anvi'o 9 "eunice"](https://anvio.org).**

`pangentro` reads the gene-clusters summary that `anvi-summarize` produces from a
pangenome (pan-db), scores every gene cluster (GC) by how *compositionally
atypical* its protein sequences are, flags an operational hypothetical-protein
(HP) status, and writes the result straight back into the pan-db so you can
sort, colour, and bin gene clusters by priority inside `anvi-display-pan`.

It is a **biology-blind triage** step: the Informational Priority Score (IPS)
does not predict function. It only tells you *where to spend a finite budget* of
downstream effort — AlphaFold/ESMFold, Foldseek, manual curation, pan-GWAS — by
ranking the proteins that are least like the rest of the proteome.

```
IPS  =  J'  ·  √JSD  ·  (1 − e^(−L/L₀))
```

| term | meaning |
|------|---------|
| `J'` | Pielou evenness of the amino-acid composition (normalised Shannon entropy). A low-complexity filter, background-independent. |
| `√JSD` | Jensen–Shannon **distance** to a background proteome — a true metric. The primary ranking driver. |
| `1 − e^(−L/L₀)` | length penalty; short ORFs (≪ `L₀`, default 100 aa) are down-weighted. |

---

## ⚠️ Migration note (v0.5 — full anvi'o transition)

`pangentro` is now **anvi'o-native and eggNOG-mapper-aware end to end**. The old
Bakta → Panaroo path (and the `preparo` helper that fed Panaroo) is
**deprecated and removed** — `pangentro` no longer reads Panaroo's
`pan_genome_reference.fa`. If you have a v0.3.x project, rebuild it through the
anvi'o workflow below. There is a single supported track now.

---

## Install

```bash
# inside your anvi'o conda environment (recommended, so anvi-* are on PATH)
pip install -e ".[stats]"      # stats extra = scipy + matplotlib for `compare`
# or just the core:
pip install -e .
```

This installs a `pangentro` console command. Core needs only `numpy` + `pandas`;
`compare` additionally needs `scipy` (+ `matplotlib` for the figure). anvi'o
itself is optional — only required if you want pangentro to call
`anvi-summarize` / import for you, or to print with anvi'o's own terminal style.

---

## The workflow

### 1 · Upstream is plain anvi'o (annotate **before** the genomes storage)

```bash
# one contigs-db per genome
anvi-gen-contigs-database -f genome_01.fa -o genome_01.db
anvi-run-hmms -c genome_01.db

# --- functional annotation -------------------------------------------------
# COG20 (gives both functions AND the single-letter COG categories pangentro
# uses to recognise category 'S' = function unknown):
anvi-run-ncbi-cogs -c genome_01.db

# eggNOG-mapper: run it on the gene-call amino-acid sequences, then import.
anvi-get-sequences-for-gene-calls -c genome_01.db --get-aa-sequences -o g01.faa
emapper.py -i g01.faa --itype proteins -o g01 --cpu 8
#   convert g01.emapper.annotations → an anvi'o functions table with columns:
#     gene_callers_id  source  accession  function  e_value
#   choosing a clear `source` name (e.g. EGGNOG) — that name becomes the
#   summary column pangentro auto-detects. (Put the eggNOG COG letters under a
#   source/column whose name ends in CATEGORY, or pass --category-source later.)
anvi-import-functions -c genome_01.db -i g01_functions.txt
```

> **Order matters.** Import functions into each contigs-db *before* you build the
> genomes storage — anvi'o snapshots functions at that point, so anything added
> later will not appear in the pangenome summary.

```bash
# genomes storage → pangenome → a default collection → summary
anvi-gen-genomes-storage -e external-genomes.txt -o MY-GENOMES.db
anvi-pan-genome -g MY-GENOMES.db -n MYPROJECT
anvi-script-add-default-collection -p MYPROJECT/MYPROJECT-PAN.db -C DEFAULT
anvi-summarize -p MYPROJECT/MYPROJECT-PAN.db -g MY-GENOMES.db \
               -C DEFAULT -o SUMMARY/
```

`anvi-summarize` writes `SUMMARY/MYPROJECT_gene_clusters_summary.txt.gz` (amino
acids by default — do **not** pass `--report-DNA-sequences`; pangentro scores in
protein space and will warn if it sees nucleotides).

### 2 · Run pangentro

```bash
pangentro integrate \
    --gene-clusters-summary SUMMARY/MYPROJECT_gene_clusters_summary.txt.gz \
    --pan-db MYPROJECT/MYPROJECT-PAN.db \
    --null-draws 200 \
    -o pangentro_out/ --project-name MYPROJECT
```

Don't have a summary yet? Give pangentro the pan-db + storage and it will run
`anvi-summarize` for you:

```bash
pangentro integrate -p MYPROJECT/MYPROJECT-PAN.db -g MY-GENOMES.db \
    --summary-dir SUMMARY/ -o pangentro_out/ --project-name MYPROJECT
```

---

## What it writes

| file | what it is |
|------|------------|
| `*_pangentro_metrics.txt` | the full per-GC table (every column below) |
| `*_items_for_anvio.txt` | **anvi'o misc-data-items** file — import to decorate the pangenome |
| `*_HP_prioritized.txt` | hypothetical GCs only, sorted by IPS (your triage worklist) |
| `*_HP_collection.txt` | an **anvi'o collection** of the top-N HP clusters, binned by HP status |
| `*_run_report.json` | provenance: version, input SHA-256, parameters, environment, counts, QC |

Key columns in the metrics table:

- `pangenome_category` — `core` / `soft_core` / `shell` / `cloud`, from GC
  occupancy across genomes (thresholds are explicit, tunable flags).
- `IPS_core_mean` / `IPS_core_max` / `IPS_self_mean` — IPS against the shared
  background vs. each gene's own genome (see *Two backgrounds* below).
- `IPS_percentile` — rank (0–100) of the ranking IPS; sort by this in anvi'o.
- `IPS_core_z`, `IPS_core_emp_p` — **null-model calibration** (see below).
- `HP_status` — `ORFan` / `uncharacterized` / `annotated`.
- `efficiency`, `mean_compression_ratio`, `sqrt_JSD_core` — diagnostics.

### 3 · Push it into the pan-db

```bash
anvi-import-misc-data pangentro_out/MYPROJECT_items_for_anvio.txt \
    -p MYPROJECT/MYPROJECT-PAN.db --target-data-table items
anvi-import-collection pangentro_out/MYPROJECT_HP_collection.txt \
    -p MYPROJECT/MYPROJECT-PAN.db -C pangentro_HP
anvi-display-pan -p MYPROJECT/MYPROJECT-PAN.db -g MY-GENOMES.db
```

In the interactive interface, **order items by `IPS` (or `IPS_percentile`) and
colour by `HP_status`** to surface the highest-priority hypothetical clusters,
then open a GC to inspect its aligned homologues.

`pangentro integrate --import --pan-db ...` will run those import commands for
you. anvi'o has no overwrite flag for items misc-data, so re-runs need
`--overwrite`, which first clears pangentro's own layers via
`anvi-delete-misc-data` before re-importing.

---

## The null model (`--null-draws N`)

Raw IPS is a magnitude; the null model turns it into a **calibrated, length-
controlled significance**. For each gene, pangentro draws `N` synthetic proteins
of the *same length* whose residues are i.i.d. from the background distribution
(sampled directly as a Multinomial over the 20 residues — no sequences are
materialised), scores them, and reports the observed IPS as:

- `IPS_core_emp_p` — one-sided empirical p-value, `(1 + #{null ≥ obs}) / (N + 1)`.
  Small ⇒ this protein is more atypical than a typical protein of its length.
- `IPS_core_z` — z-score against the null mean/sd.

Because the length penalty is identical for a gene and its length-matched null,
it cancels out — so the empirical p isolates **compositional** atypicality and
absorbs the length-dependent sampling noise of short proteins (the exact effect
the length-bias QC warns about). Start with `N = 200`; use `1000` for a figure.
The null is reproducible with `--seed` and is computed for `k = 1` (the default).

> A residue-*permutation* null would be degenerate here: shuffling a sequence
> does not change its composition, so it does not change IPS at `k = 1`. The
> background re-draw is the meaningful null.

---

## Statistics — `compare`

```bash
pangentro compare -i pangentro_out/MYPROJECT_pangentro_metrics.txt \
    -o cmp/ --group-by pangenome_category        # or: --group-by HP_status
```

Runs Kruskal–Wallis across groups, then **Holm–Bonferroni-corrected** pairwise
Mann–Whitney U tests with a rank-biserial effect size (`> 0` ⇒ the first group
is stochastically larger). Writes `group_summary.tsv`, `statistical_tests.tsv`
(with `p_adjusted`), and a seeded box/strip-plot figure (`--no-figure` to skip).
`--hypothetical-only` restricts the comparison to HP clusters.

---

## Two backgrounds

- **`IPS_core`** — √JSD against the concatenated background proteome
  (`--bg-categories`, `core` by default). One shared reference space, so values
  are comparable across all GCs. A core GC's own sequences are part of this
  background (mild, *conservative* circularity — it deflates core scores).
- **`IPS_self`** — √JSD against each gene's own genome. No cross-genome
  circularity; the cleaner primary evidence for "HPs concentrate in the
  accessory genome", at the cost of a slightly different reference per genome.

Pick which one decorates anvi'o with `--ips-column`.

---

## Hypothetical-protein status

| status | meaning |
|--------|---------|
| `ORFan` | no hit in any detected annotation source |
| `uncharacterized` | a hit, but only DUFs / "unknown function" / "hypothetical" descriptions, or only a weak COG category (`S` by default, `--weak-cog-categories`) |
| `annotated` | at least one real functional description, or a non-weak COG category |

`putative`/`probable` hits are treated as **annotated** — "putative ABC
transporter" is a real hypothesis, not an HP.

---

## Reproducibility & sensitivity

- Every run emits `*_run_report.json` with the input **SHA-256**, full parameter
  set, package versions, and the anvi'o backend — drop it next to your figures.
- Category thresholds (`--core`, `--soft-core`, `--shell`), the background
  (`--bg-categories`), `--l0`, and `--weak-cog-categories` are all explicit
  flags **made to be varied** in a sensitivity analysis.
- `integrate` reports a length-bias QC (Spearman ρ of IPS vs length) and warns
  if the length term is dominating the ranking.

---

## Development

```bash
pip install -e ".[test]"
pytest -q
```

CI runs the suite on Python 3.9 / 3.11 / 3.12 on every push (see
[`.github/workflows/ci.yml`](.github/workflows/ci.yml)).

---

## License

MIT — see [LICENSE](LICENSE). Part of **bioentro-suite**.
