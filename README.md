# bioentro suite

**Informational Zoom Suite for Biological Sequences**

A collection of command-line tools that apply information theory to biological sequences — measuring Shannon entropy, Jensen-Shannon divergence, Kolmogorov complexity, and an Informational Priority Score (IPS) to characterize and prioritize hypothetical proteins for experimental validation.

Developed as part of a thesis on the application of information theory to the functional prediction of hypothetical proteins in clinical isolates of *Pseudomonas aeruginosa*, with generalization to any bacterial system.

---

## Tools

| Tool | Description |
|---|---|
| `preparo` | Prepares genome files for the Bakta → Panaroo → bioentro pipeline: detects naming patterns (NCBI, assembler, custom), renames files, verifies Bakta outputs, and generates ready-to-run Panaroo commands |
| `bioentro` | Computes informational metrics across four biological levels: genome, genes (CDS), individual proteins, and whole proteome |
| `netentro` | Builds informational similarity networks, predicts functional class of hypothetical proteins by centroid distance, and validates predictions with leave-one-out cross-validation |
| `pangentro` | Integrates bioentro metrics with Panaroo pangenome output to analyze IPS distribution across core, shell, and cloud genomes |

---

## Installation

### From GitHub (recommended for development)

```bash
git clone https://github.com/jcmenjr/bioentro-suite.git
cd bioentro-suite
pip install -e .
```

The `-e` flag installs in *editable* mode: changes to the source files take effect immediately without reinstalling.

### From GitHub (stable, no source needed)

```bash
pip install git+https://github.com/jcmenjr/bioentro-suite.git
```

After installation, `preparo`, `bioentro`, `netentro`, and `pangentro` are available as shell commands.

---

## Full pipeline

```
Raw .fna genomes (NCBI / assembler / own sequencing)
         │
         ▼
   preparo detect   →  sample_map.tsv  (review & edit names)
   preparo rename   →  renamed .fna with clean locus tag prefixes
         │
         ▼  Bakta  (manual — see rename output for the exact command)
         │
   preparo check    →  verify Bakta completeness and ID consistency
   preparo panaroo  →  generate ready-to-run Panaroo script
         │
         ▼  Panaroo  (manual)
         │
   pangentro run-bioentro  →  run bioentro on all isolates automatically
   pangentro integrate     →  cross bioentro results with pangenome matrix
   pangentro compare       →  statistics + figure: IPS by core/shell/cloud
         │
   netentro network   →  similarity network visualization
   netentro predict   →  functional class prediction by centroid distance
   netentro validate  →  LOO cross-validation of prediction accuracy
```

---

## Quick start

### Step 0 — Prepare genomes with `preparo`

```bash
# Scan genome directory and generate naming map
preparo detect -i genomes/ -o sample_map.tsv

# Edit sample_map.tsv — set meaningful names in the 'new_name' column
# Rules: unique, ≤12 chars, alphanumeric + underscore, no spaces

# Rename files (rewrites FASTA headers for downstream traceability)
preparo rename -i genomes/ -m sample_map.tsv -o renamed/

# After running Bakta, verify all outputs are complete and IDs match
preparo check -i bakta_outputs/

# Generate ready-to-run Panaroo command
preparo panaroo -i bakta_outputs/ -o run_panaroo.sh -t 16
bash run_panaroo.sh
```

### Step 1 — Compute informational metrics with `bioentro`

```bash
# Whole genome vs uniform background
bioentro -i genome.fna -g genome.gff3 -m genome -o genome_metrics.tsv

# Per-CDS vs full-genome background (detects HGT, pathogenicity islands)
bioentro -i genome.fna -g genome.gff3 -m genes -o genes_metrics.tsv

# Per-protein vs full-proteome background (IPS for prioritization)
bioentro -i proteins.faa -m protein -o protein_metrics.tsv

# Whole proteome vs uniform background
bioentro -i proteins.faa -m proteome-global -o proteome_metrics.tsv
```

### Step 2 — Pangenomic analysis with `pangentro`

```bash
# Run bioentro on all isolates automatically
pangentro run-bioentro -f bakta_outputs/ -o bioentro_results/

# Integrate with Panaroo output (both IPS modes)
pangentro integrate \
    -p panaroo_output/gene_presence_absence.csv \
    -b bioentro_results/ \
    -r panaroo_output/pan_genome_reference.fa \
    --ips-mode both \
    -o pangenome_metrics.tsv

# Statistical comparison and figure
pangentro compare -i pangenome_metrics.tsv -o comparison/
```

### Step 3 — Network analysis with `netentro`

```bash
# Visualize informational similarity network for a target protein
netentro network -i protein_metrics.tsv -t LOCUS_TAG -o network.png

# More neighbors, higher resolution
netentro network -i protein_metrics.tsv -t LOCUS_TAG -k 30 --top-ips 15 -o network.png

# List all available protein IDs
netentro network -i protein_metrics.tsv --list-ids

# Predict functional class by informational centroid distance
netentro predict -i protein_metrics.tsv -t LOCUS_TAG -o prediction.tsv

# Validate prediction accuracy with leave-one-out cross-validation
netentro validate -i protein_metrics.tsv -o validation.tsv
```

---

## preparo subcommands

| Subcommand | When to use | Description |
|---|---|---|
| `detect` | Before Bakta | Scans genome files, auto-detects names from FASTA headers (`[organism=...]`, `[strain=...]`), generates `sample_map.tsv` draft |
| `rename` | Before Bakta | Renames `.fna` files and rewrites FASTA headers using the edited map |
| `check` | After Bakta | Verifies `.gff3`, `.faa`, and other expected files exist; checks ID prefix consistency |
| `panaroo` | After Bakta | Collects all `.gff3` files and writes a ready-to-run Panaroo shell script |

> Note: Panaroo may fail on GFF3 files containing selenocysteine codons (`transl_except=Sec`). Workaround: `grep -v "transl_except" input.gff3 > clean.gff3`. This is a known upstream bug (reported to Panaroo developers).

---

## netentro subcommands

| Subcommand | Description |
|---|---|
| `network` | Builds and saves an informational similarity network image |
| `predict` | Ranks functional classes by centroid distance; confidence is relative to intra-class dispersion |
| `validate` | Leave-one-out cross-validation: predicts each annotated protein and reports per-class precision, recall, and F1 |

**Confidence formula (predict mode):**

$$\text{Confidence}(C_k) = \max\left(0,\ 1 - \frac{d(t,\ C_k)}{\sigma_{C_k}}\right)$$

Where $\sigma_{C_k}$ is the mean distance of class members to their centroid. Confidence = 0.5 means the target is as far from the centroid as the average class member. Confidence = 0.0 means the target lies outside the typical class radius — the prediction is unreliable regardless of rank.

> `netentro network` and `netentro predict` are exploratory tools. The network visualizes the informational space of the proteome; it does not assert functional equivalence. Use `netentro validate` to quantify method accuracy for your specific dataset before drawing biological conclusions.

---

## pangentro subcommands

| Subcommand | Description |
|---|---|
| `run-bioentro` | Automates bioentro protein mode on all `.faa` files in a directory |
| `integrate` | Crosses per-isolate bioentro TSVs with Panaroo `gene_presence_absence.csv`; supports two IPS modes |
| `compare` | Kruskal-Wallis + pairwise Mann-Whitney U tests; generates publication figure |

**IPS modes in `integrate`:**

| Mode | Background | Interpretation |
|---|---|---|
| `individual` | Each protein's own isolate proteome | Atypicality within the organism's own compositional context |
| `pangenomic` | Concatenated pan-proteome (all cluster representatives) | Atypicality within the species-level compositional space |
| `both` | Both (recommended) | Enables cross-mode comparison for local vs species-level signals |

**Pangenomic categories (Panaroo defaults):**

| Category | Frequency | Biological interpretation |
|---|---|---|
| core | ≥ 99% of genomes | Essential, conserved — housekeeping functions |
| soft_core | 95–99% | Near-universal — mostly conserved |
| shell | 15–95% | Accessory — niche adaptation, HGT candidates |
| cloud | < 15% | Rare — recent acquisition, possible mobile elements |

---

## Background distributions

| Mode | Background | Biological rationale |
|---|---|---|
| `genome` | Uniform (theoretical max) | Measures absolute informational efficiency of the genome |
| `genes` | Full genome k-mer distribution | Detects CDS outliers: HGT candidates, pathogenicity islands |
| `protein` | Full proteome k-mer distribution | Detects proteins with atypical amino acid composition |
| `proteome-global` | Uniform (theoretical max) | Measures absolute informational efficiency of the proteome |

---

## Metrics

| Metric | Description | Range |
|---|---|---|
| `H_max` | Maximum theoretical Shannon entropy | bits |
| `H_real` | Observed Shannon entropy of the sequence | bits |
| `Efficiency` | H_real / H_max — compositional complexity | [0, 1] |
| `KL_bg` | KL(sequence ∥ background) — directional, asymmetric | ≥ 0 bits |
| `JSD_bg` | JSD(sequence ∥ background) — symmetric, bounded | [0, 1] |
| `sqrt_JSD` | √JSD — proper metric distance | [0, 1] |
| `Kolmogorov` | zlib compression ratio — long-range complexity proxy | [0, 1] |
| `IPS` | Efficiency × √JSD × length_penalty — priority score | [0, 1] |

> `IPS` identifies compositionally atypical sequences as candidates for further investigation — not as predictions of biological relevance. Relevance interpretation depends on the research question and biological context.

---

## Known limitations

- `netentro predict` accuracy varies by functional class. Classes with strong compositional signatures (membrane proteins, secreted proteins) predict more reliably than classes defined by positional motifs (polymerases, regulators). Always run `netentro validate` to quantify accuracy for your specific dataset.
- `IPS_pan` values (pangenomic mode) are lower in absolute terms than individual-mode IPS because the pan-proteome background is more diverse. Statistical differences between categories are meaningful but absolute values should not be compared across datasets with different pangenome sizes.
- Selenoprotein genes (`transl_except=Sec`) cause Panaroo to crash even with `--remove-invalid-genes`. Workaround: remove with `grep -v "transl_except"`. Document excluded genes as a methodological note.

---

## Requirements

- Python ≥ 3.9
- biopython ≥ 1.79
- pandas ≥ 1.5
- numpy ≥ 1.23
- networkx ≥ 3.0
- matplotlib ≥ 3.6
- scipy ≥ 1.9
- scikit-learn ≥ 1.1

All dependencies are installed automatically by `pip install`.

---

## Development

```bash
# Install with development extras
pip install -e ".[dev]"

# Run tests
pytest

# Lint and format
ruff check .
ruff format .

# Type check
mypy bioentro/
```

---

## Repository structure

```
bioentro-suite/
├── bioentro/
│   ├── __init__.py       # Package metadata and public API
│   ├── bioentro.py       # Informational metrics (genome, genes, protein, proteome)
│   ├── netentro.py       # Similarity network, functional prediction, LOO validation
│   ├── preparo.py        # Genome preparation: detect, rename, check, panaroo
│   └── pangentro.py      # Pangenomic integration: run-bioentro, integrate, compare
├── conda.recipe/
│   └── meta.yaml         # Conda/Bioconda recipe
├── workflow/
│   ├── rules/
│   │   └── bioentro.smk  # Snakemake rules for pipeline integration
│   └── config.yaml       # Example Snakemake configuration
├── tests/
│   ├── __init__.py
│   ├── test_bioentro.py  # Unit tests for metric functions
│   └── test_netentro.py  # Unit tests for network functions
├── Snakefile             # Example Snakemake pipeline
├── pyproject.toml        # Package configuration and dependencies
├── LICENSE
└── README.md
```

---

## Citation

If you use bioentro suite in your research, please cite:

> Méndez, J. (2026). *bioentro suite — Informational Zoom Suite for Biological Sequences*. Centro de Investigación en Alimentación y Desarrollo (CIAD).

```bibtex
@software{mendez2026bioentro,
  author  = {Méndez, Julio},
  title   = {bioentro suite — Informational Zoom Suite for Biological Sequences},
  url     = {https://github.com/jcmenjr/bioentro-suite},
  version = {0.1.0},
  year    = {2026},
  orcid   = {0009-0007-2468-6261}
}
```

---

## License

MIT — see [LICENSE](LICENSE) for details.

Developed by Julio Méndez · [jcmenjr@gmail.com](mailto:jcmenjr@gmail.com)