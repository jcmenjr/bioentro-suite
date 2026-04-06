# bioentro suite

**Informational Analysis Suite for Biological Sequences**

A collection of command-line tools that apply information theory to biological sequences — measuring Shannon entropy, Jensen-Shannon divergence, Kolmogorov complexity, and an Informational Priority Score (IPS) to characterize and prioritize hypothetical proteins for experimental validation.

Developed as part of a thesis on the application of information theory to the functional prediction of hypothetical proteins in clinical isolates of *Pseudomonas aeruginosa*.

---

## Tools

| Tool | Description |
|---|---|
| `preparo` | Prepares genome files for the Bakta → Panaroo → bioentro pipeline: detects naming patterns (NCBI, assembler, custom), renames files consistently, verifies Bakta outputs, and generates ready-to-run Panaroo commands |
| `bioentro` | Computes informational metrics across four biological levels: genome, genes (CDS), individual proteins, and whole proteome |
| `netentro` | Builds informational similarity networks, predicts functional class of hypothetical proteins by centroid distance, and validates predictions with leave-one-out cross-validation |

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

After installation, `preparo`, `bioentro`, and `netentro` are available as shell commands.

---

## Full pipeline

```
Raw .fna genomes (NCBI / assembler / own sequencing)
         │
         ▼
   preparo detect   →  sample_map.tsv (review & edit names)
   preparo rename   →  renamed .fna with clean locus tag prefixes
         │
         ▼  (run Bakta manually — see rename output for the exact command)
         │
   preparo check    →  verify Bakta completeness and ID consistency
   preparo panaroo  →  generate ready-to-run Panaroo script
         │
         ▼  (run Panaroo manually)
         │
   bioentro         →  informational metrics (.tsv per biological level)
   netentro network →  similarity network visualization
   netentro predict →  functional class prediction by centroid distance
   netentro validate→  LOO cross-validation of prediction accuracy
```

---

## Quick start

### Step 0 — Prepare genomes with `preparo`

```bash
# Scan your genome directory and generate a naming map
preparo detect -i genomes/ -o sample_map.tsv

# Edit sample_map.tsv — set meaningful names in the 'new_name' column
# Rules: unique, ≤12 chars, alphanumeric + underscore, no spaces

# Rename files (also rewrites FASTA headers for downstream traceability)
preparo rename -i genomes/ -m sample_map.tsv -o renamed/

# After running Bakta, verify all outputs are complete and IDs match
preparo check -i bakta_outputs/

# Generate a ready-to-run Panaroo command script
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

### Step 2 — Analyze with `netentro`

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
| `detect` | Before Bakta | Scans genome files, auto-detects names from FASTA headers, generates `sample_map.tsv` draft |
| `rename` | Before Bakta | Renames `.fna` files and rewrites FASTA headers using the edited map |
| `check` | After Bakta | Verifies `.gff3`, `.faa`, and other expected files exist and IDs are consistent |
| `panaroo` | After Bakta | Collects all `.gff3` files and writes a ready-to-run `panaroo` shell script |

`detect` auto-detects names from FASTA headers for NCBI downloads (`[organism=...]`, `[strain=...]`) and falls back to the short accession when headers have no metadata. Assembler outputs (SPAdes, Flye, Unicycler) get sequential placeholders that you edit manually.

---

## netentro subcommands

| Subcommand | Description |
|---|---|
| `network` | Builds and saves an informational similarity network image |
| `predict` | Ranks functional classes by centroid distance; reports confidence relative to intra-class dispersion |
| `validate` | Leave-one-out cross-validation: predicts each annotated protein and reports per-class precision, recall, and F1 |

**Confidence formula (predict mode):**

$$\text{Confidence}(C_k) = \max\left(0,\ 1 - \frac{d(t,\ C_k)}{\sigma_{C_k}}\right)$$

Where $\sigma_{C_k}$ is the mean distance of class members to their centroid. A confidence of 0.5 means the target is as far from the centroid as the average class member. A confidence of 0.0 means the target lies outside the typical class radius — the prediction is unreliable regardless of rank.

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

> `IPS` is only computed for `genes` and `protein` modes, where a background comparison is meaningful.

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
│   └── preparo.py        # Genome preparation: detect, rename, check, panaroo
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

> Méndez, J. (2026). *bioentro suite — Informational Analysis Suite for Biological Sequences*. Centro de Investigación en Alimentación y Desarrollo (CIAD).

```bibtex
@software{mendez2026bioentro,
  author  = {Méndez, Julio},
  title   = {bioentro suite — Informational Analysis Suite for Biological Sequences},
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