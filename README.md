# bioentro suite
 
**Informational Analysis Suite for Biological Sequences**
 
A collection of command-line tools that apply information theory to biological sequences — measuring Shannon entropy, Jensen-Shannon divergence, Kolmogorov complexity, and an Informational Priority Score (IPS) to characterize and prioritize hypothetical proteins for experimental validation.
 
Developed as part of a thesis on the application of information theory to the functional prediction of hypothetical proteins in clinical isolates of *Pseudomonas aeruginosa*.
 
---
 
## Tools
 
| Tool | Description |
|---|---|
| `bioentro` | Computes informational metrics across four biological levels: genome, genes (CDS), individual proteins, and whole proteome |
| `netentro` | Builds an informational similarity network to infer putative function of hypothetical proteins by proximity to annotated ones |
 
---
 
## Installation
 
### From GitHub (recommended for development)
 
```bash
git clone https://github.com/<your-username>/bioentro-suite.git
cd bioentro-suite
pip install -e .
```
 
The `-e` flag installs in *editable* mode: changes to the source files take effect immediately without reinstalling.
 
### From GitHub (stable, no source needed)
 
```bash
pip install git+https://github.com/<your-username>/bioentro-suite.git
```
 
After installation, both `bioentro` and `netentro` are available as shell commands.
 
---
 
## Quick start
 
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
 
### Step 2 — Build a similarity network with `netentro`
 
```bash
# Basic usage: focus on one hypothetical protein
netentro -i protein_metrics.tsv -t LOCUS_TAG -o network.png
 
# More neighbors, higher resolution
netentro -i protein_metrics.tsv -t LOCUS_TAG -k 30 --top-ips 15 -o network.png
 
# List available protein IDs first
netentro -i protein_metrics.tsv --list-ids
```
 
---
 
## Background distributions
 
The choice of background is biologically motivated:
 
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
│   ├── bioentro.py       # Core informational metrics tool
│   └── netentro.py       # Similarity network tool
├── tests/
│   ├── __init__.py
│   ├── test_bioentro.py  # Unit tests for metric functions
│   └── test_netentro.py  # Unit tests for network functions
├── pyproject.toml        # Package configuration and dependencies
├── LICENSE
└── README.md
```
 
---
 
## Citation
 
If you use bioentro suite in your research, please cite:
 
> [Méndez,J.], [2026]. *Application of Information Theory to the Functional Prediction of Hypothetical Proteins in Clinical Isolates of Pseudomonas aeruginosa*. [Center of Research for Food and Development].
 
 ```bibtex
@software{mendez2026bioentro-suite,
  author = {Méndez, Julio},
  title = {bioentro suite — Informational Analysis Suite for Biological Sequences},
  url = {[https://github.com/jcmenjr/bioentro](https://github.com/jcmenjr/bioentro-suite)},
  version = {0.1.0}},
  year = {2026},
  orcid = {0009-0007-2468-6261}
}
```
---
 
## License
 
MIT — see [LICENSE](LICENSE) for details.

This project is licensed under the MIT License. You are free to use, modify, and distribute it, provided that proper credit is given to the original author.

Developed by: Julio Méndez

Contact: jcmenjr@gmail.com