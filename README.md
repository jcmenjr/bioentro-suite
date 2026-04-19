![Logo](assets/logo.png)
# bioentro-suite
**Informational Triage Suite for Microbial Genomics**

`bioentro-suite` is a bioinformatics platform designed to prioritize hypothetical proteins (HPs) by quantifying their informational cost and pangenomic dynamics. By applying metrics from Information Theory, the suite filters computational noise and identifies evolutionary innovations with high functional potential in clinical isolates.

---

## The Problem: The "Blind Alley" of Hypothetical Proteins
In bacterial genomics, up to 40% of annotated genes are "hypothetical proteins." Validating all of them experimentally is cost-prohibitive. `bioentro-suite` solves this bottleneck by acting as an **informational triage funnel**, reducing thousands of candidates to a prioritized list based on structural complexity and species-level atypicality.

## Core Workflow
The suite integrates into standard genomic pipelines:
`Genomes (FASTA)` → `Bakta (Annotation)` → `Panaroo (Pangenome)` → **`bioentro-suite`** → `Structural Modeling (AlphaFold) / Lab Validation`.

---

## Main Tools

| Tool | Function |
| :--- | :--- |
| `preparo` | Automates GFF3/FAA preprocessing and generates ready-to-run Panaroo commands. |
| `bioentro` | Computes entropy, divergence, and IPS scores at the individual sequence level. |
| `pangentro` | Integrates informational scores with pangenome architecture (Core, Shell, Cloud). |

---

## Mathematical Foundation: IPS (Informational Priority Score)

The central metric of the suite is the **IPS**, a heuristic designed to identify sequences that represent a significant metabolic and informational investment for the bacteria.

$$IPS = \text{Efficiency} \times \sqrt{JSD} \times \text{Length Penalty}$$

### Component Rationale:
1. **Informational Efficiency ($H_{real} / H_{max}$):** Based on Shannon Entropy. It acts as a proxy for structural complexity. High values indicate rich, non-repetitive sequences with a higher probability of possessing foldable globular domains (optimal for AlphaFold modeling).
2. **Root Jensen-Shannon Divergence ($\sqrt{JSD}$):** Measures the compositional "foreignness" of the sequence by comparing it strictly against the species' **Core Proteome**. This isolates Horizontal Gene Transfer (HGT) events.
3. **Length Penalty ($1 - e^{-L/L_0}$):** - Suppresses scores for short assembly artifacts and annotation errors.
   - Models **metabolic cost**: A large, atypical sequence represents a higher energetic burden, suggesting strong selective pressure for its conservation.

---

## Pangenomic Integration (`pangentro`)

`pangentro` identifies the "sweet spot" of adaptive evolution: **Hypothetical Proteins in the Shell Genome with high IPS.**

- **Core Background Calibration:** Pangenomic IPS calculation uses the Core genome as the informational reference, ensuring that atypicality is measured against the stable identity of the species.
- **Statistical Filtering:** Enables Kruskal-Wallis and Mann-Whitney U tests to confirm if accessory sequences (Shell) carry a significantly higher informational load than the basal machinery (Core).

---

## Installation

```bash
# Clone the repository
git clone https://github.com/jcmenjr/bioentro-suite.git
cd bioentro-suite

# Install in development mode
pip install -e .

```

---

## Repository structure

```
bioentro-suite/
├── assets/
│   ├── logo.png   
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

Developed by Julio Méndez · [jcmenjr@gmail.com](mailto:jcmenjr@gmail.com) (I need help)
