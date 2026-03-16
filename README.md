# bioentro 🧬
**Genomic and Proteomic Information Complexity Analysis**

`Bioentro` is a bioinformatics tool designed to quantify informational complexity and evolutionary optimization in biological prokaryotic DNA sequences. By applying the principles of **Shannon’s Information Theory** and **Kullback–Leibler Divergence**, it characterizes the "Dark Proteome" and accessory genomes through k-mer distribution analysis.


---

## ✨ Features
- **Protein Mode:** Amino acid sequence analysis to predict structural stability and identify potentially disordered regions.
- **DNA Mode:** Calculation of genomic metrics integrating sequence data (.fna) and annotation metadata (.gff).
- **Statistical Reliability:** Automatic "Traffic Light" classification (`Reliable` / `Low_Sample`) based on sequence length relative to the k-mer search space ($20^k / 2$).
- **Candidate Prioritization:** Designed as a high-throughput pre-filter for large-scale screening prior to AlphaFold structural prediction or wet-lab validation.

## 🚀 Installation

`Bioentro` requires Python 3.7+. Since the project includes a `pyproject.toml` file, dependencies (such as `Biopython`) will be installed automatically.

```bash
# Clone the repository
git clone https://github.com/jcmenjr/bioentro.git
cd bioentro

# Install the package
pip install .

```

## 📊 Usage

Once installed, you can run `bioentro` directly from your terminal

1. For proteins (ideal for Hypothetical proteins)
   
Ideal for prioritizing hypothetical proteins. It is recommended to use k=1 for very short sequences and k=2 for structural motifs.

```bash

bioentro -a data/proteins.faa -k 1 -o results.tsv


```

2. Nucleotide-based genomic analysis (big picture)

This mode evaluates the **global informational landscape** of the entire genome or large contigs. It provides a single summary of informational efficiency and entropy for the whole sequence. This is ideal for comparative genomics between different species or to establish a baseline of "genomic texture" across different isolates.


```bash

bioentro -f data/genome.fna -g data/annotation.gff -k 6 -o results_dna.tsv

```

## 🧠 Interpretation of metrics

Bioentro quantifies the "informational signature" of biological sequences. Here is how to interpret the results:

| Metric | Definition | Biological Interpretation |
| :--- | :--- | :--- |
| **$H_{real}$** | **Shannon Entropy** | Measures the lexical variety of k-mers. **Low values** often indicate repetitive regions or low-complexity sequences. **High values** suggest complex, diverse sequences. |
| **$KL\_Div$** | **Kullback-Leibler Divergence** | Measures the "distance" between the observed sequence and a random uniform distribution. High KL suggests a **strong evolutionary signature** (the sequence is far from random). |
| **Efficiency (%)** | **Informational Efficiency** | The ratio between $H_{real}$ and the theoretical maximum ($H_{max}$). This helps identify an **optimization gradient** |
| **Status** | **Data Reliability** | `Reliable`: The sequence is long enough for the chosen $k$. <br> `Low_Sample`: The sequence is too short; results should be interpreted with caution. |

### 💡 Quick Guide for Candidate Prioritization:
* **High Efficiency + High KL:** Strongest candidates for functional, well-folded proteins (enzymes, stable structures).
* **Low Efficiency + High KL:** Potential factors of virulence or structural proteins with repetitive motifs (adhesins, etc.).
* **Low KL:** Sequences that resemble genomic "noise" or recently acquired, non-optimized elements.

Note: The biological interpretations provided by Bioentro are based on information theory models and evolutionary hypotheses. While high informational efficiency is a strong proxy for structural stability and selective pressure, these results should be treated as predictions. Users are encouraged to validate high-priority candidates using structural modeling (e.g., AlphaFold), conservation analysis, or experimental assays.

## 🎓 Citation

If you use **bioentro** in your research, please cite it as follows:

**APA Style:**
> Méndez, J. (2026). bioentro: Genomic and Proteomic Information Complexity Analysis. GitHub repository. https://github.com/jcmenjr/bioentro

**BibTeX:**
```bibtex
@software{mendez2026bioentro,
  author = {Méndez, Julio},
  title = {Bioentro: Genomic and Proteomic Information Complexity Analysis},
  url = {[https://github.com/jcmenjr/bioentro](https://github.com/jcmenjr/bioentro)},
  version = {1.0.0},
  year = {2026}
}
```
---

## 📜 License

This project is licensed under the MIT License. You are free to use, modify, and distribute it, provided that proper credit is given to the original author.

Developed by: Julio Méndez

Contact: jcmenjr@gmail.com
