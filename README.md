# TOGA-mini

Lightweight version of [TOGA](https://github.com/hillerlab/TOGA), 
published in [Science (2023)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10193443/), 
designed for rapid orthologous locus prediction from alignment chains.

Unlike full TOGA, TOGA-mini focuses solely on predicting orthologous loci for any gene type, including protein-coding genes, lncRNAs, or custom annotations. This leaves downstream analysis (e.g., realignment, ORF prediction) flexible and customizable.

## 🚀 Features

- Fast orthologous locus prediction via genomic alignment chains.
- Supports protein-coding and non-coding genes (e.g., lncRNAs).
- Suitable for rapid screening or as input for custom post-processing pipelines.
- Slimmed-down, standalone tool.

## 📦 Installation

Train classification models:

```
./chain_class_models/train_toga_chain_class_model.py
```

Build C modules:

```bash
./build_c.sh 
```

A unified `configure.sh` installation script is planned.

## 📂 Usage

To be added. For now, refer to the original TOGA documentation for general guidance.

## ⚠️ Disclaimer

The codebase retains historical artifacts from various TOGA development stages and may be messy in parts. Cleanup and modernization are ongoing.

## 🔗 References

- Original TOGA repo: https://github.com/hillerlab/TOGA
- TOGA publication: [Science, 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10193443/)
