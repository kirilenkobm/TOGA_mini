# TOGA-mini

Lightweight version of [TOGA](https://github.com/hillerlab/TOGA), 
published in [Science (2023)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10193443/), 
designed for rapid orthologous locus prediction from alignment chains.

Unlike full TOGA, TOGA-mini focuses solely on predicting orthologous loci for any gene type, including protein-coding genes, lncRNAs, or custom annotations. This leaves downstream analysis (e.g., realignment, ORF prediction) flexible and customizable.

**Note: TOGA-mini is not versioned in sync with full TOGA.**

## 🚀 Features

- Fast orthologous locus prediction via genomic alignment chains.
- Supports protein-coding and non-coding genes (e.g., lncRNAs).
- Suitable for rapid screening or as input for custom post-processing pipelines.
- Slimmed-down, standalone tool.

## 📦 Installation

```bash
conda env create -f environment.yaml
./configure.sh
```

## 📂 Usage

```
toga_mini.py [-h] chain_file transcript_file isoforms_file target_2bit query_2bit out_orthologous_regions_mapping out_classification_table

positional arguments:
  chain_file            Path to genome alignment file in chain format
  transcript_file       Path to transcripts bed12 file
  isoforms_file         Isoforms mapping
  target_2bit           Path to reference 2bit file
  query_2bit            Path to query 2bit file
  out_orthologous_regions_mapping
                        Output with orthologous regions
  out_classification_table
                        Classification table with predictions

options:
  -h, --help            show this help message and exit
 ```

## 🔗 References

- Original TOGA repo: https://github.com/hillerlab/TOGA
- TOGA publication: [Science, 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10193443/)
