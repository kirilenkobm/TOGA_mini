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

## TODO

- Isoforms graph.

## 📦 Installation

```bash
conda env create -f environment.yaml
./configure.sh
```

## 📂 Usage

```
usage: toga_mini.py [-h] chain_file transcript_file isoforms_file reference_chrom_sizes out_orthologous_regions_mapping out_classification_table

positional arguments:
  chain_file            Path to genome alignment file in chain format
  transcript_file       Path to transcripts bed12 file
  isoforms_file         Isoforms mapping TSV with columns: gene_id, transcript_id
  reference_chrom_sizes Path to reference chromosome sizes file (tab-separated: chrom_name	size)
  out_orthologous_regions_mapping
                        Output with orthologous regions
  out_classification_table
                        Classification table with predictions

options:
  -h, --help            show this help message and exit
 ```

### Isoforms file format

`isoforms_file` must be a TSV with a header row that includes `gene_id` and `transcript_id` columns (see `test_input/hg38.chr21.isoforms.tsv` for an example).

## 🔗 References

- Original TOGA repo: https://github.com/hillerlab/TOGA
- TOGA publication: [Science, 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10193443/)
