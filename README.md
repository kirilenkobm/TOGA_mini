# TOGA-mini

Please see the [original repository](https://github.com/hillerlab/TOGA).

Science [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10193443/).

## Motivation

Light-version of the TOGA pipeline that only predicts orthologous loci for a set of reference genes.
Additionally, it can work not only on protein-coding genes, but on any kind of them (like lnc-RNA).
Provides material for the further post-processing.

Codebase may be a bit messy because different parts of toga were written in different epochs.

## Installation

TODO: a single `configure.sh` script

Create classification models:

```
./chain_class_models/train_toga_chain_class_model.py
```

Build C modules:

```bash
./build_c.sh 
```

