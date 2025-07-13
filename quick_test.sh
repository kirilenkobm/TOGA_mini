#!/usr/bin/env bash
# TODO: proper testing suite
rm -rf hg38_vs_mm39_lncRNA_toga
rm -rf nextflow_stuff

 ./toga_mini.py ../CURIA/datum/hg38.mm39.allfilled.chain ../CURIA/toga_input/hg38.lncRNA.bed ../CURIA/datum/hg38.2bit ../CURIA/datum/mm39.2bit -i ../CURIA/toga_input/hg38.isoforms.tsv --project_dir hg38_vs_mm39_lncRNA_toga --nextflow_dir nextflow_stuff
