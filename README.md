## Snakemake Pipeline for Constructing *de novo* transcriptomes

This Snakemake workflow builds *de novo* transcriptomes from short read RNA-seq samples. If you want to run this pipeline, you will need to change the following files to match paths on your system
- config.yaml
- files.yaml 
- ref/tissues.txt
- ref/subtissues.txt

Additionally, a sampleTable is requried, and must be formatted in the same way the sampleTable*.tsv files are in this repo

## DAG

![alt text](https://github.com/vinay-swamy/eyeintegration_splicing/blob/master/dag.svg)
