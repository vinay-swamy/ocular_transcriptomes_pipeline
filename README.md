## Snakemake Pipeline for Constructing *de novo* transcriptomes

This Snakemake workflow builds *de novo* transcriptomes from short read RNA-seq samples. If you want to run this pipeline, you will need to change the following files to match paths on your system
- config.yaml
- files.yaml 
- ref/tissues.txt
- ref/subtissues.txt

Additionally, a sampleTable is requried, and must be formatted in the same way the sampleTable*.tsv files are in this repo

## DAG

![alt text](https://github.com/vinay-swamy/eyeintegration_splicing/blob/master/dag.svg)

## Workflow

### External libraries/tools you will need to download

#### Python

- pandas 
- pyhvgs
- pyfaidx
- Bio

#### R

- argparse
- data.table
- DBI
- glue
- matrixStats
- parallel
- RBedtools
- tidyverse
- tximport
- yaml

#### Software

- agat

### high level overview

The reference annotation has been hardcoded(links to gencode website) in `config.yaml`, and will be downloaded as part of the rule `downloadAnnotation` . All software is versioned via biowulf module system; there is a chance that the specific versions software will no longer be available. If this is the case, modify the relevant tool version in the `config` to the closest version.

Here are the major steps for the pipeline:

- download and pre-process annotation files
- Align fastqs and build initial, sample-level transcriptomes 
- merge transcriptomes to tissue level, then merge all tissue to full gtf
- calculate and cleanup ORFs for novel transcripts 
- create full annotated master GTF
- predict effect of 



### Alignment
- **Do not re-run alignment unless you absolutely have to**
- The pipeline requires on-disk fastq files to run, set by `fastq_path` in `config.yaml`. Currently its in `/data/OGVFB_BG/EiaD_2019_05/`
- the BAM files from the last run are stored in `/data/swamyvs/DNTX_STARbams`


