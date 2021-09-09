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
- predict effect of variants based on new annotation
- calculate transcript coverage 
- quantify with salmon 
- annotate ORFs with potenial function
- create sqlite db for use with shiny app

Note that several of these steps will run concurrently. The most time consuming steps are alignment and building the initial transcriptomes; these should not be re-run if possible

See the Snakefile and scripts for more in depth comments on what to do 

### download and pre-process annotation
(downloadAnnotation, liftover_intronic_variants, clean_phylop_and_snps)
- **Do not re-run unless you absolutely have to**
- we need to pull annotation from a bunch of different places, so there's a lot of stuff going on 
- some of the annotation files will need to be transformed a little to get it into the right format(see rules for more specific info)
-

### Alignment
- **Do not re-run unless you absolutely have to**
- The pipeline requires on-disk fastq files to run, set by `fastq_path` in `config.yaml`. Currently its in `/data/OGVFB_BG/EiaD_2019_05/`
- the BAM files from the last run are stored in `/data/swamyvs/DNTX_STARbams` 
- sometimes STAR hangs and the alignment will fail, so might have to try it a couple times 
- STAR does have its own option to sort output, but it was causing a lot of problems so I made sorting its own rule

### transcriptome construction
- **Do not re-run unless you absolutely have to**
- build on a per-sample wildcard level, with default parameters
- have extensively tested out some of the other modes, this is the best way to do it

## Important files
- "trackfile" - output of gffcompares, has the prefix `.tracking`. This is often refered to as `conv_tab` in scripts. Its outputed by `gffcompare`.  gffcompare is used to compare transcripts across multiple gtf files. this file contains the union of all transcripts being compared, with a unique id for each, and which others gtf(s) have transcripts that match it

