#!/bin/bash
mkdir -p ref/snps/
mkdir -p ref/phylop_20/
rsync -avz rsync://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/BED/ ref/snps
rsync -avz rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP20way/ ref/phylop_20/
