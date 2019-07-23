library(tidyverse)
source('~/scripts/write_gtf.R')
args=commandArgs(trailingOnly = T)
working_dir=args[1]
infile=args[2]
ref_gtf_file=args[3]
outfile=args[4]
setwd(working_dir)
"The gtf outputed by GFF compare has a few inconsistancies
- the xloc gene_ids don't line up exactly with gene names, with multi mapping on both sides
- gene names and other infor are only provided for the transcript, not exons which is real annoying
- this script is to fix those problems, as well as add ENS gene ids back using the comp gtf.
"

gn2gi <- rtracklayer::readGFF(ref_gtf_file) %>% select(gene_name, gene_id) %>% distinct %>% filter(!duplicated(.[,'gene_name']))

in_gtf <- rtracklayer::readGFF(infile) %>% mutate(gffc_loc=gene_id, gene_id=transcript_id)
tx2g <- in_gtf %>% 
    filter(type == 'transcript') %>% 
    select(transcript_id, gene_name) %>% 
    distinct %>% 
    mutate(gene_name=replace(gene_name, is.na(gene_name), transcript_id[is.na(gene_name)])) %>% distinct %>% 
    left_join(gn2gi) %>% mutate(gene_id=replace(gene_id, is.na(gene_id), transcript_id[is.na(gene_id)]))


out_gtf <- in_gtf %>% select(-gene_name, -gene_id) %>% left_join( tx2g)

write_gtf3(out_gtf,outfile)
