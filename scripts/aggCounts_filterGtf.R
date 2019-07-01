library(tidyverse)
library(tximport)
source('~/scripts/write_gtf.R')
args=commandArgs(trailingOnly = T)
working_dir=args[1]
gtf_file <- args[2]
qfolder <- args[3]
tx_outfile <- args[4]
filt_gtf_outfile <- args[5]
setwd(working_dir)
gtf <- rtracklayer::readGFF(gtf_file)
t2g <- gtf %>% filter(type=='transcript') %>% select(transcript_id,gene_id)
t2g$gene_id[is.na(t2g$gene_id)] <- t2g$transcript_id[is.na(t2g$gene_id)]# deal with any missing gene names
q_files <- list.files(qfolder,pattern = 'quant.sf', recursive = T,full.names = T )
sample_names <- strsplit(q_files, '/')%>% sapply(function(x)x[[4]])# ***haven't figured out a good way to deal with names, yet so have to adjust every time 
txi <- tximport(q_files, 'salmon', txOut = T, countsFromAbundance = 'lengthScaledTPM')
colnames(txi$counts) <- sample_names
tx_counts <- txi$counts %>% as.data.frame() %>% mutate(transcript_id=rownames(.)) %>% select(transcript_id, everything())
tx_counts_exp <- tx_counts %>% filter(rowSums(.[,-1]) >= (ncol(.)-1) )
gtf_exp <- filter(gtf, transcript_id %in%  tx_counts_exp$transcript_id)
write_tsv(tx_counts_exp, tx_outfile)
write_gtf3(gtf_exp, filt_gtf_outfile)
