library(tidyverse)
library(tximport)
args=commandArgs(trailingOnly = T)
working_dir=args[1]
gtf_file <- args[2]
qfolder= args[3]
tx_outfile <- args[4]
gene_outfile <- args[5]
setwd(working_dir)

t2g <- rtracklayer::readGFF(gtf_file) %>% filter(type=='transcript') %>% select(transcript_id,gene_id)
t2g$gene_id[is.na(t2g$gene_id)] <- t2g$transcript_id[is.na(t2g$gene_id)]
q_files <- list.files(qfolder,pattern = 'quant.sf', recursive = T,full.names = T )
sample_names <- strsplit(q_files, '/')%>% sapply(function(x)x[[2]])
txi <- tximport(q_files, 'salmon', txOut = T, countsFromAbundance = 'lengthScaledTPM')
colnames(txi$counts) <- sample_names
tx_counts <- txi$counts %>% as.data.frame() %>% mutate(transcript_id=rownames(.)) %>% select(transcript_id, everything())

txi <- summarizeToGene(txi = txi, tx2gene = t2g,countsFromAbundance = 'lengthScaledTPM')
colnames(txi$counts) <- sample_names
gene_counts <- txi$counts %>% as.data.frame %>% mutate(gene_id=rownames(.)) %>% select(gene_id, everything())
save(tx_counts,file=tx_outfile )
save(gene_counts,file = gene_outfile)
