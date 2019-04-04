library(tximport)
library(data.table)
library(tidyverse)
# sonnesson 2016 inspired ID of low-usage tx
# http://biorxiv.org/content/early/2015/08/24/025387
# working dir, biowulf2
args=commandArgs(trailingOnly=T)
#args <- c('/data/swamyvs/eyeintegration_splicing/',
#          'results/all_tissues.combined.gtf',
#s          'results/tx_for_removal.txt')


working_dir <- args[1]
gtf_file = args[2]
removal_file = args[3]
setwd(working_dir)
# pull in salmon files
files <-paste0('quant_files/',list.files(path='quant_files',recursive=TRUE,pattern='quant.sf'))
# Gene TX to name conversion
anno <- rtracklayer::readGFF(gtf_file)%>%dplyr::filter(type=='transcript') %>% select(transcript_id, gene_id)
#save(anno,file = 'ref/txdb.Rdata')
# pull counts
txi <- tximport(files, type = "salmon", tx2gene = anno,  txOut = T)
# get counts
tx_c <- data.frame(txi$counts)
# get gene name added
tx_c <- merge(anno,tx_c, by.x='transcript_id',by.y='row.names')
#  sum counts by gene for all samples and calculate tx usage ratio
gene_sums <- summarizeToGene(txi = txi,tx2gene = anno)$counts 
genes <- rownames(gene_sums)
#%>% rowSums %>%as.data.frame 
gene_sums <- data.frame( gene_id = genes,sums= rowSums(gene_sums))

gene_sums_tx <- left_join(tx_c[,1:2], gene_sums, by='gene_id')%>% arrange(gene_id)
tx_c <- dplyr::arrange(tx_c,gene_id)
all_ratios <- tx_c[,-(1:2)]/gene_sums_tx$sums
all_ratios[is.nan(as.matrix(all_ratios))] <- 0
# find number of samples for each transcripts which are < 5% of the total
low_usage <- which(rowSums(all_ratios)<=.05)
print(paste(length(low_usage), 'Transcripts Removed'))
tx_to_remove <- tx_c[low_usage,'transcript_id']
write(tx_to_remove,removal_file,sep='\n')
