library(tximport)
library(data.table)
library(tidyverse)
# sonnesson 2016 inspired ID of low-usage tx
# http://biorxiv.org/content/early/2015/08/24/025387
# working dir, biowulf2
args=commandArgs(trailingOnly=T)
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
gene_sums <- summarizeToGene(txi = txi,tx2gene = anno)$counts[,samples_to_keep]%>%rowSums
gene_sums_tx <- merge(tx_c[,1:3], as.data.frame(gene_sums), by.x='gene_name', all.x=T, by.y='row.names' )
tx_c <- dplyr::arrange(tx_c,gene_name)
all_ratios <- tx_c[,-(1:3)]/gene_sums_tx$gene_sums
all_ratios[is.nan(as.matrix(all_ratios))] <- 0
# find number of samples for each transcripts which are < 5% of the total
low_usage <- which(rowSums(all_ratios)<=.05)
print(paste(length(low_usage), 'Transcripts Removed'))
tx_to_remove <- tx_c[low_usage,'transcript_id']
write(tx_to_remove,removal_file,sep='\n')
