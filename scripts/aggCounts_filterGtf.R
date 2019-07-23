library(tidyverse)
library(tximport)
library(matrixStats)
source('~/scripts/write_gtf.R')
#args <- c('~/NIH/dev_eyeintegration_splicing/', 'data/gtfs/raw_tissue_gtfs/RPE_Fetal.Tissue_st.gtf', 'data/rawST_tx_quant_files/RPE_Fetal.Tissue')
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
proc_boostrap <- function(file){
    BS_raw <- read_tsv(file, col_names = F)
    name <- str_split(file,'/')[[1]][4]
    res <- tibble(transcript_id=BS_raw$X1, !!name := BS_raw[,-1] %>% as.matrix %>% {matrixStats::rowVars(.)} )
    return(res)
}
which_var_fail <- function(all_bs){
    novel_bs <- filter(all_bs, grepl('MSTRG', transcript_id))
    ref_bs <- filter(all_bs, !grepl('MSTRG', transcript_id) )
    novel_medvar <- novel_bs[,-1] %>% as.matrix() %>% rowMedians()
    ref_medvar <- ref_bs[,-1] %>% as.matrix %>% rowMedians()
    refvar_95 <- quantile(ref_medvar, .95)
    print(refvar_95)
    failed_novel <- filter(novel_bs, novel_medvar > refvar_95) 
    print( paste(nrow(failed_novel)/nrow(novel_bs), ' percent novel tanscripts removed'))
    return(failed_novel)
}

bs_files <- list.files(qfolder, pattern='quant_bootstraps.tsv', recursive = T, full.names = T)
quant_bs <- lapply(bs_files, proc_boostrap) %>% reduce(left_join) %>% filter(transcript_id %in% gtf_exp$transcript_id)
high_var_tx <- which_var_fail(quant_bs)
gtf_exp_var_pass <- filter(gtf_exp, !transcript_id %in% high_var_tx$transcript_id)
tx_counts_exp_var_pass <- filter(tx_counts_exp, !transcript_id %in% high_var_tx$transcript_id)
write_tsv(tx_counts_exp_var_pass, tx_outfile)
write_gtf3(gtf_exp_var_pass, filt_gtf_outfile)














