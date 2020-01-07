library(tidyverse)
library(matrixStats)
library(parallel)
args <- commandArgs(trailingOnly = T)
#save(args, file='testing/fixgtfargs.Rdata')
wd <- args[1]
path_to_quant_files <- args[2]
gtf_file<- args[3]
out_gtf_file <- args[4]
setwd(wd)
source('~/scripts/write_gtf.R')
source('~/scripts/read_salmon.R')

#Read in salmon bootstrap variance and identify transcripts to remove.

proc_boostraps <- function(file){
    BS_raw <- suppressMessages(read_tsv(file, col_names = F) )
    name_idx <- str_split(file, '/')[[1]] %>% grep('quant_bootstrap', .) %>% {. - 1}
    name <- str_split(file,'/')[[1]][name_idx]
    res <- tibble(transcript_id=BS_raw$X1, !!name := BS_raw[,-1] %>% as.matrix %>% {matrixStats::rowVars(.)} )
    return(res)
}

var_sum <- function(path_to_quant_files){
    bs_files <- list.files(path = path_to_quant_files,
                           pattern = 'quant_bootstraps.tsv.gz', full.names = T, recursive = T)
  

    all_bs <- lapply(bs_files, proc_boostraps) %>% reduce(left_join)#%>%
        #rename(!!s_tissue:=transcript_id) %>% inner_join(converter_tab) %>% select(transcript_id, everything())
    novel_bs <- filter(all_bs, !grepl('ENST', transcript_id))
    ref_bs <- filter(all_bs, grepl('ENST', transcript_id) )
    novel_medvar <- novel_bs %>% select(-transcript_id) %>% as.matrix() %>% rowMedians()
    ref_medvar <- ref_bs %>% select(-transcript_id) %>% as.matrix %>% rowMedians()

    refvar_95 <- quantile(ref_medvar, .95)
    paste('~removed', sum(novel_medvar > refvar_95), 'transcripts based on var' )
    return(filter(novel_bs, novel_medvar > refvar_95) %>% pull(transcript_id))
}

#TODO
quant_sum <- function(path_to_quant_files){
  quant_tab <- read_salmon(path_to_quant_files)
  remove <- rowSums(quant_tab[,-1]) == 0 # remove transcripts with all 0 counts 
  print(paste('~removed', sum(removed), 'transcripts based on exp' ))
  return(quant_tab %>% .[remove, ] %>% pull(transcript_id))
}

gtf<- rtracklayer::readGFF(gtf_file)
var_transcripts_to_remove <- var_sum(path_to_quant_files)
quant_tx_to_remove <- quant_sum(path_to_quant_files)
#write(x=transcripts_to_remove, file = paste0('data/misc/transcripts_with_high_var',{},'.txt', sep = '\n')
out_gtf = filter(gtf, !transcript_id %in% var_transcripts_to_remove, !transcript_id %in% quant_tx_to_remove)
write_gtf3(out_gtf, out_gtf_file)

