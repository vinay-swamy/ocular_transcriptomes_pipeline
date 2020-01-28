library(tidyverse)
library(matrixStats)
library(parallel)
library(argparse)

parser <- ArgumentParser()
parser$add_argument('--workingDir', action = 'store', dest = 'wd')
parser$add_argument('--pathToQuantFiles',action = 'store', dest = 'path_to_quant_files')
parser$add_argument('--gtfFile', action = 'store', dest = 'gtf_file')
parser$add_argument('--outGtfFile', action = 'store', dest = 'out_gtf_file')
list2env(parser$parse_args(), .GlobalEnv)

#args <- commandArgs(trailingOnly = T)
#save(args, file='testing/fixgtfargs.Rdata')
# wd <- args[1]
# path_to_quant_files <- args[2]
# gtf_file<- args[3]
# out_gtf_file <- args[4]
# setwd(wd)
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
    print(paste('~removed', sum(novel_medvar > refvar_95), 'transcripts based on var' ) )
    return(filter(novel_bs, novel_medvar > refvar_95) %>% pull(transcript_id))
}

# remove transcripts that have 0 counts across all samples 
quant_sum <- function(path_to_quant_files){
  quant_tab <- read_salmon(path_to_quant_files)
  remove <- rowSums(quant_tab[,-1]) == 0 # remove transcripts with all 0 counts 
  print(paste('~removed', sum(remove), 'transcripts based on exp' ))
  remove_tx <- quant_tab %>% .[remove, ] %>% pull(transcript_id)
  quan_tab_filt <- quant_tab[!remove,]
  return(list(quant_tab=quan_tab_filt, remove_tx=remove_tx))
}

single_tx_removal <- function(quant_tab, gtf){
    avg_quant <- tibble(transcript_id = quant_tab$transcript_id, 
                        avg_quant = rowMeans(quant_tab[,-1]))
    #only remove single exon tx that are part of genes with more than one transcript, unless they are a novel loci 

    single_exon_tx <- gtf %>%
      filter(type == 'exon') %>% 
      group_by(transcript_id) %>% 
      summarise(n_exon=n()) %>% 
      filter(n_exon == 1) %>% pull(transcript_id)
    multi_tx_gene <- gtf %>% 
      filter(type == 'transcript') %>% 
      group_by(gene_name) %>% 
      summarise(n_tx=n()) %>% 
      filter(n_tx>1) %>% pull(gene_name)
    novel_loci <- gtf %>% filter(class_code == 'u') %>% 
      pull(transcript_id)
    #keep transcripts that are single exon and belonging to a gene with multiple transcripts, or transcripts that are both
    # a single exon trancript and a novel loci 
  valid_single_exon_tx <-  gtf %>% filter((  ((transcript_id %in% single_exon_tx ) & (gene_name %in% multi_tx_gene)) | 
                      ((transcript_id %in% single_exon_tx) & (transcript_id %in% novel_loci))
                   )
                  ) %>% pull(transcript_id) %>% unique()
  avg_quant_single_tx <- filter(avg_quant, transcript_id %in% valid_single_exon_tx)
  avg_quant_multi_tx <- filter(avg_quant, !transcript_id %in% valid_single_exon_tx)
  outlier_mark <- quantile(avg_quant_multi_tx$avg_quant, .75) * 1.5 # outliers are defined as 1.5* 75th percentile 
  single_tx_below_outlier_mark <- avg_quant_single_tx %>% filter( avg_quant <= outlier_mark) %>% pull(transcript_id)
  print(paste0('~removed ', length(single_tx_below_outlier_mark), ' single exon txs'))
  return(single_tx_below_outlier_mark)
}

gtf<- rtracklayer::readGFF(gtf_file)
var_transcripts_to_remove <- var_sum(path_to_quant_files)
res <- quant_sum(path_to_quant_files)
quant_tx_to_remove <- res[['remove_tx']]
single_tx_to_remove <- single_tx_removal(quant_tab = res[['quant_tab']], gtf = gtf)


#write(x=transcripts_to_remove, file = paste0('data/misc/transcripts_with_high_var',{},'.txt', sep = '\n')
out_gtf = gtf %>% filter(!transcript_id %in% var_transcripts_to_remove, 
                         !transcript_id %in% quant_tx_to_remove,
                         !transcript_id %in% single_tx_to_remove)
write_gtf3(out_gtf, out_gtf_file)

