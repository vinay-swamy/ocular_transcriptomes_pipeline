library(tidyverse)
library(tximport)
library(parallel)
source('~/scripts/read_salmon.R')
args <- c('/data/swamyvs/eyeintegration_splicing/', 'sampleTableV6.tsv', 'data/rawST_tx_quant_files/',
          'data/misc/TCONS2MSTRG.tsv', 'data/exp_files/all_tissue_quant.tsv.gz')
args <- commandArgs(trailingOnly = F)
wd <- args[1]
sample_file <- args[2]
path_to_quant <- args[3]
tcons2mstrg_file <- args[4]
out_exp_file <- args[5]



read_quant <- function(path, s_tissue){
    quant <- read_salmon(path) %>% rename(!!s_tissue := transcript_id)
}

sample_table <- read_tsv(sample_file)
subtissues <- unique(sample_table$subtissue)
quant_paths <- paste0(path_to_quant, subtissues)
tc2ms <- read_tsv(tcons2mstrg_file)
all_quant <- mclapply(seq_along(subtissues), function(i) read_quant(quant_paths[i], subtissues[i]), mc.cores = 12) %>% 
    reduce(left_join, .init =tc2ms)
cols <- colnames(all_quant)
cols_to_keep <-  c('transcript_id', cols[cols %in% sample_table$sample])
all_quant <- all_quant[,cols_to_keep]















