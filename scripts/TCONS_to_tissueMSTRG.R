library(tidyverse)
library(tximport)
#args <- c('~/NIH/eyeintegration_splicing/', 'data/quant_files', 'sampleTableDev.tsv','testing/all_tissues.tracking', 'testing/qout.tsv', 'testing/tcontab.tsv')
args <- commandArgs(trailingOnly = T)
wd <- args[1]
quant_path <- args[2]
sample_table_file <- args[3]
track_file <- args[4]
out_quant_table <- args[5]
TCONS_to_MSTRG_file <- args[6]

nm_col <- function(col){
    keep <- sapply(tissues, function(x) any(grepl(x, col)))
    return(tissues[keep])
}

process_columns <- function(tab,col_name){
    col <- tab %>% pull(!!col_name) 
    name_col <- tab %>% pull(transcript_id)
    det <- suppressWarnings(str_split(col, ':|\\|'))
    tcons2oid <- lapply(1:length(det), function(i)  det[[i]][-1] %>%
                            .[grepl('MSTRG\\.\\d+\\.\\d+', .) | grepl('ENST', .)] %>%
                            {tibble(transcript_id=rep(name_col[i], length(.)),oId= . )}) %>% 
        bind_rows() %>% rename(!!col_name:=oId) %>% distinct 
    return(tcons2oid)
}

read_salmon <- function(qdir,tissue){
    qfiles <- list.files(paste(qdir, tissue, sep = '/'), pattern = 'quant.sf', recursive = T, full.names = T)
    names <- str_split(qfiles, '/') %>% sapply(function(x) x[4])
    txi <- tximport::tximport(files=qfiles, type='salmon', txOut = T, countsFromAbundance = 'lengthScaledTPM')
    colnames(txi$counts) <- names 
    counts <- txi$counts %>% as.data.frame %>% mutate(!!tissue:=rownames(.)) %>% select(!!tissue, everything())
}

setwd(wd)
sample_table <- read_tsv(sample_table_file, col_names = c('sample', 'run', 'paired', 'tissue', 'subtissue', 'origin'))
track_tab <- read_tsv(track_file, col_names = F) %>% select(-X3, -X4)
tissues <- unique(sample_table$subtissue)
names <- c('transcript_id', 'gene_id',  apply(track_tab[,-(1:2)], 2, nm_col))
colnames(track_tab) <- names
cn <- colnames(track_tab)[-(1:2)]
tcons2mstrg <- lapply(cn, function(col) process_columns(track_tab,col)) %>% reduce(full_join)
all_quant <- lapply(tissues, function(x) suppressMessages(read_salmon(quant_path,x)))
tis_cols_idx <- ncol(tcons2mstrg)
complete_quant <- reduce(all_quant, left_join, .init = tcons2mstrg) %>% .[,-(2:tis_cols_idx)]
save(complete_quant, file = out_quant_table)
write_tsv(tcons2mstrg, TCONS_to_MSTRG_file)









