library(tidyverse)
library(tximport)
library(parallel)
#args <- c('~/NIH/eyeintegration_splicing/', 'data/quant_files', 'sampleTableDev.tsv','testing/all_tissues.tracking', 'testing/qout.tsv', 'testing/tcontab.tsv')
args <- commandArgs(trailingOnly = T)
wd <- args[1]
quant_path <- args[2]
sample_table_file <- args[3]
track_file <- args[4]
out_quant_table <- args[5]
TCONS_to_MSTRG_file <- args[6]

nm_col <- function(col){
    col=col[col!='-']
    name=str_split(col[1], ':|\\|')[[1]][2] %>% str_split('_MSTRG') %>% .[[1]] %>% .[1]
    return(name)
}

process_columns <- function(tab,col_name){
    print(col_name)
    tab <- tab %>% filter( (!! rlang::sym(col_name)) != '-')
    col <- tab %>% pull(!!col_name)
    name_col <- tab %>% pull(transcript_id)
    det <- suppressWarnings(str_split(col, ':|\\|'))
    z <- '
    most of the oIds have a 1:1 mapping to tx ids, so split into 2 chunks for faster run time, like a lot faster
    _s denotes simple case, _c dentotes complex case
    '

    d_lengths <- sapply(det, length)
    base <- min(d_lengths)
    simple <- d_lengths == base
    det_s <- det[simple]
    name_col_s <- name_col[simple]
    tx_simple <-  lapply(1:length(det_s), function(i)  det_s[[i]][3] %>%
                c(name_col_s[i], .)) %>%
         do.call(rbind, .) %>% as.data.frame(stringsAsFactors=F)
    colnames(tx_simple) <- c('transcript_id', col_name)
    #%>% rename(!!col_name:=oId) %>% distinct
    det_c <- det[!simple]
    name_col_c <- name_col[!simple]
    tx_comp <- lapply(1:length(det_c), function(i)  det_c[[i]][-1] %>%
                            .[grepl('MSTRG\\.\\d+\\.\\d+', .) | grepl('ENST', .)] %>%
                            {tibble(transcript_id=rep(name_col_c[i], length(.)),oId= . )}) %>%
        bind_rows() %>% rename(!!col_name:=oId) %>% distinct

    return(list(simple=tx_simple, comp=tx_comp))
}

read_salmon <- function(qdir,tissue){
    qfiles <- list.files(paste(qdir, tissue, sep = '/'), pattern = 'quant.sf', recursive = T, full.names = T)
    name_idx <- str_split(qfiles[1], '/')[[1]] %>% grep('quant.sf', .) %>% {. -1}
    names <- str_split(qfiles,'/') %>% sapply(function(x) x[name_idx])
    txi <- tximport::tximport(files=qfiles, type='salmon', txOut = T, countsFromAbundance = 'lengthScaledTPM')
    colnames(txi$counts) <- names
    counts <- txi$counts %>% as.data.frame %>% mutate(!!tissue:=rownames(.)) %>% select(!!tissue, everything())
}

process_salmon_complex <- function(path, tissue){
    quant <- suppressMessages(read_salmon(quant_path,tissue))
    quant_simple <- inner_join(tc2mstrg_simple %>% select(transcript_id, !!tissue),quant )
    multi_map <- tc2mstrg_complex[[tissue]]
    which_max_exp <- multi_map %>% inner_join(quant) %>% mutate(avg_exp=rowSums(.[,-(1:2)])) %>% group_by(transcript_id) %>%
        do(.[which.max(.$avg_exp), tissue])
    quant_comp <- inner_join(which_max_exp, quant)
    return(list(quant=bind_rows(quant_simple, quant_comp) %>% select(- !!tissue), comp=which_max_exp ))
}



setwd(wd)
sample_table <- read_tsv(sample_table_file)
track_tab <- read_tsv(track_file, col_names = F) %>% select(-X3, -X4)
tissues <- unique(sample_table$subtissue)
names <- c('transcript_id', 'gene_id',  apply(track_tab[,-(1:2)], 2, nm_col))
colnames(track_tab) <- names
cn <- colnames(track_tab)[-(1:2)]
tcons2mstrg <- mclapply(cn, function(col) process_columns(track_tab,col), mc.cores = 8)
tc2mstrg_simple <- lapply(tcons2mstrg, function(x) x[['simple']]) %>% reduce(full_join)
tc2mstrg_complex <-lapply(tcons2mstrg, function(x) x[['comp']])

all_quant <- lapply(tissues, function(x) suppressMessages(read_salmon(quant_path, x)))


tis_cols_idx <- ncol(tc2mstrg_simple)
complete_quant <- reduce(all_quant, left_join, .init = tc2mstrg_simple) %>% .[,-(2:tis_cols_idx)]
save(complete_quant, file = out_quant_table)
save(tc2mstrg_complex, file='data/misc/multi_mapping_transcripts.Rdata')
write_tsv(tc2mstrg_simple, TCONS_to_MSTRG_file)
