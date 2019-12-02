library(tidyverse)
library(tximport)
library(parallel)
#args <- c('~/NIH/eyeintegration_splicing/', 'data/quant_files', 'sampleTableDev.tsv','testing/all_tissues.tracking', 'testing/qout.tsv', 'testing/tcontab.tsv')
args <- commandArgs(trailingOnly = T)
wd <- args[1]
sample_table_file <- args[2]
track_file <- args[3]
raw_gtf_file <- args[4]
ref_gtf_file <- args[5]
TCONS_to_MSTRG_file <- args[6]
out_gtf_file <- args[7]
#save(args, file='/tmp/smoobargs.rdata')
source('~/scripts/write_gtf.R')
# nm_col <- function(col){
#     col=col[col!='-']
#     name=str_split(col[1], '_0|\\|')[[1]][3] %>% str_split('_MSTRG') %>% .[[1]] %>% .[1]
#     return(name)
# }

nm_col_clean <- function(col){
    raw_name <- col %>% .[.!='-'] %>% .[!grepl('ENST', .)] %>% head
    name <- str_split(raw_name, '\\|')[[1]][2] %>% str_split('_0') %>% .[[1]] %>% .[1]
}#slightly slower version of nm_col, but works for all cases


# get_col_slow <- function(tissue){
#     k <- apply(track_tab[,-(1:4)], 2, function(x) x[x!='-'] %>% any(grepl(tissue, .)))
#     stopifnot(length(k) > 1)
#     names[which(k)] <- tissue
# }

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
    tx_comp=data.frame()
    if(sum(!simple) > 0){
        det_c <- det[!simple]
        name_col_c <- name_col[!simple]

        tx_comp <- lapply(1:length(det_c), function(i)  det_c[[i]][-1] %>%
                                .[grepl('MSTRG\\.\\d+\\.\\d+', .) | grepl('ENST', .)] %>%
                                {tibble(transcript_id=rep(name_col_c[i], length(.)),oId= . )}) %>%
            bind_rows() %>% rename(!!col_name:=oId) %>% distinct
        }
    return(list(simple=tx_simple, comp=tx_comp))
}


setwd(wd)
sample_table <- read_tsv(sample_table_file)
track_tab <- read_tsv(track_file, col_names = F) 
tissues <- unique(sample_table$subtissue)
names <- c('transcript_id', 'gene_id','refid', 'class_code',  apply(track_tab[,-(1:4)], 2, nm_col_clean))

stopifnot(sum(!tissues %in% names) == 0)#check

colnames(track_tab) <- names
cn <- colnames(track_tab)[-(1:4)]
tcons2mstrg <- mclapply(cn, function(col) process_columns(track_tab,col), mc.cores = 8)
tc2mstrg_simple <- lapply(tcons2mstrg, function(x) x[['simple']]) %>% reduce(full_join)

gfc_gtf <- rtracklayer::readGFF(raw_gtf_file)
ref_gtf <- rtracklayer::readGFF(ref_gtf_file)
ref_gtf_tx <- ref_gtf %>% filter(type == 'transcript') %>% select(seqid, strand, start, end)
gfcgtf_reftx_absmatch <- gfc_gtf %>% filter(type  == 'transcript', class_code == '=') %>% 
    select(seqid, strand, start, end, transcript_id) %>% inner_join(ref_gtf_tx) %>% pull(transcript_id)
targ <- gfc_gtf$transcript_id == 'transcript' & gfc_gtf$class_code == '='
gfc_gtf$class_code[targ] <- '*'
modify <- gfc_gtf$transcript_id %in% gfcgtf_reftx_absmatch & gfc_gtf$class_code == '*'
gfc_gtf$class_code[modify] <- '='

tc2oid <- gfc_gtf %>% filter(type == 'transcript') %>% 
    mutate(new_id=replace(transcript_id, class_code == '=', cmp_ref[class_code == '='])) %>% 
    select(transcript_id, new_id, class_code, gene_name, oId)
final_gtf <- gfc_gtf %>% select(-class_code, -gene_name, -oId) %>% left_join(tc2oid) %>% select(-transcript_id) %>% 
    rename(transcript_id=new_id)

tcons2mstrg_complete <- tc2oid %>% select(transcript_id, new_id) %>% left_join(tc2mstrg_simple, .) %>% 
    select(-transcript_id) %>% rename(transcript_id=new_id) %>% select(transcript_id, everything())
write_tsv(tcons2mstrg_complete, TCONS_to_MSTRG_file)
write_gtf3(final_gtf, out_gtf_file)


