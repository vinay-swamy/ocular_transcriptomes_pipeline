library(tidyverse)
library(data.table)
library(glue)
library(yaml)
library(parallel)
library(argparse)

parser <- ArgumentParser()
parser$add_argument('--workingDir', action = 'store', dest='working_dir')
parser$add_argument('--filesYaml', action = 'store', dest = 'files_yaml')

#files_yaml <- '/data/swamyvs/ocular_transcriptomes_pipeline/files.yaml'
list2env(parser$parse_args(), .GlobalEnv)
files <- read_yaml(files_yaml)
path_to_merged_gtf <- files$combined_gtf_path
path_to_filt_gtfs <- files$raw_gtf_path
path_to_final_tissue_gtfs <- files$final_gtf_path

setwd(working_dir)
source('scripts/write_gtf.R')
nm_col_clean <- function(col){
    raw_name <- col %>% .[.!='-'] %>% .[!grepl('ENST', .)] %>% head
    name <- str_split(raw_name, '_\\d|\\|')[[1]][3] 
}#slightly slower version of nm_col, but works for all cases


process_columns <- function(tab,col_name){
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

read_track_file <- function(track_file){
    track_tab <- fread(track_file, sep = '\t', header = F) %>% as_tibble
    names <- c('transcript_id', 'gene_id','refid', 'class_code',  apply(track_tab[,-(1:4)], 2, nm_col_clean))
    colnames(track_tab) <- names
    cn <- colnames(track_tab)[-(1:4)]
    tcons2mstrg <- mclapply(cn, function(col) process_columns(track_tab,col), mc.cores = min(length(cn), parallel::detectCores() - 10))
    conv_tab <-  lapply(tcons2mstrg, function(x) x[['simple']]) %>% reduce(full_join) %>% inner_join(track_tab[, 1:4],.)
    num_comp <- sapply(tcons2mstrg, function(x) nrow(x[['comp']])) %>% sum
    message(glue('{num_comp} transcripts that map to multiple reference transcripts'))
    return(conv_tab)
}

conv_tab <- read_track_file(files$all2tissue_tracking) %>% mutate(refid = str_split(refid, '\\|') %>% sapply(function(x) x[2]))
sample_table <- read_tsv(files$sample_table_file)
subtissues <- unique(sample_table$subtissue)
full_gtf <- rtracklayer::readGFF(files$base_all_tissue_gtf)
#### there are some cases where the same reference transcript maps to multiple DNTX transcripts 
conv_tab_ref_bad_match <- filter(conv_tab, class_code == '=') %>% 
    filter(duplicated(refid)) %>% 
    pull(refid) %>% 
    {filter(conv_tab, refid %in% ., class_code == '=') }
ref_gtf <- rtracklayer::readGFF(files$ref_GTF)
full_gtf_tx_bad_match <- filter(full_gtf, transcript_id %in% conv_tab_ref_bad_match$transcript_id, type == 'exon') %>% 
    select(seqid, strand, start, end, transcript_id)
pre_num_exons_tx <- full_gtf %>% group_by(transcript_id) %>% summarise(nex=n())

refgtf_tx_bad_match <- filter(ref_gtf, transcript_id %in% conv_tab_ref_bad_match$refid, type == 'exon' ) %>%
    select(seqid, strand, start, end)

tx_match_fail <- full_gtf_tx_bad_match %>% anti_join(refgtf_tx_bad_match) %>% pull(transcript_id) %>% unique

conv_tab <- conv_tab %>% mutate(class_code = replace(class_code, transcript_id %in% tx_match_fail & class_code == '=', '!='))
####
fwrite(conv_tab, files$all2tissue_convtab, sep = '\t')

####now make the pan-eye gtfs
path_to_eye_gtfs <- 'data/gtfs/pan_eye'
eye_tissues <- c('Retina_Adult.Tissue', 'Retina_Fetal.Tissue', 'RPE_Adult.Tissue', 
                 'RPE_Fetal.Tissue', 'Cornea_Adult.Tissue', 'Cornea_Fetal.Tissue')
exp_in_eye <- conv_tab %>% select(eye_tissues) %>% apply(2, function(x) !is.na(x)) %>% {rowSums(.) >0}
eye_conv_tab <- conv_tab %>% filter(exp_in_eye) %>%  select(colnames(conv_tab)[1:4],eye_tissues)
eye_gtf <- filter(full_gtf, transcript_id %in% eye_conv_tab$transcript_id)
fwrite(eye_conv_tab,glue('{path_to_eye_gtfs}.convtab'), sep = '\t')
write_gtf3(eye_gtf, glue('{path_to_eye_gtfs}.gtf'))






