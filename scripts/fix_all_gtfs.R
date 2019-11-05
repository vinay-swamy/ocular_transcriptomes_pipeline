library(tidyverse)
library(matrixStats)
library(parallel)
args <- c('/data/swamyvs/eyeintegration_splicing/', 'sampleTableV6.tsv', 'ref/gencode_comp_ano.gtf',
          'data/rawST_tx_quant_files/', 'data/gtfs/raw_tissue_gtfs/', 'data/gtfs/filtered_tissue/',
          'data/gffcomp_dir/all_tissues.combined.gtf', 'data/gffcomp_dir/all_tissues.tracking' ,
          'data/gtfs/all_tissues.combined.gtf', 'data/misc/TCONS2MSTRG.tsv')
args <- commandArgs(trailingOnly = T)
save(args, file='testing/fixgtfargs.Rdata')
wd <- args[1]
sample_table_file <- args[2]
ref_gtf_file <- args[3]
path_to_quant_files <- args[4]
path_to_raw_gtfs <- args[5]
path_to_filtered_gtfs <- args[6]
gfc_gtf_with_st_file <- args[7]
track_file <- args[8]
out_gtf_file <- args[9]
final_track_file <- args[10]

setwd(wd)

sample_table <- read_tsv(sample_table_file)
subtissues <- unique(sample_table$subtissue)
source('~/scripts/write_gtf.R')

#process the tracking file from the merges
#----
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
gfc_gtf_with_st <- rtracklayer::readGFF(gfc_gtf_with_st_file )
track_tab <- read_tsv(track_file, col_names = F) %>% select(-X3, -X4)
tissues <- unique(sample_table$subtissue)
names <- c('transcript_id', 'gene_id', 'stringtie',  apply(track_tab[,-(1:3)], 2, nm_col))
colnames(track_tab) <- names
cn <- colnames(track_tab)[-(1:2)]
tcons2mstrg <- mclapply(cn, function(col) process_columns(track_tab,col), mc.cores = 12)
tc2mstrg_simple <- lapply(tcons2mstrg, function(x) x[['simple']]) %>% reduce(full_join)
tcm_mat <- apply(tc2mstrg_simple[,-1],2, function(col) !is.na(col))
det_st <- tcm_mat[,1]
det_in_tissues <- tcm_mat[,-1] %>% {rowSums(.) > 0 } # if a transcript was deteced in at least one tissue, we want to keep it
tc2mstrg_det <- filter(tc2mstrg_simple, det_in_tissues) 

#----

#Read in salmon bootstrap variance and identify transcripts to remove.
#----
proc_boostraps <- function(file){
    BS_raw <- suppressMessages(read_tsv(file, col_names = F) )
    name_idx <- str_split(file, '/')[[1]] %>% grep('quant_bootstrap', .) %>% {. - 1}
    name <- str_split(file,'/')[[1]][name_idx]
    res <- tibble(transcript_id=BS_raw$X1, !!name := BS_raw[,-1] %>% as.matrix %>% {matrixStats::rowVars(.)} )
    return(res)
}

var_sum <- function(s_tissue){
    bs_files <- list.files(path = paste0(path_to_quant_files,s_tissue),
                           pattern = 'quant_bootstraps.tsv.gz', full.names = T, recursive = T)
    converter_tab <- tc2mstrg_det %>% select(transcript_id, stringtie, !!s_tissue) %>% filter(!is.na(.[,s_tissue])) %>%
        mutate(transcript_id=replace(transcript_id, grepl('ENST', stringtie),stringtie[grepl('ENST', stringtie)] )) %>%
        select(-stringtie)

    all_bs <- lapply(bs_files, proc_boostraps) %>% reduce(left_join) %>%
        rename(!!s_tissue:=transcript_id) %>% inner_join(converter_tab) %>% select(transcript_id, everything())
    novel_bs <- filter(all_bs, !grepl('ENST', transcript_id))
    ref_bs <- filter(all_bs, grepl('ENST', transcript_id) )
    novel_medvar <- novel_bs %>% select(-transcript_id, -(!!s_tissue)) %>% as.matrix() %>% rowMedians()
    ref_medvar <- ref_bs %>% select(-transcript_id, -(!!s_tissue)) %>% as.matrix %>% rowMedians()

    refvar_95 <- quantile(ref_medvar, .95)

    return(filter(novel_bs, novel_medvar > refvar_95) %>% pull(transcript_id))
}

transcripts_to_remove <- mclapply(subtissues, var_sum, mc.cores = 12) %>% combine
write(x=transcripts_to_remove, file = 'data/misc/transcripts_with_high_var.txt', sep = '\n')
#----



#write the master gtf.
#----
gfc_gtf <- filter(gfc_gtf_with_st, transcript_id %in% tc2mstrg_det$transcript_id, #only keep gffc transcripts,
                  !transcript_id %in% transcripts_to_remove) # remove high var tx.
conv_tab <- filter(gfc_gtf, type == "transcript") %>% select(transcript_id,gene_id, gene_name) %>%
    inner_join(tc2mstrg_det %>% select(transcript_id, stringtie)) %>%
    mutate(master_transcript_id=replace(transcript_id, grepl('ENST', stringtie),stringtie[grepl('ENST', stringtie)] )) %>%
    select(-stringtie)

gfc_gtf_full <- gfc_gtf %>% select(-gene_id, -oId, -class_code, -tss_id, -contained_in, -gene_name, -cmp_ref) %>%
    left_join(conv_tab)

tc2mstrg_final <- gfc_gtf_full %>% filter(type == "transcript") %>%
    select(transcript_id, master_transcript_id,gene_name, gene_id) %>%
    inner_join(tc2mstrg_det )



gn2gi <- rtracklayer::readGFF(ref_gtf_file) %>% select(gene_name, gene_id) %>% distinct %>% filter(!duplicated(.[,'gene_name']))

in_gtf <- gfc_gtf_full %>% select(-transcript_id) %>% rename(transcript_id=master_transcript_id) %>%
    mutate(gffc_loc=gene_id, gene_id=transcript_id)
tx2g <- in_gtf %>%
    filter(type == 'transcript') %>%
    select(transcript_id, gene_name) %>%
    distinct %>%
    mutate(gene_name=replace(gene_name, is.na(gene_name), transcript_id[is.na(gene_name)])) %>% distinct %>%
    left_join(gn2gi) %>% mutate(gene_id=replace(gene_id, is.na(gene_id), transcript_id[is.na(gene_id)]))


out_gtf <- in_gtf %>% select(-gene_name, -gene_id) %>% left_join( tx2g)
write_gtf3(out_gtf, out_gtf_file)
write_tsv(tc2mstrg_final, final_track_file)
#----
#write the tissue gtfs
fix_tissue_specific_gtf <- function(s_tissue){
    gtf <- rtracklayer::readGFF(paste0(path_to_raw_gtfs, s_tissue, '_st.gtf'))
    id_tab <- select(tc2mstrg_final, master_transcript_id=transcript_id, !!s_tissue) %>% filter(!is.na(.[,s_tissue])) %>%
        rename(transcript_id= !!s_tissue) %>% filter()
    res <- gtf %>% inner_join(id_tab) %>% select(-transcript_id) %>% rename(transcript_id=master_transcript_id)
    print(c(s_tissue, nrow(gtf)-nrow(res) ) )
    write_gtf3(res, paste0(path_to_filtered_gtfs, s_tissue, '.gtf'))
}

mclapply(subtissues, fix_tissue_specific_gtf, mc.cores = 12)

save.image(file = 'testing/fix_all_gtfs_image.Rdata')












