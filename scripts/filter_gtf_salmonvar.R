library(tidyverse)
library(matrixStats)
library(parallel)
args <- commandArgs(trailingOnly = T)
save(args, file='testing/fixgtfargs.Rdata')
wd <- args[1]
path_to_quant_files <- args[3]
gtf_file<- args[4]
out_gtf_file <- args[5]
setwd(wd)
source('~/scripts/write_gtf.R')


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
  

    all_bs <- lapply(bs_files, proc_boostraps) %>% reduce(left_join) %>%
        rename(!!s_tissue:=transcript_id) %>% inner_join(converter_tab) %>% select(transcript_id, everything())
    novel_bs <- filter(all_bs, !grepl('ENST', transcript_id))
    ref_bs <- filter(all_bs, grepl('ENST', transcript_id) )
    novel_medvar <- novel_bs %>% select(-transcript_id, -(!!s_tissue)) %>% as.matrix() %>% rowMedians()
    ref_medvar <- ref_bs %>% select(-transcript_id, -(!!s_tissue)) %>% as.matrix %>% rowMedians()

    refvar_95 <- quantile(ref_medvar, .95)

    return(filter(novel_bs, novel_medvar > refvar_95) %>% pull(transcript_id))
}
gtf<- rtracklayer::readGFF(gtf_file)
transcripts_to_remove <- var_sum(path_to_quant_files)
#write(x=transcripts_to_remove, file = 'data/misc/transcripts_with_high_var.txt', sep = '\n')
out_gtf = filter(out_gtf, !transcript_id %in% transcripts_to_remove)
write_gtf3(out_gtf, out_gtf_file)

