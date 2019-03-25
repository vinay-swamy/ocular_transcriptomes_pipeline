#setwd('~/NIH/eyeintegration_splicing/')
library(tidyverse)
args <-c('~/NIH/eyeintegration_splicing/',
         'PSI',
         'results_b38/all_tissues.incCts.tsv')
args <- commandArgs(trailingOnly = T)
working_dir=args[1]
type=args[2]
outfile=args[3]
setwd(working_dir)

combine_allTissues <- function(event, type){
    #sample_table <- read_tsv('sampleTableV3.tsv', col_names = c('sample','run','paired','tissue','subtissue','origin'))
    files0 <- list.files('rmats_clean', pattern = paste0(type,'.',event), recursive = T, full.names = T)
    tissues <- strsplit(files0,'/') %>% sapply(function(x) x[[2]])
    event_df_list <- lapply(files0,function(x)suppressMessages(read_tsv(x)))
    keep <-  lapply(event_df_list,nrow) > 0 & tissues != 'synth'
    event_df_list <- event_df_list[keep]
    all_tissue_df <- reduce(event_df_list, full_join)
    all_tissue_df[is.na(all_tissue_df)] <- 0
    event_t <- strsplit(event, '\\.') %>% sapply(function(x)x[[1]])
    all_tissue_df <- all_tissue_df %>% mutate(event=event_t)
    all_tissue_df
}

all_events <- c('SE.MATS.JC.txt', 'RI.MATS.JC.txt', 'MXE.MATS.JC.txt', 'A5SS.MATS.JC.txt','A3SS.MATS.JC.txt')
all_events_df <- lapply(all_events, function(x) combine_allTissues(event = x, type=type)) %>% bind_rows 
#  a single exon can be involved in multiple events ie in SE and RI; these can be controadictory as well, ie the the same exon
# might be a SE but also called as an MXE, involving the same exons ; a lot of exons are called twice like this, but   
# from looking at a few of these, it seems like the SE events have the most reads; so for now I'm going to just drop the duplicated ones
# most likely going to use some kind of read threshold when classifying these events
all_events_df_f <- all_events_df[!duplicated(all_events_df[,c('seqid','strand','start','end')]),]
write_tsv(all_events_df_f, outfile)
                         


