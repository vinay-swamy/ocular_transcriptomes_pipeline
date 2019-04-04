library(tidyverse)
args <- commandArgs(trailingOnly = T)
working_dir <- args[1]
out_med_file <- args[2]
out_inc_file <- args[3]
setwd(working_dir)

med_files <- list.files('results/complete_rmats_output', pattern ='medCounts.tsv' , full.names = T )[-3]
inc_files <- list.files('results/complete_rmats_output', pattern ='incLevel.tsv' , full.names = T )[-3]
evs <- c('A3SS','A5SS','RI','SE')
med_tabs <- lapply(med_files, function(x) read_tsv(x))
inc_tabs <- lapply(inc_files, function(x) read_tsv(x))
med_tab_cn <- c(c('GeneID', 'seqid'	,'strand',	'start', 'end',	'upstreamES',	'upstreamEE',	'downstreamES',	'downstreamEE'), 
        colnames(med_tabs[[1]])[-c(1:9)])
inc_tab_cn <- c(c('GeneID', 'seqid'	,'strand',	'start', 'end',	'upstreamES',	'upstreamEE',	'downstreamES',	'downstreamEE'), 
                colnames(inc_tabs[[1]])[-c(1:9)])
med_tabs <- lapply(evs, function(x) med_files[grepl(x,med_files)] %>% 
                                    read_tsv(col_names = med_tab_cn, skip = 1) %>% 
                                    mutate(event_type=x)) %>% do.call(rbind,.)
inc_tabs <- lapply(evs, function(x) inc_files[grepl(x,med_files)] %>% read_tsv(col_names = inc_tab_cn, skip = 1) %>% 
                                    mutate(event_type=x)) %>% 
                                    do.call(rbind,.)

dups <- filter(med_tabs[2:5], duplicated(med_tabs[2:5])) %>% split(1:nrow(.))
all <- split(med_tabs[2:5], 1:nrow(med_tabs))
distinct <- !all%in%dups
med_tabs <- med_tabs[distinct,]
inc_tabs <- inc_tabs[distinct,]
write_tsv(med_tabs, out_med_file)
write_tsv(inc_tabs, out_inc_file)









