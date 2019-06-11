library(tidyverse)
args <- commandArgs(trailingOnly = T)
wd <- args[1]
sample_file <- args[2]
ref_gtf_file <- args[3]
wsize <- as.numeric(args[4])
tx_quant <- args[5]
grow_full_end_tab <- args[6]
grow_full_start_tab <- args[7]
ref_full_end_tab <- args[8]
ref_full_start_tab <- args[9]
grow_end_longer <- args[10]
grow_start_longer <- args[11]
ref_end_longer <- args[12]
ref_start_longer <- args[13]
save(args, file = 'testing/args.rdata')
setwd(wd)
ref_gtf <- rtracklayer::readGFF(ref_gtf_file)
sample_table <- read_tsv(sample_file, col_names = c('sample', 'run', 'paired','tissue', 'subtissue', 'origin'))
med_sample_count = sample_table %>% group_by(subtissue) %>% summarise(n=n()) %>% pull(n) %>% median
load(tx_quant)
counts= tx_counts
rm(tx_counts)
#lets try min 5 counts per sample
counts <- counts %>% rename(ID=transcript_id) %>% {filter(., rowSums(.[,-1]) > med_sample_count/2 )}
#variables for data collection - windowsize,number of isoforms for a given tx( fixing at 2 rn)
#lets fix window size, and only keep exons that have 2 isoforms.

end_longer_all <- filter(ref_gtf, type == 'exon',transcript_id %in% counts$ID) %>% group_by(seqid, strand, start ) %>%
    summarise(n=length(unique(end)),max_end=max(end), min_end=min(end)) %>%
    mutate(short_length=min_end-start, long_length=max_end-min_end) %>% filter(n>1)
end_longer <- end_longer_all %>% filter(n==2, short_length>=wsize, long_length>=wsize) %>% ungroup %>%
    mutate(wstart=min_end-wsize, wend=min_end+wsize, name=paste0('EL_', 1:nrow(.)), score=1000)
start_longer_all <- filter(ref_gtf, type == 'exon', transcript_id %in% counts$ID) %>% group_by(seqid, strand, end ) %>%
    summarise(n=length(unique(start)),max_start=max(start), min_start=min(start)) %>%
    mutate(short_length=end-max_start, long_length=max_start-min_start) %>% filter(n>1)
start_longer <- start_longer_all %>% filter(n==2,short_length>=wsize, long_length>=wsize ) %>% ungroup %>%
    mutate(wstart=max_start-wsize, wend=max_start+wsize, name=paste0('SL_', 1:nrow(.)), score=1000 )
exons_no_change <- ref_gtf %>% filter(type =='exon', transcript_id %in% counts$ID) %>% anti_join(end_longer_all %>% select(seqid, strand, start) ) %>%
    anti_join(start_longer_all %>% select(seqid, strand, end )) %>% mutate(length=end-start) %>% filter(length>wsize)
write_tsv(end_longer_all, grow_full_end_tab)
write_tsv(start_longer_all, grow_full_start_tab)
end_longer %>% select(seqid, wstart, wend, name, score, strand) %>%  write_tsv(grow_end_longer, col_names = F)
start_longer %>% select(seqid, wstart, wend, name, score, strand) %>%  write_tsv(grow_start_longer, col_names = F)


set.seed(98789)
index <- sample(1:nrow(exons_no_change), nrow(exons_no_change))
exons_no_change_end_longer <- index[1:(length(index)/2)] %>% exons_no_change[.,] %>%
    mutate(wstart=end-wsize, wend=end+wsize, id=paste0('ref_EL_', 1:nrow(.)), score=1000)
exons_no_change_start_longer <- index[(length(index)/2 +1):length(index)] %>% exons_no_change[.,] %>%
    mutate(wstart=start-wsize, wend=start+wsize, id=paste0('ref_SL_', 1:nrow(.)), score=1000)
write_tsv(exons_no_change_end_longer, ref_full_end_tab)
write_tsv(exons_no_change_start_longer, ref_full_start_tab)
exons_no_change_end_longer %>% select(seqid, wstart, wend, id, score, strand) %>%  write_tsv(ref_end_longer, col_names = F)
exons_no_change_start_longer %>% select(seqid, wstart, wend, id, score, strand) %>%  write_tsv(ref_start_longer, col_names = F)
