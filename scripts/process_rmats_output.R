library(tidyverse)

#going to need location_files
#args <- c('~/NIH/eyeintegration_splicing/', 'sampleTableV5.tsv', 'rmats_out/', 'ref/rmats_locs/', 'testing/smo.df', 'testing/dfdf.df' )
args <- commandArgs(trailingOnly = T)
wd=args[1]
sample_table_file <- args[2]
rmats_out_dir <- args[3]
rm_loc_dir <- args[4]
all_tissue_psi_outfile <- args[5]
all_tissue_incCounts_outfile <- args[6] 

sample_table <- read_tsv(sample_table_file, col_names = c('sample', 'run', 'paired', 'tissue', 'subtissue', 'origin'))
subtissues <- sample_table %>% filter(paired == 'y') %>% pull(subtissue) %>% unique()
rm_locs <- paste0(rm_loc_dir, subtissues,'.rmats.txt' ) %>% .[!grepl('synth',.)]

process_events_for_tissue <- function(loc){
    z<- "
        if you check here https://www.pnas.org/content/pnas/111/51/E5593.full.pdf?with-ds=yes you can see how the different event types are defined.
        Its possible for an alternative splice site exon to be involved in a skipped exon event ie ig given A-B-C, where B is has ASS variant B*m,
        you ca have A-C, A-B-C, A-B*-C as possible variants, whereas in rMATs the ASS events have only 2 possibilities A-B-C, A-B*-C. I believe that 
        the skipped exon case is a better measure of PSI,because it counts read to more exons(so does MXE), so what I'm going to do is use the SE as 
        the default case, and thne only pull exons from the other files that are not present in the skipped exon file 
        "
    
    samples <- scan(loc, character(),  sep = ',') %>% str_split('/') %>% .[.!='']  %>% sapply(function(x) x[grep('.bam', x, fixed = T) -1]) 
    print(samples)
    tissue <- str_split(loc, '/|\\.rmats')[[1]][3]
    event_files <- paste0(rmats_out_dir, tissue,'/', c('SE', 'A3SS', 'A5SS', 'MXE', 'RI'), '.MATS.JC.txt')
    
    process_count_column <- function(df, count_col, len_col, loc_cols){
        df <- df[,c(loc_cols, count_col)] %>% distinct # we dont want to double count the same counts
        counts <-  df %>% 
            pull(!!count_col)  %>% str_split(',') %>% #split into a list and length scale counts 
            lapply(as.numeric) %>%  do.call(rbind, .) %>%  { . / df %>% pull(!!len_col)} %>% cbind(df[1:4],.) 
        colnames(counts) <- c('seqid', 'strand', 'start','end', samples)
        counts %>% group_by(seqid, strand, start, end) %>% summarise_all(sum) %>% ungroup
    }
    replace_nan <- function(x){
        x <- as.matrix(x)
        x[is.nan(x)] <- 0
        x
    }
    proc_SE_RI_files <- function(file, lcols){
        tab <- read_tsv(file)
        ex_counts <- process_count_column(tab, 'SJC_SAMPLE_1', 'SkipFormLen', lcols)
        inc_counts <-  process_count_column(tab, 'IJC_SAMPLE_1', 'IncFormLen', lcols)
        psi <- (inc_counts[,-(1:4)]/( inc_counts[,-(1:4)] + ex_counts[,-(1:4)] )) %>% replace_nan %>% {cbind(inc_counts[,1:4], .)}
        cn_psi<- c('seqid', 'strand', 'start', 'end', paste0(colnames(psi)[-(1:4)], '_psi'))
        cn_incCounts <- c('seqid', 'strand', 'start', 'end', paste0(colnames(inc_counts)[-(1:4)], '_inc'))
        colnames(inc_counts) <- cn_incCounts
        colnames(psi) <- cn_psi
        return(list(incCounts=inc_counts,
                    psi=psi))
    }
    proc_ASS_MXE_files <- function(file, long_name, short_name ){
        z <-"
        It's a little confusing, but for the ASS files, the long isoform is the 'inlcluded' isoform, and the short is the 'excluded'.
        Therefore psi for the long form is inc/(inc+ex), and the psi for the short is ex/(inc+ex)
        "
        tab <- read_tsv(file)
        long <- tab %>% select(chr, strand, contains(long_name), contains('Len'), contains('SAMPLE'))
        long_cols <- long %>% select(chr, strand, contains(long_name), contains('Len')) %>% colnames 
        long_inc_counts <- process_count_column(long, 'IJC_SAMPLE_1', 'IncFormLen',long_cols)
        long_ex_counts <- process_count_column(long, 'SJC_SAMPLE_1','SkipFormLen', long_cols)
        long_psi <-  (long_inc_counts[,-(1:4)]/(long_inc_counts[,-(1:4)]+ long_ex_counts[,-(1:4)])) %>% replace_nan %>% {cbind(long_ex_counts[,1:4],.)}
        short <- tab %>% select(chr, strand, contains(short_name), contains('Len'), contains('SAMPLE'))
        short_cols <- short %>% select(chr, strand, contains(short_name), contains('Len')) %>% colnames 
        short_inc_counts <- process_count_column(short, 'IJC_SAMPLE_1', 'IncFormLen', short_cols)
        short_ex_counts <- process_count_column(short, 'SJC_SAMPLE_1', 'SkipFormLen', short_cols)
        short_psi <- (short_ex_counts[,-(1:4)]/ (short_ex_counts[,-(1:4)] + short_inc_counts[,-(1:4)])) %>% replace_nan %>% cbind(short_ex_counts[,1:4],.)
        inc_counts <- bind_rows(long_inc_counts, short_ex_counts)
        psi <-  bind_rows(long_psi, short_psi)
        cn_psi<- c('seqid', 'strand', 'start', 'end', paste0(colnames(psi)[-(1:4)], '_psi'))
        cn_incCounts <- c('seqid', 'strand', 'start', 'end', paste0(colnames(inc_counts)[-(1:4)], '_inc'))
        colnames(inc_counts) <- cn_incCounts
        colnames(psi) <- cn_psi
        return(list(incCounts=inc_counts,
                    psi=psi))
    } 

    se_cols <- c('chr','strand', "exonStart_0base", "exonEnd","IncFormLen", "SkipFormLen" )
    ri_cols <- c("chr", "strand", "riExonStart_0base", "riExonEnd","IncFormLen", "SkipFormLen")
    SE_all <- proc_SE_RI_files(event_files[[1]], se_cols)
    
    A3_all <- proc_ASS_MXE_files(event_files[[2]], 'long', 'short')
    A5_all <- proc_ASS_MXE_files(event_files[[3]],'long', 'short')
    #it looks like the same exon can sometimes be called as long and short smh going to just remove it being called as the short version
    MXE_all <- proc_ASS_MXE_files(event_files[[4]], '1st', '2nd')
    RI_all <- proc_SE_RI_files(event_files[[5]], ri_cols)
    join_cols <- c('seqid', 'strand', 'start', 'end')
    #join the MXE first, because its got reads spanning 3 exons 
    join_lvl <- function(lvl){
        all_lvl <- SE_all[[lvl]] %>% mutate(event_file='SE')
        all_lvl <- all_lvl %>% anti_join(MXE_all[[lvl]],., by=join_cols) %>% mutate(event_file='MXE') %>% bind_rows(all_lvl, .) 
        all_lvl <- all_lvl %>% anti_join(A3_all[[lvl]],., by=join_cols) %>% mutate(event_file='A3SS') %>% bind_rows(all_lvl, .)
        all_lvl <- all_lvl %>% anti_join(A5_all[[lvl]],., by=join_cols) %>% mutate(event_file='A5SS') %>% bind_rows(all_lvl, .) 
        all_lvl <- all_lvl %>% anti_join(RI_all[[lvl]],., by=join_cols) %>% mutate(event_file='RI') %>% bind_rows(all_lvl, .) 
    }
    
    all_inc <- join_lvl('incCounts') %>% filter(!duplicated(.[,join_cols]))
    all_psi <- join_lvl('psi')  %>% filter(!duplicated(.[,join_cols]))
    return(list(incCounts=all_inc, psi=all_psi))
}
#currently no support for empty files will figure that out later 
all_tissues_results <- lapply(rm_locs, function(x) suppressMessages(process_events_for_tissue(x)))

reduce_lvl <- function(lvl){
    lapply(all_tissues_results,function(x) x[[lvl]]) %>% 
        reduce(full_join, by=c('seqid', 'strand', 'start' ,'end'))
}
all_tissue_psi <- reduce_lvl('psi')
psi_tissue_which.event <- all_tissue_psi %>% select(contains('event_file')) %>% replace(is.na(.), 'nd')
all_tissue_psi <- all_tissue_psi %>% select(-contains('event_file'))
all_tissue_incCounts <- reduce_lvl('incCounts') %>% select(-contains('event_file'))
write_tsv(all_tissue_psi, all_tissue_psi_outfile)
write_tsv(all_tissue_incCounts, all_tissue_incCounts_outfile)



