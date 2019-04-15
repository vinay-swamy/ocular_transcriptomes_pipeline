library(tidyverse)
library(matrixStats)


working_dir <- args[1]
novel_exon_exp_tabs <- args[2]
ref_exon_tab <- args[3]
exon_ca_bed <- args[4]
intron_tab <- args[5]
fusion_gene_files <- args[6]
t_tissue <- args[7]
sample_tab <- args[8]
gtf_file <- args[9]
outfile <- args[10]

setwd('~/NIH/eyeintegration_splicing/')
load('/Volumes/data/eyeintegration_splicing/results/novel_exon_expression_tables.Rdata')
ref_exon_table <- read_tsv('/Volumes/data/eyeintegration_splicing/results/ref_exon_table.tsv')
ref_bed <- read_tsv('/Volumes/data/eyeintegration_splicing/results/exons_for_coverage_analysis.bed', col_names = c('seqid', 'start', 'end'))
intron_table <- read_tsv('/Volumes/data/eyeintegration_splicing/results/intron_info_tab.tsv') %>%  
     select(-tx, ljid=id) %>% distinct


t_tissue <- 'RPE_Fetal.Tissue'
sample_table <- read_tsv('sampleTableV4.tsv', col_names = c('sample', 'run', 'paired', 'tissue', 'subtissue', 'origin'))
samples <- filter(sample_table, subtissue == t_tissue, paired =='y') %>% pull(sample)
cov_paths <- paste0('~/NIH/eyeintegration_splicing/coverage_files/', t_tissue,'/', samples, '.regions.bed.gz') %>% .[ file.exists(.)]
t2g <- rtracklayer::readGFF('/Volumes/data/eyeintegration_splicing/results/all_tissues.combined.gtf') %>% filter(type =='transcript') %>% select(transcript_id, gene_name) %>% 
    mutate(gene_name=replace(gene_name, is.na(gene_name),transcript_id[is.na(gene_name)]))
novel_exon_locs <- nx_inc_cts %>% filter(reclassified_event == 'novel_exon', is.not_long) %>%   select(seqid,strand, start, end, reclassified_event, gene_name, ljid) %>% distinct
#novel_exon_locs %>% select(seqid, strand, start, end, gene_name) %>% distinct %>% dim
#'TCONS_00133047'
#load('/Volumes/data/eyeintegration_splicing/')
all_cov <- lapply(cov_paths, function(epath)read_tsv(epath, col_names = c('seqid', 'start', 'end', 'cov')) %>% 
                      mutate(length=end-start, ls_cov=cov/length, cov=NULL))


mean_based_zscore <- function( cov_df, sample_name, novel_exon_locs, ref_exon_table){
    novel_exon_cov <- novel_exon_locs %>% inner_join(cov_df) 
    
    ref_exon_cov <- ref_exon_table %>% select(seqid, strand, start, end, gene_name) %>% distinct %>% inner_join(cov_df)
    
    ref_gene_dists <- ref_exon_cov %>% group_by(gene_name) %>% summarise(g_mean=mean(ls_cov), g_med=median(ls_cov), g_sd=sd(ls_cov))
    
    all_cov <- inner_join(novel_exon_cov, ref_gene_dists) %>% filter(ls_cov != 0) %>%  mutate(zscore=(ls_cov - g_mean)/g_sd) %>%
        select(-length, -contains('g_'),-ls_cov, -reclassified_event, !!sample_name := zscore)
    
    all_cov
}

intron_cov <- function(cov_df ,sample_name,  intron_bed){
    icov <- cov_df  %>% inner_join(intron_bed, .) %>%  distinct %>% 
        group_by(ljid) %>% summarise(intron_cov=mean(ls_cov))
    novel_exon_cov <- novel_exon_locs %>% inner_join(cov_df) 
    
    a_cov <- icov %>% 
        left_join(novel_exon_cov, .) %>% mutate( !!sample_name := ls_cov -intron_cov) %>% select(-length, -ls_cov, -intron_cov)
    
    a_cov
    
}

all_zscores <- lapply(1:length(all_cov), 
                      function(i) suppressMessages(mean_based_zscore(all_cov[[i]], samples[i],novel_exon_locs, ref_exon_table))) %>% 
    reduce( full_join)

all_icov <- lapply(1:length(all_cov), function(i) 
    suppressMessages(intron_cov(all_cov[[i]], samples[i], intron_table) %>% distinct)) %>% reduce(left_join) 
common_ids <- intersect(all_icov$ljid, all_zscores$ljid)

all_zscores <- all_zscores %>% filter(ljid %in%common_ids)
all_icov <- all_icov %>% filter(ljid %in%common_ids)


rank_exons <-  function(zsc, intron_diff, scores, t_name){
    per_sample_score <- function(id, zsc, intron_diff, scores){
        zcut=scores[[1]]
        icut=scores[[2]]
        samp <- inner_join(zsc %>% select(ljid, zscore=!!id),
                           intron_diff %>% select(ljid,iscore=!!id)) %>%
            mutate(zscore= replace_na(zscore, -Inf), iscore= replace_na(iscore, -Inf),
                   detected=(zscore >=zcut) & (iscore >=icut))
        return(samp %>% select(ljid, !!id:=detected) %>% distinct)
        
    }
    cols <- colnames(select(all_icov ,contains('H')))#hacky for now
    res <- lapply(cols, function(x) per_sample_score(x, zsc, intron_diff, scores))%>%  
    reduce(left_join) %>% left_join(select(zsc, ljid)) %>% 
        select(ljid, everything()) 
    return(res )
    
}
load('rdata/possible_fusion_genes.Rdata')
#fusions <- filter(possible_fusion_gens, ljid %in% exon_scoring$ljid)
exon_scoring <- rank_exons(all_zscores, all_icov, c(-1, .04), 'num_samps_detected') %>% 
    filter(!ljid %in% possible_fusion_gens$ljid) #%>% arrange(desc(num_samps_detected))
num_detected_per_samp <- colSums(exon_scoring[,-1])
med <- median(num_detected_per_samp)
stdv <- sd(num_detected_per_samp)
passed <- (num_detected_per_samp -med)/stdv > -2
keep <- which(passed) +1
msg <- paste0('Warning: ', sum(!passed), ' samples dropped from ', t_tissue)
print(msg)
exon_scoring <- exon_scoring[,c(1,keep)]
detect_co <- round((ncol(exon_scoring) -1)/3)
keep_exons <- rowSums(exon_scoring[,-1]) > detect_co
detected_exons <- exon_scoring[keep_exons,]

write_tsv(outfile)
