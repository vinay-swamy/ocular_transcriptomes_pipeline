library(tidyverse)
library(matrixStats)
setwd('~/NIH/eyeintegration_splicing/')
load('/Volumes/data/eyeintegration_splicing/testing/novel_exon_expression_tables.Rdata')
ref_exon_table <- read_tsv('/Volumes/data/eyeintegration_splicing/results/ref_exon_table.tsv')
ref_bed <- read_tsv('/Volumes/data/eyeintegration_splicing/results/exons_for_coverage_analysis.bed', col_names = c('seqid', 'start', 'end'))
intron_bed <- read_tsv('/Volumes/data/eyeintegration_splicing/testing/intron_info_tab.tsv')


t_tissue <- 'RPE_Fetal.Tissue'
sample_table <- read_tsv('sampleTableV4.tsv', col_names = c('sample', 'run', 'paired', 'tissue', 'subtissue', 'origin'))
samples <- filter(sample_table, subtissue == t_tissue, paired =='y') %>% pull(sample)
exon_cov_paths <- paste0('~/NIH/eyeintegration_splicing/coverage_files/', t_tissue,'/', samples, '.regions.bed.gz') %>% .[ file.exists(.)]
#intron_cov_paths <- paste0('~/NIH/eyeintegration_splicing/coverage_files/', t_tissue,'/', samples, '.intron_per_gene.bed.gz')  %>% .[ file.exists(.)] 
intron_cov_paths <- paste0('~/NIH/eyeintegration_splicing/rpe_fetal_coverage/',samples, '.regions.bed.gz')
t2g <- rtracklayer::readGFF('/Volumes/data/eyeintegration_splicing/results/all_tissues.combined.gtf') %>% filter(type =='transcript') %>% select(transcript_id, gene_name) %>% 
    mutate(gene_name=replace(gene_name, is.na(gene_name),transcript_id[is.na(gene_name)]))
novel_exon_locs <- nx_inc_cts %>% filter(reclassified_event == 'novel_exon', is.not_long) %>%   select(seqid,strand, start, end, reclassified_event, gene_name, ljid) %>% distinct
novel_exon_locs %>% select(seqid, strand, start, end, gene_name) %>% distinct %>% dim


mean_based_zscore <- function( epath, sample_name, novel_exon_locs, ref_exon_table){
    ecov <- read_tsv(epath, col_names = c('seqid', 'start', 'end', 'cov')) %>% 
        mutate(length=end-start, ls_cov=cov/length, cov=NULL)
    
    
    novel_exon_cov <- novel_exon_locs %>% inner_join(ecov) 
    
    ref_exon_cov <- ref_exon_table %>% select(seqid, strand, start, end, gene_name) %>% distinct %>% inner_join(ecov)
    
    ref_gene_dists <- ref_exon_cov %>% group_by(gene_name) %>% summarise(g_mean=mean(ls_cov), g_med=median(ls_cov), g_sd=sd(ls_cov))
    
    all_cov <- inner_join(novel_exon_cov, ref_gene_dists) %>% filter(ls_cov != 0) %>%  mutate(zscore=(ls_cov - g_mean)/g_sd) %>%
        select(-length, -contains('g_'),-ls_cov, -reclassified_event, !!sample_name := zscore)
    
    all_cov
}

intron_cov <- function(ipath, epath ,sample_name,  intron_bed){
    ecov <- read_tsv(epath, col_names = c('seqid', 'start', 'end', 'cov')) %>% 
        mutate(length=end-start, ls_cov=cov/length, cov=NULL)
    icov <- read_tsv(ipath, col_names = c('seqid', 'start', 'end', 'cov')) %>% 
        mutate(length=end-start, ls_cov=cov/length, cov=NULL) %>% left_join(intron_bed) %>% select(-tx) %>% distinct #%>% 
    #rename(!!sample_name := ls_cov)
    novel_exon_cov <- novel_exon_locs %>% inner_join(ecov) 
    
    a_cov <- icov %>% select( ljid=id, ls_cov) %>% group_by(ljid) %>% summarise(intron_cov= mean(ls_cov)) %>% 
        left_join(novel_exon_cov, .) %>% mutate( !!sample_name := ls_cov -intron_cov) %>% select(-length, -ls_cov, -intron_cov)
    
    a_cov
    
}
all_zscores <- lapply(3:length(exon_cov_paths), 
                      function(i) suppressMessages(mean_based_zscore(exon_cov_paths[i], samples[i],novel_exon_locs, ref_exon_table))) %>% 
    reduce( full_join)

all_icov <- lapply(3:length(intron_cov_paths), function(i) 
    suppressMessages(intron_cov(intron_cov_paths[i], exon_cov_paths[i], samples[i], intron_bed) %>% distinct)) %>% reduce(left_join) %>% 
    filter(!is.na(.[,8]))
common_ids <- intersect(all_icov$ljid, all_zscores$ljid)

all_zscores <- all_zscores %>% filter(ljid %in%common_ids)
all_icov <- all_zscores %>% filter(ljid %in%common_ids)


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
        select(ljid, everything()) %>% 
        mutate(total_detected=rowSums(.[,-1]))
    
    return(res %>% select(ljid, !!t_name:= total_detected) %>% distinct)
    
}

score_list <- list( c(0,.5),c(-.25,.25), c(-.5,.1),c(-.75,.08),c(-1,.06),c(-1.25,.04), c(-1.4,.01), c(-1.6,0), )
names(score_list) <- c('a','b','c','d', 'e', 'f','g','h')
exon_scoring <- lapply(1:length(score_list), 
            function(i) rank_exons(all_zscores, all_icov, score_list[[i]], names( score_list)[i]) %>% filter(!duplicated(ljid))) %>% 
    reduce(left_join) %>% mutate(total= rowSums(.[,-1])) %>% filter(total >=0) %>% arrange(desc(total)) 
load('rdata/possible_fusion_genes.Rdata')
filter(possible_fusion_gens, ljid %in% exon_scoring$ljid) -> fusions

keep <- filter(exon_scoring, !ljid %in% fusions$ljid) %>% arrange(desc(a)) %>% select(-total) %>% 
    inner_join(select(novel_exon_locs, gene_name, ljid))


i='rm_123631'
tx=filter(nx_inc_cts, ljid ==i) %>% pull(transcript_id) %>% .[1]
makeMultiBeds(gtf = gtf_plotting, tx=tx, nv=i, tissue=t_tissue)

