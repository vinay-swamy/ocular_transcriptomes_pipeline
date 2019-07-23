library(tidyverse)
library(matrixStats)

# args <- c( '/Volumes/data/eyeintegration_splicing/',
#            'results/novel_exon_workspace.Rdata',
#            'results/possible_fusion_genes.Rdata',
#            'RPE_Fetal.Tissue',
#            'sampleTableV4.tsv',
#            'results/all_tissues.combined.gtf',
#            'testing/xhmooby.ranked.tsv'
#          )

args <- commandArgs(trailingOnly = T)
working_dir <- args[1]
exon_info_workspace <- args[2]
fusion_gene_file <- args[3]
t_tissue <- args[4]
sample_tab <- args[5]
gtf_file <- args[6]
outfile <- args[7]

setwd(working_dir)
load(exon_info_workspace)
load('results/novel_exon_expression_tables.Rdata')
ref_exon_table <- ref_exon_tab
intron_table <- intron_info_tab %>%  
     select(-tx, ljid=id) %>% distinct
rm(intron_info_tab, ref_exon_tab)
sample_table <- read_tsv(sample_tab, col_names = c('sample', 'run', 'paired', 'tissue', 'subtissue', 'origin'))
samples <- filter(sample_table, subtissue == t_tissue, paired =='y') %>% pull(sample)
cov_paths <- paste0('coverage_files/', t_tissue,'/', samples, '.regions.bed.gz')#
exists <-  file.exists(cov_paths)
samples <- samples[exists]
cov_paths <- cov_paths[exists]
gfc_gtf <- rtracklayer::readGFF(gtf_file) %>% mutate(start=start-1) 
t2g <- gfc_gtf %>% filter(type =='transcript') %>% select(transcript_id, gene_name) %>% 
    mutate(gene_name=replace(gene_name, is.na(gene_name),transcript_id[is.na(gene_name)]))
l2g <- nx_inc_cts %>% select(ljid, gene_name) %>% distinct 

novel_exon_locs <- nx_inc_cts %>% filter(reclassified_event == 'novel_exon', is.not_long) %>%   
    select(seqid,strand, start, end, reclassified_event, gene_name, ljid) %>% distinct
long_locs <- ass_tab %>% filter(class =='long') %>% left_join(l2g) %>% 
    select(seqid,strand, start=l_start, end=l_end, reclassified_event, gene_name, ljid)
ref_exon_locs <- incCounts %>% filter(!ljid %in% nx_inc_cts$ljid) %>% mutate(reclassified_event= 'ref') %>% 
    inner_join(gfc_gtf %>% filter(type =='exon') %>% select(seqid, strand, start, end, transcript_id)) %>% 
    left_join(t2g) %>% 
    select(seqid,strand, start, end, reclassified_event, gene_name, ljid) %>% distinct

#novel_loci_locs <- novel_loci_tab %>% mutate(reclassified_event='novel_loci', gene_name=transcript_id) %>% 
#    select(seqid,strand, start, end, reclassified_event, gene_name, ljid) %>% distinct
    
all_exon_locs <- rbind(novel_exon_locs, long_locs, ref_exon_locs)#, novel_loci_locs)

#intron_table <- rbind(intron_table, novel_loci_intron_tab %>% select(seqid, start, end,ljid=id, strand)) %>% distinct

#long_exon_locs <- nx_inc_cts %>% filter(ljid %in% long_locs$ljid) %>% select(seqid, strand, start, end, ljid, gene_name) %>% distinct



all_cov <- lapply(cov_paths, function(epath)read_tsv(epath, col_names = c('seqid', 'start', 'end', 'cov')) %>% 
                      mutate(length=end-start, ls_cov=cov/length, cov=NULL))


mean_based_zscore <- function( cov_df, sample_name, novel_exon_locs, ref_exon_table){
    #note - ref tab is generated from eyeIntegration data, which is missing some of the genes, we need to make a tpms table
    # on our own from the ref, lossing about 8k stuff it looks like
    novel_exon_cov <- novel_exon_locs %>% inner_join(cov_df)
    
    ref_exon_cov <- ref_exon_table %>% select(seqid, strand, start, end, gene_name) %>% distinct %>% inner_join(cov_df)
    
    ref_gene_dists <- ref_exon_cov %>% group_by(gene_name) %>% summarise(g_mean=mean(ls_cov), g_med=median(ls_cov), g_sd=sd(ls_cov))
    
    all_cov <- inner_join(novel_exon_cov, ref_gene_dists) %>%  # #filter(ls_cov != 0) %>%
        mutate(zscore=(ls_cov - g_mean)/g_sd) %>%
        select(-length, -contains('g_'),-ls_cov, -reclassified_event, !!sample_name := zscore)
    
    all_cov
}



intron_cov <- function(cov_df ,sample_name, novel_exon_locs,  intron_bed){
    icov <- cov_df  %>% inner_join(intron_bed, .) %>%  distinct %>% 
        group_by(ljid) %>% summarise(intron_cov=mean(ls_cov))
    novel_exon_cov <- novel_exon_locs %>% inner_join(cov_df) 
    
    a_cov <- icov %>% 
        left_join(novel_exon_cov, .) %>% mutate( !!sample_name := ls_cov -intron_cov) %>% select(-length, -ls_cov, -intron_cov)
    
    a_cov
    
}


all_zscores <- lapply(1:length(all_cov), 
                      function(i) suppressMessages(mean_based_zscore(all_cov[[i]], samples[i],all_exon_locs, ref_exon_table))) %>% 
    reduce( full_join)



all_icov <- lapply(1:length(all_cov), function(i) 
    suppressMessages(intron_cov(all_cov[[i]], samples[i], all_exon_locs, intron_table) %>% distinct)) %>% reduce(left_join) 

####CHANGE ME BACK************************************************************
##########FUUUUUUUUUUUUUUK some how the ljids got boofed up, hopefully will get fixed in new run

all_zscores <- all_zscores %>% filter(!duplicated(.[,c('seqid', 'strand', 'start','end')]))
all_icov <- all_icov %>% filter(!duplicated(.[,c('seqid', 'strand', 'start','end')]))

all_zscores <-  all_icov %>% select(seqid, strand, start, end) %>% inner_join(all_zscores,.)
all_icov <- all_zscores %>% select(seqid, strand, start, end) %>% inner_join(all_icov,.)
####CHANGE ME BACK************************************************************
rank_exons <-  function(zsc, intron_diff, scores, t_name){
    per_sample_score <- function(id, zsc, intron_diff, scores){
        zcut=scores[[1]]
        icut=scores[[2]]
        samp <- inner_join(zsc %>% select( seqid, strand, start, end,ljid, gene_name, zscore=!!id),
                           intron_diff %>% select( seqid, strand, start, end, ljid, gene_name, iscore=!!id)) %>%
            mutate(zscore= replace_na(zscore, -Inf), iscore= replace_na(iscore, -Inf),
                   detected=(zscore >=zcut) & (iscore >=icut))
        return(samp %>% select(seqid, strand,start, end, ljid, gene_name,  !!id:=detected) %>% distinct)
        
    }
    cols <- colnames(select(all_icov ,contains('H'), contains ('SRS')))#hacky for now
    res <- lapply(cols, function(x) per_sample_score(x, zsc, intron_diff, scores))%>%  
    reduce(left_join) %>% #left_join(select(zsc, ljid)) %>% 
        select(ljid, everything()) 
    return(res)
    
}
load(fusion_gene_file)
#fusions <- filter(possible_fusion_gens, ljid %in% exon_scoring$ljid)
exon_scoring <- rank_exons(all_zscores, all_icov, c(-1, .04), 'num_samps_detected') %>% 
    filter(!ljid %in% possible_fusion_gens$ljid) %>% distinct #%>% arrange(desc(num_samps_detected))
nc_cols <- c('seqid','strand', 'start', 'end', 'ljid', 'gene_name')
num_detected_per_samp <- exon_scoring %>% select(-nc_cols) %>%  colSums() 
med <- median(num_detected_per_samp)
stdv <- sd(num_detected_per_samp)
passed <- (num_detected_per_samp -med)/stdv > -2
keep <- which(passed) +length(nc_cols)
msg <- paste0('Warning: ', sum(!passed), ' samples dropped from ', t_tissue)
print(msg)
exon_scoring <- exon_scoring %>% select(nc_cols, keep)
detect_co <- trunc((ncol(exon_scoring %>% select(-nc_cols)))/3)
keep_exons <- exon_scoring %>% select(-nc_cols) %>%  {rowSums(.) > detect_co}
detected_exons <- exon_scoring[keep_exons,] %>% inner_join(select(all_exon_locs, ljid, reclassified_event)) %>% 
    select(nc_cols, reclassified_event, everything()) %>% distinct %>% mutate(total= select(., -nc_cols, -reclassified_event) %>% rowSums())



write_tsv(detected_exons, outfile)

