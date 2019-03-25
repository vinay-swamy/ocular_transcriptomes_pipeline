library(tidyverse)
setwd('~/NIH/eyeintegration_splicing/')

tissue <- 'RPE_Fetal.Tissue'
load('results_b38/all_novel_exon_infoV2.Rdata')
load('results_b38/as_event_ls_class.Rdata')
bed <- read_tsv('/Volumes/data/eyeintegration_splicing/results/exons_for_cov_analysis_rpe.bed')
ascot_table <- read_tsv('ascot_exon_locations_formatted.tsv') %>% mutate(in.ascot=T)
retnet <- '/Volumes/data/eyeintegration_splicing/ref/retnet_hgncIDs_2017-03-28.txt'
load('limma_deg_RPE_fetal.Rdata')
ret_net_genes <- scan(retnet, character(), '\n')
mgog_gene <- scan('ref/mgog_genes.txt', character(), '\n') %>% unique
makeMultiBeds<- function(gtf,tx, nv, tissue){
    filter(gtf,transcript_id==tx, type=='exon') %>% 
        mutate(rgb=ifelse(reclassified_event != 'ref', '250,0,0', '0,0,250'),score=1000, ts=start, te=end ) %>% 
        select(seqid, start, end, transcript_id, score, strand,ts, te, rgb) -> nov_bed
    ref_txs <- filter(gtf, transcript_id==tx) %>% pull(gene_id) %>% unique %>% {filter(gtf, gene_id ==. , grepl('ENST', oId))} %>% 
        pull(transcript_id) %>% unique
    gene_name= filter(gtf,transcript_id==tx) %>% pull(gene_name) %>% .[1] 
    name <- paste0('beds2/',gene_name, '_',tissue,'_exons_bed.bed')
    #rbind(filter(bed, id == ref), filter(bed, id == nv)) %>% mutate(score=1000, strand=nov_bed$strand[1], ts=start, te=end, rgb=c( '0,0,250', '0,250,0')) %>% write_tsv('regions.bed', col_names = F)
    write_tsv(nov_bed,name, col_names = F)
    
    for( tx in ref_txs ){
        r_tx <- filter(gtf, transcript_id == tx) %>% pull(oId) %>% .[1]
        bed <- filter(gtf, transcript_id ==tx, type=='exon') %>% mutate(transcript_id=r_tx,score=1000, ts=start, te=end, rgb='0,0,250') %>% 
            select(seqid, start, end, transcript_id, score, strand,ts, te, rgb)
        name=paste0('beds2/', gene_name,'_', r_tx, '.bed')
        write_tsv(bed, name , col_names = F)
    }
}
files=list.files(path = '/Volumes/data/eyeintegration_splicing/rpe_fetal_coverage/', recursive = T, pattern = 'regions.bed.gz', full.names = T) %>% 
    .[grep('.csi', . , invert = T)]
names <- strsplit(files,'/|\\.') %>% sapply(function(x) x[[7]])
gtf_full <- rtracklayer::readGFF('results_b38/all_tissues.combined.gtf') 
gtf_full_anno <- gtf_full %>% filter(grepl('ENST', oId)) %>% pull(transcript_id) %>% 
    {filter(gtf_full, transcript_id %in% .) } %>%  mutate(reclassified_event='ref') %>% rbind( gtf_novel_exon_anno)


all_coverage_tabs <- lapply(1:length(files), function(i) 
                            read_tsv(files[i],col_names =  c('seqid', 'start','end','id',paste0(names[i],'_cov'))))
#there's dups in the bed data frame, i should really figure out
all_coverage_info<- reduce(all_coverage_tabs, .f = function(x,y) left_join(x,y, by= c('seqid', 'start','end','id')) %>% distinct) %>% 
    mutate(length=end-start) %>% select(c('seqid', 'start','end','id'), length, everything())
ls_cov_info <- select(all_coverage_info, -seqid, -start,-end ,-id) %>% { . / all_coverage_info$length} %>% as.matrix %>% # that ./ is on purpoise
    { matrixStats::rowMedians(.)} %>% {data.frame(select(all_coverage_info, seqid, start, end, id), med_cov= .)}

ljid2genename <- select(nx_med_cts_eye, ljid, transcript_id) %>% inner_join(select(gtf_novel_exon_anno, transcript_id, gene_name)) %>% 
    filter(!is.na(gene_name))


novel_exon_cov <-ls_cov_info %>% mutate(length=end-start) %>% 
    filter( id %in% nx_all$ljid, length <1000) %>% 
    rename(ljid=id) %>%  
    inner_join(ljid2genename) %>% 
    select(-med_cov, med_cov) %>% arrange(desc(med_cov)) %>% 
    left_join(ascot_table %>% select(seqid, start,end, in.ascot)) %>% mutate(in.ascot=replace_na(in.ascot, F)) %>% 
    left_join(select(limma_deg_RPE_fetal, transcript_id,adj.P.Val )) %>% left_join(select(nx_tx_exp_eye, transcript_id, tpm=RPE_Fetal.Tissue))

salmon_cut_off_lvl <- 15
novel_exon_cov_mgog <- filter(novel_exon_cov, gene_name %in% mgog_gene, tpm >salmon_cut_off_lvl)
write_tsv(novel_exon_cov_mgog, 'novel_exons_mgog_list.tsv')
#novel_exon_cov_ret_net <- filter(novel_exon_cov, gene_name %in% ret_net_genes)
for( id in novel_exon_cov_mgog$ljid){
    tx=filter(nx_med_cts_eye, ljid==id) %>% pull(transcript_id) %>% .[1]
    makeMultiBeds(gtf_full_anno, tx, id, tissue) 
}
remove<- novel_exon_cov_ret_net$ljid[!novel_exon_cov_ret_net$ljid %in% 'rm_6739' ]
novel_exon_cov <- filter(novel_exon_cov, !ljid %in% remove)

genes_to_show <- clipr::read_clip()
exons_to_show <- filter(novel_exon_cov, gene_name %in% genes_to_show ) %>% select(seqid, start, end, ljid, gene_name) %>% distinct


for( id in exons_to_show$ljid){
    tx=filter(nx_med_cts_eye, ljid==id) %>% pull(transcript_id) %>% .[1]
    makeMultiBeds(gtf_full_anno, tx, id, tissue)
}




