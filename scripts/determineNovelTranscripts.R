library(GenomicRanges)
library(tidyverse)
library(matrixStats)

# args <- c('/Volumes/data/eyeintegration_splicing/',
#           'results/all_tissues.combined.gtf',
#           'ref/gencodeAno_comp.gtf',
#           'sampleTableV4.tsv',
#           'results/all_tissues.PSI.tsv',
#           'results/all_tissues.incCts.tsv',
#           'results/salmon_gene_quant.Rdata',
#           'results/salmon_tx_quant.Rdata',
#           'testing/salmon_tissue_level_counts.Rdata',
#           'testing/novel_exon_expression_tables.Rdata')

args <- commandArgs(trailingOnly = T)
working_dir <- args[1]
gfc_gtf_file <- args[2]
ref_gtf <- args[3]
sample_table_file <- args[4]
psiTab_file <- args[5]
incCts_file <- args[6]
salmon_gene_quant<- args[7]
salmon_tx_quant <- args[8]
tissue_level_counts <- args[9]
exon_info_file <- args[10]

#save(args, file = 'testing/args_detnxtx.Rdata')


setwd(working_dir)
event_header <- c('ljid', 'seqid'	,'strand',	'start', 'end',	'event', 'length')
gfc_gtf <- rtracklayer::readGFF(gfc_gtf_file) %>% mutate(start=start-1)
txid_2_oid <- filter(gfc_gtf, type=='transcript') %>% select(transcript_id, oId, gene_id, gene_name)
txid_2_gname <- filter(gfc_gtf, type =='transcript') %>% select(transcript_id, gene_name) %>% 
    mutate(gene_name=replace(gene_name, is.na(gene_name),transcript_id[is.na(gene_name)]))
sample_table <- read_tsv(sample_table_file, col_names =c('sample','run','paired','tissue','subtissue','origin'))
psiTab <- read_tsv(psiTab_file) %>% mutate(ljid =paste('rm', 1:nrow(.), sep = '_'), length=end-start) %>% 
    select(ljid,seqid, strand, start, end,event,length, everything())


incCounts <- read_tsv(incCts_file) %>% mutate(ljid =paste('rm', 1:nrow(.),sep = '_'),length=end-start) %>% 
    select(ljid,seqid, strand, start, end,event,length, everything())
#old_inc_counts <-read_tsv('~/NIH/eyeintegration_splicing/results_b38/all_tissues.incCts.tsv')

load(salmon_tx_quant)
load(salmon_gene_quant)

make_tissue_level_counts <- function(st, exp_df){
    subtissues <- pull(st, subtissue) %>% unique
    med_cts <-  lapply(subtissues, function(x) filter(st, subtissue == x) %>% pull(sample) %>% 
                           exp_df[,.] %>% as.matrix %>% rowMedians) %>%
        do.call(cbind,.)
    colnames(med_cts) <- subtissues
    med_cts
}

med_tx_counts <- make_tissue_level_counts(sample_table, tx_counts) %>% as.data.frame %>% mutate(transcript_id=tx_counts$transcript_id) %>%
    select(transcript_id, everything(), -grep('synth', colnames(.)))
med_gene_counts <- make_tissue_level_counts(sample_table, gene_counts) %>% as.data.frame %>% mutate(GeneID=gene_counts$gene_id) %>% 
    select(GeneID, everything(), -grep('synth', colnames(.)))
save(med_tx_counts, med_gene_counts, file = tissue_level_counts)
###########################################################################
all_events_count_distributions <- incCounts %>% select(-event_header) %>% apply(2, function(x) quantile(x, seq(0,1,.1))) %>% t
#old_counts_distribution <-  old_inc_counts %>% select(contains('_incC')) %>% apply(2, function(x) quantile(x, seq(0,1,.1))) %>% t
pheatmap::pheatmap(log2(all_events_count_distributions[,-(1:6)]), cluster_rows = F, cluster_cols = F)

# based on this paper https://link.springer.com/article/10.1007%2Fs10709-007-9139-4 exons are not super big, so i think 
# removing exons longer than 1500 bp should be more than enough

rmats_minimally_expressed <- incCounts



#all_ref_exons
gencode_ref<-  rtracklayer::readGFF(ref_gtf) %>% mutate(start=start - 1)  %>% filter(type =='exon') %>% select(seqid, strand,start,end) %>% distinct %>% 
    mutate(seqid=as.character(seqid))
ensembl_ref <- rtracklayer::readGFF('ref/ensembl_r95.gtf') %>% mutate(start=start - 1)  %>% filter(type =='exon') %>% select(seqid, strand,start,end) %>% distinct %>% 
    mutate(seqid=as.character(seqid), seqid= paste0('chr', seqid))
refseq_ucsc <- rtracklayer::readGFF('ref/ucsc_refseq.gtf') %>% mutate(start=start-1, seqid=as.character(seqid)) #%>% 
gene2seq <- refseq_ucsc %>%  filter(type =='transcript') %>%  select(chroms=seqid, gene=gene_name) %>% distinct
refseq_ucsc <- refseq_ucsc %>% filter(type=='exon') %>% select(seqid, strand, start, end) %>% distinct
refseq_ncbi <- rtracklayer::readGFF('ref/refseq_r95.gff') %>% as.data.frame() %>% 
    left_join(gene2seq, by='gene') %>% filter(type=='exon') %>%  select(-seqid) %>% select(seqid=chroms, strand, start, end) %>% 
    filter(!is.na(seqid)) %>% distinct %>% mutate(start=start-1)



all_ref_exons <- rbind(gencode_ref, ensembl_ref, refseq_ncbi,refseq_ucsc ) %>% 
    mutate(seqid=as.character(seqid), strand=as.character(strand)) %>% distinct
save(all_ref_exons, file='rdata/all_ref_exons.Rdata')
#ref_exon_full <- split(all_ref_exons, 1:nrow(all_ref_exons))
# next, create a set of exons with some level of novelty
all_novel_exons <- filter(gfc_gtf,  type=='exon')  %>% 
    select(seqid, strand,start,end) %>% mutate(seqid=as.character(seqid)) %>% distinct

novel_string_tie_exons <- anti_join(all_novel_exons, all_ref_exons)

#this is all "novel" exons that were detected in rmats
novel_st_exons_in_rmats <- inner_join(novel_string_tie_exons, incCounts) %>%  select(seqid, strand, start, end, ljid) %>%  distinct


# now lets break these down a little more
ref_exon_starts <- all_ref_exons %>% select(seqid, strand, start) %>% distinct

ref_exon_ends <- all_ref_exons %>% select(seqid, strand, end) %>% distinct

nv_start_in_ref <- inner_join(novel_st_exons_in_rmats, ref_exon_starts) %>% {novel_st_exons_in_rmats$ljid %in% .$ljid }
nv_end_in_ref <- inner_join(novel_st_exons_in_rmats, ref_exon_ends) %>% {novel_st_exons_in_rmats$ljid %in% .$ljid}


fully_novel_exons <- (!nv_start_in_ref) & (!nv_end_in_ref)# truely novel exon has both start and end not in ref
a3ss_novel_exons <-  (nv_start_in_ref) & (!nv_end_in_ref)# a3ss novel exon has a known start, new end
a5ss_novel_exons <- (!nv_start_in_ref) & (nv_end_in_ref)# a5ss novel eoxn has new start, known end
ri_novel_exons <- (nv_start_in_ref) & (nv_end_in_ref)# not a known exon, but starts and ends on known exons
determine_overlap_Exons <-  function(df, ref){
    ref_range <- GRanges(seqnames = ref$seqid, ranges = IRanges(start = ref$start,end= ref$end), strand = ref$strand)
    nx_range <-  GRanges(seqnames = df$seqid, ranges = IRanges(start = df$start,end= df$end), strand = df$strand, 
                         names=df$ljid)
    k <- findOverlaps(nx_range, ref_range, select = 'first')# %>% as.data.frame %>% select(seqid=seqnames,strand, start, end, -width)
    return(tibble(ljid=df$ljid,ss_exon=!is.na(k) ))
    
}


novel_st_exons_in_rmats <- novel_st_exons_in_rmats %>%
    mutate(reclassified_event= case_when(fully_novel_exons ~ "novel_exon",
                                         a3ss_novel_exons ~ "A3SS",
                                         a5ss_novel_exons ~ "A5SS",
                                         ri_novel_exons ~ "RI")) %>% 
    #mutate(length=end -start) %>% filter(length <=1500)
    mutate(length= end -start, is.not_long=length <1500) 

novel_st_exons_in_rmats <- determine_overlap_Exons(novel_st_exons_in_rmats %>% filter(reclassified_event=='novel_exon') %>% distinct,all_ref_exons) %>% 
    left_join(novel_st_exons_in_rmats,.) %>%
    mutate(ss_exon=replace_na(ss_exon, F), reclassified_event = replace(reclassified_event, ss_exon, 'subset_exon'), ss_exon=NULL) %>% 
    inner_join( incCounts) %>% 
    inner_join(select(gfc_gtf, transcript_id, seqid, strand, start, end)) %>% 
    select(ljid, seqid, start, end, transcript_id, reclassified_event, everything())

###figure out skipped exon transcripts

gtf_se<- left_join(gfc_gtf, select(novel_string_tie_exons, seqid, strand, start, end)%>%distinct %>% mutate(is.novel=T), 
                   by=c('seqid','strand' ,'start','end')) %>%
    mutate(is.novel=replace(is.novel, is.na(is.novel), F))
skipped_exon_tx <- filter(gtf_se, grepl('MSTRG', oId)) %>% pull(transcript_id) %>% 
    {filter(gtf_se, transcript_id %in% .)} %>% group_by(transcript_id) %>% summarise(no_new_exon=all(!is.novel)) %>%
    filter(no_new_exon)
txid2gene <- filter(gfc_gtf, type=='transcript') %>% select(transcript_id, gene_name) %>% distinct %>% 
    mutate(gene_name=replace(gene_name, is.na(gene_name), transcript_id[is.na(gene_name)]))

nx_inc_cts <- novel_st_exons_in_rmats %>% left_join(txid2gene) %>% select(-contains('incC'), contains('incC'))
nx_psi <- inner_join( nx_inc_cts%>%select(transcript_id, reclassified_event, ljid, gene_name ), psiTab) %>% select(-contains('PSI'), contains('PSI'))
nx_skipped_exon <- skipped_exon_tx %>% left_join(txid2gene)

save(nx_inc_cts, nx_psi, nx_skipped_exon, file = exon_info_file )

# 
# gfc_gtf %>% filter(type == 'transcript') %>% mutate(sscore=999) %>%  select(seqid, start, end, transcript_id) %>% 
#      write_tsv('testing/st_tx.bed', col_names = F)
# gencode_comp_gtf <- rtracklayer::readGFF('ref/gencodeAno_comp.gtf')
# gencode_comp_gtf %>% filter(type =='gene',!grepl('RP\\d', gene_name), !grepl('MIR', gene_name), !grepl('RF', gene_name)) %>% 
#     mutate(score =999) %>% 
#     select(seqid, start, end, gene_name) %>% 
#     write_tsv('testing/gc_genes.txt', col_names = F)
#com='bedtools intersect -s -loj -a testing/st_tx.bed  -b testing/gc_genes.txt  > testing/overlap.bed'
#system2(command = com, wait = T)

# res <- read_tsv('testing/overlap.bed',col_names = F) %>% dplyr::rename(transcript_id=X4, gene_name =X8) %>% 
#     select(transcript_id, gene_name) %>% distinct 
# k <- res %>% group_by(transcript_id) %>% summarise(n=n()) %>% 
#     filter(n>1, transcript_id %in% nx_inc_cts$transcript_id)
# 
# filter(res, transcript_id %in% k$transcript_id) %>% inner_join(nx_inc_cts) -> possible_fusion_gens
# save(possible_fusion_gens, file ='~/NIH/eyeintegration_splicing/rdata/possible_fusion_genes.Rdata')

#GRanges(seqnames = all_tx$seqid, ranges = IRanges(start = all_tx$start, end = all_tx$end), strand = )




