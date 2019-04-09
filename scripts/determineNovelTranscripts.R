library(tidyverse)
library(matrixStats)

args <- c('/Volumes/data/eyeintegration_splicing/',
          'results/all_tissues.combined.gtf',
          'ref/gencodeAno_comp.gtf',
          'sampleTableV4.tsv',
          'results/all_tissues.PSI.tsv',
          'results/all_tissues.incCts.tsv',
          'results/salmon_gene_quant.Rdata',
          'results/salmon_tx_quant.Rdata',
          'testing/salmon_tissue_level_counts.Rdata',
          'testing/novel_exon_expression_tables.Rdata')

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
process_refSeq_ano <- function(df){
    acc <- df %>% select(seqid, chromosome) %>% filter(!is.na(chromosome), chromosome != 'Unknown') %>% 
        mutate(chr=paste0('chr',chromosome))
    bed <- df %>% filter(type == 'exon') %>% select(seqid, strand, start, end) %>% inner_join(acc) %>%
        select(-chromosome,-seqid ) %>% 
        select(seqid=chr, strand, start, end) %>% distinct %>% mutate(start =start-1)
    bed
}



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
ensembl_ref <- rtracklayer::readGFF('~/NIH/eyeintegration_splicing/ref/ensembl_r95.gtf') %>% mutate(start=start - 1)  %>% filter(type =='exon') %>% select(seqid, strand,start,end) %>% distinct %>% 
    mutate(seqid=as.character(seqid), seqid= paste0('chr', seqid))
refseq_ncbi_ref <- process_refSeq_ano(rtracklayer::readGFF('~/NIH/eyeintegration_splicing/ref/refseq_r95.gff') %>% as.data.frame)
refseq_ucsc <- rtracklayer::readGFF('ref/ucsc_refseq.gtf') %>% mutate(start=start-1, seqid=as.character(seqid)) %>% select(seqid, strand, start, end) %>% distinct

all_ref_exons <- rbind(gencode_ref, ensembl_ref, refseq_ncbi_ref,refseq_ucsc ) %>% 
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

#?why so many retained introns? are the intronic from all the samples getting amplified reads from each sample 
novel_st_exons_in_rmats <- novel_st_exons_in_rmats %>%
    mutate(reclassified_event= case_when(fully_novel_exons ~ "novel_exon",
                                         a3ss_novel_exons ~ "A3SS",
                                         a5ss_novel_exons ~ "A5SS",
                                         ri_novel_exons ~ "RI"))  %>% 
    inner_join( incCounts) %>%  mutate(length= end -start, is.not_long=length <1500) %>% 
    inner_join(select(gfc_gtf, transcript_id, seqid, strand, start, end)) %>% select(ljid, seqid, start, end, transcript_id, reclassified_event, everything())

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
