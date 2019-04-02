library(tidyverse)
library(matrixStats)


args <- c('/Volumes/data/eyeintegration_splicing/',
         'results/all_tissues.combined.gtf',
         'sampleTableV4.tsv',
         'results/all_tissues.PSI.tsv',
         'results/all_tissues.incCts.tsv',
         'results/salmon_gene_quant.Rdata',
         'results/salmon_tx_quant.Rdata',
         'ref/gencodeAno_comp.gtf',
         'results/salmon_tissue_level_counts.Rdata',
         'results/as_event_ls_classV2.Rdata')

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
#all_events_count_distributions <- incCounts %>% select(-event_header) %>% apply(2, function(x) quantile(x, seq(0,1,.1))) %>% t
#old_counts_distribution <-  old_inc_counts %>% select(contains('_incC')) %>% apply(2, function(x) quantile(x, seq(0,1,.1))) %>% t
#k <- rbind(all_events_count_distributions, old_counts_distribution)
#pheatmap::pheatmap(k, cluster_rows = F, clustering_distance_cols = F)
#Heatmap(log2(k+1), cluster_rows = F, cluster_columns = F)
#length_distribution <- incCounts %>% pull(length) %>% quantile( seq(0,1,.1))
# seems like a threshold of counts >.25 might be good enough to remove out noise - .25 is ~top 25% of expression for each sample
# which is about what it was at before I length scaled the counts
# based on this paper https://link.springer.com/article/10.1007%2Fs10709-007-9139-4 exons are not super big, so i think 
# removing exons longer than 1500 bp should be more than enough
incCts_cut_off_lvl <- .25
salmon_cut_off_lvl <- 15

expressed_in_any_tissue <-  rowSums(med_tx_counts[,-1] > salmon_cut_off_lvl) >0
stringtie_min_exp <- med_tx_counts[expressed_in_any_tissue,]
stringtie_novel_min_exp <- filter(txid_2_oid,grepl('MSTRG',oId)) %>% pull(transcript_id) %>%
    {filter(stringtie_min_exp, transcript_id %in% .)}

mat <- select(incCounts, -event_header)
keep_cts <- rowSums(mat > incCts_cut_off_lvl) > 0
incCounts_filter <- incCounts[keep_cts,]
keep_len <- incCounts_filter$length >1500
rmats_minimally_expressed <- incCounts_filter[!keep_len,]

# rm <- select(rmats_minimally_expressed, seqid, strand, start, end, ljid)
# st <- filter(gfc_gtf, transcript_id %in% stringtie_min_exp$transcript_id, type == 'exon') %>%  select(seqid, strand, start, end) %>% distinct
# k <- anti_join(rm, st) %>% pull(ljid) %>% {filter(rmats_minimally_expressed, ljid %in% .)}

#####Second, determine which exons have been identified as novel; break these down further to reclassify events 


# create a set of reference exons by selecting all exons from reference transcripts.    


all_ref_exons <-  rtracklayer::readGFF(ref_gtf) %>% mutate(start=start - 1)  %>% filter(type =='exon') %>% select(seqid, strand,start,end) %>% distinct %>% 
    mutate(seqid=as.character(seqid))
#all_ref_exons <- rbind(all_ref_exons_gc, all_ref_exons_st) %>% mutate(seqid=as.character(seqid), strand=as.character(strand)) %>% distinct
ref_exon_full <- split(all_ref_exons, 1:nrow(all_ref_exons))

# next, create a set of exons with some level of novelty
all_novel_exons <- filter(gfc_gtf,  type=='exon')  %>% 
    select(seqid, strand,start,end) %>% mutate(seqid=as.character(seqid))
all_novel_exon_list <- split(all_novel_exons, 1:nrow(all_novel_exons))
novel <- !all_novel_exon_list%in%ref_exon_full
novel_string_tie_exons <- all_novel_exons[novel,] %>% distinct %>% mutate(is.novel=T)
#this is all "novel" exons that were detected in rmats

###It would be a good idea to record the  tx/exon that are in string tie and not in rmats, and look have a separate evaluation with
### different criteria, doing the same for rmats.
##############
# novel_st_exons_NOT_in_rmats <- anti_join(novel_string_tie_exons, rmats_minimally_expressed) %>%
#     left_join(select(gfc_gtf, transcript_id, seqid, strand, start, end))
# st <- filter(gfc_gtf, transcript_id %in% stringtie_min_exp$transcript_id, type == 'exon') %>%  select(seqid, strand, start, end) %>% distinct
# rmats_exons_NOT_in_st <- anti_join(rmats_minimally_expressed,st)
# save(desc,novel_st_exons_NOT_in_rmats, rmats_exons_NOT_in_st, file = 'results/novel_exons_not_agreeing.Rdata')
# rm(desc,novel_st_exons_NOT_in_rmats, rmats_exons_NOT_in_st)
##############


novel_st_exons_in_rmats <- inner_join(novel_string_tie_exons, rmats_minimally_expressed, by=c('seqid','strand', 'start','end'))

# now lets break these down a little more
ref_exon_starts <- split(all_ref_exons%>%select(seqid,strand, start), 1:nrow(all_ref_exons))
ref_exon_ends <- split(all_ref_exons%>%select(seqid,strand, end),1:nrow(all_ref_exons)) 

nv_ex_all <- split(novel_st_exons_in_rmats%>%select(seqid, strand, start, end), 1:nrow(novel_st_exons_in_rmats))
nv_ex_start <- split(novel_st_exons_in_rmats%>%select( seqid, strand, start), 1:nrow(novel_st_exons_in_rmats))
nv_ex_end <- split(novel_st_exons_in_rmats%>%select( seqid, strand,end), 1:nrow(novel_st_exons_in_rmats))

fully_novel_exons <- (!nv_ex_start%in%ref_exon_starts) & (!nv_ex_end%in%ref_exon_ends)# truely novel exon has both start and end not in ref
a3ss_novel_exons <-  (nv_ex_start%in%ref_exon_starts) & (!nv_ex_end%in%ref_exon_ends)# a3ss novel exon has a known start, new end
a5ss_novel_exons <- (!nv_ex_start%in%ref_exon_starts) & (nv_ex_end%in%ref_exon_ends)# a5ss novel eoxn has new start, known end
ri_novel_exons <- (nv_ex_start%in%ref_exon_starts) & (nv_ex_end%in%ref_exon_ends)# not a known exon, but starts and ends on known exons
k <- fully_novel_exons+a3ss_novel_exons+a5ss_novel_exons+ri_novel_exons
#?why so many retained introns? are the intronic from all the samples getting amplified reads from each sample 
novel_st_exons_in_rmats <- novel_st_exons_in_rmats %>%
    mutate(reclassified_event= case_when(fully_novel_exons ~ "novel_exon",
                                         a3ss_novel_exons ~ "A3SS",
                                         a5ss_novel_exons ~ "A5SS",
                                         ri_novel_exons ~ "RI"))%>%
    inner_join(select(gfc_gtf, transcript_id, seqid, strand, start, end)) %>% select(event_header, transcript_id, reclassified_event, everything())

###figure out skipped exon transcripts

gtf_se<- left_join(gfc_gtf, select(novel_string_tie_exons, seqid, strand, start, end)%>%distinct %>% mutate(is.novel=T), 
                   by=c('seqid','strand' ,'start','end')) %>%
    mutate(is.novel=replace(is.novel, is.na(is.novel), F))
skipped_exon_tx <- filter(gtf_se, grepl('MSTRG', oId)) %>% pull(transcript_id) %>% 
{filter(gtf_se, transcript_id %in% .)} %>% group_by(transcript_id) %>% summarise(no_new_exon=all(!is.novel)) %>%
    filter(no_new_exon)
nx_inc_cts <- novel_st_exons_in_rmats
nx_psi <- inner_join( nx_inc_cts%>%select(transcript_id, reclassified_event, is.novel, ljid ), psiTab)
nx_skipped_exon <- skipped_exon_tx
save(nx_inc_cts, nx_psi, nx_skipped_exon, file = exon_info_file )

# #determine novel transcript expression within the eye
# eye_tissues <- grepl('Retina|RPE|Cornea', colnames(stringtie_min_exp))
# expressed_in_eye <-  rowSums(stringtie_min_exp[,eye_tissues] > salmon_cut_off_lvl) >=1
# stringtie_expressed_in_eye <- stringtie_min_exp[expressed_in_eye,]
# stringtie_novel_eye <- filter(txid_2_oid, grepl('MSTRG',oId)) %>% pull (transcript_id) %>% 
#     {filter(stringtie_expressed_in_eye, transcript_id %in% .)}
# 
# #determine exon expression in the eye
# eye_cols <- grepl('Retina|RPE|Cornea', colnames(rmats_minimally_expressed))
# expressed_in_eye <- rowSums(rmats_minimally_expressed[,eye_cols] > incCts_cut_off_lvl) >=1
# rmats_expressed_eye <- rmats_minimally_expressed[expressed_in_eye,]
# #keep novel tx/exons that are expressed in the eye
# novel_st_exons_in_rmats_eye <- filter(novel_st_exons_in_rmats, transcript_id %in% stringtie_expressed_in_eye$transcript_id, 
#                                       ljid%in% rmats_expressed_eye$ljid)
# 
# 
# 
# 
# 
# nc_col_order=colnames(novel_st_exons_in_rmats_eye)[1:10]
# nx_inc_cts_eye <- novel_st_exons_in_rmats_eye %>% select(nc_col_order, everything())
# nx_psi_eye <- inner_join( nx_inc_cts_eye%>%select(transcript_id, reclassified_event, is.novel, ljid ), incTab) %>% 
#     select(nc_col_order, everything())
# nx_tx_exp_eye <- filter(med_tx_counts, transcript_id%in%nx_inc_cts_eye$transcript_id)
# nx_tx_skipped_exon_tx_eye <- filter(stringtie_expressed_in_eye, transcript_id%in% skipped_exon_tx$transcript_id)
# 
# gtf_annotated <- left_join(gfc_gtf, select(novel_st_exons_in_rmats_eye, seqid, strand, start, end, reclassified_event)%>%distinct, 
#                            by=c('seqid','strand' ,'start','end')) %>%
#     mutate(reclassified_event=replace(reclassified_event, is.na(reclassified_event), 'ref'))
# 
# 
# gtf_novel_exon_anno <- filter(gtf_annotated, transcript_id %in% nx_tx_skipped_exon_tx_eye$transcript_id | 
#                                              transcript_id%in% nx_inc_cts_eye$transcript_id)
# write_tsv(gtf_novel_exon_anno, 'novel_exon_gtf.tsv')
# save(nx_psi_eye, nx_inc_cts_eye, nx_tx_exp_eye, nx_tx_skipped_exon_tx_eye, gtf_novel_exon_anno, gtf_annotated, file = exon_info_file)
# 
# #rtracklayer::readGFF('ref/gencodeGFF3.gff') %>% as.data.frame %>% write_tsv('ref/gencodeGFF3.gff.tsv')
