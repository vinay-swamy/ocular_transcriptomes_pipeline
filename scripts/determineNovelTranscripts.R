library(tidyverse)
library(matrixStats)

# args <- c('~/NIH/eyeintegration_splicing/',
#           'results_b38/all_tissues.combined.gtf',
#           'sampleTableV4.tsv',
#           'results_b38/all_rmats_events_tissues.incLevels.tsv',
#           'results_b38/all_rmats_events_tissues.medCounts.tsv',
#           'rdata/salmon_no_ssn_gene_quant.Rdata',
#           'rdata/salmon_no_ssn_tx_quant.Rdata',
#           'results_b38/salmon_tissue_level_counts.Rdata',
#           'results_b38/as_event_ls_class.Rdata')

args <- commandArgs(trailingOnly = T)

working_dir <- args[1]
gfc_gtf_file <- args[2]
sample_table_file <- args[3]
incTab_file <- args[4]
medCts_file <- args[5]
salmon_gene_quant<- args[6]
salmon_tx_quant <- args[7]
tissue_level_counts <- args[8]
exon_info_file <- args[9]



event_header <- c('ljid', 'GeneID', 'seqid'	,'strand',	'start', 'end',	'upstreamES',	'upstreamEE',	'downstreamES',	'downstreamEE')
gfc_gtf <- rtracklayer::readGFF(gfc_gtf_file) %>% mutate(start=start-1)
txid_2_oid <- filter(gfc_gtf, type=='transcript') %>% select(transcript_id, oId)
sample_table <- read_tsv(sample_table_file, col_names =c('sample','run','paired','tissue','subtissue','origin'))
incTab <- read_tsv(incTab_file) %>% mutate(ljid =paste('rm', 1:nrow(.), sep = '_')) %>% select(ljid, everything(), -grep('synth', colnames(.)))
medCounts <- read_tsv(medCts_file) %>% mutate(ljid =paste('rm', 1:nrow(.),sep = '_')) %>% select(ljid, everything(), -grep('synth', colnames(.)))
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
#####first, filter out salmon counts and rmats counts based on expression
#set medct levels at 25, seems like its a reasonable level, and salmon level to 15,based on looking at expression quantiles,
# also going to reduce the inc levle cuttoff to .5, at least for this list
txid_2_oid <- filter(gfc_gtf, type=='transcript') %>% select(transcript_id, oId)
salmon_cut_off_lvl <- 15
medcts_cut_off_lvl <- 25
#inc_cut_off_lvl=.5 lets just ignore inc levels entirely for now
# stringtie expressed in any tissue
expressed_in_any_tissue <-  rowSums(med_tx_counts[,-1] > salmon_cut_off_lvl) >0
stringtie_min_exp <- med_tx_counts[expressed_in_any_tissue,]
stringtie_novel_min_exp <- filter(txid_2_oid,grepl('MSTRG',oId)) %>% pull(transcript_id) %>% 
    {filter(stringtie_min_exp, transcript_id %in% .)}
#rmats rminimally expressed
nc_cols <- c('ljid',event_header,'event_type')
mat=medCounts%>%select(-nc_cols)
highly_expressed_exons <- rowSums(mat > medcts_cut_off_lvl) >0
rmats_minimally_expressed<- medCounts[highly_expressed_exons,]


#####Second, determine which exons have been identified as novel; break these down further to reclassify events 


# create a set of reference exons by selecting all exons from reference transcripts.    
all_ref_exons <- filter(gfc_gtf, grepl('ENST', oId)) %>% pull(transcript_id) %>% 
                    {filter(gfc_gtf, transcript_id %in% . , type=='exon')} %>% select(seqid, strand,start,end) %>% 
    mutate(seqid=as.character(seqid))
ref_exon_full <- split(all_ref_exons[,-1], 1:nrow(all_ref_exons))

# next, create a set of exons with some level of novelty
all_novel_exons <- filter(gfc_gtf, transcript_id %in% stringtie_novel_min_exp$transcript_id, type=='exon')  %>% 
    select( seqid, strand,start,end)
all_novel_exon_list <- split(all_novel_exons[,-1], 1:nrow(all_novel_exons))
novel <- !all_novel_exon_list%in%ref_exon_full

#this table is all novel exons
novel_string_tie_exons <- all_novel_exons[novel,] %>% distinct %>% mutate(is.novel=T)
#this is all "novel" exons that were detected in rmats
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

#?why so many retained introns? are the intronic from all the samples getting amplified reads from each sample 
novel_st_exons_in_rmats <- novel_st_exons_in_rmats %>%
    mutate(reclassified_event= case_when(fully_novel_exons ~ "novel_exon",
                                         a3ss_novel_exons ~ "A3SS",
                                         a5ss_novel_exons ~ "A5SS",
                                         ri_novel_exons ~ "RI")) %>% 
    inner_join(select(gfc_gtf, transcript_id, seqid, strand, start, end))




#determine novel transcript expression within the eye
eye_tissues <- grepl('Retina|RPE|Cornea', colnames(stringtie_min_exp))
expressed_in_eye <-  rowSums(stringtie_min_exp[,eye_tissues] > salmon_cut_off_lvl) >=1
stringtie_expressed_in_eye <- stringtie_min_exp[expressed_in_eye,]
stringtie_novel_eye <- filter(txid_2_oid, grepl('MSTRG',oId)) %>% pull (transcript_id) %>% 
    {filter(stringtie_expressed_in_eye, transcript_id %in% .)}

#determine exon expression in the eye
eye_cols <- grepl('Retina|RPE|Cornea', colnames(rmats_minimally_expressed))
expressed_in_eye <- rowSums(rmats_minimally_expressed[,eye_cols] > medcts_cut_off_lvl) >=1
rmats_expressed_eye <- rmats_minimally_expressed[expressed_in_eye,]
#keep novel tx/exons that are expressed in the eye
novel_st_exons_in_rmats_eye <- filter(novel_st_exons_in_rmats, transcript_id %in% stringtie_expressed_in_eye$transcript_id, 
                                      ljid%in% rmats_expressed_eye$ljid)

###figure out skipped exon transcripts

gtf_se<- left_join(gfc_gtf, select(novel_string_tie_exons, seqid, strand, start, end)%>%distinct %>% mutate(is.novel=T), 
                   by=c('seqid','strand' ,'start','end')) %>%
    mutate(is.novel=replace(is.novel, is.na(is.novel), F))
skipped_exon_tx <- filter(gtf_se, grepl('MSTRG', oId)) %>% pull(transcript_id) %>% 
{filter(gtf_se, transcript_id %in% .)} %>% group_by(transcript_id) %>% summarise(no_new_exon=all(!is.novel)) %>%
    filter(no_new_exon)



nc_col_order=c(colnames(novel_st_exons_in_rmats_eye)[1:12],'event_type', 'reclassified_event')
nx_med_cts_eye <- novel_st_exons_in_rmats_eye %>% select(nc_col_order, everything())
nx_incLvl_eye <- inner_join( nx_med_cts_eye%>%select(transcript_id, reclassified_event, is.novel, ljid ), incTab) %>% 
    select(nc_col_order, everything())
nx_tx_exp_eye <- filter(med_tx_counts, transcript_id%in%nx_med_cts_eye$transcript_id)
nx_tx_skipped_exon_tx_eye <- filter(stringtie_expressed_in_eye, transcript_id%in% skipped_exon_tx$transcript_id)

gtf_annotated <- left_join(gfc_gtf, select(novel_st_exons_in_rmats_eye, seqid, strand, start, end, reclassified_event)%>%distinct, 
                           by=c('seqid','strand' ,'start','end')) %>%
    mutate(reclassified_event=replace(reclassified_event, is.na(reclassified_event), 'ref'))


gtf_novel_exon_anno <- filter(gtf_annotated, transcript_id %in% nx_tx_skipped_exon_tx_eye$transcript_id | 
                                             transcript_id%in% nx_med_cts_eye$transcript_id)
write_tsv(gtf_novel_exon_anno, 'novel_exon_gtf.tsv')
save(nx_incLvl_eye, nx_med_cts_eye, nx_tx_exp_eye, nx_tx_skipped_exon_tx_eye, gtf_novel_exon_anno, file = exon_info_file)

#rtracklayer::readGFF('ref/gencodeGFF3.gff') %>% as.data.frame %>% write_tsv('ref/gencodeGFF3.gff.tsv')
