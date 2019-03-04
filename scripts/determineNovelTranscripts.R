library(tidyverse)
library(matrixStats)

args <- c('~/NIH/eyeintegration_splicing/',
          'RPE',
          'results_b38/all_tissues.combined.gtf',
          'results_b38/all_rmats_events_tissues.incLevels.tsv',
          'results_b38/all_rmats_events_tissues.medCounts.tsv',
          'sampleTableV4.tsv',
          'results_b38/salmon_tx_quant.Rdata',
          'results_b38/salmon_gene_quant.Rdata')

args <- commandArgs(trailingOnly = T)

working_dir <- args[1]
gfc_gtf_file <- args[2]
incTab_file <- args[3]
medCts_file <- args[4]
sample_table_file <- args[5]
salmon_tx_quant <- args[6]
salmon_gene_quant <- args[7]



event_header <- c('ljid', 'GeneID', 'seqid'	,'strand',	'start', 'end',	'upstreamES',	'upstreamEE',	'downstreamES',	'downstreamEE')


gfc_gtf <- rtracklayer::readGFF(gfc_gtf_file) %>% rename(TxID=transcript_id) %>% mutate(start=start-1)
txid_2_oid <- filter(gfc_gtf, type=='transcript') %>% select(TxID, oId)
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

med_tx_counts <- make_tissue_level_counts(sample_table, tx_counts) %>% as.data.frame %>% mutate(TxID=tx_counts$transcript_id) %>%
    select(TxID, everything(), -grep('synth', colnames(.)))
med_gene_counts <- make_tissue_level_counts(sample_table, gene_counts) %>% as.data.frame %>% mutate(GeneID=gene_counts$gene_id) %>% 
    select(GeneID, everything(), -grep('synth', colnames(.)))

#set cuttoff levels at 25, seems like its a reasonable level, based on looking at expression quantiles
txid_2_oid <- filter(gfc_gtf, type=='transcript') %>% select(TxID, oId)
salmon_cut_off_lvl <- 25
rmats_cut_off_lvl <- 25
# stringtie expressed in any tissue
expressed_in_any_tissue <-  rowSums(med_tx_counts[,-1] > salmon_cut_off_lvl) >=1
stringtie_all_expressed_tx <- med_tx_counts[expressed_in_any_tissue,]

k <- apply(tx_counts[,-1], 2, summary) %>%  t()


#rmats rminimally expressed
lvl=.95
nc_cols <- c('ljid',event_header,'event_type')
# remove exons with low counts 
# remove exons that are expressed aboce cutoff level
mat=medCounts%>%select(-nc_cols)
highly_expressed_exons <- rowSums(mat > rmats_cut_off_lvl) >0
high_in_all_tissues <- rowSums(mat > lvl) == (ncol(mat))
incTab_filt <- incTab[highly_expressed_exons & !high_in_all_tissues,]
ctsTab_filt <- medCounts[highly_expressed_exons & ! high_in_all_tissues,]

#determine overlap between the rmats and stringtie
gtf_expressed <- filter(gfc_gtf, TxID %in% stringtie_all_expressed_tx$TxID)
gtf_exp_exons <- filter(gtf_expressed,type=='exon') %>% select(TxID, seqid,strand, start,end)
rmats_exon_exp_tab <- select(incTab_filt, ljid, seqid,strand, start,end)
rmats_id_to_tx_id <- inner_join(gtf_exp_exons, rmats_exon_exp_tab, by=c('seqid','strand', 'start','end'))
rmats_st_baseline <- filter(ctsTab_filt, ljid %in% rmats_id_to_tx_id$ljid)
stringtie_rm_baseline <- filter(med_tx_counts, TxID %in% rmats_id_to_tx_id$TxID)
stringtie_novel_expressed_all <-  filter(txid_2_oid, grepl('MSTRG',oId)) %>% pull (TxID) %>% 
                                    {filter(stringtie_rm_baseline, TxID %in% .)} 

#determine novel transcript expression within the eye
eye_tissues <- grepl('Retina|RPE|Cornea', colnames(stringtie_rm_baseline))
expressed_in_eye <-  rowSums(stringtie_rm_baseline[,eye_tissues] > salmon_cut_off_lvl) >=1
stringtie_expressed_in_eye <- stringtie_rm_baseline[expressed_in_eye,]
stringtie_novel_eye <- filter(txid_2_oid, grepl('MSTRG',oId)) %>% pull (TxID) %>% 
    {filter(stringtie_expressed_in_eye, TxID %in% .)}

#determine exon expression in the eye
eye_cols <- grepl('Retina|RPE|Cornea', colnames(rmats_st_baseline))
expressed_in_eye <- rowSums(rmats_st_baseline[,eye_cols] >rmats_cut_off_lvl) >=1
rmats_expressed_eye <- rmats_st_baseline[expressed_in_eye,] %>% inner_join(rmats_id_to_tx_id, by='ljid') 

#determinr which stringtie novel exons are detected in rmats; 
##first determine which st exons are novel

gtf.genes_with_novel_tx <- pull(stringtie_novel_expressed_all, TxID) %>% {filter(gfc_gtf, TxID %in% . )} %>% pull(GeneID) %>%
{filter(gfc_gtf, GeneID %in% .)}
# create a set of reference exons by selecting all exons from reference transcripts.    
all_ref_exons <- filter(gfc_gtf, grepl('ENST', oId)) %>% pull(TxID) %>% 
                    {filter(gfc_gtf, TxID %in% . , type=='exon')} %>% select(TxID, seqid, strand,start,end)
all_ref_exon_list <- split(all_ref_exons[,-1], 1:nrow(all_ref_exons))


all_novel_exons <- filter(gfc_gtf, TxID %in% stringtie_novel_expressed_all$TxID, type=='exon')  %>% 
    select(TxID, seqid, strand,start,end)
all_novel_exon_list <- split(all_novel_exons[,-1], 1:nrow(all_novel_exons))


novel <- !all_novel_exon_list%in%all_ref_exon_list
novel_string_tie_exons <- all_novel_exons[novel,] %>% mutate(is.novel=T)
novel_st_exons_in_rmats <- inner_join(novel_string_tie_exons, rmats_st_baseline, by=c('seqid','strand', 'start','end'))

gtf_novel_exon_anno <-  filter(gfc_gtf,TxID %in% stringtie_novel_expressed_all$TxID) %>% 
    left_join(novel_st_exons_in_rmats[,-1], by=c('seqid','strand', 'start','end')) %>% 
    mutate(is.novel=replace(is.novel, is.na(is.novel),F))
n_tx_with_novel_exon <- gtf_novel_exon_anno  %>% group_by(TxID) %>% summarise(any(is.novel))


novel_st_exons_in_rmats_eye <- filter(novel_st_exons_in_rmats, TxID %in% stringtie_expressed_in_eye$TxID, 
                                      ljid%in% rmats_expressed_eye$ljid)
novel_exon_med_cts <- filter(rmats_st_baseline, ljid%in%novel_st_exons_in_rmats_eye$ljid)  %>% inner_join(rmats_id_to_tx_id[,c(1,6)],. ,by='ljid')
novel_exon_inc_lvl <- filter(incTab_filt, ljid%in% novel_st_exons_in_rmats_eye$ljid) %>% inner_join(rmats_id_to_tx_id[,c(1,6)],. ,by='ljid')
novel_exon_tx_counts <- filter(med_tx_counts, TxID%in% novel_exon_med_cts$TxID )
novel_exon_gtf <- filter(gtf_novel_exon_anno, TxID%in% novel_st_exons_in_rmats_eye$TxID)
tx_to_write <- filter(gfc_gtf, grepl('ENST', oId)) %>% pull(TxID) %>% c(novel_st_exons_in_rmats_eye$TxID)
write(tx_to_write, 'tx_to_use_for_comparison.txt',sep = '\n' )
gtf_novel_wide <- filter(gfc_gtf, TxID%in% novel_st_exons_in_rmats_eye$TxID)
write_tsv(gtf_novel_wide, 'novel_exon_gtf.tsv')
save(novel_exon_med_cts, novel_exon_inc_lvl, novel_exon_tx_counts, novel_exon_gtf, gtf_novel_exon_anno, n_tx_with_novel_exon, file = 'all_novel_exon_info.Rdata')
gdata::keep(novel_exon_med_cts, novel_exon_inc_lvl, novel_exon_tx_counts, novel_exon_gtf, gtf_novel_exon_anno, n_tx_with_novel_exon, sure=T)

















