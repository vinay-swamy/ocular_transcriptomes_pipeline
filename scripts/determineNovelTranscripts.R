library(tidyverse)
library(matrixStats)

args <- c('~/NIH/eyeintegration_splicing/',
          'results_b38/all_tissues.combined.gtf',
          'results_b38/all_rmats_events_tissues.incLevels.tsv',
          'results_b38/all_rmats_events_tissues.medCounts.tsv',
          'sampleTableV4.tsv',
          'rdata/salmon_no_ssn_gene_quant.Rdata',
          'rdata/salmon_no_ssn_tx_quant.Rdata')

args <- commandArgs(trailingOnly = T)

working_dir <- args[1]
gfc_gtf_file <- args[2]
incTab_file <- args[3]
medCts_file <- args[4]
sample_table_file <- args[5]
salmon_gene_quant<- args[6]
salmon_tx_quant <- args[7]



event_header <- c('ljid', 'GeneID', 'seqid'	,'strand',	'start', 'end',	'upstreamES',	'upstreamEE',	'downstreamES',	'downstreamEE')


gfc_gtf <- rtracklayer::readGFF(gfc_gtf_file) %>% dplyr::rename(TxID=transcript_id) %>% mutate(start=start-1)
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
# these are rmats exons found in stringtie
rmats_st_baseline <- filter(ctsTab_filt, ljid %in% rmats_id_to_tx_id$ljid)
# these are stringtie transcripts with a detected rmats exon
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

# next, create a set of novel exonss
all_novel_exons <- filter(gfc_gtf, TxID %in% stringtie_novel_expressed_all$TxID, type=='exon')  %>% 
    select(TxID, seqid, strand,start,end)
all_novel_exon_list <- split(all_novel_exons[,-1], 1:nrow(all_novel_exons))
novel <- !all_novel_exon_list%in%all_ref_exon_list
#this table is all novel exons
novel_string_tie_exons <- all_novel_exons[novel,] %>% mutate(is.novel=T)
#this is all novel exons that were detected in rmats
novel_st_exons_in_rmats <- inner_join(novel_string_tie_exons, rmats_st_baseline, by=c('seqid','strand', 'start','end'))

novel_st_exons_in_rmats_eye <- filter(novel_st_exons_in_rmats, TxID %in% stringtie_expressed_in_eye$TxID, 
                                      ljid%in% rmats_expressed_eye$ljid)
nx_med_cts_eye <- filter(medCounts, ljid %in% novel_st_exons_in_rmats_eye$ljid)
nx_incLvl_eye <- filter(incTab, ljid %in% novel_st_exons_in_rmats_eye$ljid)


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
rtracklayer::readGFF('ref/gencodeGFF3.gff') %>% as.data.frame %>% write_tsv('ref/gencodeGFF3.gff.tsv')



library(IRanges)
#scripts/determine_novel_exon_impact

#yet another episode of WIWT
determine_novel_exon_impact_but_not_like_a_dumbass <- function(gtf,ref, nx){
    contained_within <- function(feat_df, nx){
        # create an interval for a given feature; don't need to worry about strand
        feat_df <- arrange(feat_df, start)
        n <- nrow(feat_df)
        if(n==0) return(F)
        ref_rang <- IRanges(start = feat_df$start[1], end = feat_df$end[n])
        nx_rang <- IRanges(start=nx$start, end=nx$end)
        k <- findOverlaps(ref_rang, nx_rang)
        return(length(k) >0) 
    }
    
    identify_contained_feature <- function(feature_list, nx){
        res <- sapply(feature_list, function(y) contained_within(y,nx))
        #print(res)
        if(res[[2]]) return('CDS')
        if(res[[1]]) return('5UTR')
        if(res[[3]]) return('3UTR')
        
        return('OUTSIDE')
        
    }
    gtf.tx <- filter(gtf, transcript_id==nx$transcript_id[1])
    ref_tx_id <- gtf.tx$cmp_ref[1] 
    if(is.na(ref_tx_id)) return('TOSS:no_ref_id')
    ref.tx <- filter(ref,transcript_id==ref_tx_id)
    if(nrow(ref.tx)==0) return('TOSS:ref_not_found')
    if(ref.tx$transcript_type[1]!='protein_coding') return('NO_PC')
    ref_start <- filter(ref.tx, type=='transcript') %>% pull(start) 
    feature_list <- list(
        fputr = filter(ref.tx, type=='five_prime_UTR'),
        cds = filter(ref.tx, type=='CDS'),
        tputr = filter(ref.tx, type=='three_prime_UTR'))
    res <- identify_contained_feature(feature_list, nx)
    if(res=='OUTSIDE'){
        if(ref.tx$strand[1]=='+') return(ifelse(nx$start<ref_start, 'US', 'DS'))
        if(ref.tx$strand[1]=='-') return(ifelse(nx$start>ref_start, 'US', 'DS'))
    }
    return(res)
}


nx_l <- dplyr::rename(novel_exon_inc_lvl, transcript_id=TxID) %>% select(transcript_id, seqid, strand, start, end) %>% split(1:nrow(.))
ref <- rtracklayer::readGFF('~/NIH/eyeintegration_splicing/ref/gencodeGFF3.gff') %>% as.data.frame
gtf <- dplyr::rename(gfc_gtf, transcript_id=TxID)
determine_novel_exon_impact_but_not_like_a_dumbass(gtf, ref, nx_l[[1]])
simple_consequence <-sapply(nx_l, function(x) determine_novel_exon_impact_but_not_like_a_dumbass(gtf, ref, x))
tab <- table(simple_consequence)
barplot(tab[!grepl('TOSS', names(tab))], ylim = c(0,3000))
k= simple_consequence=='5UTR'

novel_exon_5putr <- nx_l[k] %>% do.call(rbind,.)












