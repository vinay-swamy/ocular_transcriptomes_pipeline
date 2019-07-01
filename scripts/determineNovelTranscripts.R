library(GenomicRanges)
library(tidyverse)
library(matrixStats)

args <- c('/Volumes/data/eyeintegration_splicing/',
          'data/gtfs/all_tissues.combined.gtf',
          'ref/gencode_comp_ano.gtf',
          'ref/ensembl_ano.gtf',
          'ref/refseq_ncbi.gff3',
          'ref/ucsc.gtf',
          'sampleTableDev.tsv',
          'data/rmats/all_tissues.PSI.tsv',
          'data/rmats/all_tissues.incCts.tsv',
          'testing/novel_exon_expression_tables.Rdata')

args <- commandArgs(trailingOnly = T)
working_dir <- args[1]
gfc_gtf_file <- args[2]
ref_gtf <- args[3]
ensbl_gtf <- args[4] 
refseq_gtf <- args[5]
ucsc_gtf <- args[6]
sample_table_file <- args[7]
psiTab_file <- args[8]
incCts_file <- args[9]
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

rmats_minimally_expressed <- incCounts

k <- anti_join(incCounts, gencode_ref)


#all_ref_exons
gencode_ref<-  rtracklayer::readGFF(ref_gtf) %>% mutate(start=start - 1)  %>% filter(type =='exon') %>% select(seqid, strand,start,end) %>% distinct %>% 
    mutate(seqid=as.character(seqid))
ensembl_ref <- rtracklayer::readGFF(ensbl_gtf) %>% mutate(start=start - 1)  %>% filter(type =='exon') %>% select(seqid, strand,start,end) %>% distinct %>% 
    mutate(seqid=as.character(seqid), seqid= paste0('chr', seqid))
refseq_ucsc <- rtracklayer::readGFF(ucsc_gtf ) %>% mutate(start=start-1, seqid=as.character(seqid)) #%>% 
gene2seq <- refseq_ucsc %>%  filter(type =='transcript') %>%  select(chroms=seqid, gene=gene_name) %>% distinct
refseq_ucsc <- refseq_ucsc %>% filter(type=='exon') %>% select(seqid, strand, start, end) %>% distinct
refseq_ncbi <- rtracklayer::readGFF(refseq_gtf) %>% as.data.frame() %>% 
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
    mutate(length= end -start, is.not_long=length <1500) %>% filter(is.not_long)

novel_st_exons_in_rmats <- determine_overlap_Exons(novel_st_exons_in_rmats %>% filter(reclassified_event=='novel_exon') %>% distinct,all_ref_exons) %>% 
    left_join(novel_st_exons_in_rmats,.) %>%
    mutate(ss_exon=replace_na(ss_exon, F), reclassified_event = replace(reclassified_event, ss_exon, 'subset_exon'), ss_exon=NULL) %>% 
    inner_join( incCounts) %>% 
    inner_join(select(gfc_gtf, transcript_id, seqid, strand, start, end)) %>% 
    select(ljid, seqid, start, end, transcript_id, reclassified_event, everything())


wsize <- function(l) return(c(l[2]-l[1], l[3]-l[2]))

k <- a3s[,4:6]
m <- apply(k, 1, function(x) wsize(sort(x))) %>% t() %>% as.data.frame 




a3s <-  novel_st_exons_in_rmats %>% filter(reclassified_event == 'A3SS') %>% select(ljid, seqid, strand, start, alt_end=end) %>% 
    inner_join(gencode_ref) %>% mutate(delta = alt_end - end) %>% cbind(., apply(.[,4:6], 1, function(x) wsize(sort(x))) %>% t() %>% as.data.frame ) %>% 
    dplyr::rename(us_size=V1, ds_size=V2) %>% filter(us_size>=40, ds_size>=40) %>% distinct 

a3s %>% group_by(ljid) %>% summarise(n=n()) %>% pull(n) %>% table 
exon <- a3s[1,2:5] %>% dplyr::rename(end=alt_end) %>% mutate(is.novel=T)

process_a3s_events <- function(exon){
    refs <- exon %>% select(everything(),alt_end=end) %>% inner_join(gencode_ref) %>% mutate(is.novel=F) %>% select(colnames(exon))
    rbind(refs, exon) %>% arrange(desc(end))
    if(nrow(refs)==2){
        
        
        
    }else 
    rbind(refs, exon) %>% arrange(desc(end))
    
    
    
}







###figure out skipped exon transcripts
# gtf_se<- left_join(gfc_gtf, select(novel_string_tie_exons, seqid, strand, start, end)%>%distinct %>% mutate(is.novel=T), 
#                    by=c('seqid','strand' ,'start','end')) %>%
#     mutate(is.novel=replace(is.novel, is.na(is.novel), F))
# skipped_exon_tx <- filter(gtf_se, grepl('MSTRG', oId)) %>% pull(transcript_id) %>% 
#     {filter(gtf_se, transcript_id %in% .)} %>% group_by(transcript_id) %>% summarise(no_new_exon=all(!is.novel)) %>%
#     filter(no_new_exon)
# txid2gene <- filter(gfc_gtf, type=='transcript') %>% select(transcript_id, gene_name) %>% distinct %>% 
#     mutate(gene_name=replace(gene_name, is.na(gene_name), transcript_id[is.na(gene_name)]))








nx_inc_c ts <- novel_st_exons_in_rmats %>% left_join(txid2gene) %>% select(-contains('incC'), contains('incC'))
nx_psi <- inner_join( nx_inc_cts%>%select(transcript_id, reclassified_event, ljid, gene_name ), psiTab) %>% select(-contains('PSI'), contains('PSI'))
nx_skipped_exon <- skipped_exon_tx %>% left_join(txid2gene)

save(nx_inc_cts, nx_psi, nx_skipped_exon, file = exon_info_file )

