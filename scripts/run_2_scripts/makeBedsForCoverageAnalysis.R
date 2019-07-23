library(tidyverse)

args=c('/Volumes/data/eyeintegration_splicing/',
       'results/all_tissues.combined.gtf',
       'ref/gencodeAno_comp.gtf',
       'sampleTableV4.tsv',
       'results/salmon_tissue_level_counts.Rdata',
       'results/novel_exon_expression_tables.Rdata',
       'results/novel_exon_workspace.Rdata',
       'results/exons_for_coverage_analysis.bed'
       )
args <- commandArgs(trailingOnly = T)
working_dir <- args[1]
gfc_gtf_file <- args[2]
ref_gtf_file <- args[3]
sample_file <- args[4]
salmon_count_file <- args[5]
exon_info_file <- args[6]
#outputs start here----------/
exon_info_workspace <- args[7]
exon_bed_file <- args[8]

setwd(working_dir)
###part 1 make reference exon set using eyeintegration tpms
load(exon_info_file)
ref_gtf <- rtracklayer::readGFF(ref_gtf_file) %>% mutate(start =start -1)
t2g <- ref_gtf %>%filter(type =='transcript') %>% select(gene_name, transcript_id) %>% distinct
sample_design <- read_tsv('sampleTableV4.tsv', col_names = c('sample', 'run', 'paired', 'tissue', 'subtissue', 'origin'))
tx_quant <- read_tsv('ref/2019_tx_TPM_03.tsv.gz') %>% select(ID, which(colnames(.) %in% sample_design$sample))  
tx_quant <- tx_quant %>% mutate(transcript_id=str_split(ID, '\\(|\\)') %>% sapply(function(x) x[[2]]),  ID=NULL) %>% inner_join(t2g, .)
sample_design <- sample_design %>% .[sample_design$sample %in% colnames(tx_quant),]
subtissues <- sample_design %>% pull(subtissue) %>% unique

top_tx_by_tissue_2 <-  subtissues %>% lapply(function(x) filter(sample_design, subtissue ==x) %>% 
                                                 pull( sample) %>% 
                                                 {select(tx_quant, gene_name, transcript_id, .)} %>%  
                                                 mutate(med_ct= matrixStats::rowMedians(.[,-(1:2)] %>% as.matrix)) %>% 
                                                 filter(med_ct > 1) %>% select(gene_name, transcript_id) %>% mutate(exp=T)
                                             # roughly defining expression as median tpm of 1, so pretty basic might have to tighten this criteria to improve results.
) %>%  reduce(full_join, by=c('gene_name', 'transcript_id'))

colnames(top_tx_by_tissue_2) <- c('gene_name', 'transcript_id',paste0('exp.', subtissues))
top_tx_by_tissue_2[is.na(top_tx_by_tissue_2)] <- 0
ref_exon_tab <- filter(ref_gtf, type =='exon') %>% select(seqid,strand, start,end, transcript_id) %>% inner_join(top_tx_by_tissue_2) %>% distinct
#write_tsv(ref_exon_tab, ref_exon_table)
nx_bed <- nx_inc_cts %>% select(seqid, start, end) %>% distinct
ref_exon_bed <-ref_exon_tab %>% select(seqid, start, end) %>% distinct %>% rbind(nx_bed)

##############################################
#part2

load(salmon_count_file)
salmon_cut_off_lvl <- 15
rmats_cut_off_lvl <- .25
gfc_gtf <- rtracklayer::readGFF(gfc_gtf_file) %>% mutate(start=start-1)
ref_gtf <- filter(gfc_gtf, grepl('ENST', oId)) %>% pull(transcript_id) %>% 
    {filter(gfc_gtf, transcript_id %in% . )} 
transcripts_bed <- gfc_gtf %>%  select(seqid, start, end, transcript_id)
ref_exon_bed <- rbind(ref_exon_bed, filter(gfc_gtf, type =='exon') %>% select(seqid, start, end)) %>% distinct


#write_tsv(transcripts_bed, transcript_loc_bed, col_names = F)

all_ref_exons <- filter(ref_gtf, type=='exon') %>% select(seqid, strand,start,end) %>% 
    mutate(seqid=as.character(seqid)) %>% distinct
############################################################################################
#local intron cov
gtf_ano <- incCounts %>% select(seqid, strand, start, end, ljid) %>% distinct %>% mutate(is.rmats=T) %>% 
    left_join(gfc_gtf, .) %>% filter(type=='exon') %>% mutate(is.rmats=replace_na(is.rmats, F)) 
calc_intron <- function(df){
    locs <- which(df$is.rmats)
    k <- lapply(locs, function(l){
        if((l-1)>0 &((l+1)<=nrow(df))){
            us_e <- c(df[(l-1),'end'],df[l,'end'])
            ds_e <- c( df[l,'start'], df[(l+1),'start'])
            data_frame(seqid=df$seqid[1:2], start=us_e, end=ds_e, id=rep(df$ljid[l],2), tx=df$transcript_id[1:2], strand=df$strand[1:2])
        }else if((l-1)==0){
            us_e <- df[l,'end']
            ds_e <- df[(l+1),'start']
            data_frame(seqid=df$seqid[1], start=us_e, end=ds_e, id=df$ljid[l], tx=df$transcript_id[1], strand=df$strand[1])
        }else if((l+1)==(nrow(df)+1)){
            us_e <- c(df[(l-1),'end'])
            ds_e <- c( df[l,'start'])
            data_frame(seqid=df$seqid[1], start=us_e, end=ds_e, id=df$ljid[l], tx=df$transcript_id[1], strand=df$strand[1]) 
        }else {
            data.frame(seqid=NA, start=NA, end=NA, id=NA, tx=NA,strand=NA)
        }
    }) %>% bind_rows
    return(k)
}



intron_info_tab <- gtf_ano %>% filter(is.rmats) %>% pull(transcript_id) %>% {filter(gtf_ano, transcript_id %in% . , type=='exon')} %>%
    split(.$transcript_id) %>%  lapply(calc_intron) %>% bind_rows %>% filter(!is.na(seqid), !is.na(end)) 
#write_tsv(intron_info_tab,'results/intron_info_tab.tsv')
intron_bed <-  intron_info_tab %>% select(seqid, start, end) %>% distinct
#write_tsv(intron_bed,'testing/introns.bed', col_names = F)

ref_exon_bed <- intron_bed  %>%  rbind(ref_exon_bed) %>% distinct

########
#novel_tx's
novel_loci_tab <- filter(gfc_gtf, type =='transcript') %>% filter(is.na(gene_name)) %>% pull(transcript_id) %>% 
    {filter(gfc_gtf, transcript_id %in% ., type =='exon')} %>% unite(ljid, transcript_id, exon_number,remove = F) %>% 
    select(seqid, strand, start, end, ljid, transcript_id) %>% mutate(is.rmats=T)

gtf_loci_ano <- novel_loci_tab %>% select(seqid, strand, start, end, ljid) %>% distinct %>% mutate(is.rmats=T) %>% 
    left_join(gfc_gtf, .) %>% filter(type=='exon') %>% mutate(is.rmats=replace_na(is.rmats, F)) 
    
novel_loci_intron_tab <- gtf_loci_ano %>% filter(is.rmats) %>% pull(transcript_id) %>% {filter(gtf_ano, transcript_id %in% . , type=='exon')} %>%
    split(.$transcript_id) %>%  lapply(calc_intron) %>% bind_rows %>% filter(!is.na(seqid), !is.na(end)) 

ref_exon_bed <- rbind( ref_exon_bed, novel_loci_tab %>% select(seqid, start, end),
                       novel_loci_intron_tab %>%  select(seqid, start,end))
    

#######
#? is it better to use all ref exons, (even ones not expresed in eye) and then take the exon that creates the longest coverage between
###A3SS and RI's :  we can treat RI's as a3ss in practice as they stradle known exons, and we assume that knownn exons will be covered the sameish
# okay this looks really confusing, but its not that bad. both A3/A5 events can be broken into 

RI_long <- filter(nx_inc_cts, reclassified_event =='RI', is.not_long) %>% 
    select(ljid,reclassified_event,seqid, strand, start, new_end=end) %>% 
    distinct %>% inner_join( all_ref_exons %>% rename(ref_end=end)) %>% rename(new_start=start, end=new_end) %>% 
    inner_join( all_ref_exons %>% rename(ref_start=start)) %>% rename(s_end=end, s_start=new_start, l_start=ref_start,l_end=ref_end ) %>%
    mutate(class='long', delta=NA) %>% select(ljid, reclassified_event, seqid, strand, class, delta, contains('s_'), contains('l_'))



df <- filter(nx_inc_cts, reclassified_event =='A3SS', is.not_long) %>% 
    select(ljid,reclassified_event,seqid, strand, start, new_end=end) %>% 
    distinct %>% inner_join( all_ref_exons %>% rename(ref_end=end)) %>% 
    mutate(delta=new_end -ref_end) 
A3_long <- df %>% 
    filter(delta >0) %>% arrange(desc(delta))  %>% 
    mutate(class='long', l_start=ref_end+1) %>% rename(s_start=start, s_end=ref_end, l_end=new_end) %>% 
    select(ljid, reclassified_event, seqid, strand, class, delta, contains('s_'), contains('l_'))

A3_short <- df %>%
    filter(delta <0) %>%  arrange(delta) %>% 
    mutate(class='short', l_start=new_end+1) %>% rename(s_start=start, l_end=ref_end, s_end=new_end) %>%
    select(ljid, reclassified_event, seqid, strand, class, delta, contains('s_'), contains('l_'))

#write A3 all
###A5SS
df <- filter(nx_inc_cts, reclassified_event=='A5SS', is.not_long) %>% select(ljid, reclassified_event, seqid, strand, new_start=start, end) %>% 
    distinct %>%  inner_join( all_ref_exons %>% rename(ref_start=start)) %>% 
    mutate(delta=new_start -ref_start)

a5_long <- df %>%  filter( delta <0) %>% 
    arrange(delta)  %>% mutate(class='long', l_end= ref_start -1) %>% 
    rename(l_start=new_start, s_start=ref_start, s_end=end )%>%
    select(ljid, reclassified_event, seqid, strand, class, delta, contains('s_'), contains('l_'))

a5_short <- df %>% filter( delta >0) %>% 
    arrange(desc(delta)) %>% 
    mutate(class='short', l_end=new_start -1) %>%
    rename( l_start=ref_start, s_start=new_start, s_end=end) %>%
    select(ljid, reclassified_event, seqid, strand, class, delta, contains('s_'), contains('l_'))

#### novel exons
df <- filter(nx_inc_cts, reclassified_event=='novel_exon') %>% select(ljid, seqid, strand,start,end, reclassified_event) %>% 
    distinct %>% rbind(., all_ref_exons%>%mutate(ljid='.', reclassified_event='ref'))     %>% 
    mutate(reclassified_event=replace_na(reclassified_event, 'ref'),ljid=replace_na(ljid,'.')) %>% arrange(start)

ljids <- filter(df, ljid!='.') %>% pull(ljid)
nx_all <-  lapply(ljids, function(x){i <- which(df$ljid == x)[1] # find first occurence of novel exon and
while(df[i-1,'reclassified_event']=='novel_exon'){ i <- i-1} # find the index of closest reference exon
data.frame(df[i,], ref_start=df[i-1,'start'], ref_end=df[i-1,'end'])
}) %>% # create output df;
    bind_rows %>% mutate(class='nx', ref_id=paste('nx', 1:nrow(.), sep = '_'))


ass_tab <- rbind(a5_long, a5_short, A3_long, A3_short, RI_long)

complete_bed <- rbind(ass_tab %>% select(seqid, start=s_start, end=s_end),
                      ass_tab %>% select(seqid, start=l_start, end=l_end),
                      nx_all %>% select(seqid, start, end),
                      nx_all %>% select(seqid, start=ref_start,end=ref_end),
                      stringsAsFactors=F) %>% rbind(ref_exon_bed) %>% distinct %>% filter(!is.na(seqid), !is.na(start),!is.na(end))

true_novel_exons <- nx_all
ref_bed <- complete_bed
save(ass_tab, true_novel_exons ,ref_exon_tab, intron_info_tab, nx_inc_cts, nx_psi, nx_skipped_exon,ref_bed, 
     novel_loci_intron_tab, novel_loci_tab, file = exon_info_workspace)
write_tsv(complete_bed, exon_bed_file, col_names = F)

