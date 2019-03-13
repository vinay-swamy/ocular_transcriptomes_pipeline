library(tidyverse)
args=c('~/NIH/eyeintegration_splicing/', 'results_b38/all_tissues.combined.gtf', 
       'all_novel_exon_info.Rdata', 'results_b38/as_event_ls_class.Rdata',
       '/Volumes/data/eyeintegration_splicing/results/exons_for_cov_alys.bed')
args <- commandArgs(trailingOnly = T)
working_dir <- args[1]
gfc_gtf_file <- args[2]
exon_info_file <- args[3]
event_ls_file <- args[4]
exon_bed_file <- args[5]


setwd(working_dir)
load(exon_info_file)
gfc_gtf <- rtracklayer::readGFF(gfc_gtf_file) %>% mutate(start=start-1)
ref_gtf <- filter(gfc_gtf, grepl('ENST', oId)) %>% pull(transcript_id) %>% 
{filter(gfc_gtf, transcript_id %in% . )} 
all_ref_exons <- filter(ref_gtf, type=='exon') %>% select(seqid, strand,start,end) %>% 
    mutate(seqid=as.character(seqid)) %>% distinct


#? is it better to use all ref exons, (even ones not expresed in eye) and then take the exon that creates the longest coverage between
###A3SS and RI's :  we can treat RI's as a3ss in practice as they stradle known exons, and we assume that knownn exons will be covered the sameish
df <- filter(nx_med_cts_eye, reclassified_event=='A3SS' | reclassified_event == 'RI') %>% 
    select(ljid,reclassified_event,seqid, strand, start, new_end=end) %>% 
    distinct %>% unite(spl,seqid,strand,start, remove = F)
A3_new_longer_than_ref <- inner_join(df, all_ref_exons %>% rename(ref_end=end)) %>% 
    mutate(delta=new_end -ref_end)  %>% 
    filter(delta >0) %>% arrange(desc(delta)) %>% filter(!duplicated(spl)) %>% mutate(class='long')

A3_All <- filter(df, !spl%in%A3_new_longer_than_ref$spl) %>% 
    inner_join(all_ref_exons %>% rename(ref_end=end))%>% 
    mutate(delta=new_end -ref_end)%>%
    filter(delta <0) %>%  arrange(delta) %>% filter(!duplicated(spl)) %>% 
    mutate(class='short') %>% rbind(A3_new_longer_than_ref) %>% mutate(ref_id=paste('A3SS', 1:nrow(.), sep = '_'))
#write A3 all
###A5SS
df <- filter(nx_med_cts_eye, reclassified_event=='A5SS') %>% select(ljid, reclassified_event, seqid, strand, new_start=start, end) %>% 
    distinct %>% unite(spl,seqid,strand,end, remove = F)

a5_new_longer_than_ref <- inner_join(df, all_ref_exons %>% rename(ref_start=start)) %>% 
    mutate(delta=new_start -ref_start) %>%  filter( delta <0) %>% 
    arrange(delta) %>% filter(!duplicated(spl)) %>% mutate(class='long')

a5_all <-filter(df, !spl%in%a5_new_longer_than_ref$spl) %>% 
    inner_join(all_ref_exons %>% rename(ref_start=start)) %>% 
    mutate(delta=new_start -ref_start) %>% filter( delta >0) %>% 
    arrange(desc(delta)) %>% filter(!duplicated(spl)) %>%
    mutate(class='short') %>% rbind(a5_new_longer_than_ref) %>% mutate(ref_id=paste('A5SS', 1:nrow(.), sep = '_'))
#### novel exons
df <- filter(nx_med_cts_eye, reclassified_event=='novel_exon') %>% select(ljid, seqid, strand,start,end, reclassified_event) %>% 
    distinct %>% rbind(., all_ref_exons%>%mutate(ljid='.', reclassified_event='ref'))     %>% 
    mutate(reclassified_event=replace_na(reclassified_event, 'ref'),ljid=replace_na(ljid,'.')) %>% arrange(start)

ljids <- filter(df, ljid!='.') %>% pull(ljid)
nx_all <-  lapply(ljids, function(x){i <- which(df$ljid == x)[1] # find first occurence of novel exon and
                                     while(df[i-1,'reclassified_event']=='novel_exon'){ i <- i-1} # find the index of closest reference exon
                                     data.frame(df[i,], ref_start=df[i-1,'start'], ref_end=df[i-1,'end'])
                                     }) %>% # create output df;
    bind_rows %>% mutate(class='nx', ref_id=paste('nx', 1:nrow(.), sep = '_'))

complete_bed <- rbind(A3_All %>% filter(class=='long') %>% select(seqid, start=ref_end, end=new_end, id =ljid),#long novel region
                      A3_All %>% filter(class=='long') %>% select(seqid, start, end=ref_end, id= ref_id),# long ref region
                      A3_All %>% filter(class=='short') %>%select(seqid, start, end=new_end,id=ljid), # short novel region
                      A3_All %>% filter(class=='short') %>%select(seqid, start=new_end,end=ref_end,id=ref_id), # short ref region
                      a5_all %>% filter(class=='long') %>% select(seqid, start=new_start, end=ref_start,id= ljid), # long novel
                      a5_all %>% filter(class=='long') %>% select(seqid, start=ref_start, end,id=ref_id),# long ref
                      a5_all %>% filter(class=='short') %>%select(seqid, start=new_start, end,id= ljid), #short novel
                      a5_all %>% filter(class=='short') %>%select(seqid, start=ref_start, end=new_start,id=ref_id),#
                      nx_all %>% select(seqid, start, end,id= ljid),
                      nx_all %>% select(seqid, start=ref_start,end=ref_end, id=ref_id),
                      stringsAsFactors=F)
save(A3_All,a5_new_longer_than_ref, nx_all, file = event_ls_file)
write_tsv(complete_bed, exon_bed_file, col_names = F)
    

    
    
    
    









