library(tidyverse)
library(pheatmap)


stringtie_gtf_file <- 'results/all_tissues.combined.gtf'
incLvlTab_file <- 'results/rmats_all_tissue_incLevels.tsv'
medCounts_files <-  'results/rmats_all_tissue_medCounts.tsv'

st_gtf <- rtracklayer::readGFF(stringtie_gtf_file) %>% mutate(start=start-1, seqid=as.character(seqid))

incTab <- read_tsv(incLvlTab_file) %>% 
  rename(seqid=chr,start=exonStart_0base,end=exonEnd) %>% select( -synth_incLevel)
medCounts <- read_tsv(medCounts_files) %>% 
  rename(seqid=chr,start=exonStart_0base,end=exonEnd) %>% select( -synth_incLevel)
{
st_gtf_list <- select(st_gtf, seqid,strand,start,end) %>% split( 1:nrow(st_gtf))
incTab_list <- select(incTab,seqid,strand,start,end) %>% split(1:nrow(incTab))
st_gtf$in_rmats <- st_gtf_list%in%incTab_list
rmats_in_st <- incTab_list%in%st_gtf_list
sum(rmats_in_st) #rn they are all in! woot 
novel_tx_ids <- filter(st_gtf, grepl('MSTRG',oId)) %>% pull(transcript_id)
rmats_tx_ids <- filter(st_gtf,in_rmats) %>% pull(transcript_id) %>%unique
rmats_gtf <- filter(st_gtf, transcript_id %in% rmats_tx_ids)}
# lets pick exons that have a median count 50 in at least 1 tissue
highly_expressed_exons <- rowSums(medCounts[,10:36] > 50) >0
incTab_filt <- incTab[highly_expressed_exons,]
itf <- incTab_filt[,10:35]
#pheatmap(incTab_filt[,10:35], treeheight_row = 0)
lvl <-.95
high_inc_one_tissue <- rowSums(itf > lvl) == 1
low_inc_other_tissue <- rowSums(itf < (1-lvl)) == 25 
tissue_specfic <- high_inc_one_tissue & low_inc_other_tissue
sum(tissue_specfic)
inctab_tissue_specifc <- incTab_filt[tissue_specfic,]
View(filter(inctab_tissue_specifc,Lens_Stem.Cell.Line_incLevel > .5 ))

pheatmap(inctab_tissue_specifc[,10:35], treeheight_row = 0,show_rownames = F)
