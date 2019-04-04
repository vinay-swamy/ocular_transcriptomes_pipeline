library(tidyverse)
gff_file <- 'results/stringtie_alltissues_cds.gff3'
gencode_gtf_file <- 'ref/gencodeAno_bsc.gtf'
gff <- rtracklayer::readGFF(gff_file) %>% as.data.frame 
gff_txdb <-  filter(gff, type=='mRNA') %>% select(transcript_id,oId,gene_id,cmp_ref)
novel_txdb <- filter(gff_txdb, grepl('MSTRG',oId), !is.na(cmp_ref))
to_blast_gff3 <- filter(gff, oId%in% novel_txdb$cmp_ref) %>%left_join(novel_txdb, .,  by='cmp_ref')%>% 
  select(transcript_id.x,transcript_id.y)%>%unique
write_tsv(to_blast_gff3, path='ref/noveltx_2_ref_tx.tsv', col_names = F)

gencodeGTF <- rtracklayer::readGFF(gencode_gtf_file) %>% filter(type=='transcript')
gencodeGTF_pc <- filter(gencodeGTF, gene_type=='protein_coding') %>% select(start, end, gene_id, transcript_id) %>% 
  mutate( transcript_id= strsplit(transcript_id,'\\.')%>%sapply(function(x) x[[1]]))

df <- mutate(gencodeGTF_pc, len=end-start) %>% split(.['gene_id']) %>% lapply(function(x) arrange(x,len) %>% .[c(1,nrow(x)),])%>% 
  lapply(function(x) cbind(x[1,], x[2,]) ) %>% do.call(rbind,.) 
dfi <- df[,c(4,9)]
same <- dfi[,1]==dfi[,2]
dfi <- dfi[!same,]

znf <- filter(gff,gene_name=='ZNF469')
