library(tidyverse)
library(RBedtools)
source('~/scripts/write_gtf.R')
args <- c('~/NIH/dev_eyeintegration_splicing/', '~/NIH/occular_transcriptomes_paper/all_tissues.combined.gtf', 
          '/Volumes/data/eyeintegration_splicing/sampleTableV6.tsv', 
          '~/NIH/eyeintegration_splicing/dl_data/all_tissues.combined_transdecoderCDS.gff3.gz',
          '~/NIH/occular_transcriptomes_paper/data/all_tissues.combined_V1.Rdata',
          '~/NIH/occular_transcriptome_shiny/all_tissues.combined_NovelAno.gtf')
args <- commandArgs(trailingOnly = T)
working_dir <- args[1]
gfc_gtf_file <- args[2]
sample_table_file <- args[3]
gff3_file <- args[4]
outfile <- args[5]
gtf_ano_outfile <- args[6]
setwd(working_dir)

if(!file.exists('rdata/all_ref_tx_exons.rdata')){
  base_ref <- rtracklayer::readGFF('ref/gencode_comp_ano.gtf') %>% mutate(origin='gencode', seqid=as.character(seqid))
  ensembl_ref <- rtracklayer::readGFF('ref/ensembl_ano.gtf') %>% mutate(seqid=as.character(seqid), origin='ensembl')
  chr <-  paste0('chr', ensembl_ref$seqid %>% as.character)
  non_chr <-  ensembl_ref$seqid %>% as.character %>% as.numeric %>% is.na
  ensembl_ref <- ensembl_ref %>% mutate(seqid=replace(seqid %>% as.character,!non_chr, chr[!non_chr]))
  
  refseq_ucsc <- rtracklayer::readGFF('ref/ucsc.gtf')  %>% mutate(origin='ucsc', seqid=as.character(seqid))
  
  refseq_ncbi <- rtracklayer::readGFF('ref/refseq_ncbi.gff3') %>% as.tibble
  seq2reg <- refseq_ncbi %>% filter(type == 'region', !is.na(chromosome), chromosome!='Unknown') %>% select(seqid, chromosome) %>%
    distinct %>%
    mutate(seqid=as.character(seqid), rep_seq = paste0('chr', chromosome), rep_seq=case_when(rep_seq=='chrX' ~ 'X',
                                                                                             rep_seq=='chrY' ~ 'Y',
                                                                                             TRUE ~ rep_seq))
  refseq_ncbi <- inner_join(refseq_ncbi, seq2reg[,-2]) %>% rename(region_id=seqid, seqid=rep_seq) %>%
    select(seqid, region_id, everything()) %>% mutate(origin='ncbi')
  
  all_transcripts <- rbind(
    base_ref %>% filter(type == 'transcript') %>% select(seqid, strand, start, end, origin),
    ensembl_ref %>% filter(type == 'transcript') %>% select(seqid, strand, start, end, origin),
    refseq_ucsc %>% filter(type == 'transcript') %>% select(seqid, strand, start, end, origin),
    refseq_ncbi %>% filter(type == 'transcript') %>% select(seqid, strand, start, end, origin)
  ) %>% mutate(seqid=as.character(seqid)) %>% filter(!duplicated(.[,-5]))
  
  table(all_transcripts$origin)
  all_exons <- rbind(
    base_ref %>% filter(type == 'exon') %>% select(seqid, strand, start, end, origin),
    ensembl_ref %>% filter(type == 'exon') %>% select(seqid, strand, start, end, origin),
    refseq_ucsc %>% filter(type == 'exon') %>% select(seqid, strand, start, end, origin),
    refseq_ncbi %>% filter(type == 'exon') %>% select(seqid, strand, start, end, origin)
  ) %>%
    mutate(seqid=as.character(seqid)) %>% filter(!duplicated(.[,-5]))
  save(all_transcripts, all_exons, file='rdata/all_ref_tx_exons.rdata')
} else {load('rdata/all_ref_tx_exons.rdata')}
gfc_gtf <- rtracklayer::readGFF(gfc_gtf_file)

novel_transcripts <- anti_join(gfc_gtf %>% filter(type == 'transcript'), all_transcripts) %>% filter(!grepl('TCONS', gene_name))
novel_loci <- anti_join(gfc_gtf %>% filter(type == 'transcript'), all_transcripts) %>% filter(grepl('TCONS', gene_name))
novel_single_exon_tx <- novel_transcripts$transcript_id %>% {filter(gfc_gtf, transcript_id  %in% .)} %>% group_by(transcript_id) %>%
  summarise(count=n()) %>% filter(count == 2) %>% pull(transcript_id) %>%  {filter(novel_transcripts, transcript_id %in% .)}
novel_transcripts <- filter(novel_transcripts, !transcript_id %in% novel_single_exon_tx$transcript_id)

# remove novel loci that overlap with known genes
novel_loci_bed <- novel_loci %>% filter(type == 'transcript') %>%  mutate(score=999) %>% 
  select(seqid, start, end, transcript_id, score, strand) %>% from_data_frame %>% RBedtools('sort',i=.)

intersect <- all_transcripts %>% mutate(score=999) %>% select(seqid, start, end, origin, score, strand) %>% 
  from_data_frame %>% 
  RBedtools('sort', output = 'stdout', i=.) %>% 
  RBedtools('intersect',options = '-loj -s',a=novel_loci_bed, b=.  ) %>% 
  to_data_frame



novel_loci_distinct <- filter(intersect, X8 == -1) %>% pull(X4) %>% {filter(novel_loci, transcript_id %in% .)} 


novel_exons <- gfc_gtf %>% 
  filter(type == 'exon', !transcript_id %in% novel_loci$transcript_id, !transcript_id %in% novel_single_exon_tx$transcript_id ) %>%
  select(seqid, strand, start, end) %>% 
  anti_join( all_exons) %>% distinct %>% 
  mutate(id=paste0('nvl_exon', 1:nrow(.)))

nvl_start <- anti_join(novel_exons %>% select(seqid, strand, start, id), all_exons) %>% { novel_exons$id  %in% .$id}
nvl_end <- anti_join(novel_exons %>% select(seqid, strand, end, id), all_exons) %>% {novel_exons$id  %in% .$id }

novel_exons <- novel_exons %>% mutate(nv_type=case_when(nvl_start & nvl_end ~ 'novel_exon',
                                                        !nvl_start & nvl_end ~ 'A3SS',
                                                        nvl_start & !nvl_end ~ 'A5SS',
                                                        !nvl_start & !nvl_end ~ 'RI'))




gfc_gtf_ano <- filter(gfc_gtf, !transcript_id %in% novel_loci$transcript_id)
gfc_gtf_ref <- filter(gfc_gtf_ano, !transcript_id %in% novel_transcripts$transcript_id)

gfc_gtf_full <-  gfc_gtf_ano %>% filter(transcript_id  %in% novel_transcripts$transcript_id) %>% select(seqid, strand, start, end) %>%
  distinct %>% anti_join(gfc_gtf_ref) %>% anti_join(all_exons) %>% anti_join(all_transcripts) %>%  
  mutate(is.novel=T) %>% left_join(gfc_gtf_ano, .) %>% mutate(is.novel=replace_na(is.novel, F))


# 
# uniq_tss <-  gfc_gtf_full %>% filter(exon_number==1) %>% 
#     select(seqid, strand, start) %>% distinct %>% 
#     mutate(tss_id=paste0('TSS_', 1:nrow(.)))

same_start <-  gfc_gtf_full %>% filter(exon_number==1) %>%
  select(seqid, strand, start, end, gene_name) %>% distinct %>%  
  group_by(seqid, strand ,start) %>% 
  summarise(count=n(), max_end=max(end), s_gene_name=first(gene_name)) %>% filter( count >1)

uniq_starts <- gfc_gtf_full %>% filter(exon_number==1) %>%
  select(seqid, strand, start, end, gene_name) %>% distinct %>% 
  anti_join(same_start) %>% 
  bind_rows(., same_start %>% select(seqid, strand, start, end=max_end, gene_name=s_gene_name))

multi_start_genes <- uniq_starts %>% group_by(gene_name) %>% summarise(count=n()) %>% filter(count >1) %>% pull(gene_name)
uniq_start_multi_gene <- novel_exons %>% mutate(novel_start=nvl_start) %>% filter(novel_start) %>% 
  select(seqid, strand, start, novel_start) %>% distinct %>%  
  left_join(uniq_starts, .) %>% mutate(novel_start=replace_na(novel_start, F)) %>% filter(gene_name %in% multi_start_genes)



terminal_exons <-  gfc_gtf_full %>% filter(type == 'exon') %>% group_by(transcript_id) %>% 
  summarise(seqid=last(seqid), strand=last(strand), start=last(start), end=last(end), gene_name=last(gene_name)) %>% 
  select(-transcript_id) %>%  distinct

same_ends <- terminal_exons %>% group_by(seqid, strand, end) %>% 
  summarise(min_start=min(start), count=n(), s_gn=first(gene_name)) %>% filter(count>1)
uniq_ends <- terminal_exons %>% anti_join(same_ends) %>% 
  bind_rows(same_ends %>% select(seqid, strand, start=min_start, end, gene_name =s_gn))
multi_end_genes <- uniq_ends %>% group_by(gene_name) %>% summarise(count=n()) %>% filter(count>1) %>% pull(gene_name)

uniq_ends_multi_gene <- novel_exons %>% mutate(novel_end=nvl_end) %>% select(seqid, strand, end, novel_end) %>% 
  filter(novel_end) %>% distinct %>% 
  left_join(uniq_ends, .) %>% mutate(novel_end=replace_na(novel_end, F)) %>% filter(gene_name %in% multi_end_genes)

novel_exons_TSES <- novel_exons %>% left_join( uniq_start_multi_gene %>% select(seqid, strand, start, novel_start) %>% distinct) %>% 
  left_join(uniq_ends_multi_gene %>% select(seqid, strand, end, novel_end)) %>% rename(novelTSS=novel_start, novelTES=novel_end)
novel_exons_TSES[is.na(novel_exons_TSES)] <- F

novel_exons_TSES <- novel_exons_TSES %>% mutate(nv_type_rc = case_when(novelTSS ~ 'novel_TSS',
                                                                       novelTES ~ 'novel_TES',
                                                                       TRUE ~ nv_type))
save(uniq_start_multi_gene,all_exons, all_transcripts, novel_exons_TSES,  
     uniq_ends_multi_gene, novel_loci_distinct, novel_transcripts, file = outfile)
#now lets make a formated_gtf that we can use for everything
##this should have: transcript_type:pc, or not pc., exon_coding as well as is.novel exon, 
##novel exon type, first or last exon, single exon

### is the novel exon in the cds part 3 -  only check if the novel regions of novel exons are actually in the CDS
gff3 <- rtracklayer::readGFF(gff3_file) %>% as_tibble %>%
  mutate(ID=str_extract(ID,'TCONS_[0-9]+|ENSG[0-9]+'), type=as.character(type))
cds_bed <- gff3 %>% filter(type == 'CDS') %>% 
  group_by(ID) %>% 
  summarise(seqid=first(seqid), strand=first(strand), start=min(start), end=max(end)) %>% 
  mutate(score=999) %>% 
  select(seqid, start, end, transcript_id=ID, score, strand) %>% 
  from_data_frame %>% 
  RBedtools('sort', i=.)
exon_bed <- all_exons %>% mutate(score=888) %>% 
  select(seqid, start, end, origin, score, strand) %>% 
  from_data_frame %>% 
  RBedtools('sort', i=.)

novel_exons_tx <- novel_exons_TSES %>% 
  inner_join(gfc_gtf %>% filter(type == 'exon') %>% select( seqid,strand, start, end, transcript_id ))
novel_seq_bed <- novel_exons_tx  %>% 
  mutate(score= 777, cp_id=paste(id, transcript_id, sep = ';')) %>% # join novel exon id and txid
  select(seqid, start, end, cp_id, score, strand) %>% 
  from_data_frame %>% 
  RBedtools('sort',output = 'stdout', i=.) %>% 
  RBedtools('subtract', options = '-s',output = 'stdout',  a=., b=exon_bed) %>% 
  RBedtools('intersect', options = '-s -wo', a=., b=cds_bed) %>%
  to_data_frame


exon_locations <- novel_seq_bed %>% 
  mutate(transcript_id= str_split(X4, ';') %>% sapply(function(x) x[2]),
         id=str_split(X4, ';') %>% sapply(function(x) x[1])) %>% 
  filter(transcript_id == X10) %>% 
  select(id, transcript_id) %>% 
  mutate(exon_location='CDS') %>% 
  left_join(novel_exons_tx,.) %>% 
  mutate(exon_location= case_when(is.na(exon_location) & (!transcript_id %in% gff3$ID) ~ 'NC', 
                              is.na(exon_location) & (transcript_id %in% gff3$ID)  ~ 'UTR',
                              TRUE ~ exon_location)
         ) %>% distinct

tcons2mstrg <- gfc_gtf_full %>% filter(type == 'transcript') %>%  select(transcript_id, oId) %>% distinct
exon_info <- novel_exons_TSES %>% select(seqid, strand, start, end, novel_exon_id=id, novel_exon_type=nv_type_rc)
last_exons <- gfc_gtf %>% filter(type == 'exon') %>% select(seqid, strand, start, end, transcript_id) %>% 
  group_by(transcript_id) %>% 
  summarise(seqid=last(seqid), strand=last(strand), start=last(start), end=last(end)) %>%
  mutate(is.last='TRUE')
single_exons <- gfc_gtf %>% filter(type == 'exon') %>% group_by(transcript_id) %>% summarise(count=n()) %>% 
  filter(count == 1) %>% pull(transcript_id)

complete_gtf <- gfc_gtf %>% 
  select(-oId, -class_code, -tss_id, -contained_in, -cmp_ref, -cmp_ref_gene) %>% #remove junk
  left_join(tcons2mstrg) %>% #add fixed oId
  mutate(transcript_type= ifelse(transcript_id %in% gff3$ID, 'protein_coding', 'noncoding'), 
         type=as.character(type)) %>% # add transcript type
  left_join(exon_locations) %>% # add exon type - pc/nc/utr/  
  left_join(last_exons) %>% # mark which exons are last exons
  left_join(exon_info) %>% #novel exon type - TSS/TES, single exons
  mutate(is.last=replace_na(is.last, 'FALSE'),is.singleExon=ifelse(transcript_id %in% single_exons, 'TRUE', 'FALSE') )

write_gtf3(complete_gtf, gtf_ano_outfile)                                  







