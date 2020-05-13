library(tidyverse)
library(argparse)
library(data.table)
library(RBedtools)
library(parallel)
library(yaml)
library(glue)

parser <- ArgumentParser()
parser$add_argument('--workingDir', action = 'store',dest = 'working_dir')
parser$add_argument('--fileYaml', action = 'store', dest = 'file_yaml')
parser$add_argument('--agatGff', action  = 'store', dest = 'agat_gff3_file')
list2env(parser$parse_args(), .GlobalEnv)

source('~/scripts/write_gtf.R')
setwd(working_dir)
files <- read_yaml(file_yaml)
filt_gtf_path <- files$raw_gtf_path
final_gtf_path <-files$final_gtf_path
######
# working_dir <- '/data/swamyvs/ocular_transcriptomes_pipeline/old_data/'
# file_yaml <- '/data/swamyvs/ocular_transcriptomes_pipeline/files.yaml'
# agat_gff3_file <- 'data/vep/gtf_startstop_added.gff'
###### 
#read in agat gff, compare its annotated start stops to actual cds,  remove cds that shirnk/grow too much and then 
# appropriately reclassify everything 
agat_gff3 <- rtracklayer::readGFF(agat_gff3_file) %>% as_tibble %>% 
    unnest(Parent) %>%  
    mutate(Parent = as.character(Parent), 
           ID = replace(ID, type%in%c('start_codon', 'stop_codon'),Parent[type%in%c('start_codon', 'stop_codon')]),
           type = as.character(type), 
           transcript_id = str_extract_all(ID, 'DNTX_\\d+') %>% sapply(function(x) x[1]) ,
           exon_id = str_remove_all(ID,'p\\d+\\.'), 
           exon_id = replace(exon_id, !grepl('exon', exon_id), NA))
    

gff3_has_valid_start_stop <- agat_gff3  %>% 
    group_by(transcript_id) %>% 
    summarise(has_start_stop = sum(type %in% c('stop_codon', 'start_codon'))) %>% 
    filter(has_start_stop == 2)
# calcualte cds lenghhs
base_dntx_gtf <- rtracklayer::readGFF(files$base_all_tissue_gtf) %>% 
    mutate(exon_id = paste0(transcript_id,'.',type, exon_number), 
           exon_id = replace(exon_id, !grepl('exon', exon_id), NA))  

ref_gtf <- rtracklayer::readGFF(files$ref_GTF)

ref_cds_lengths <- ref_gtf %>% 
    filter(transcript_type == 'protein_coding', type == 'CDS') %>% 
    mutate(length = end-start) %>% 
    group_by(gene_name, transcript_id) %>% 
    summarise(length =  (sum(length )+n()) / 3) %>%
    group_by(gene_name) %>% 
    summarise(max_length = max(length), min_length = min(length) )
tx2gene <- base_dntx_gtf %>% filter(type == 'transcript') %>% select(transcript_id, gene_id, gene_name) %>% distinct %>% 
    mutate(gene_name = replace(gene_name, is.na(gene_name), transcript_id[is.na(gene_name)]))
co <- 200
dntx_lengthcomp_fail <- agat_gff3 %>% 
    filter(transcript_id %in% gff3_has_valid_start_stop$transcript_id, type == 'CDS') %>% 
    mutate(length = end-start) %>% 
    group_by(transcript_id) %>% 
    summarise(length = (sum(length )+n()) / 3)%>% 
    inner_join(tx2gene,.) %>% 
    inner_join(ref_cds_lengths) %>% 
    mutate(max_diff = length - max_length, min_diff = length - min_length) %>% 
    filter(max_diff >co | min_diff < -co)

gff3_valid_clean <- agat_gff3 %>% 
    filter(type !='gene', 
           transcript_id %in% gff3_has_valid_start_stop$transcript_id, 
          !transcript_id %in% dntx_lengthcomp_fail$transcript_id) %>% 
    mutate(transcript_type = 'protein_coding', 
           transcript_biotype = 'protein_coding',
           type = replace(type, type == 'mRNA', 'transcript')) %>% 
    select(-Name, -Parent, -ID)
nc_gtf <- base_dntx_gtf %>% 
    filter(!transcript_id %in% gff3_valid_clean$transcript_id) %>%
    mutate(transcript_type = 'noncoding',
           transcript_biotype = 'noncoding') %>% 
    select(colnames(gff3_valid_clean))


tx2code <- base_dntx_gtf %>% filter(type == 'transcript') %>% select(transcript_id, class_code, oId) %>% distinct 
clean_gtf <- bind_rows(gff3_valid_clean, nc_gtf) %>% 
    inner_join(tx2gene) %>% 
    inner_join(tx2code) %>% 
    mutate(exon_number = str_split(exon_id, 'exon') %>% sapply(function(x) x[2]))


gene_gtf <- clean_gtf %>% filter(type == 'transcript') %>%  group_by(gene_name) %>% 
    summarise(seqid = first(seqid), source = first(source), type='gene', start = min(start), end  = max(end), 
              score = first(score), strand  = first(strand), phase = first(phase),gene_id = first(gene_id))
final_gtf <- bind_rows(clean_gtf, gene_gtf) %>% 
    mutate(seqid = factor(seqid, levels = levels(base_dntx_gtf$seqid))) %>% 
    arrange(seqid, start) 

gene_order <- unique(final_gtf$gene_name)
transcript_order <- unique(final_gtf$transcript_id)
type_order <- c('gene', 'transcript', 'exon', 'five_prime_UTR', 'start_codon', 'CDS', 'stop_codon', 'three_prime_UTR')
final_gtf_sorted <- final_gtf %>% 
    mutate(gene_name = factor(gene_name, levels = gene_order), 
           transcript_id = replace_na(transcript_id, '^'),
           transcript_id = factor(transcript_id, levels = c('^', transcript_order)),
           type = factor(type, levels = type_order)
           ) %>% 
    arrange(gene_name, transcript_id, type) %>% 
    select(-exon_id) %>% 
    mutate(transcript_id = as.character(transcript_id), 
           transcript_id = replace(transcript_id, transcript_id == '^', NA))

format_tissue_specific_info <- function(ctab, s_subtissue, t_gtf){
    master_ctab <- ctab %>% 
        select(transcript_id, !!s_subtissue) %>% 
        filter(!(is.na(.[,s_subtissue])))
    tissue_conv <- glue('{filt_gtf_path}{s_subtissue}.convtab') %>%
        fread(sep ='\t') %>% 
        as_tibble %>% rename(!!s_subtissue := transcript_id) %>% 
        inner_join(master_ctab, .) %>% 
        select(-(!!s_subtissue))
    tissue_gtf <- t_gtf %>% filter(transcript_id %in% tissue_conv$transcript_id)
    final_gtf_file <- glue('{final_gtf_path}{s_subtissue}.gtf')
    write_gtf3(df = tissue_gtf, final_gtf_file)
    det_df <- apply(tissue_conv[,-(1:4)],2, function(x) !is.na(x)) %>% as_tibble %>%  bind_cols(tissue_conv[,1], .)
    fwrite(det_df,glue('{final_gtf_path}{s_subtissue}.detdf'), sep = '\t')
    fwrite(tissue_conv,glue('{final_gtf_path}{s_subtissue}.convtab'), sep = '\t')
}

sample_table <- read_tsv(files$sample_table_file)
conv_tab <- fread(files$all2tissue_convtab, sep = '\t') %>% as_tibble
subtissues <- unique(sample_table$subtissue)
mclapply(subtissues, function(x) format_tissue_specific_info(conv_tab, x, final_gtf_sorted),
         mc.cores = min(length(subtissues), parallel::detectCores() - 10) ) 


############################################################
####### now run the analysis to make the master gtf ########
############################################################

load(files$ref_tx_exon_rdata)
gfc_gtf <- final_gtf_sorted

novel_transcripts <- anti_join(gfc_gtf %>% filter(type == 'exon'), all_exons) %>% 
    filter(!grepl('DNTX', gene_name)) %>% 
    pull(transcript_id) %>% {filter(gfc_gtf, type == 'transcript', transcript_id %in% .)}
built_ref_tx <- filter(gfc_gtf, type == 'transcript', class_code == '=') %>% inner_join(all_transcripts) 
novel_loci <- anti_join(gfc_gtf %>% filter(type == 'transcript'), all_transcripts) %>% filter(class_code =='u')

novel_single_exon_tx <- novel_transcripts$transcript_id %>% {filter(gfc_gtf, transcript_id  %in% .)} %>% group_by(transcript_id) %>%
    summarise(count=n()) %>% filter(count == 2) %>% pull(transcript_id) %>%  {filter(novel_transcripts, transcript_id %in% .)}
novel_transcripts <- filter(novel_transcripts, !transcript_id %in% novel_single_exon_tx$transcript_id)

# remove novel loci that overlap with known genes, padded 5kb upstream and downstreqam, and overlap known repeat regions
## 
pad=2000
all_transcript_bed <- all_transcripts %>%
    mutate(score=999, start=start - pad, end= end+pad,start=replace(start, start<0, 0)) %>% 
    select(seqid, start, end, origin, score, strand) %>% 
    from_data_frame %>% 
    RBedtools('sort',  i=.) 

no_intersect <- novel_loci %>%
    mutate(score=999) %>% 
    select(seqid, start, end, transcript_id, score, strand) %>% 
    from_data_frame %>% 
    RBedtools('sort', output = 'stdout', i=.) %>% 
    RBedtools('intersect',options = '-v -s',  a=., b=all_transcript_bed) %>% 
    #RBedtools('intersect', options = '-v -s', a=.,b=repeat_bed_file) %>% 
    to_data_frame



novel_loci_distinct <- novel_loci %>% filter(transcript_id %in% no_intersect$X4)
novel_loci_failed <- novel_loci %>% filter(!transcript_id %in% no_intersect$X4)

write(novel_loci_distinct$transcript_id, file = files$novel_loci_txids, sep = '\n')


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


# 
# 
gfc_gtf_ano <- filter(gfc_gtf, !transcript_id %in% novel_loci$transcript_id)
gfc_gtf_ref <- filter(gfc_gtf_ano, !transcript_id %in% novel_transcripts$transcript_id)

gfc_gtf_full <-  gfc_gtf_ano %>% filter(transcript_id  %in% novel_transcripts$transcript_id) %>% select(seqid, strand, start, end) %>%
    distinct %>% anti_join(gfc_gtf_ref) %>% anti_join(all_exons) %>% anti_join(all_transcripts) %>%
    mutate(is.novel=T) %>% left_join(gfc_gtf_ano, .) %>% mutate(is.novel=replace_na(is.novel, F))
# 

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



terminal_exons <-  gfc_gtf_full %>% 
    filter(type == 'exon') %>% 
    group_by(transcript_id) %>% 
    summarise(seqid=last(seqid), strand=last(strand), start=last(start), end=last(end), gene_name=last(gene_name)) %>% 
    select(-transcript_id) %>% distinct

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

save(uniq_start_multi_gene, 
     all_exons,
     all_transcripts, 
     novel_exons_TSES,  
     uniq_ends_multi_gene, novel_loci_distinct, novel_transcripts, file = files$exon_class_rdata)
#now lets make a formated_gtf that we can use for everything
##this should have: transcript_type:pc, or not pc., exon_coding as well as is.novel exon, 
##novel exon type, first or last exon, single exon


### is the novel exon in the cds part 3 -  only check if the novel regions of novel exons are actually in the CDS
# gff3 <- rtracklayer::readGFF(gff3_file) %>% as_tibble %>%
#     mutate(ID=str_extract(ID,'DNTX_[0-9]+|ENSG[0-9]+'), type=as.character(type))
# cds_bed <- gff3 %>% filter(type == 'CDS') %>% 
#     group_by(ID) %>% 
#     summarise(seqid=first(seqid), strand=first(strand), start=min(start), end=max(end)) %>% 
#     mutate(score=999) %>% 
#     select(seqid, start, end, transcript_id=ID, score, strand) %>% 
#     from_data_frame %>% 
#     RBedtools('sort', i=.)
# exon_bed <- all_exons %>% mutate(score=888) %>% 
#     select(seqid, start, end, origin, score, strand) %>% 
#     from_data_frame %>% 
#     RBedtools('sort', i=.)
# 
# novel_exons_tx <- novel_exons_TSES %>% 
#     inner_join(gfc_gtf %>% filter(type == 'exon') %>% select( seqid,strand, start, end, transcript_id ))
# novel_seq_bed <- novel_exons_tx  %>% 
#     mutate(score= 777, cp_id=paste(id, transcript_id, sep = ';')) %>% # join novel exon id and txid
#     select(seqid, start, end, cp_id, score, strand) %>% 
#     from_data_frame %>% 
#     RBedtools('sort',output = 'stdout', i=.) %>% 
#     RBedtools('subtract', options = '-s',output = 'stdout',  a=., b=exon_bed) %>% 
#     RBedtools('intersect', options = '-s -wo', a=., b=cds_bed) %>%
#     to_data_frame
# 
# 
# exon_locations <- novel_seq_bed %>% 
#     mutate(transcript_id= str_split(X4, ';') %>% sapply(function(x) x[2]),
#            id=str_split(X4, ';') %>% sapply(function(x) x[1])) %>% 
#     filter(transcript_id == X10) %>% 
#     select(id, transcript_id) %>% 
#     mutate(exon_location='CDS') %>% 
#     left_join(novel_exons_tx,.) %>% 
#     mutate(exon_location= case_when(is.na(exon_location) & (!transcript_id %in% gff3$ID) ~ 'NC', 
#                                     is.na(exon_location) & (transcript_id %in% gff3$ID)  ~ 'UTR',
#                                     TRUE ~ exon_location)
#     ) %>% distinct

tcons2mstrg <- gfc_gtf_full %>% filter(type == 'transcript') %>%  select(transcript_id, oId) %>% distinct
exon_info <- novel_exons_TSES %>% select(seqid, strand, start, end, novel_exon_id=id, novel_exon_type=nv_type_rc)
last_exons <- gfc_gtf %>% filter(type == 'exon') %>% 
    select(seqid, strand, start, end, transcript_id) %>% 
    group_by(transcript_id) %>% 
    summarise(seqid=last(seqid), strand=last(strand), start=last(start), end=last(end)) %>%
    mutate(is.last='TRUE')
single_exons <- gfc_gtf %>% filter(type == 'exon') %>% group_by(transcript_id) %>% summarise(count=n()) %>% 
    filter(count == 1) %>% pull(transcript_id)

first_exon_df <- gfc_gtf %>% 
    filter(type == 'exon', exon_number == 1) %>% 
    select(seqid, strand, start, end, transcript_id) %>% 
    mutate(is.first = T)


complete_gtf <- gfc_gtf %>% 
    left_join(tcons2mstrg) %>% #add fixed oId
    mutate(type=as.character(type)) %>% # add transcript type
    #left_join(exon_locations) %>% # add exon type - pc/nc/utr/  
    left_join(last_exons) %>% # mark which exons are last exons
    left_join(exon_info) %>% #novel exon type - TSS/TES, single exons
    left_join(first_exon_df) %>% 
    mutate(is.last=replace_na(is.last, 'FALSE'),is.singleExon=ifelse(transcript_id %in% single_exons, 'TRUE', 'FALSE'), 
           is.first = replace_na(is.first, 'FALSE'))

write_gtf3(complete_gtf, files$anno_all_tissue_gtf)                                  







