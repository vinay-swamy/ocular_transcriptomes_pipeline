library(tidyverse)
library(RBedtools)

# args <- c('/Volumes/data/eyeintegration_splicing/',
#           '/Volumes/data/eyeintegration_splicing/data/gtfs/raw_tissue_gtfs/Retina_Adult.Tissue.combined.gtf',
#           '/Volumes/data/eyeintegration_splicing/ref/gencode_comp_ano.gtf',
#           '/Volumes/data/eyeintegration_splicing/data/gtfs/raw_tissue_gtfs/Retina_Adult.Tissue.tracking',
#           'Retina_Adult.Tissue',
#           '/Volumes/data/eyeintegration_splicing/sampleTableFull.tsv',
#           '/Volumes/data/eyeintegration_splicing/ref/UCSC_grch38_repats.bed',
#           'tout.gtf')

source('~/scripts/write_gtf.R')

args <- commandArgs(trailingOnly = T)

wd <- args[1]
gtf_file <- args[2]
ref_gtf_file <- args[3]
tracking_file <- args[4]
t_tissue <- args[5]
sample_table_file <- args[6]
repeats_file <- args[7]
out_gtf_file <- args[8]
setwd(wd)
#save(args, file='/tmp/args.Rdata')
# nm_col <- function(col){
#     col=col[col!='-']
#     name=str_split(col[1], '\\.\\d|:|\\|')[[1]][2] %>% str_split('_MSTRG') %>% .[[1]] %>% .[1]
#     return(name)
# }
nm_col_clean <- function(col){
    raw_name <- col %>% .[.!='-'] %>% .[!grepl('ENST', .)] %>% head
    name=str_split(raw_name, '\\.\\d|:|\\|')[[1]][2] %>% str_split('_MSTRG') %>% .[[1]] %>% .[1]
    return(name)
}
if(!file.exists('rdata/all_ref_tx_exons.rdata')){
    source('scripts/make_allexonsrdata.R')
} else {load('rdata/all_ref_tx_exons.rdata')}



gtf <- rtracklayer::readGFF(gtf_file)
ref_gtf <- rtracklayer::readGFF(ref_gtf_file)
track_tab <- read_tsv(tracking_file, col_names=F)
names <- c('transcript_id', 'gene_id','refid','code', apply(track_tab[,-(1:4)], 2, nm_col_clean))


if(any(is.na(names)) ){
   stop()
}
colnames(track_tab) <- names

#first, identify de novo transcripts that absolutely match the refernece annotation
ref_gtf_tx <- ref_gtf %>% filter(type == 'transcript')%>% select(seqid, strand, start, end, refid=transcript_id)
gtf_tx <- gtf %>% filter(type == 'transcript') 
gffc_ref_absmatch <- gtf_tx %>% filter(class_code == '=') %>% inner_join(ref_gtf_tx) %>% pull(transcript_id)

# next remove transcripts that dont meet detection req
det_df <- apply(track_tab[,-(1:4)],2, function(x) x!='-') %>% as.data.frame %>%  bind_cols(track_tab[,1:4],.)
num_det_df <-det_df %>% mutate(num_det=rowSums(det_df[,-(1:4)])) %>%   
    select(transcript_id, gene_id, refid, code,num_det) 
sample_table <- read_tsv(sample_table_file)

load('ref/core_tight.Rdata')
#load('/Volumes/data/eyeintegration_splicing/ref/core_tight.Rdata')
sample_Table_studies <- core_tight %>% 
    select(sample=sample_accession, study_accession) %>% distinct %>% left_join(sample_table,.) %>%
    mutate(study_accession= case_when(is.na(study_accession) & tissue == 'RPE' ~ 'OGVFB',
                                      is.na(study_accession) & tissue == 'Retina' ~ 'EMTAB',
                                      T ~ study_accession
    )) %>% filter(subtissue == t_tissue, sample %in% colnames(track_tab))
studies <- unique(sample_Table_studies$study_accession)
co=3
nsamp <- sample_table %>% filter(subtissue == t_tissue) %>% nrow 
nsamp_co <- max(trunc(.1*nsamp), 3)
if(length(studies) >=co){
    # if there are 3 or more studies, keep tx present in at least 3 of them  
    det_in_study <- function(study){
        filter(sample_Table_studies, study_accession %in% study) %>% pull(sample) %>% 
            {select(det_df, transcript_id, .)}  %>% { rowSums(.[,-1])}
    }
    det_by_study <- lapply(studies,det_in_study ) %>% bind_cols() 
    keep_study <- det_by_study %>% {rowSums(. >=1 )  >= co}
    keep_nsamp <- det_by_study %>% {rowSums(.) >=nsamp_co }
    keep_tx <- filter(det_df, keep_study, keep_nsamp) %>% pull(transcript_id)
}else{
    # these are 2 or less studies -  most of these are gtex samples 
    nstudy=length(studies)
    det_in_study <- function(study){
        filter(sample_Table_studies, study_accession %in% study) %>% pull(sample) %>% 
            {select(det_df, transcript_id, .)}  %>% { rowSums(.[,-1])}
    }
    det_by_study <- lapply(studies,det_in_study ) %>% bind_cols() 
    keep_co <- rowSums(det_by_study) >=co
    keep_study <- rowSums(det_by_study >=1) == nstudy
    keep_tx <- filter(det_df, keep_co, keep_study) %>% pull(transcript_id)
    
}










keep_codes <- c('=','+','c','k','m','n','j', 'u')
num_det_df_chess_kc <- filter(num_det_df, transcript_id %in% keep_tx, code %in% keep_codes) %>%
    mutate(code = case_when(code == '=' & transcript_id %in% gffc_ref_absmatch ~ '=',
                            code == '=' & !transcript_id %in% gffc_ref_absmatch ~ '~',
                            T ~ code) )

table(num_det_df_chess_kc$code)


filt_gtf <- num_det_df_chess_kc %>% select(transcript_id,class_code=code, num_det) %>% 
    inner_join(gtf %>% select(-class_code), .)
# remove novel loci that dont meet criteria (with 5kb of a known gene, overlapping a repeat region)
novel_loci <- anti_join(filt_gtf %>% filter(type == 'transcript'), all_transcripts) %>% filter(class_code =='u')
pad=5000
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
    RBedtools('intersect',options = '-v -s', output = 'stdout', a=., b=all_transcript_bed) %>% 
    RBedtools('intersect', options = '-v -s', a=.,b=repeats_file) %>% 
    to_data_frame
novel_loci_failed <- novel_loci %>% filter(!transcript_id %in% no_intersect$X4) %>% pull(transcript_id)

filt_gtf <- filt_gtf %>% filter(!transcript_id %in% novel_loci_failed)


tc2oid <- filt_gtf %>% filter(type == 'transcript') %>% 
    mutate(new_id=replace(transcript_id, class_code == '=', cmp_ref[class_code == '='])) %>% 
    select(transcript_id, new_id, class_code, gene_name, oId)
final_gtf <- filt_gtf %>% select(-class_code, -gene_name, -oId) %>% left_join(tc2oid) %>% select(-transcript_id) %>% 
    rename(transcript_id=new_id)

write_gtf3(final_gtf, out_gtf_file)







