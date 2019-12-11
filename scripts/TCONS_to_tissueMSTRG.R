library(tidyverse)
library(tximport)
library(parallel)
# args <- c('/Volumes/data/eyeintegration_splicing/','sampleTableFull.tsv' ,'data/gffcomp_dir/all_tissues.tracking',
#           'data/gffcomp_dir/all_tissues.combined.gtf', 'ref/gencode_comp_ano.gtf', 'testing/t2m.tsv', 'testing/ogtf.gtf')

args <- commandArgs(trailingOnly = T)
wd <- args[1]
sample_table_file <- args[2]
track_file <- args[3]
raw_gtf_file <- args[4]
ref_gtf_file <- args[5]
TCONS_to_MSTRG_file <- args[6]
out_gtf_file <- args[7]
#save(args, file='/tmp/smoobargs.rdata')
source('~/scripts/write_gtf.R')
# nm_col <- function(col){
#     col=col[col!='-']
#     name=str_split(col[1], '_0|\\|')[[1]][3] %>% str_split('_MSTRG') %>% .[[1]] %>% .[1]
#     return(name)
# }

nm_col_clean <- function(col){
    raw_name <- col %>% .[.!='-'] %>% .[!grepl('ENST', .)] %>% head
    name <- str_split(raw_name, '\\|')[[1]][2] %>% str_split('_0') %>% .[[1]] %>% .[1]
}#slightly slower version of nm_col, but works for all cases



process_columns <- function(tab,col_name){
    print(col_name)
    tab <- tab %>% filter( (!! rlang::sym(col_name)) != '-')
    col <- tab %>% pull(!!col_name)
    name_col <- tab %>% pull(transcript_id)
    det <- suppressWarnings(str_split(col, ':|\\|'))
    z <- '
    most of the oIds have a 1:1 mapping to tx ids, so split into 2 chunks for faster run time, like a lot faster
    _s denotes simple case, _c dentotes complex case
    '

    d_lengths <- sapply(det, length)
    base <- min(d_lengths)
    simple <- d_lengths == base
    det_s <- det[simple]
    name_col_s <- name_col[simple]
    tx_simple <-  lapply(1:length(det_s), function(i)  det_s[[i]][3] %>%
                c(name_col_s[i], .)) %>%
         do.call(rbind, .) %>% as.data.frame(stringsAsFactors=F)
    colnames(tx_simple) <- c('transcript_id', col_name)
    #%>% rename(!!col_name:=oId) %>% distinct
    tx_comp=data.frame()
    if(sum(!simple) > 0){
        det_c <- det[!simple]
        name_col_c <- name_col[!simple]

        tx_comp <- lapply(1:length(det_c), function(i)  det_c[[i]][-1] %>%
                                .[grepl('MSTRG\\.\\d+\\.\\d+', .) | grepl('ENST', .)] %>%
                                {tibble(transcript_id=rep(name_col_c[i], length(.)),oId= . )}) %>%
            bind_rows() %>% rename(!!col_name:=oId) %>% distinct
        }
    return(list(simple=tx_simple, comp=tx_comp))
}


setwd(wd)
sample_table <- read_tsv(sample_table_file) %>% filter(subtissue != 'synth')
track_tab <- read_tsv(track_file, col_names = F) 
tissues <- unique(sample_table$subtissue)
names <- c('transcript_id', 'gene_id','refid', 'class_code',  apply(track_tab[,-(1:4)], 2, nm_col_clean))

stopifnot(sum(!tissues %in% names) == 0)#check

colnames(track_tab) <- names
cn <- colnames(track_tab)[-(1:4)]
tcons2mstrg <- mclapply(cn, function(col) process_columns(track_tab,col), mc.cores = 6)
tc2mstrg_simple <- lapply(tcons2mstrg, function(x) x[['simple']]) %>% reduce(full_join)

gfc_gtf <- rtracklayer::readGFF(raw_gtf_file)
ref_gtf <- rtracklayer::readGFF(ref_gtf_file)
ref_gtf_tx <- ref_gtf %>% filter(type == 'transcript') %>% select(seqid, strand, start, end)
gfcgtf_reftx_absmatch <- gfc_gtf %>% filter(type  == 'transcript', class_code == '=') %>% 
    select(seqid, strand, start, end, transcript_id) %>% inner_join(ref_gtf_tx) %>% pull(transcript_id)



'

NOTE:
gffcompare marks exact intron chain matches as `=`. This allows for mutliple TSS/TES, some of which may be novel. 
Therefore, the true exact match tx are those that have `=` as the class code, and start and end at the same location as
a reference transcript. So I have to extract transcts that have `=` as the code, and changed them to `*`, but the `*`
changed to `+` as part of a later step.

'



max_number <- gfc_gtf$transcript_id %>% 
    unique %>% 
    str_split('_') %>% 
    sapply(function(x) x[2]) %>% 
    as.numeric() %>% 
    max

tx2code <- gfc_gtf %>% 
    filter(type == 'transcript') %>% 
    select(transcript_id,cmp_ref,class_code, oId, gene_name) %>%
    distinct %>% 
    mutate(new_cmp_ref=replace(cmp_ref, is.na(cmp_ref), transcript_id[is.na(cmp_ref)]), 
           new_class_code=case_when(class_code == '=' & transcript_id %in% gfcgtf_reftx_absmatch ~ '=',
                                    class_code == '=' & !transcript_id %in% gfcgtf_reftx_absmatch ~ '*',
                                    TRUE ~ class_code
                                    ),
           new_transcript_id= replace(transcript_id, class_code =='*', 
                                      paste0('DNTX_0', seq(max_number+1,max_number+1+sum(class_code =='*')) )   
                                      ),
           new_gene_name= replace(gene_name, is.na(gene_name), transcript_id[is.na(gene_name)])
           ) %>% 
    select(-cmp_ref, -class_code, -gene_name)

column_order <- colnames(gfc_gtf)

final_gtf <- gfc_gtf %>% 
    select(-oId, -gene_name) %>%  
    left_join(tx2code) %>% 
    select(-class_code, -transcript_id, -cmp_ref) %>% 
    rename(class_code=new_class_code, cmp_ref=new_cmp_ref, 
           transcript_id=new_transcript_id, gene_name=new_gene_name ) %>% 
    .[,column_order]

tcons2mstrg_complete <- tx2code %>% select(transcript_id, new_transcript_id, new_class_code) %>% left_join(tc2mstrg_simple,.) %>% 
    select(-transcript_id) %>% rename(transcript_id=new_transcript_id,class_code= new_class_code ) %>% 
    select(transcript_id,class_code, everything())



write_tsv(tcons2mstrg_complete, TCONS_to_MSTRG_file)
write_gtf3(final_gtf, out_gtf_file)


