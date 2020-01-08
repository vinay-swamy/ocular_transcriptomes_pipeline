library(tidyverse)
library(matrixStats)
args <- c('/Volumes/data/ocular_transcriptomes_pipeline/',
          'data/gtfs/all_tissues.combined_NovelAno.gtf',
          'data/misc/TCONS2MSTRG.tsv',
          'data/all_tissue_quant.Rdata',
          'sampleTableFull.tsv',
          '/Volumes/data/ocular_transcriptomes_pipeline/data/misc/final_dd/REPLACE.dd.tsv')


working_dir <- args[1]
gtf_file <- args[2]
tc2m_file <- args[3]
quant_file <- args[4]
sample_file <- args[5]
dd_stem <- args[6]
outfile <- args[7]
setwd(working_dir)

process_det_files <- function(det_file, tissue){
    df <- read_tsv(det_file) %>% select(-refid, -gene_id, -code) %>% select(transcript_id, everything())
    num_det <- df[,-1] %>% rowSums() %>% {tibble(transcript_id=df$transcript_id, !!tissue :=. )} 
    return(num_det)
}

calc_isoform_percentage <- function(t_tissue){
    df <- counts_by_tissue %>% select(transcript_id, gene_name, new_tx_id, !!t_tissue)
    tt_sym <- as.symbol(t_tissue)
    df_gene_sums <- df %>% 
        select(-transcript_id) %>% 
        group_by(gene_name) %>%  
        summarise(gene_sums:= sum(!!tt_sym)) %>% 
        left_join(df, .) %>% 
        mutate(piu = .[[t_tissue]] / .[['gene_sums']], !!t_tissue :=NULL ) %>% 
        select(transcript_id, gene_name, new_tx_id, !!t_tissue:=piu)
    return(df_gene_sums)
    
}

replace_nan <- function(df) {
    df_fixed <- lapply(colnames(df),function(col) pull(df, col) %>%  
                           {replace(., is.nan(.), 0)}) %>% bind_cols %>% as_tibble
    colnames(df_fixed) <- colnames(df)
    return(df_fixed)
    
}


gtf <- rtracklayer::readGFF(gtf_file)
gtf <- gtf %>% mutate(new_tx_id= replace(oId, !grepl('ENST', oId), transcript_id[!grepl('ENST', oId)]))
tc2m <- read_tsv(tc2m_file)
sample_table <- read_tsv(sample_file) %>% filter(subtissue != 'synth')
subtissues <- unique(sample_table$subtissue)[1:5]#ONLY FOR TESTING PURPOSES
t2g <- gtf %>% filter(type == 'transcript') %>% select(transcript_id, gene_name,new_tx_id) %>% distinct

all_det <- lapply(subtissues, function(tis) 
                    gsub(pattern = 'REPLACE', replacement = tis, x = dd_stem) %>% {process_det_files(.,tis)}) %>% 
    reduce(full_join) %>% inner_join(t2g, .)
all_det[is.na(all_det)] <- 0    


samples_per_tissue <-  sample_table %>% group_by(subtissue) %>% summarise(count=n()) %>% filter(subtissue %in% subtissues)# REMOVE ME
samples_per_tissue_list <- 1/samples_per_tissue$count
names(samples_per_tissue_list) <- samples_per_tissue$subtissue

frac_samp_det <- as.matrix(all_det[,names(samples_per_tissue_list)]) %*% diag(samples_per_tissue_list) 
colnames(frac_samp_det) <- names(samples_per_tissue_list)
frac_samp_det <- frac_samp_det %>% as_tibble() %>% bind_cols(all_det[,c('transcript_id','gene_name', 'new_tx_id')], .)



load(quant_file)
all_quant[is.na(all_quant)] <- 0
counts_by_tissue <- lapply(subtissues,
                               function(tis) filter(sample_table, subtissue == tis) %>% pull(sample) %>%
                                   {all_quant[,c('transcript_id', .)]} %>%
                                   mutate(!!tis := rowMedians(.[,-1] %>% as.matrix)) %>%
                                   select(transcript_id, !!tis)
                           ) %>%
    reduce(left_join) %>% inner_join(t2g, .)

piu_raw <- lapply(colnames(counts_by_tissue)[-(1:3)], calc_isoform_percentage) %>% reduce(left_join)
piu <-replace_nan(piu_raw)



# fix transcript id column
gtf <- gtf %>% select(-transcript_id) %>% rename(transcript_id = new_tx_id)
frac_samp_det <- frac_samp_det %>% select(-transcript_id) %>% rename(transcript_id = new_tx_id)
all_det <-  all_det %>% select(-transcript_id) %>% rename(transcript_id = new_tx_id)
piu <- piu %>% select(-transcript_id) %>% rename(transcript_id = new_tx_id)
tc2m <- tc2m %>% inner_join(t2g, .) %>% select(-transcript_id) %>% rename(transcript_id = new_tx_id)
tissue_det <- tc2m %>% select(subtissues) %>% apply(2, function(x) !is.na(x)) %>% as_tibble %>% bind_cols(tc2m[,1:3],.)


save(gtf, all_det, frac_samp_det, piu, tc2m,tissue_det, file=outfile)

# 
# keep_tx <- tc2m %>% select(subtissues) %>% apply(2, function(x) !is.na(x)) %>% rowSums() %>% {. >0} %>% {tc2m[.,'transcript_id']}
# dev_gtf <- gtf %>% filter(transcript_id %in% keep_tx)
# dev_frac_det <- frac_samp_det %>% filter(transcript_id %in% keep_tx)
# dev_all_det <- all_det %>% filter(transcript_id %in% keep_tx)
# dev_piu <- piu %>% filter(transcript_id %in% keep_tx)
# dev_tissue_det <- tissue_det %>% filter(transcript_id %in% keep_tx)
# dev_subtissues <- subtissues
# dev_gene_names <- unique(dev_gtf$gene_name)
# save(dev_gtf,dev_frac_det, dev_all_det, dev_piu, dev_tissue_det,dev_gene_names, dev_subtissues, file = '~/NIH/ocular_transcriptomes_shiny/data/numdet.Rdata')
# 






