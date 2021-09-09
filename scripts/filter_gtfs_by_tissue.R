library(tidyverse)
library(data.table)
library(parallel)
library(glue)
library(argparse)
parser <- ArgumentParser()
parser$add_argument('--workingDir', action = 'store', dest = 'working_dir')
parser$add_argument('--sampleTable', action = 'store', dest = 'sample_table_file')
parser$add_argument('--subtissue', action = 'store', dest = 't_subtissue')
parser$add_argument('--gtfDir', action = 'store', dest = 'prefix_for_gtfs')
parser$add_argument('--stringtieQuant', action = 'store', dest = 'path_to_stringtie_gtfs')
##########
# working_dir='/data/swamyvs/ocular_transcriptomes_pipeline/'
# sample_table_file="sampleTableFullV3.tsv"             
# t_subtissue="Thyroid"             
# prefix_for_gtfs="data/gtfs/raw_tissue_gtfs/Thyroid"
# path_to_stringtie_gtfs="st_out/"
#############3    

list2env(parser$parse_args(), .GlobalEnv)


source('scripts/write_gtf.R')
source('scripts/shared_funcs.R')

## create input and output filenames
track_file <- glue('{prefix_for_gtfs}.tracking') ## generated with gffcompare
raw_merged_gtf <- glue('{prefix_for_gtfs}.combined.gtf') ## generated with gffcompare
filtered_convtab_file <- glue('{prefix_for_gtfs}.convtab')
filtered_gtf_file <-  glue('{prefix_for_gtfs}.combined.filtered.gtf')



setwd(working_dir)


## read in gtf output of stringtie, to pull out a matrix of transcript_id, TPM
## TPM column will be sample name, ie SRSXXXX
read_stringtie <- function(file, name ){
    df <- data.table::fread(file, sep='\t', header = F) %>% as_tibble %>% filter(V3 == "transcript")
    if (nrow(df) == 0) {return(tibble())}
    tibble(
        !!name := str_extract(df$V9, 'transcript_id \".+\";') %>% str_split('\"') %>% sapply(function(x)x[2]),
        TPM=str_extract(df$V9, 'TPM \".+\";') %>% str_split('\"') %>% sapply(function(x)x[2]) %>% as.numeric 
    )
    
}

#### process track file into an id conversion table 

conv_tab <- read_track_file(track_file) ## sourced from scripts/shared_funcs.R

#### read stringtie quantification and output a transcripts x per-sample TPM matrix
### use conv tab to map to common identifiers
## NAs will be generated for transcripts NOT detected in a sample
samples <- read_tsv(sample_table_file) %>% filter(subtissue  == t_subtissue) %>% pull(sample)
st_files <- paste0(path_to_stringtie_gtfs, samples, '.gtf')
all(file.exists(st_files))
stringtie_quant <- lapply(seq_along(st_files), function(i) read_stringtie(st_files[i], samples[i]) %>% 
                              inner_join(conv_tab[c('transcript_id', samples[i])]) %>% 
                              select(-(!!samples[i]))  %>% rename(!!samples[i]:=TPM)) %>% 
    reduce(full_join) %>% 
    select(transcript_id, everything())
## NAs are transcripts not detected in a sample, so set to 0
stringtie_quant[is.na(stringtie_quant)] <- 0

## keep transcripts that have an averge TPM of 1 ( remember this is on a per-tissue basis, so we should see consistency here)
quant_filtered <- stringtie_quant %>% filter(rowMeans(.[,samples])>=1)
conv_tab_filtered <- conv_tab %>% filter(transcript_id %in% quant_filtered$transcript_id)
gtf <- rtracklayer::readGFF(raw_merged_gtf) ## this is from the gffcompare call previously in this rule
filtered_gtf <- gtf %>% filter(transcript_id %in% conv_tab_filtered$transcript_id)
write_gtf3(filtered_gtf, filtered_gtf_file)
fwrite(conv_tab_filtered, file = filtered_convtab_file, sep = '\t')














