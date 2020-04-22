library(tidyverse)
library(tximport)
library(glue)
library(argparse)
# args <- commandArgs(trailingOnly = T)

parser <- ArgumentParser()
parser$add_argument('--workingDir', action = 'store', dest = 'wd')
parser$add_argument('--pathToQuant', action = 'store', dest = 'path_to_quant')
parser$add_argument('--sampleTable', action = 'store', dest = 'sample_table_file')
parser$add_argument('--outDir', action = 'store', dest = 'out_dir')
list2env(parser$parse_args(), .GlobalEnv)
# path_to_quant <- 'data/salmon_quant/dntx/'
# sample_table_file <- 'sampleTableFullV3.tsv'
# out_dir <- 'data/'
source('scripts/read_salmon.R')
setwd(wd)

sample_table <- read_tsv(sample_table_file)

all_quant <- lapply(unique(sample_table$subtissue), function(x) read_salmon(glue('{path_to_quant}{x}/'),
                                                                            which_return = 'transcript', 
                                                                            quant_type = 'abundance')) #%>% 
all_quant <- all_quant %>% reduce(full_join)
save(all_quant, file = glue('{out_dir}/all_tissue_quant.Rdata'))
eye_samples <- sample_table %>% filter(tissue %in% c('Retina', 'RPE', 'Cornea')) %>% pull(sample)
eye_quant <- all_quant %>% select(transcript_id, eye_samples) %>% 
    filter(apply(.[,-1], 2,function(x) !is.na(x)) %>% {rowSums(.) > 0})
save(eye_quant, file = glue('{out_dir}/pan_eye_quant.Rdata'))













