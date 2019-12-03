library(tidyverse)
library(tximport)
args <- commandArgs(trailingOnly = T)
wd <- args[1]
path_to_quant <- args[2]
out_exp_file <- args[3]
setwd(wd)
exp_files <- list.files(path_to_quant, pattern = '.RDS', recursive = T, full.names = T)
all_quant <- lapply(exp_files, readRDS) %>% reduce(full_join)
save(all_quant, file = out_exp_file)














