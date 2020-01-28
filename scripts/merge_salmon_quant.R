library(tidyverse)
library(tximport)
library(argparse)
# args <- commandArgs(trailingOnly = T)

parser <- ArgumentParser()
parser$add_argument('--workingDir', action = 'store', dest = 'wd')
parser$add_argument('--pathToQuant', action = 'store', dest = 'path_to_quant')
parser$add_argument('--outExpFile', action = 'store', dest = 'all_quant')
list2env(parser$parse_args(), .GlobalEnv)
# wd <- args[1]
# path_to_quant <- args[2]
# out_exp_file <- args[3]
setwd(wd)
exp_files <- list.files(path_to_quant, pattern = '.RDS', recursive = T, full.names = T)
all_quant <- lapply(exp_files, readRDS) %>% reduce(full_join)
save(all_quant, file = out_exp_file)














