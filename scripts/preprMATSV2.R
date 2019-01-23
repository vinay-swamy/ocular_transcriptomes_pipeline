library(tidyverse)
setwd('/data/swamyvs/eyeintegration_splicing')
args = commandArgs(trailingOnly=TRUE)
sample_design <- read.table(args[1],stringsAsFactors = F, header = F, sep = '\t')
colnames(sample_design) <- c('sample_accession', 'run_accession', 'paired','tissue','subtissue','origin')
sample_design <-sample_design %>%
    mutate(path=sapply(sample_design$sample_accession,function(x)paste0('STARbams_realigned/',x,'/Aligned.out.bam')),stringsAsFactors = F)

#ignore single ended files for now
sample_design <- sample_design[!is.na(sample_design$path),]
ttypes=unique(sample_design$subtissue)
for(type in ttypes){
  df <- filter(sample_design, subtissue==type ,paired=='y')
  writeLines(df$path,paste0('ref/rmats_locs/',gsub(' ', '.' ,type),'.rmats.txt'),sep = ',')
 
}

st_full <- read_tsv('sampleTableFull.tsv', 
                    col_names = c('sample_accession', 'run_accession', 'paired','tissue','subtissue','origin')) %>%
    filter(subtissue=="RPE_Stem.Cell.Line")
