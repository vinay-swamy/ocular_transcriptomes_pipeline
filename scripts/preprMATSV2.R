library(tidyverse)
args = commandArgs(trailingOnly=TRUE)
working_dir=args[1]
sample_file=args[2]
bam_dir=args[3]
setwd(working_dir)

sample_design <- read.table(sample_file, stringsAsFactors = F, header = F, sep = '\t')
colnames(sample_design) <- c('sample_accession', 'run_accession', 'paired','tissue','subtissue','origin')
sample_design <-sample_design %>%
    mutate(path=sapply(sample_design$sample_accession,function(x)paste0(bam_dir, x,'/Sorted.out.bam')),stringsAsFactors = F)

#ignore single ended files for now
sample_design <- sample_design[!is.na(sample_design$path),]
ttypes=unique(sample_design$subtissue)
for(type in ttypes){
  df <- filter(sample_design, subtissue==type ,paired=='y')
  writeLines(df$path,paste0('ref/rmats_locs/',gsub(' ', '.' ,type),'.rmats.txt'),sep = ',')

}
