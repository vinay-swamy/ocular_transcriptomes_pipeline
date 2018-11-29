setwd('~/NIH/eyeintegration_splicing/')
library(tidyverse)
good_samples <- read_csv('~/NIH/autoRNAseq/smoothed_filtered_tpms.csv')%>%colnames
sample_table <- read.table('sampleTable1128_tissues.tab',header = F,stringsAsFactors = F,sep='\t')
sample_table_filtered <- filter(sample_table, V1%in%good_samples)
eye_tissues <- c('Retina','RPE','Cornea',"EyeLid",'Lens','ESC')
samples_by_body_tissue <- sample_table_filtered%>%
  filter(!V4%in%eye_subtissues)%>%split(.[,4])
set.seed(3526)
synth_body <- lapply(samples_by_body_tissue,function(x)sample(x[,1],7))%>%do.call(c,.)
samples_eye <- filter(sample_table_filtered, V4%in%c('Retina','RPE','Cornea'))
samples_body <- filter(sample_table_filtered,V1%in%synth_body)
full_samp_tab <- rbind(samples_eye,samples_body,stringsAsFactors=F)
write.table(full_samp_tab,'sampleTableV2_tissues.tab',row.names = F,col.names = F,sep = '\t',quote = F)
full_samp_tab[!full_samp_tab$V4%in%c('Retina','RPE','Cornea'),4:5] <- 'body'
write.table(full_samp_tab,'sampleTableV2.tab',row.names = F,col.names = F,sep = '\t',quote = F)



