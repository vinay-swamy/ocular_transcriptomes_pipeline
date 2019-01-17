library(tidyverse)
set.seed(123321)
samples <- read_tsv('../autoRNAseq/sampleTable1218_tissues.tab', col_names = F)
qcd_samples <- read_csv('ref/QCdsapmles1218.txt') %>% colnames(.)
samples <- filter(samples, X1%in%qcd_samples)
eye <- filter(samples,X5%in%c("Retina_Adult.Tissue", "RPE_Adult.Tissue",'Cornea_Adult.Tissue','Lens_Stem.Cell.Line'))
noteye <- filter(samples,!X4%in%c('Retina','RPE','Cornea','Retinal.Endothelium','Lens','eyelid','ESC', 'Brain'))
brain <- filter(samples,X4=='Brain') %>% sample_n(65)
noteye <- mutate(noteye,X6=X5, X5=X4) 
brain <- mutate(brain,X6=X5, X5=X4)
res <- rbind(eye,noteye,brain, stringsAsFactors=F)
p_sample <- function(df){
    k=5
    if (nrow(df) < 10 ) return(data.frame()) 
    sample_n(df,k)
}

k<- split(res, res$X4) %>% lapply(function(x) p_sample(x)) %>% do.call(rbind,.)
k <- mutate(k, X6=X4, X5='synth', X4='synth')
newsamp <- filter(res, !X1%in%k$X1)
fin <- rbind(newsamp,k, stringsAsFactors=F)

write_tsv(fin,'/Volumes/data/eyeintegration_splicing/sampleTableV3.tsv',col_names = F)

