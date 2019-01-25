library(tidyverse)
set.seed(123321)
samples <- read_tsv('sampleTableFull.tsv', col_names = c('sample','run','paired','tissue','subtissue','origin'))
new_samples <- 'HM7FMBBXX_16496358_S54
HM7FMBBXX_16496358_S54
HM7FMBBXX_16496360_S56
HM7FMBBXX_16496360_S56
HM7FMBBXX_16496356_S52
HM7FMBBXX_16496356_S52
HN3CMBBXX_17042835_S49
HN3CMBBXX_17042835_S49
HN3CMBBXX_17042837_S51
HN3CMBBXX_17042837_S51
HN3CMBBXX_17042836_S50
HN3CMBBXX_17042836_S50
HM7HWBBXX_16424741_S52
HM7HWBBXX_16424741_S52
HM7MVBBXX_16424741_S63
HM7MVBBXX_16424741_S63
HM7FMBBXX_16424741_S61
HM7FMBBXX_16424741_S61
HM7HWBBXX_16424750_S61
HM7HWBBXX_16424750_S61
HM7FMBBXX_16424750_S70
HM7FMBBXX_16424750_S70
HM7MVBBXX_16424750_S72
HM7MVBBXX_16424750_S72' %>% strsplit('\n') %>% unlist %>% unique
in_house_sample_df <- data.frame(sample=new_samples, run='.', paired='y',tissue='RPE',
                                 subtissue='RPE_Stem.Cell.Line', origin='Stem.Cell.Line', stringsAsFactors = F)
qcd_samples <- read_csv('ref/QCdsapmles1218.txt') %>% colnames(.)
samples <- filter(samples,sample%in%qcd_samples)
samples <- rbind(samples, in_house_sample_df, stringsAsFactors=F)
# relabel samples
eye_tissues <- c("Retina", "RPE",'Cornea')
relabel_samps <- function(df, tis){
    t_df <- filter(df, tissue==tis)
    adult_tissue <- filter(t_df, subtissue==paste0(tis,'_Adult.Tissue'))
    adult_tissue$tissue <- adult_tissue$subtissue <- paste0(tis,'_Adult.Tissue')
    fetal_like_tissue <- filter(t_df, !subtissue==paste0(tis,'_Adult.Tissue') & 
                                    !subtissue==paste0(tis,'_Cell.Line') )
    fetal_like_tissue$origin <- fetal_like_tissue$subtissue
    fetal_like_tissue$tissue <- fetal_like_tissue$subtissue <- paste0(tis,'_Fetal.Tissue')
    return(rbind(adult_tissue, fetal_like_tissue))
}
lens <- filter(samples, tissue=='Lens') %>% mutate(origin=tissue, tissue='Lens_Stem.Cell.Line',
                                                   subtissue='Lens_Stem.Cell.Line')


eye <- lapply(eye_tissues, function(x) relabel_samps(samples,x)) %>% do.call(rbind,.) %>% rbind(lens)
eye_tissues <- c(unique(eye$subtissue),'synth')
noteye <- filter(samples,!tissue%in%c('Retina','RPE','Cornea','Lens')) %>% 
    mutate(subtissue=gsub('\\.','',subtissue) %>% gsub('\\(|\\)','-',.)) %>% mutate(origin=tissue, tissue=subtissue)
set.seed(123321)
synth <- split(noteye, noteye$subtissue) %>% lapply(function(x) sample_n(x,2)) %>% do.call(rbind,. ) %>% 
    mutate(origin=tissue, tissue='synth', subtissue='synth')
noteye <- noteye %>% filter(!sample%in%synth$sample)
res <- rbind(eye, noteye, synth)

tissues <- pull(res, tissue) %>% unique
subtissues <- filter(res, paired=='y') %>% pull(tissue) %>% unique
write(tissues,'/Volumes/data/eyeintegration_splicing/ref/tissues.txt',sep = '\n')
write(subtissues,'/Volumes/data/eyeintegration_splicing/ref/subtissues.txt',sep = '\n')
write_tsv(res,'/Volumes/data/eyeintegration_splicing/sampleTableV4.tsv',col_names = F)

