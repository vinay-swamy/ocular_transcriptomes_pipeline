#setwd('~/NIH/eyeintegration_splicing/')
library(tidyverse)
setwd('/data/swamyvs/eyeintegration_splicing')
k <- '
an exon can be involved in 2 seperate events across different tissues 
ie cn be binary in tissue a and in b  but are2 didffernt events; right now will remvoe and examine later
make a wide table with all counts

'
#args <-c('SE.MATS.JC.txt','test_inc.tsv','test_med.tsv')
args <- commandArgs()
event=args[1]
outfile_inclvl=args[2]
outfile_medcts=args[3]

combine_allTissues <- function(event, outfile_inclvl,outfile_medcts){
    event_header <- list(SE.MATS.JC.txt=c('GeneID', 'chr'	,'strand',	'exonStart_0base',	'exonEnd',	'upstreamES',	'upstreamEE',	'downstreamES',	'downstreamEE'),
                         RI.MATS.JC.txt=c('GeneID', 'chr'	,'strand',	'riExonStart_0base',	'riExonEnd'	,'upstreamES'	,'upstreamEE'	,'downstreamES'	,'downstreamEE'),
                         MXE.MATS.JC.txt=c('GeneID', 'chr',	'strand',	'1stExonStart_0base',	'1stExonEnd',	'2ndExonStart_0base',	'2ndExonEnd'	,'upstreamES',	'upstreamEE',	'downstreamES',	'downstreamEE'),
                         A5SS.MATS.JC.txt=c('GeneID', 'chr',	'strand',	'longExonStart_0base',	'longExonEnd',	'shortES',	'shortEE',	'flankingES',	'flankingEE'),
                         A3SS.MATS.JC.txt=c('GeneID', 'chr',	'strand',	'longExonStart_0base',	'longExonEnd'	,'shortES',	'shortEE'	,'flankingES',	'flankingEE')
    )
    event='SE.MATS.JC.txt'
    
    
    #sample_table <- read_tsv('sampleTableV3.tsv', col_names = c('sample','run','paired','tissue','subtissue','origin'))
    files0 <- list.files('rmats_clean', pattern = paste0('bin.', event), recursive = T, full.names = T)
    tissues <- strsplit(files0,'/') %>% sapply(function(x) x[[2]])
    event_df_list <- lapply(files0,function(x)suppressMessages(read_tsv(x)))
    names(event_df_list) <- tissues
    keep <-  lapply(event_df_list,nrow) > 0
    event_df_list <- event_df_list[keep]
    all_event_df <- do.call(rbind,event_df_list)
    synthcol <- select(all_event_df, 
                        which(colnames(all_event_df)%in%event_header[[event]] | grepl('2',colnames(all_event_df))))
                        #keep cols with location info, and synth count info.
    keep_1 <- !duplicated(synthcol)# remove duplicates from cbind
    tmp <- synthcol[keep_1,]
    dup_2 <- duplicated(tmp[,1:4])# exon can be part of 2 separate events in different tissues, so remove those
    dup_list <- tmp[dup_2,1:4] %>% split(.,1:nrow(.))
    tmp_list <- split(tmp[,1:4], 1:nrow(tmp))
    keep2 <- !tmp_list%in%dup_list
    synthcol_distinct <- tmp[keep2,]
    
    make_full_table <- function(df_l, col, exon_df){
        t_cols <- c(event_header[[event]],col)
        exon_df_full <- select(exon_df,c(event_header[[event]],
                                         gsub('1','2',col)))
        exon_df <- select(exon_df,event_header[[event]])# get all events 
        df_l <- lapply(df_l, function(x) select(x, t_cols))
        full_tab<- reduce(.x = df_l,.f = left_join,by=event_header[[event]], .init = exon_df) %>%
            left_join(exon_df_full,by=event_header[[event]])
    
        colnames(full_tab) <-  c(event_header[[event]],paste0(names(df_l), '_incLevel'),'synth_incLevel')
        full_tab[is.na(full_tab)] <- 0
        return(full_tab)
    }
    
    full_tab_incLevel <- make_full_table(df_l = event_df_list,col ="IncLevel1_avg",exon_df = synthcol_distinct )
    full_tab_med_counts <- make_full_table(df_l = event_df_list,col ="IJC_SAMPLE_1_med",exon_df = synthcol_distinct )
    # make_full_wide_table <- function(){
    # files1 <- list.files('rmats_clean', pattern = paste0('wide.', event), recursive = T, full.names = T)
    # wide_df_list <- lapply(files1,function(x)suppressMessages(read_tsv(x)))
    # 
    # synthcol_wide <- lapply(wide_df_list, function(x) select(x, which(colnames(x)%in%event_header[[event]] | grepl('_2',colnames(x))))) %>%
    #     do.call(rbind,.)
    # 
    # }
    write_tsv(full_tab_incLevel,  outfile_inclvl)
    write_tsv(full_tab_med_counts, outfile_medcts)
}

combine_allTissues(event = event,
                   outfile_inclvl = outfile_inclvl, 
                   outfile_medcts = outfile_medcts)





