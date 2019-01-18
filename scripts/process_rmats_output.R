
#setwd('~/NIH/eyeintegration_splicing/')
setwd('/data/swamyvs/eyeintegration_splicing')
library(tidyverse)
#args=c('rmats_out/RPE_Adult.Tissue/A3SS.MATS.JC.txt','A3SS.MATS.JC.txt', 'test_raw.tsv ','test_bin.tsv',
#       'test_multi.tsv')
args= commandArgs( trailingOnly = T)
file= args[1]
event= args[2]
outfile_raw <-args[3]
outfile_bin <- args[4]
outfile_multi <- args[5]
process_rmats_output <- function(file,event,outfile_raw,outfile_bin,outfile_multi){
    event_header <- list(SE.MATS.JC.txt=c('chr'	,'strand',	'exonStart_0base',	'exonEnd',	'upstreamES',	'upstreamEE',	'downstreamES',	'downstreamEE'),
                           RI.MATS.JC.txt=c('chr'	,'strand',	'riExonStart_0base',	'riExonEnd'	,'upstreamES'	,'upstreamEE'	,'downstreamES'	,'downstreamEE'),
                           MXE.MATS.JC.txt=c('chr',	'strand',	'1stExonStart_0base',	'1stExonEnd',	'2ndExonStart_0base',	'2ndExonEnd'	,'upstreamES',	'upstreamEE',	'downstreamES',	'downstreamEE'),
                           A5SS.MATS.JC.txt=c('chr',	'strand',	'longExonStart_0base',	'longExonEnd',	'shortES',	'shortEE',	'flankingES',	'flankingEE'),
                           A3SS.MATS.JC.txt=c('chr',	'strand',	'longExonStart_0base',	'longExonEnd'	,'shortES',	'shortEE'	,'flankingES',	'flankingEE')
    )

    dup_header  <- list(SE.MATS.JC.txt=c('chr'	,'strand',	'exonStart_0base',	'exonEnd'),
                      RI.MATS.JC.txt=c('chr'	,'strand',	'riExonStart_0base',	'riExonEnd'	),
                      MXE.MATS.JC.txt=c('chr',	'strand',	'X1stExonStart_0base',	'X1stExonEnd',	'X2ndExonStart_0base',	'X2ndExonEnd'),
                      A5SS.MATS.JC.txt=c('chr',	'strand',	'longExonStart_0base',	'longExonEnd'),
                      A3SS.MATS.JC.txt=c('chr',	'strand',	'longExonStart_0base',	'longExonEnd')
    )


    parse_count_info <- function(df,col){
        t_cols <- pull(df, col) %>% lapply(function(x) strsplit(x,',')%>%unlist%>%as.numeric) #%>% unlist # %>%as.numeric
        flattened_counts <- lapply(t_cols, function(x){x[is.na(x)] <- 0; x})
        t_mean <- sapply(flattened_counts, mean)
        t_median <- sapply(flattened_counts,median)
        t_sd <-  sapply(flattened_counts,sd)
        final <- data.frame(t_mean,t_median,t_sd)
        colnames(final) <- paste(col,c('avg','med','sd'),sep = '_')
        return(final)
    }
    target_event <- read_tsv(file = file)
    if( nrow(target_event)==0){
      print('empty event file')
      k=colnames(target_event)
      writeLines(k, outfile_raw, sep='\t')
      writeLines(k, outfile_bin, sep='\t')
      writeLines(k, outfile_multi, sep='\t')
      return(0)
    }

    countsCol<-c('IJC_SAMPLE_1','SJC_SAMPLE_1','IJC_SAMPLE_2','SJC_SAMPLE_2',"IncLevel1","IncLevel2")
    procd_cols <- lapply(countsCol,function(x)parse_count_info(target_event,x)) %>%do.call(cbind,.)
    pvals <- mutate(target_event, PValue=replace(PValue, PValue==0,2.2e-16), FDR=replace(FDR,FDR==0,2.2e-16)) %>%
        select(PValue, FDR)
    out_df <- data.frame(select(target_event,c("GeneID","geneSymbol",event_header[[event]])),
                         procd_cols,
                         pvals)
    write_tsv(out_df, outfile_raw)

    out_df <- filter(out_df, IJC_SAMPLE_1_med > 10) # sample 1 is our target sample, keep ex-juncs with more than 10 counts
    k <- select(out_df, dup_header[[event]])
    dups <- k[duplicated(k),]
    dup_list <- split(dups,1:nrow(dups))
    exon_list <- split(k, 1:nrow(k))
    bool <- exon_list%in%dup_list
    uniq_exons <- out_df[!bool,]
    duped_exons <- out_df[bool,]
    write_tsv(uniq_exons,outfile_bin)
    write_tsv(duped_exons,outfile_multi)
}

process_rmats_output(file = file,
                     event = event,
                     outfile_raw = outfile_raw,
                     outfile_bin = outfile_bin,
                     outfile_multi = outfile_multi)
