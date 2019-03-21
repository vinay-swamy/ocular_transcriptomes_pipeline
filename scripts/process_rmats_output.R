library(tidyverse)

# args=c('~/NIH/eyeintegration_splicing/',
#        'rmats_out/RPE_Adult.Tissue/RI.MATS.JC.txt',
#        'RI.MATS.JC.txt',
#        'sampleTableV4.tsv',
#        'RPE_Adult_tissue',
#        'rmats_clean/RPE_Adult.Tissue/psi.RI.MATS.JC.txt',
#        'rmats_clean/RPE_Adult.Tissue/incCts.RI.MATS.JC.txt')
args= commandArgs( trailingOnly = T)
working_dir=args[1]
event_file=args[2]
event= args[3]
sample_file <- args[4]
tissue <- args[5]
outfile_incCts <- args[6]
outfile_psi <- args[7]
setwd(working_dir)


process_rmats_output <- function(file,event,sample_file,outfile_incCts, outfile_psi){
    event_header <- list(SE.MATS.JC.txt=c('chr'	,'strand',	'exonStart_0base',	'exonEnd',	'upstreamES',	'upstreamEE',	'downstreamES',	'downstreamEE'),
                         RI.MATS.JC.txt=c('chr'	,'strand',	'riExonStart_0base',	'riExonEnd'	,'upstreamES'	,'upstreamEE'	,'downstreamES'	,'downstreamEE'),
                        MXE.MATS.JC.txt=c('chr',	'strand',	'1stExonStart_0base',	'1stExonEnd',	'2ndExonStart_0base',	'2ndExonEnd'	,'upstreamES',	'upstreamEE',	'downstreamES',	'downstreamEE'),
                       A5SS.MATS.JC.txt=c('chr',	'strand',	'longExonStart_0base',	'longExonEnd',	'shortES',	'shortEE',	'flankingES',	'flankingEE'),
                       A3SS.MATS.JC.txt=c('chr',	'strand',	'longExonStart_0base',	'longExonEnd'	,'shortES',	'shortEE'	,'flankingES',	'flankingEE')
    )

    dup_header  <- list(SE.MATS.JC.txt=c('chr'	,'strand',	'exonStart_0base',	'exonEnd'),
                        RI.MATS.JC.txt=c('chr'	,'strand',	'riExonStart_0base',	'riExonEnd'	),
                       MXE.MATS.JC.txt=c('chr',	'strand',	'X1stExonStart_0base',	'X1stExonEnd'),
                      A5SS.MATS.JC.txt=c('chr',	'strand',	'longExonStart_0base',	'longExonEnd'),
                      A3SS.MATS.JC.txt=c('chr',	'strand',	'longExonStart_0base',	'longExonEnd')
    )


    parse_count_info <- function(df,col){
        
        len <- ifelse(grepl('IJC', col),'IncFormLen','SkipFormLen') %>% {pull(df, .)}
        t_cols <- pull(df, col) %>% lapply(function(x) strsplit(x,',')%>%unlist%>%as.numeric) #%>% unlist # %>%as.numeric
        flattened_counts <- lapply(t_cols, function(x){x[is.na(x)] <- 0; x})
        t_median <- sapply(flattened_counts,median) / len
        final <- data.frame(t_median)
        colnames(final) <- paste('med', col , sep = '_')
        return(final)
    }



    target_event<- read_tsv(file = event_file)
    sample_table <- read_tsv(file = sample_file, col_names = c('sample','run','paired','tissue','subtissue','origin'))
    if( nrow(target_event)==0){
        print('empty event file')
        k=c('seqid','strand', 'start', 'end', 'Z')
        writeLines(k, outfile_psi, sep='\t')
        writeLines(k, outfile_incCts, sep='\t')
        return(0)
    }
    if(event == 'MXE.MATS.JC.txt'){

        # split into similar columns; the ijc sample1 is inclusion cts for 1st, sjc sample1 is for sample 2
        exon_1 <- select(target_event, chr, strand,`1stExonStart_0base`, `1stExonEnd`, IJC_SAMPLE_1, IncFormLen) %>% 
            {mutate(.,med_IJC_SAMPLE_1=parse_count_info(.,'IJC_SAMPLE_1')[,1] / IncFormLen, 
                        IJC_SAMPLE_1=NULL,IncFormLen=NULL, med_SJC_SAMPLE_1=0)}
        #yeah this looks real fuckin weird but its right - going to say that each mxe exon is used 100% of the time
        # and since the SJC_sample 1 counts are for the 2nd exon, 
        exon_2 <- select(target_event, chr, strand,`2ndExonStart_0base`, `2ndExonEnd`, SJC_SAMPLE_1, SkipFormLen) %>% 
            {mutate(.,med_IJC_SAMPLE_1=parse_count_info(.,'SJC_SAMPLE_1')[,1] / SkipFormLen, 
                        SJC_SAMPLE_1=NULL,SkipFormLen=NULL, med_SJC_SAMPLE_1=0 )}
        colnames(exon_1) <- colnames(exon_2) <- c('seqid', 'strand','start','end','med_IJC_SAMPLE_1', 'med_SJC_SAMPLE_1')
        
        exon_cts <- rbind(exon_1, exon_2, stringsAsFactors=F) %>% group_by(seqid,strand, start, end) %>%  
            summarise(incCts=sum(med_IJC_SAMPLE_1)) %>% mutate(PSI=1)
        incLvl <- select(exon_cts, -incCts)
        colnames(incLvl) <- c('seqid','strand', 'start', 'end', paste(tissue, 'PSI',sep = '_'))
        medCts <- select(exon_cts, -PSI)
        colnames(medCts) <- c('seqid','strand', 'start', 'end', paste(tissue, 'incCts',sep = '_'))
        write_tsv(incLvl, outfile_psi)
        write_tsv(medCts, outfile_incCts)
        return(0)

    }

    countsCol<-c('IJC_SAMPLE_1','SJC_SAMPLE_1')


    exon_cts <- lapply(countsCol,function(x)parse_count_info(target_event,x)) %>%do.call(cbind,.) %>% 
        {cbind(target_event[,event_header[[event]]], .) }%>%  
        group_by_(.dots = dup_header[[event]]) %>% 
        summarise(incCts=sum(med_IJC_SAMPLE_1),exclCts=sum(med_SJC_SAMPLE_1)) %>% 
        mutate(PSI=incCts/(incCts + exclCts)) %>% mutate(PSI=replace(PSI, is.nan(PSI),0))

    incLvl <- select(exon_cts, dup_header[[event]], PSI)
    colnames(incLvl) <- c('seqid','strand', 'start', 'end', paste(tissue, 'PSI',sep = '_'))
    medCts <- select(exon_cts, dup_header[[event]], incCts)
    colnames(medCts) <- c('seqid','strand', 'start', 'end', paste(tissue, 'incCts',sep = '_'))
    
    write_tsv(incLvl, outfile_psi)
    write_tsv(medCts, outfile_incCts)

}

process_rmats_output(file = file,
                     event = event,
                     sample_file = sample_file,
                     outfile_incCts = outfile_incCts,
                     outfile_psi = outfile_psi)
