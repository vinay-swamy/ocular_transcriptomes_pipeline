library(tidyverse)
library(pheatmap)


gfc_gtf_file <- 'results/all_tissues.combined.gtf'
stringtie_gtf_file <- 'results/all_tissues_st.gtf'
incLvlTab_file <- 'results/rmats_all_tissue_incLevels.tsv'
medCounts_files <-  'results/rmats_all_tissue_medCounts.tsv'
event <- 'SE.MATS.JC.txt'


analyze_ouput <- function(gfc_gtf_file,stringtie_gtf_file,incLvlTab_file,medCounts_files, event) {
    event_header <- list(SE.MATS.JC.txt=c('GeneID', 'seqid'	,'strand',	'start', 'end',	'upstreamES',	'upstreamEE',	'downstreamES',	'downstreamEE'),
                         RI.MATS.JC.txt=c('GeneID', 'seqid'	,'strand',	'start', 'end'	,'upstreamES'	,'upstreamEE'	,'downstreamES'	,'downstreamEE'),
                         MXE.MATS.JC.txt=c('GeneID', 'seqid',	'strand',	'1stExonStart_0base',	'1stExonEnd',	'2ndExonStart_0base',	'2ndExonEnd'	,'upstreamES',	'upstreamEE',	'downstreamES',	'downstreamEE'),
                         A5SS.MATS.JC.txt=c('GeneID', 'seqid',	'strand',	'start', 'end',	'shortES',	'shortEE',	'flankingES',	'flankingEE'),
                         A3SS.MATS.JC.txt=c('GeneID','seqid',	'strand',	'start', 'end'	,'shortES',	'shortEE'	,'flankingES',	'flankingEE')
    )
    
    writeLines('\n')
    print(event)
    incTab <- suppressMessages( read_tsv(incLvlTab_file) %>% select( -synth_incLevel))
    os <- length(event_header[[event]])
    colnames(incTab) <- c(event_header[[event]], colnames(incTab)[(os+1) : ncol(incTab)])
    
    
    medCounts <- suppressMessages( read_tsv(medCounts_files)%>% select( -synth_incLevel))
    colnames(medCounts) <- c(event_header[[event]], colnames(medCounts)[(os+1) : ncol(medCounts)])
    
    incTab <- data.frame(ljid=paste0('rm_',1:nrow(incTab)), incTab, stringsAsFactors = F)
    
    
    
    # lets pick exons that have a median count 50 in at least 1 tissue
    nc_cols <- event_header[[event]]
    highly_expressed_exons <- rowSums(medCounts%>%select(-nc_cols) > 50) >0
    incTab_filt <- incTab[highly_expressed_exons,]
    ctsTab_filt <- medCounts[highly_expressed_exons,]
    print('number of highly expressed events')
    print(nrow(incTab_filt))
    
    itf <- select(incTab_filt, -ljid, -nc_cols)
    lvl <-.95
    high_inc_one_tissue <- rowSums(itf > lvl) == 1
    low_inc_other_tissue <- rowSums(itf < (1-lvl)) == (ncol(itf) - 1)  
    tissue_specfic <- high_inc_one_tissue & low_inc_other_tissue
    inctab_tissue_specifc <- incTab_filt[tissue_specfic,]
    ctstab_tissue_specifc <- ctsTab_filt[tissue_specfic,]
    print('number of tissue specfic events')
    print(nrow(inctab_tissue_specifc))
    save(inctab_tissue_specifc,ctstab_tissue_specifc, file = paste0('results/',strsplit(event,'\\.')[[1]][1],'_tissue_specfic.Rdata'))
    writeLines('\n')
    
    eye_tissues <- paste0(c('Retina_Adult.Tissue','RPE_Adult.Tissue','Cornea_Adult.Tissue','Lens_Stem.Cell.Line'),
                          '_incLevel')
    eye_specfic <- lapply(eye_tissues, function(x) filter(inctab_tissue_specifc, inctab_tissue_specifc[,x] >= lvl ))
    names(eye_specfic) <- c('Retina','RPE','Cornea','Lens')
    lapply(eye_specfic, nrow)
    
    
    
    
    
    gfc_gtf <- rtracklayer::readGFF(gfc_gtf_file) %>% mutate(start=start-1, seqid=as.character(seqid))
    gfc_gtf_rmats <- left_join(gfc_gtf, select(incTab, ljid, seqid, strand, start,end),
                              by=c('seqid','strand','start','end'))
    filter(gfc_gtf_rmats, !is.na(ljid)) %>% pull(ljid) %>% unique %>%length
    
    
    summarize_tissue <- function(tissue, df_list, gfc_gtf_rmats){
        print(tissue)
        df <- df_list[[tissue]]
        print('number of specfic exons')
        print(nrow(df))
        ret_ids <- df %>% pull(ljid)
        gfc_tis <- filter(gfc_gtf_rmats, ljid%in%ret_ids) %>% pull(transcript_id) %>% 
        {filter(gfc_gtf_rmats, transcript_id %in% . , type=='transcript') }
        print('number of transcripts with specific exons ')
        print(nrow(gfc_tis))
        print('number of novel transcripts with specfic exons')
        print(sum(grepl('MSTRG',gfc_tis$oId)))
        print('number of known trasncritpts with specifc exons')
        print(sum(!grepl('MSTRG',gfc_tis$oId)))
        writeLines('\n')
    }
    
    eye_tissues <- c('Retina','RPE','Cornea','Lens')
    k <- lapply(eye_tissues, function(x) summarize_tissue(x,eye_specfic,gfc_gtf_rmats))
}
analyze_ouput(gfc_gtf_file = 'results/all_tissues.combined.gtf',
              stringtie_gtf_file = 'results/all_tissues_st.gtf',
              incLvlTab_file = 'results_bw/all_tissues.A5SS.incLevel.tsv',
              medCounts_files =  'results_bw/all_tissues.A5SS.medCounts.tsv',
              event = 'A5SS.MATS.JC.txt')
analyze_ouput(gfc_gtf_file = 'results/all_tissues.combined.gtf',
             stringtie_gtf_file = 'results/all_tissues_st.gtf',
             incLvlTab_file = 'results_bw/all_tissues.A3SS.incLevel.tsv',
             medCounts_files =  'results_bw/all_tissues.A3SS.medCounts.tsv',
             event = 'A3SS.MATS.JC.txt')
analyze_ouput(gfc_gtf_file = 'results/all_tissues.combined.gtf',
              stringtie_gtf_file = 'results/all_tissues_st.gtf',
              incLvlTab_file = 'results_bw/all_tissues.RI.incLevel.tsv',
              medCounts_files =  'results_bw/all_tissues.RI.medCounts.tsv',
              event = 'RI.MATS.JC.txt')
