setwd('~/NIH/eyeintegration_splicing/')
library(tidyverse)
library(IRanges)
library(parallel)
args= commandArgs()
gtf_file=args[1]
cds_file=args[2]
cores=args[4]
outfile=args[3]
gtf <- rtracklayer::readGFF(gtf_file)
td_cds <- rtracklayer::readGFF(cds_file) %>% as.data.frame()


scale_tx_length <- function(gtf){

    gtf <- gtf %>% mutate(scaled_start=0,scaled_end=0)
    offset <- 0
    for(i in 1:nrow(gtf)){
        start=gtf$start[i]
        gtf[i,'scaled_start'] <- gtf[i,'start']-start + offset
        gtf[i,'scaled_end'] <- gtf[i,'end']-start + offset
        offset <- gtf[i,'scaled_end']
    }
    return(gtf %>% mutate(scaled_length= scaled_end-scaled_start))    
}

get_cds_genomic_coords <- function(t_tx,cds_loc){
    tx_rang <- IRanges(start = t_tx$scaled_start,end = t_tx$scaled_end)
    cds_rang<- IRanges(start = cds_loc, end = cds_loc)
    hit <- IRanges::findOverlaps(query = cds_rang,subject = tx_rang)
    idx <-  subjectHits(hit)
    return(idx)
    
}

merge_gtf_cds <- function(tx,gtf,td_cds){
    a=Sys.time()
    t_tx <- filter(gtf, transcript_id==tx) %>% mutate(ID=NA,parent=NA, start=start-1,length=end-start )
    tx_line <- t_tx[1,] %>% mutate(type='mRNA',scaled_start=NA,scaled_end=NA, scaled_length=NA, 
                                   ID=tx, parent = NA)
    t_tx <- t_tx[-1,] %>% mutate(ID=paste(tx,'exon',exon_number,sep = '_'), parent=tx)
    td_cds_tx <- filter(td_cds,seqid==tx,type=='gene') %>% mutate(start=start-1,length=end-start) 
    td_cds_tx$length
    if(sum(t_tx$length) != td_cds_tx$length[1]) return('spliced tx does not match cds ')
    cds <- filter(td_cds,seqid==tx,type=='CDS')
    cds_start <- cds$start[1]
    cds_end <- cds$end[1]
    t_tx <- scale_tx_length(t_tx)
    start_idx <- get_cds_genomic_coords(t_tx,cds_start)
    end_idx <- get_cds_genomic_coords(t_tx,cds_end)
    cds_df <- t_tx[start_idx:end_idx, ] %>% mutate(type ='CDS', ID=paste(tx,'CDS',exon_number,sep = '_'),
                                                   parent=tx, phase = cds$phase[1] )
    TSS <- cds_df[1,'start'] + (cds_start -cds_df[1,'scaled_start'])
    cds_df[1,'start'] <- TSS
    n = nrow(cds_df)
    TES <- cds_df[n,'start'] + (cds_end- cds_df[n,'scaled_start'])
    cds_df[n,'end'] <- TES
    fputr <- t_tx[1:start_idx,] %>% mutate(type='five_prime_UTR',ID=paste(tx,'5putr',exon_number,sep = '_'),
                                            parent=tx)
    k=nrow(fputr)
    fputr[k,'end'] <- TSS
    tputr <- t_tx[end_idx:nrow(t_tx),] %>% mutate(type ='three_prime_UTR',ID=paste(tx,'3putr',exon_number,sep = '_'), 
                                                  parent=tx)
    tputr[1,'start'] <-TES
    res <-rbind(tx_line,t_tx,fputr,cds_df,tputr) %>%
        mutate(length=end-start, gene_id= .[1,'gene_id'], gene_name= .[1,'gene_name'], oId=.[1,'oId'], 
               cmp_ref= .[1,'cmp_ref']) %>%
        dplyr::select(-c('class_code','tss_id','contained_in','cmp_ref_gene','scaled_start','scaled_end','scaled_length'))
    b=Sys.time()
    print(b-a)
    return(res)
}

write_gff3 <- function(df, file){
    to_write <- lapply(colnames(df[9:ncol(df)]), function(x) pull(df,x) %>% paste(x,., sep = '=')) %>% 
        do.call(cbind,.) %>% split(.,1:nrow(df)) %>%lapply( function(x) grep('=NA',x,invert = T) %>% x[.]) %>%
        sapply( function(x) paste(x,collapse = ';')) 
    res <- cbind(df[,1:8],to_write,stringsAsFactors=F)
    res[is.na(res)] <- '.'
    write('##gff-version 3',file = file)
    write_tsv(res,file, append = T)
}

tx_list <- pull(td_cds,seqid) %>% unique
fin <- mclapply(tx_list[1:10],function(x) merge_gtf_cds(tx = x,gtf = gtf ,td_cds = td_cds), mc.cores = cores) %>% 
    do.call(rbind,.) %>% dplyr::select(-c('transcript_id','gene_id','gene_name','oId','cmp_ref','exon_number'),
                          c('transcript_id', 'gene_id','gene_name','oId','cmp_ref','exon_number'))

write_gff3(df = fin,file = outfile)


