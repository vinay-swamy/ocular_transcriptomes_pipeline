setwd('/data/swamyvs/eyeintegration_splicing')
library(tidyverse)
library(IRanges)
library(parallel)
args= commandArgs(trailingOnly = T)
#args=c('results/all_tissues.combined.gtf',
#       'results/combined_stringtie_tx.fa.transdecoder.gff3',
#       'testout.gff3',
#       1)
gtf_file=args[1]
cds_file=args[2]
outfile=args[3]
cores=args[4]


print(gtf_file)
print(cds_file)
gtf <- rtracklayer::readGFF(gtf_file)

td_cds <- rtracklayer::readGFF(cds_file) %>% as.data.frame()


scale_tx_length <- function(gtf){
    # add releative lengths to exons ie coordinates in spliced tx
    gtf <- gtf %>% mutate(scaled_start=1,scaled_end=1)
    gtf[1,] <- gtf[1,]%>% mutate(scaled_end=end-start) 
    offset <- 0
    if(nrow(gtf)==1) return(gtf %>% mutate(scaled_length= scaled_end-scaled_start))
    for(i in 2:nrow(gtf)){
        start=gtf$scaled_end[i-1] + 1
        gtf[i,'scaled_start'] <-start
        gtf[i,'scaled_end'] <- start+  gtf[i,'end']-gtf[i,'start']
    }
    return(gtf %>% mutate(scaled_length= scaled_end-scaled_start))
}

get_cds_genomic_coords <- function(t_tx,cds_loc){
    tx_rang <- IRanges(start = t_tx$scaled_start,end = t_tx$scaled_end)
    cds_rang<- IRanges(start = cds_loc, end = cds_loc)
    hit <- IRanges::findOverlaps(query = cds_rang,subject = tx_rang)
    idx <-  subjectHits(hit)
    return(idx[1])

}

merge_gtf_cds <- function(tx,gtf,td_cds,gene){
    # The main challenge here is that the cds file has coordinates relative to each individual gene
    # so we have to map the features from cds to gtf. We do this by first c=finding relative coordinates for each tx 
    # in the gtf, then overlapping cds coordinates, and then create the cds, 5putr, 3putr types and merge to form 
    # a complete entry in gff3 format.
    
    t_tx <- filter(gtf, transcript_id==tx) %>% mutate(ID=NA,parent=NA,length=end-start )
    tx_line <- t_tx[1,] %>% mutate(type='mRNA',scaled_start=NA,scaled_end=NA, scaled_length=NA,
                                   ID=tx, parent = gene)
    t_tx <- t_tx[-1,] %>% mutate(ID=paste(tx,'exon',exon_number,sep = '_'), parent=tx)
    td_cds_tx <- filter(td_cds,seqid==tx,type=='gene') %>% mutate(length=end-start)


    if(nrow(td_cds_tx)==0) return(data.frame())
    cds <- filter(td_cds,seqid==tx,type=='CDS')
    cds_start <- cds$start[1]
    cds_end <- cds$end[nrow(cds)]
    t_tx <- scale_tx_length(t_tx)
    start_idx <- get_cds_genomic_coords(t_tx,cds_start)
    end_idx <- get_cds_genomic_coords(t_tx,cds_end)
    cds_df <- t_tx[start_idx:end_idx, ] %>% mutate(type ='CDS', ID=paste(tx,'CDS',exon_number,sep = '_'),
                                                   parent=tx, phase = cds$phase[1] )
    TSS <- cds_df[1,'start'] + (cds_start -cds_df[1,'scaled_start'])
    n = nrow(cds_df)
    TES <- cds_df[n,'start'] + (cds_end- cds_df[n,'scaled_start'])
    # DONT CHANGE VARIABLES UNTIL THEY ARE DONE BEING USED!!!!!!!
    cds_df[1,'start'] <- TSS
    cds_df[n,'end'] <- TES

    #now that we have made the right cds entry, we make the rest
    fuckWithPhase <- function(cdf){
      cdf <- mutate(cdf,length=end-start +1)
      sphase <- cdf$length[1]%%3
      if(nrow(cdf)==1) return(cdf)
      for( i in 2:nrow(cdf)){
        sphase=cdf[i,'phase']=sphase
        sphase <- sphase+ cdf[i,'length']
        sphase <- sphase%%3
      }
      cdf
    }
    cds_df <- fuckWithPhase(cds_df)
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
    # b=Sys.time()
    # print(b-a)
    return(res)
}

addGeneToTx <- function(gene,gtf,td_cds,tx2g){
    gene_gtf=filter(gtf, gene_id==gene) %>% arrange(start)
    outdf=gene_gtf[1,] %>% mutate(seqid=as.character(seqid),type='gene', source=as.character(source) , start=min(gene_gtf$start),end=max(gene_gtf$end), transcript_id=NA, ID=gene_id,
                                  parent=NA, tss_id=NULL,class_code=NULL, contained_in=NULL, cmp_ref_gene=NULL, length=end-start)
    t_list <-  filter(tx2g, gene_id==gene)%>%pull(transcript_id) %>% unique
    out_list=lapply(t_list,function(x)merge_gtf_cds(tx = x,gtf = gtf ,td_cds = td_cds, gene = gene)) %>% do.call(rbind,.) %>%
      mutate(seqid=as.character(seqid),type=as.character(type), source=as.character(source))
    outdf=rbind(outdf,out_list,stringsAsFactors=F)
    return(outdf)
}
wrapper <- function(x,gtf,td_cds,tx2g){
    out = tryCatch(
        {
            addGeneToTx(gene = x,gtf = gtf ,td_cds = td_cds, tx2g = tx2g)
        },
        error=function(cond)
        {
            message(paste('error',x))
            return(NA)
        },
        warning=function(cond)
        {
            message(paste('warning',x))
            return(NULL)
        },
        finally={}
    )
}


write_gff3 <- function(df, file){
    #formatparent and format id are to mimic ensembl gff3
    formatParent <- function(type, parent) {
      ol=character( length = length(type))
      ol[type!='mRNA'] <- 'transcript:'
      ol[type=='mRNA'] <- 'gene:'
      ol[type=='gene'] <- ''
      return(paste0(ol,parent))
      }
    formatID <- function(type,ID){
      ol=character( length = length(type))
      ol[type=='mRNA'] <- 'transcript:'
      ol[type=='gene'] <- 'gene:'
      return(paste0(ol,ID))
    }

    df <- df %>%mutate(parent=formatParent(type,parent)) %>% dplyr::rename(Parent=parent) %>%
      mutate(ID=formatID(type,ID), biotype=ifelse(type=='mRNA' | type=='gene' , 'protein_coding' , NA))

    to_write <- lapply(colnames(df[9:ncol(df)]), function(x) pull(df,x) %>% paste(x,., sep = '=')) %>%
        do.call(cbind,.) %>% split(.,1:nrow(df)) %>%lapply( function(x) grep('=NA',x,invert = T) %>% x[.]) %>%
        sapply( function(x) paste(x,collapse = ';'))
    res <- cbind(df[,1:8],to_write,stringsAsFactors=F)
    res[is.na(res)] <- '.'
    write('##gff-version 3',file = file)
    write_tsv(res,file, append = T)
}
tx_list <- pull(td_cds,seqid) %>% unique %>% as.character()
tx2g <- filter(gtf, type=='transcript', transcript_id %in% tx_list)%>% select(gene_id, transcript_id)
gene_list <- pull(tx2g,gene_id) %>% unique
print(length(tx_list))
fin <- mclapply(gene_list,function(x) wrapper(x = x,gtf = gtf ,td_cds = td_cds, tx2g = tx2g), mc.cores = cores) %>%
    do.call(rbind,.) %>% dplyr::select(-c('transcript_id','gene_id','gene_name', 'oId','cmp_ref','exon_number'),
                                       c('transcript_id', 'gene_id','gene_name', 'oId','cmp_ref','exon_number'))

fin=fin %>% mutate(strand=gsub('*','+', strand, fixed=T))
write_gff3(df = fin,file = outfile)



