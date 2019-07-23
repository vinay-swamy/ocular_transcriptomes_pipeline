library(tidyverse)
library(IRanges)
library(parallel)
args= commandArgs(trailingOnly = T)
# args=c('~/NIH/eyeintegration_splicing/',
#         'results_b38/all_tissues.combined.gtf',
#         'results_b38/transdecoder_results/combined_stringtie_tx.fa.transdecoder.gff3',
#         'results_b38/corrected_lengths.tsv',
#         'testout.gff3',
#         2)
working_dir= args[1]
gtf_file=args[2]
cds_file=args[3]
length_correction_file=args[4]
outfile=args[5]
cores=args[6]
setwd(working_dir)

correct_CDS_lengths <- function(cor_file, cds_tab){
    cor_tab <- read_tsv(cor_file)
    #this table contains transcriptids and the length of their adjusted tx(80% of time no change)
    colnames(cor_tab) <- c('transcript_id', 'cor_len')
    new_starts <- filter(cds_tab,type=='CDS') %>% select(seqid, start, end) %>% dplyr::rename(transcript_id=seqid) %>%
        left_join(cor_tab, by='transcript_id') %>%
        mutate(aalen=(end-start+1)/3 -1, aa_removed= aalen- cor_len, new_start=start+aa_removed*3) %>% pull(new_start)
        #mutate(new_end=(aalen+1)*3+new_start-1, diff=new_end-end) %>% select(start,new_start,end,new_end, diff)
        #first, calulate the length of the original protein sequence (aalen), then figure out how many aa's
        #were removed (aa removed), then conver this back to nt and add to start, and then fix the end because thats
        #fucked up too.

    new_starts[new_starts<0] <- NA
    cds_tab$start[cds_tab$type=='CDS'] <- new_starts
    cds_tab <- filter(cds_tab, is.na(start)) %>% pull(seqid) %>% {filter(cds_tab, !seqid%in% .)}
    return(cds_tab)

}

scale_tx_length <- function(gtf){
    # add releative lengths to exons ie coordinates in spliced tx
    gtf <- gtf %>% mutate(scaled_start=1,scaled_end=1)
    if(gtf$strand[1]=='+'){
        gtf[1,] <- gtf[1,]%>% mutate(scaled_end=end-start)
        offset <- 0
        if(nrow(gtf)==1) return(gtf %>% mutate(scaled_length= scaled_end-scaled_start))
        for(i in 2:nrow(gtf)){
            start=gtf$scaled_end[i-1] + 1
            gtf[i,'scaled_start'] <-start
            gtf[i,'scaled_end'] <- start+  gtf[i,'end']-gtf[i,'start']
        }
        return(gtf %>% mutate(scaled_length= scaled_end-scaled_start))
    }else{
        n=nrow(gtf)
        gtf[n,] <- gtf[n,]%>% mutate(scaled_end=end-start)
        offset <- 0
        if(nrow(gtf)==1) return(gtf %>% mutate(scaled_length= scaled_end-scaled_start))
        for(i in (n-1):1){
            start=gtf$scaled_end[i+1] + 1
            gtf[i,'scaled_start'] <-start
            gtf[i,'scaled_end'] <- start+  gtf[i,'end']-gtf[i,'start']
        }
        return(gtf %>% mutate(scaled_length= scaled_end-scaled_start))
    }
}

get_cds_genomic_coords <- function(t_tx,cds_loc){
    tx_rang <- IRanges(start = t_tx$scaled_start,end = t_tx$scaled_end)
    cds_rang<- IRanges(start = cds_loc, end = cds_loc)
    hit <- IRanges::findOverlaps(query = cds_rang,subject = tx_rang, maxgap = 0)
    idx <-  subjectHits(hit)
    return(idx[1])

}

#fix +/- early

##ITS THE STRAND!!!!!!!
merge_gtf_cds <- function(tx,gtf,td_cds,gene){
    # The main challenge here is that the cds file has coordinates relative to each individual gene
    # so we have to map the features from cds to gtf. We do this by first c=finding relative coordinates for each tx
    # in the gtf, then overlapping cds coordinates, and then create the cds, 5putr, 3putr types and merge to form
    # a complete entry in gff3 format.
    # need to subtract 1 from gtf start to allow overlap with a single
    t_tx <- filter(gtf, transcript_id==tx) %>% mutate(ID=NA,parent=NA, start=start, length=end-start )
    tx_line <- t_tx[1,] %>% mutate(type='mRNA',scaled_start=NA,scaled_end=NA, scaled_length=NA,
                                   ID=tx, parent = gene)
    t_tx <- t_tx[-1,] %>% mutate(ID=paste(tx,'exon',exon_number,sep = '_'), parent=tx, length=end-start)
    td_cds_tx <- filter(td_cds,seqid==tx,type=='gene') %>% mutate(length=end-start)
    if(nrow(td_cds_tx)==0) return(data.frame())
    if(t_tx$strand[1]=='+'){
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
        cds_df[n,'end'] <- TES -1

        #now that we have made the right cds entry, we make the rest
        fputr <- t_tx[1:start_idx,] %>% mutate(type='five_prime_UTR',ID=paste(tx,'5putr',exon_number,sep = '_'), parent=tx)
        k=nrow(fputr)
        # cds coordinaes cannot overlap with utr coordinates.
        fputr[k,'end'] <- TSS-1
        tputr <- t_tx[end_idx:nrow(t_tx),] %>% mutate(type ='three_prime_UTR',ID=paste(tx,'3putr',exon_number,sep = '_'), parent=tx)
        tputr[1,'start'] <-TES
        res <-rbind(tx_line,t_tx,fputr,cds_df,tputr) %>%
            mutate(length=end-start, gene_id= .[1,'gene_id'], gene_name= .[1,'gene_name'], oId=.[1,'oId'],
                   cmp_ref= .[1,'cmp_ref']) %>%
            dplyr::select(-c('class_code','tss_id','contained_in','cmp_ref_gene','scaled_start','scaled_end','scaled_length'))
        res <-res %>%  filter(!start>=end)

        # b=Sys.time()
        # print(b-a)
        return(res)
    }else{
        # the negative strand been ruining my day since 2016
        # since we are using the negative strand, the genomic coodinatses are
        cds <- filter(td_cds,seqid==tx,type=='CDS')
        cds_start <- cds$start[1]
        cds_end <- cds$end[nrow(cds)]
        t_tx <- scale_tx_length(t_tx)
        start_idx <- get_cds_genomic_coords(t_tx,cds_start)
        end_idx <- get_cds_genomic_coords(t_tx,cds_end)
        cds_df <- t_tx[end_idx:start_idx, ] %>% mutate(type ='CDS', ID=paste(tx,'CDS',exon_number,sep = '_'),
                                                       parent=tx, phase = cds$phase[1] )
        n = nrow(cds_df)
        TSS <- cds_df[n,'end'] - (cds_start -cds_df[n,'scaled_start'])
        TES <- cds_df[1,'end'] - (cds_end- cds_df[1,'scaled_start']) +1 #vonk
        # DONT CHANGE VARIABLES UNTIL THEY ARE DONE BEING USED!!!!!!!
        cds_df[n,'end'] <- TSS
        cds_df[1,'start'] <- TES
        cds_df <-  mutate(cds_df, length=end-start)

        #now that we have made the right cds entry, we make the rest
        fputr <- t_tx[start_idx:nrow(t_tx),] %>% mutate(type='five_prime_UTR',ID=paste(tx,'5putr',exon_number,sep = '_'), parent=tx)

        fputr[1,'start'] <- TSS+1
        tputr <- t_tx[1:end_idx,] %>% mutate(type ='three_prime_UTR',ID=paste(tx,'3putr',exon_number,sep = '_'), parent=tx)
        k=nrow(tputr)
        tputr[k,'end'] <-TES-1
        res <-rbind(tx_line,t_tx,fputr,cds_df,tputr) %>%
            mutate(length=end-start, gene_id= .[1,'gene_id'], gene_name= .[1,'gene_name'], oId=.[1,'oId'],
                   cmp_ref= .[1,'cmp_ref']) %>%
            dplyr::select(-c('class_code','tss_id','contained_in','cmp_ref_gene','scaled_start','scaled_end','scaled_length'))
        res <-res %>%  filter(start<end)

        # b=Sys.time()
        # print(b-a)
        return(res)



    }
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

evaluate_gff_accuracy <- function(ref_gff, new_gff){
    types <- c('exon','CDS')

    check_accuracy <- function(ref,new,t){
        ref <- filter(ref,type==t) %>% select(start,end)
        new <- filter(new,type==t) %>% select(start,end)
        k <- all.equal(ref, new, tolerance=0)
        if(!isTRUE(k))return( k %>% strsplit(':') %>% sapply(function(x) x[[3]]))
        return(k)

    }

    eval <- lapply(types , function(x) check_accuracy(ref_gff,new_gff,x))
    names(eval) <- types
    eval

}


gtf <- rtracklayer::readGFF(gtf_file)
td_cds <- rtracklayer::readGFF(cds_file) %>% as.data.frame()
td_cds <- correct_CDS_lengths(length_correction_file, td_cds)
tx_list <- pull(td_cds,seqid) %>% unique %>% as.character()
tx2g <- filter(gtf, type=='transcript', transcript_id %in% tx_list)%>% select(gene_id, transcript_id)
gene_list <- pull(tx2g,gene_id) %>% unique
print(length(tx_list))
fin <- mclapply(gene_list,function(x) wrapper(x = x,gtf = gtf ,td_cds = td_cds, tx2g = tx2g), mc.cores = cores) %>%
    do.call(rbind,.) %>% dplyr::select(-c('transcript_id','gene_id','gene_name', 'oId','cmp_ref','exon_number'),
                                       c('transcript_id', 'gene_id','gene_name', 'oId','cmp_ref','exon_number'))

fin=fin %>% mutate(strand=gsub('*','+', strand, fixed=T))
save.image(file='mergeCDS_Rimg.Rdata')
write_gff3(df = fin,file = outfile)





#**************************************
# diff <- fin$end-fin$start
# set.seed(420024)
# ref_gff <- rtracklayer::readGFF('~/NIH/eyeintegration_splicing/enst.gff') %>% as.data.frame
# ref_gff <- filter(ref_gff, transcript_type=='protein_coding')
# tx_to_test<- filter(gtf, grepl('ENST',oId), strand=='+') %>% sample_n(25) %>% filter(transcript_id%in%td_cds$seqid,) %>%
#     select(transcript_id,oId) %>% {split( .,1:nrow(.))}
#
# i=1
# dfl=list()
# for( txs in tx_to_test){
#     res <- merge_gtf_cds(tx =txs[1,1] , gtf = gtf,gene = 'SMOOBY',td_cds = td_cds)
#     ref <- filter(ref_gff, transcript_id==txs[1,2], !type%in%c('start_codon','stop_codon'))
#     dfl <-c(dfl, evaluate_gff_accuracy(ref,res))
#
# }
#
#
#
#
# corner_cases <- fin[diff<=0,] %>% select(transcript_id, oId) %>% distinct
# corner_cases_ref <- filter(corner_cases, grepl('ENST', oId))
# targ_tx <- filter(ref_gff, transcript_id=='ENST00000335137.4', !type%in%c('stop_codon','start_codon'))
# tx='TCONS_00000065'
# gene='AGRN'
# reftx='ENST00000379370.6'
# ref_tx_gf <- filter(ref_gff, transcript_id==reftx, !type%in%c('start_codon','stop_codon')) %>% mutate(length=end-start)
# gene='MMMM'
# k <- td_cds %>% filter(type=='CDS') %>% mutate(len=end-start+1, frame=len%%3) %>% pull(frame)
# k <-
# new_lengths <- read_tsv('/Volumes/data/eyeintegration_splicing/results/corrected_lengths.tsv')
# colnames(new_lengths) <- c('transcript_id', 'cor_len')
#
#
#
# table(tab$aa_removed) %>%sort(decreasing = T)
# oneoff <- filter(tab, aa_removed==1) %>% pull(transcript_id) %>% {filter(gtf, transcript_id%in%., type=='transcript')}
# twooff <- filter(tab, aa_removed==0) %>% pull(transcript_id) %>% {filter(gtf, transcript_id%in%., type=='transcript')}
# fixed <- read_tsv('fixed_tx.txt', col_names = F) %>% dplyr::rename(transcript_id=X1)%>%pull(transcript_id) %>% { filter(gtf_tx, transcript_id%in% . ) }
# filter(ref_gff, transcript_id%in%fixed$oId) %>% View
# filter(ref_gff, transcript_id%in%fixed$oId)
# bdpc <- pull(twooff, oId) %>% {filter( ref_gff, transcript_id%in% .  ,type=='transcript', transcript_type=='protein_coding') }
tx2g<- rbind(filter(gtf, type=='transcript', transcript_id %in% tx_list, strand=='+') %>% select(gene_id, transcript_id) %>% sample_n(10),
      filter(gtf, type=='transcript', transcript_id %in% tx_list, strand=='-') %>% select(gene_id, transcript_id)  %>% sample_n(10))
#*************************************
