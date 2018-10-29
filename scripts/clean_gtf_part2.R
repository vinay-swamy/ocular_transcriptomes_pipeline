library(rtracklayer)
library(dplyr)
library(biomaRt)
a <- Sys.time()
#load("/data/swamyvs/eyeintegration_splicing/tmp.RData")
load('tmp.RData')
# load and format bedtools output
genenames <- read.table('testout.bed',stringsAsFactors = F,header = F)%>%arrange(desc(V13))%>%
    dplyr::rename(transcript_id='V4',ext_gene_name='V10')
genenames <- genenames[!duplicated(genenames[,1:3]),]

genenames$type <- 'transcript'
gtf$type <- as.character(gtf$type)
out_gtf <- left_join(gtf.tmp,genenames,by=c('type','transcript_id'))

#add genenames back in
targ_tx <- genenames$transcript_id
targ_genes <- genenames$ext_gene_name
#exons are important again, could be a lot faster
for(k in 1:length(targ_tx)){
    out_gtf[out_gtf$transcript_id==targ_tx[k],'gene_name'] <- targ_genes[k]
    
}
#add remaining gene names,
out_gtf[is.na(out_gtf$gene_name),'gene_name'] <- out_gtf[is.na(out_gtf$gene_name),'gene_id']
out_gtf <- out_gtf[,1:17]

add_novel_exon_flag <- function(novel,gene){
    ref <- filter(ref_gtf,gene_name==gene, type=='exon')
    novel_ex <- filter(novel,type=='exon')
    #working for a specifc gene rn
    novel_list <- split(novel_ex[,c('start','end')],seq(nrow(novel_ex)))%>%lapply(as.numeric)
    ref_list <- split(ref[,c('start','end')],seq(nrow(ref)))%>%lapply(as.numeric)
    novel[novel$type=='exon','is_novel_exon'] <- !novel_list%in%ref_list
    return(novel)
}

out_gtf <- filter(out_gtf,gene_name%in%events.filtered$geneSymbol)
nov_knowngene <- filter(out_gtf, grepl('MSTR',transcript_id),!grepl('MSTR',gene_name))
out_gtf$is_novel_exon <- FALSE 
#SLOW 
out_gtf<- lapply(unique(out_gtf$gene_name),
           function(x) add_novel_exon_flag(filter(out_gtf,gene_name==x),gene=x))%>%do.call(rbind,.)

# nkg_samd1l_ex_tab<- filter(nov_knowngene,gene_name==i,type=='exon')
# nkg_samd1l_ex <- lapply(1:nrow(nkg_samd1l_ex), function(x) as.list(nkg_samd1l_ex[x,c('start','end')]))
# rgtf_samd11_ex_tab <- filter(ref_gtf,gene_name==i, type=='exon')
# rgtf_samd11_ex <- lapply(1:nrow(rgtf_samd11_ex),function(x) as.list(rgtf_samd11_ex[x,c('start','end')]))

#out_gtf[is.na(out_gtf)] <- '.'

# now add back the remaining gene features , not sure if this is important at all but might as well
 # ref_genes <- filter(ref_gtf,type=='gene')
 # ref_genes <- ref_genes[ref_genes$gene_name%in%gtf$gene_name,  ]
 # ref_genes <- ref_genes[,c(1:9,14,11,)]
 # cn <- c("seqid", "source", "type", "start", "end" ,"score", "strand", "phase" ,"gene_id", "transcript_id" ,"gene_name"    )
 # colnames(out_gtf) <- colnames(ref_genes) <- cn
#gtf_final <- rbind(out_gtf,ref_genes)
gtf_final <- out_gtf
# identify what transcripts are novel
# v  <- !grepl('ENST',gtf_final$transcript_id)
# gtf_final$added_by_stringtie <- ifelse(v,'y','n')
gtf_final$is_novel_exon <- ifelse(gtf_final$is_novel_exon,'y','n')

save(gtf_final,file = 'ref/gtf_unformatted.Rdata')
#gtf_final[is.na(gtf_final)] <- '.'
to_write <- apply(gtf_final[,9:18],1,as.list)%>%lapply(unlist)%>%
    lapply(na.omit)%>%lapply( function(y)lapply(names(y), function(x) paste0(x, ' "',y[x],'" '))%>%unlist%>%
                                  paste(collapse = ';'))%>%as.character()

#lapply(names(k), function(x) paste0(x, ' "',k[x],'" '))%>%unlist%>%paste(collapse = ' ')


# ls fto_write <- paste( paste0('gene_id ',gtf_final$gene_id),
#                    paste0(' transcript_id ','"',gtf_final$transcript_id,'"'),
#                    paste0(' exon_numer ','"',gtf_final$exon_number,'"'),
#                    paste0(' gene_name ','"',gtf_final$gene_name,'"'),
# 		   paste0(' added_by_stringtie ','"',gtf_final$added_by_stringtie,'"'),sep = ';')
gtf_final$to_write <- to_write
gtf_final <- arrange(gtf_final, seqid, start, desc(type))
gtf_final <- arrange(gtf_final, seqid,start, gene_id, transcript_id,type)
write.table(gtf_final[,c(1:8,19)],'ref/combined_final.gtf',sep = '\t',quote = F, col.names = F,row.names = F)
b <- Sys.time()
b-a