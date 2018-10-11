library(rtracklayer)
library(dplyr)
library(biomaRt)
load("~/NIH/autoRNAseq/tmp.RData")
genenames <- read.table('testout.bed',stringsAsFactors = F,header = F)%>%arrange(desc(V13))%>%
    rename(transcript_id='V4',ext_gene_name='V10')
genenames <- genenames[!duplicated(genenames[,1:3]),]
genenames$type <- 'transcript'
gtf$type <- as.character(gtf$type)
out_gtf <- left_join(gtf.tmp,genenames,by=c('type','transcript_id'))

#add genenames back in
targ_tx <- genenames$transcript_id
targ_genes <- genenames$ext_gene_name
for(k in 1:length(targ_tx)){
    out_gtf[out_gtf$transcript_id==targ_tx[k],'gene_name'] <- targ_genes[k]
    
}
out_gtf[is.na(out_gtf$gene_name),'gene_name'] <- out_gtf[is.na(out_gtf$gene_name),'gene_id']
out_gtf <- out_gtf[,1:12]
out_gtf[is.na(out_gtf)] <- '.'


ref_genes <- filter(ref_gtf,type=='gene')
ref_genes <- ref_genes[ref_genes$gene_name%in%gtf$gene_name,  ]
ref_genes <- ref_genes[,c(1:9,14,11,20)]
cn <- c("seqid", "source", "type", "start", "end" ,"score", "strand", "phase" ,"gene_id", "transcript_id" ,"gene_name", "exon_number"    )
colnames(out_gtf) <- colnames(ref_genes) <- cn

gtf_final <- rbind(out_gtf,ref_genes)
gtf_final[is.na(gtf_final)] <- '.'
to_write <- paste( paste0('gene_id ',gtf_final$gene_id),
                   paste0(' transcript_id ','"',gtf_final$transcript_id,'"'),
                   paste0(' exon_numer ','"',gtf_final$exon_number,'"'),
                   paste0(' gene_name ','"',gtf_final$gene_name,'"'),sep = ';')
gtf_final$to_write <- to_write  
gtf_final <- arrange(gtf_final, seqid, start, desc(type))
gtf_final <- arrange(gtf_final, seqid,start, gene_id, transcript_id,type)
write.table(gtf_final[,c(1:8,13)],'ref/combined_final.gtf',sep = '\t',quote = F, col.names = F,row.names = F)
