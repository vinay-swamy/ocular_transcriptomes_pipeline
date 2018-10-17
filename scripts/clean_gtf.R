setwd('/data/swamyvs/eyeintegration_splicing')
library(rtracklayer)
library(dplyr)
library(biomaRt)
gtf <- readGFF('ref/comb.gtf')
dups <- filter(gtf, type=='transcript')%>%.[duplicated(.[,1:7]),'transcript_id']
problem_tx <-filter(gtf,transcript_id%in%dups)
gtf <- filter(gtf,!transcript_id%in%dups)
###this script is bad and should feel bad
# for the sake of simplicity, lets remove duplicate transcripts from these each of the GTFs, then join together by transcript
# exons don't wont have tissue info, but really who cares about them anyway
read_GTF <- function(name){
    gff<- readGFF(name)%>%dplyr::select(c("seqid","type","start","end","strand","transcript_id"))
    dups <- filter(gff, type=='transcript')%>%.[duplicated(.[,1:5]),'transcript_id']
    gff <- filter(gff,!transcript_id%in%dups)
    return(gff)
}

gtf.tx <- filter(gtf,type=='transcript')
ret <- read_GTF('ref/Retina_st.gtf')
rpe <- read_GTF('ref/RPE_st.gtf')
cor <- read_GTF('ref/Cornea_st.gtf')
esc <- read_GTF('ref/ESC_st.gtf')
pod <- read_GTF('ref/body_st.gtf')
col <- colnames(gtf)
gtf.tx_tissue_specific<- left_join(filter(gtf,type=='transcript'),filter(ret,type=='transcript'),by=c("seqid",'type',"start","end","strand"))%>%
    left_join(filter(rpe,type=='transcript'),by=c("seqid",'type',"start","end","strand"))%>%
    left_join(filter(cor,type=='transcript'),by=c("seqid",'type',"start","end","strand"))%>%
    left_join(filter(esc ,type=='transcript'),by=c("seqid",'type',"start","end","strand"))#%>%
    # left_join(filter(pod,type=='transcript'),by=c("seqid",'type',"start","end","strand"))
# add uncooperative tx's back in
problem_tx <- data.frame(problem_tx,matrix(data = NA,nrow = nrow(problem_tx),ncol = 4))
colnames(gtf.tx_tissue_specific) <- colnames(problem_tx) <- c(col,'Retina','RPE','Cornea','ESC')
gtf <- left_join(gtf,gtf.tx_tissue_specific,by=col)

colnames(gtf) <- colnames(problem_tx) <- c(col,'Retina','RPE','Cornea','ESC')
gtf <- rbind(gtf,problem_tx,stringsAsFactors=F)
gdata::keep(gtf,sure=T)
#remove single exon transcripts
tx_counts <- table(gtf$transcript_id)
keep <- names(tx_counts[tx_counts>2])
gtf <-gtf[gtf$transcript_id%in%keep,]
gtf$seqid <- as.character(gtf$seqid)
gtf$type <- as.character(gtf$type)
#gtf <-gtf[grep('chr',gtf$seqid),]
#fix some formatting
ref_gtf <- readGFF('ref/gencodeAno_bsc.gtf')%>%filter(gene_type=="protein_coding" )
ref_gtf$seqid <- as.character(ref_gtf$seqid)
gtf.tmp <- gtf
ref_gtf.tmp <- ref_gtf
ref_gtf$seqid[ref_gtf$seqid=='chrX'] <- gtf$seqid[gtf$seqid=='chrX'] <- 'chr23'
ref_gtf$seqid[ref_gtf$seqid=='chrY'] <- gtf$seqid[gtf$seqid=='chrY'] <- 'chr24'
ref_gtf$seqid[ref_gtf$seqid=='chrM'] <- gtf$seqid[gtf$seqid=='chrM'] <- 'chr25'
ref_bed <- data.frame(seqid=strsplit(ref_gtf$seqid,'r')%>%sapply(function(x)x[2]),ref_gtf$start,ref_gtf$end,ref_gtf$gene_name,
                      scor='.',strand=ref_gtf$strand, stringsAsFactors = F)
write.table(ref_bed,'refgtf.bed',col.names = F,row.names = F,quote = F,sep = '\t',na='.')
gtf.missing_genename <- gtf[is.na(gtf$ref_gene_id),]%>%filter(type=='transcript')
bed <- data.frame(seqid=strsplit(gtf.missing_genename$seqid,'r')%>%sapply(function(x)x[2]),gtf.missing_genename$start,
                  gtf.missing_genename$end,gtf.missing_genename$transcript_id,score='.',strand=gtf.missing_genename$strand,
                  stringsAsFactors = F)
write.table(bed,'missing.bed',col.names = F,row.names = F,quote = F,sep ='\t',na='.')
colnames(bed) <- c('V1','V2','V3','transcript_id')
bed$V1 <- as.integer(bed$V1)
#read in bedtools result
gtf <- gtf.tmp
ref <- ref_gtf.tmp
save.image("tmp.RData")
#system2('bedtools intersect -f .5  -wo -s -a missing.bed -b refgtf.bed > testout.bed')




