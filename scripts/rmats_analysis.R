#setwd("/gpfs/gsfs9/users/swamyvs/autoRNAseq/rmats_final/Retina_Adult.Tissue")
setwd('rmats_final')  
library(dplyr)
library(rtracklayer)
library('SNPlocs.Hsapiens.dbSNP149.GRCh38',lib.loc = '~/R/3.5/library/')
library(biomaRt)
library(UpSetR)
library(ggplot2)
library(gridExtra)
# skipped_exon <- read.table("SE.MATS.JC.txt",sep = '\t', header = T,stringsAsFactors = F)
# counts <- c("IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2")
# skipped_exon[is.na(skipped_exon)] <- 1
# skipped_exon[,counts] <- apply(skipped_exon[,counts],2,function(y) strsplit(y,',')%>%sapply(function(x) as.numeric(x)%>%sum))
# total_event_counts <- rowSums(skipped_exon[,counts])
# # mean <- mean(total_event_counts)
# # dev <- sd(total_event_counts)
# # hist(total_event_counts)
# # max(total_event_counts)
# # keep events with an average of 50 counts per group, still a lot left might increase it 
# keep <- which(total_event_counts>=100)#trimmed out about half
# skipped_exon.filtered <- skipped_exon[keep,]
# retnet <- read.csv('../../ref/RetNet_genes.csv',stringsAsFactors = F,header = F)$V1
# SE_retnet <- skipped_exon[skipped_exon$geneSymbol%in%retnet,]
# unique(SE_retnet$geneSymbol)%>%length
# 
# mart <- useMart('ENSEMBL_MART_SNP','hsapiens_snp')
# snps.chr15 <- getSNPlocs('15')
##PART 1 looking at gtf
source('plot_testing V5.broken.R')
tissues <- c('Retina','RPE','Cornea','ESC')
gtf <- load('ref/gtf_unformatted.Rdata')
gtf.tx <-dplyr:: filter(gtf_final,type=='transcript')
sample_table <- read.table('samplesrun_1010.tissue',header = F, stringsAsFactors = F,sep = '\t')
tissue_count <- table(sample_table$V4)%>%.[c('Retina','RPE','Cornea','ESC')]
stringTie_results <- data.frame(sample_number=tissue_count,stringsAsFactors = F)
novel <- grepl('MSTR', gtf.tx$transcript_id)%>%filter(gtf.tx,.)
novel_count <- nrow(novel)
ref <- nrow(gtf.tx)-novel_count
gtfs_per_tissues <- lapply(tissues, function(x) paste0('ref/', x,'_st.gtf')%>%readGFF%>%filter(type=='transcript'))
for (x in tissues){
    gtf.tx[!is.na(gtf.tx[,x]),x] <- gtf.tx[!is.na(gtf.tx[,x]),'transcript_id']

}
novel_tx_by_tissue_gtf <- lapply(gtfs_per_tissues,function(x) is.na(x[,'reference_id'])%>%sum)%>%unlist
names(novel_tx_by_tissue) <- tissues
novel_by_tissues <- lapply(tissues, function(x) (grepl('MSTR',gtf.tx[,x]))%>%sum)                       
stringTie_results$novel_transcipts <- novel_tx_by_tissue_gtf
tx_by_tissue <- lapply(tissues, function(x) gtf.tx[,x]%>%na.omit)
names(tx_by_tissue) <- tissues
novel_tx_by_tissue <- lapply(tissues, function(x) novel[,x]%>%na.omit)
names(novel_tx_by_tissue) <- tissues
k <- lapply(1:nrow(gtf.tx), function(x) is.na(gtf.tx[x,tissues]))%>%lapply(all)%>%unlist
missinall <- gtf.tx[k,]%>%filter( is.na(ref_gene_id))
upset(fromList(tx_by_tissue),order.by = 'freq') #all tx
upset(fromList(novel_tx_by_tissue),order.by = 'freq')



sample_table <- read.table('samplesrun_1010.tissue',stringsAsFactors = F, header = F,sep = '\t')
colnames(sample_table) <- c('sample','run','paired','tissue','subtissue','origin')
samples_eye <- filter(sample_table,tissue%in%c('Retina','RPE','Cornea','ESC'))
body_synth <- read.table('ref/synth_body.txt',stringsAsFactors = F, header = F)
sample_table <- rbind(samples_eye,filter(sample_table, sample%in%body_synth$V1))
sample_table[!sample_table$tissue%in%c('Retina','RPE','Cornea','ESC'),c('tissue','subtissue')] <- 'body'
files <- paste0('quant_files/',sample_table$sample,'/quant.sf')
load('ref/gtf_unformatted_new.Rdata')
txdb <- filter(gtf_final,type=='transcript')%>%select(transcript_id,gene_name)
txdb$comb_name <- paste0(txdb$gene_name,'-',txdb$transcript_id)
files <- paste0('quant_files/',sample_table$sample,'/quant.sf')
txi <- tximport::tximport(files = files,type = 'salmon',tx2gene = txdb,countsFromAbundance = 'lengthScaledTPM',txOut=T)
colnames(txi$counts) <- sample_table$sample
counts <- txi$counts
txi.counts_by_tissue <- lapply(tissues, function(x) filter(sample_table,tissue==x)[,'sample']%>%txi$counts[,.])%>%lapply(as.data.frame)
counts_novel_by_tissues <-lapply(txi.counts_by_tissue,function(x) x[rowSums(x) > 50*ncol(x),])%>%lapply( function(x) x[grepl('MSTR',rownames(x)),])
tx_novel_by_tis <- lapply(counts_novel_by_tissues,rownames)
names(tx_novel_by_tis) <- tissues
i <- lapply(tissues, function(x) setdiff(tx_novel_by_tis[[x]],tx_novel_by_tis[!tissues%in%x])%>%unlist%>%length)
#have to manuallly save upsets, so annoying
tx_novel_by_tis_count <- lapply(tx_novel_by_tis, length)
stringTie_results$novel_expressed_transcripts <- tx_novel_by_tis_count
load('tables_for_poster.Rdata')
mytheme <- gridExtra::ttheme_default(
    core = list(fg_params=list(cex = 3.0)),
    colhead = list(fg_params=list(cex = 3.0)),
    rowhead = list(fg_params=list(cex = 3.0)))

grid.table(ret_signoveltx,theme=mytheme)
grid.table(rpe_signovelret,theme=mytheme)


upset(fromList(tx_novel_by_tis[-4]),order.by = 'freq',text.scale = 2.5)
ggsave('upset_major_tissues.tiff',device = 'tiff',width = 9, height = 6, units = 'in',res = 300)
subtissue <- unique(sample_table$subtissue)[1:10]
txi.counts_by_subtissue <- lapply(subtissue, function(x) filter(sample_table,subtissue==x)[,'sample']%>%txi$counts[,.])%>%lapply(as.data.frame)
counts_by_subtissue <-lapply(txi.counts_by_subtissue,function(x) x[rowSums(x) > 50*ncol(x),])%>%lapply( function(x) x[grepl('MSTR',rownames(x)),])
novel_tx_by_subtissue <- lapply(counts_by_subtissue,rownames)
names(novel_tx_by_subtissue) <- subtissue
upset(fromList(novel_tx_by_subtissue[c('Retina_Adult.Tissue','Retina_Stem.Cell.Line')]),order.by = 'freq',text.scale = 2.5)
ggsave('upset_Retina.tiff',device = 'tiff',width = 9, height = 6, units = 'in',dpi = 300)
upset(fromList(novel_tx_by_subtissue[c('RPE_Adult.Tissue','RPE_Fetal.Tissue','RPE_Cell.Line','RPE_Stem.Cell.Line')]),order.by = 'freq',text.scale = 2.5)
ggsave('upset_RPE.tiff',device = 'tiff',width = 9, height = 6, units = 'in',dpi = 300)
upset(fromList(novel_tx_by_subtissue[c('Cornea_Adult.Tissue','Cornea_Fetal.Tissue','Cornea_Cell.Line')]),order.by = 'freq',text.scale = 2.5)
ggsave('upset_Cornea.tiff',device = 'tiff',width = 9, height = 6, units = 'in',dpi = 300)
#look at novel tx's
novel_tx <- filter(txdb,grepl('MSTR',transcript_id))
tx.counts.novel <- txi.counts[rownames(txi.counts)%in%novel_tx$transcript_id,]%>%as.data.frame()
keep <- which(rowSums(tx.counts.novel)>50*ncol(tx.counts.novel)) #~avg 50 counts per sample, still keeping a lot of stuff 
tx.counts.novel <- tx.counts.novel[keep,]
colnames(tx.counts.novel) <-colnames(txi.counts) <-  sample_table$sample
counts_by_tissue <- lapply(unique(sample_table$subtissue), function(x) filter(sample_table,subtissue==x)[,'sample']%>%tx.counts.novel[,.]%>%rowMeans)%>%do.call(cbind,.)%>%as.data.frame
colnames(counts_by_tissue) <- unique(sample_table$subtissue)  
new_rows <- filter(txdb, transcript_id%in%rownames(counts_by_tissue))%>%select(comb_name)
rownames(counts_by_tissue) <- new_rows$comb_name
n_genes <- rownames(counts_by_tissue)%>%lapply(function(x) strsplit(x,'-')[[1]][1])%>%unlist%>%unique
rpe_genes[rpe_genes%in%n_genes]# see below for RPE genes
RPE_at <- read.table("rmats_analysis/RPE_Adult.Tissue_PE_VS_body_PE/SE.MATS.JC.txt" ,header = T, sep = '\t',stringsAsFactors = F)
splice_graph(all_events = RPE_at,gtf = gtf_final,gene='RLBP1',sub_tissue = 'RPE.Adult.Tissue')
#didnt find anything good
ret_at <- read.table('rmats_analysis/Retina_Adult.Tissue_VS_body/SE.MATS.JC.txt',header = T, sep = '\t',stringsAsFactors = F)
splice_graph(gene="TEAD1",all_events = ret_at,gtf = gtf_final,sub_tissue = 'Retina_adult')


#mak- got something but its a little big
#trnt - maybe a retained intron
#RAX2 - maybe
#AIPL1
#DHX1
all_events <- ret_at
events.filtered <- filter(all_events,new_fdr<=.05)
sig_ret <- filter(events.filtered, geneSymbol%in%ret_genes)
#all.equal(x%>%as.numeric,rmats_exon_coords[j,c('start','end')]%>%as.numeric,tolerance = n,check.names=F

ret_at <- read.table('rmats_analysis/Retina_Adult.Tissue_VS_body/SE.MATS.JC.txt',header = T, sep = '\t',stringsAsFactors = F)
RPE_at <- read.table('rmats_analysis/RPE_Adult.Tissue_PE_VS_body_PE/SE.MATS.JC.txt', header = T, sep = '\t',stringsAsFactors = F)
cornea_at <- read.table('rmats_analysis/Cornea_Adult.Tissue_SE_VS_body_SE/SE.MATS.JC.txt', header = T, sep = '\t',stringsAsFactors = F)
sig_rpe_retnet <- filter(RPE_at,new_fdr<=.05)%>%filter(geneSymbol%in%ret_genes)
sig_rpre_rpe <- filter(RPE_at,new_fdr<=.05)%>%filter(geneSymbol%in%rpe_genes)# none
sig_cornea_retnet <-  filter(cornea_at,new_fdr<=.05)%>%filter(geneSymbol%in%ret_genes)#none

findSigNovelTx <- function(gtf_final,all_events){
    final=list()
    list_dex=1
    for(i in unique(all_events$geneSymbol)){
        gtf.gene <- filter(gtf_final,gene_name==i,type=='exon')
        gtf.gene[,c('start','end')] <- gtf.gene[,c('start','end')]-1 
        events <- filter(all_events,geneSymbol==i)
        exon_list <- split(gtf.gene[,c('start','end')],seq(nrow(gtf.gene)))%>%lapply(as.numeric)
        for(event in 1:nrow(events) ){
            k <- events[event,]
            p <- lapply(exon_list,function(x) all.equal(x,k[,5:6]%>%as.numeric,check.names=F,tolerance=1.0e-6)%>%isTRUE)%>%unlist
            l <- gtf.gene[p,]%>%filter(grepl('MSTR',transcript_id))
            if(nrow(l)>0){
                res <- list(gtf.gene[p,],k)
                m<- data.frame(gene=i,tx =as.character(l$transcript_id)%>%paste(collapse = ','),qval=k[,'new_fdr'],k[,c("SJC_SAMPLE_1_avg","IJC_SAMPLE_1_avg","SJC_SAMPLE_2_avg","IJC_SAMPLE_2_avg")])
                final[[list_dex]] <- m
                list_dex <- list_dex+1
                }   
        }

    }
    final <- do.call(rbind,final)
    return(final)
}

ret_signoveltx <- final
rpe_signovelret <- findSigNovelTx(gtf_final = gtf_final ,all_events = sig_rpe_retnet)

save(ret_signoveltx,rpe_signovelret,file = 'tables_for_poster.Rdata')

################################################################################
## part2 looking at rmats out_put

all_SE_novel_events <- read.table('rmats_analysis/all.SE.novelevents.txt',stringsAsFactors = F, header = T,sep = '\t')%>%
    .[!duplicated(all_SE_novel_events[,c("chr","strand","exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE" )]),]
all_SE_novel_events$novel <- TRUE
retAdult_SE <- read.table('rmats_analysis/Retina_Adult.Tissue_VS_body/SE.MATS.JC.txt',header=T,stringsAsFactors = F,sep = '\t')
countsCol<-c('IJC_SAMPLE_1','SJC_SAMPLE_1','IJC_SAMPLE_2','SJC_SAMPLE_2')
k=list()
for(i in countsCol) k <- c(k,paste(i,c('counts','avg','med','sd'),sep = '_'))
k <- unlist(k)
colnames(retAdult_SE) <- c(colnames(retAdult_SE)[1:10],k,'pval','qval')
events.filtered <- filter(retAdult_SE,qval<0.05)
events_with_novel <- left_join(retAdult_SE,all_SE_novel_events[,4:12],by=c("chr","strand","exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE" ))
events_with_novel
k <- duplicated(all_SE_novel_events[,c("chr","strand","exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE" )])

####
sample_table <- read.table('samplesrun_1010.tissue',stringsAsFactors = F, header = F,sep = '\t')
colnames(sample_table) <- c('sample','run','paired','tissue','subtissue','origin')
samples_eye <- filter(sample_table,tissue%in%c('Retina','RPE','Cornea','ESC'))
body_synth <- read.table('ref/synth_body.txt',stringsAsFactors = F, header = F)
sample_table <- rbind(samples_eye,filter(sample_table, sample%in%body_synth$V1))
sample_table[!sample_table$tissue%in%c('Retina','RPE','Cornea','ESC'),c('tissue','subtissue')] <- 'body'
files <- paste0('quant_files/',sample_table$sample,'/quant.sf')
load('ref/gtf_unformatted_new.Rdata')
txdb <- filter(gtf_final,type=='transcript')%>%select(transcript_id,gene_name)
txdb$comb_name <- paste0(txdb$gene_name,'-',txdb$transcript_id)
txi <- tximport::tximport(files = files,type = 'salmon',tx2gene = txdb,countsFromAbundance = 'lengthScaledTPM',txOut=T)
txi.counts <- txi$counts
keep <- which(rowSums(txi.counts)>1000)
txi.counts <- txi.counts[keep,]
colnames(txi.counts) <- sample_table$sample
txi.counts_by_tissue <- lapply(unique(sample_table$subtissue), function(x) filter(sample_table,subtissue==x)[,'sample']%>%txi.counts[,.]%>%rowMeans)%>%do.call(cbind,.)%>%as.data.frame
colnames(txi.counts_by_tissue) <- unique(sample_table$subtissue)
new_rows <- filter(txdb, transcript_id%in%rownames(txi.counts_by_tissue))%>%select(comb_name)
rownames(txi.counts_by_tissue) <- new_rows$comb_name
#look at RPE genes
rpe_genes <-  c('POU5F1','KLF4','SOX2','MITF','TFEC','PAX6','OTX2','SOX9','SIX3','ROR','LHX2','CRX','GAS1','GPNMB','MYRIP','RAB27A','OCA2','DCT','TYROSINASE','TYRP1','PMEL',
                'ALDH1A3','RPE65','RDH5','RLBP1','CDH1','CDH3','BEST1','COL11A1','TRPM1','CLDN16')
rpe_gene_counts <- lapply(rpe_genes,function(x)txi.counts_by_tissue[grepl(x,rownames(txi.counts_by_tissue)),])%>%do.call(rbind,.)
files <- list.files(path='rmats_analysis/', pattern = 'SE.MATS.JC.txt',recursive = T)%>%paste0('rmats_analysis/',.)
RPE_at <- read.table(files[10],header = T, sep = '\t',stringsAsFactors = F)
RPE_cL <- read.table(files[8],header = T, sep = '\t',stringsAsFactors = F)
RPE_scl <- read.table(files[6],header = T, sep = '\t',stringsAsFactors = F)

g='CRX'
splice_graph(all_events = RPE_at,gtf = gtf_final,gene=g,sub_tissue = 'RPE.Adult.Tissue')
splice_graph(all_events = RPE_scl,gtf = gtf_final,gene=g,sub_tissue = 'RPE.Stem.Cell.Line')
splice_graph(all_events = RPE_cL,gtf = gtf_final,gene = g, sub_tissue = 'RPE.Cell.Line')

#look at retina genes
ret_genes <- read.csv('ref/RetNet_genes.csv',header = F, stringsAsFactors = FALSE)$V1
ret_gene_counts <- lapply(ret_genes$V1,function(x)txi.counts_by_tissue[grepl(x,rownames(txi.counts_by_tissue)),])%>%do.call(rbind,.)
ret_at <- read.table('rmats_analysis/Retina_Adult.Tissue_VS_body/SE.MATS.JC.txt',header = T, sep = '\t',stringsAsFactors = F)
splice_graph(all_events = ret_at,gtf = gtf_final,gene = 'ZC3H14',sub_tissue = 'RetinaAdultTissue')



files <- list.files(path='rmats_analysis/', pattern = 'SE.MATS.JC.txt',recursive = T)%>%paste0('rmats_analysis/',.)
for(i in files){
    save(se,gtf_final,file='prlbm.Rdata')
    se <- read.table(i,header = T, sep = '\t',stringsAsFactors = F)
    print(i)
    try(splice_graph(gene='BEST1', gtf= gtf_final,all_events = se),outFile = stdout())
    }

splice_graph(all_events = read.table(files[1],header = T, sep = '\t',stringsAsFactors = F),gtf = gtf_final,gene='BEST1')

