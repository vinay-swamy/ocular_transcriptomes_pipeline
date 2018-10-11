#setwd("/gpfs/gsfs9/users/swamyvs/autoRNAseq/rmats_final/Retina_Adult.Tissue")
setwd('rmats_final')  
library(dplyr)
library('SNPlocs.Hsapiens.dbSNP149.GRCh38',lib.loc = '~/R/3.5/library/')
library(biomaRt)
skipped_exon <- read.table("SE.MATS.JC.txt",sep = '\t', header = T,stringsAsFactors = F)
counts <- c("IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2")
skipped_exon[is.na(skipped_exon)] <- 1
skipped_exon[,counts] <- apply(skipped_exon[,counts],2,function(y) strsplit(y,',')%>%sapply(function(x) as.numeric(x)%>%sum))
total_event_counts <- rowSums(skipped_exon[,counts])
# mean <- mean(total_event_counts)
# dev <- sd(total_event_counts)
# hist(total_event_counts)
# max(total_event_counts)
# keep events with an average of 50 counts per group, still a lot left might increase it 
keep <- which(total_event_counts>=100)#trimmed out about half
skipped_exon.filtered <- skipped_exon[keep,]
retnet <- read.csv('../../ref/RetNet_genes.csv',stringsAsFactors = F,header = F)$V1
SE_retnet <- skipped_exon[skipped_exon$geneSymbol%in%retnet,]
unique(SE_retnet$geneSymbol)%>%length

mart <- useMart('ENSEMBL_MART_SNP','hsapiens_snp')
snps.chr15 <- getSNPlocs('15')
