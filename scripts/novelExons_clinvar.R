library(biomaRt)
library(tidyverse)
library(GenomicRanges)
library(limma)
# sel_vars <- read.csv('nc+pot_var.csv')[,1:14]
# novel_exon_tab$ljid <- 1:nrow(novel_exon_tab)
# #sel_vars$nx_id <- novel_exon_tab$ [sel_vars$nx_id]
# sel_vars <- dplyr::rename(sel_vars, ljid=nx_id)
# sel_var_exp <- left_join(sel_vars,novel_exon_tab[,c(1:8,12,16:21)],by=c('ljid'))
# load('diff_exp_by_tissue.Rdata')
# ttf <- topTableF(efit,number = 30000000,adjust.method = 'fdr')
# ttf$transcript_id <- rownames(ttf)
# sel_var_exp_de <- left_join(sel_var_exp,ttf, by=c('transcript_id'))
# save(sel_var_exp_de,file='retNet_inClinvar_complete_info.Rdata') 
# ^ to make first variatn table



mart <- useMart('ENSEMBL_MART_SNP',dataset = 'hsapiens_snp')
listAttributes(mart = mart)
#clinvar_full <- read.table('ref/variant_summary_11.05.18.txt', stringsAsFactors = F, sep = '\t', header = T, comment.char = '', fill = T)

at <- c('refsnp_id',  'validated' ,'clinical_significance' ,'pmid', 'consequence_type_tv')
qu <- c('chr_name','start','end')
val <- list(chr_name=1,start=gtf_final$start[1],end=gtf_final$end[1])
p <- getBM(attributes = at,filters = qu,values = val,mart = mart)
chr_name
chrom_start
chrom_end
chrom_strand

#rs1635544 has reveal intresting results, but we gotta take it further
# k <- clinvar_in_rn[!grepl('p.',clinvar_in_rn$Name),]
# k_bm <- getBM(attributes = at,filters = 'snp_filter',values = list(paste0('rs',k$RS...dbSNP.)), mart = mart)
# k <- rename(k,RS...dbSNP.='refsnp_id')
# k_bm$refsnp_id <-  strsplit(k_bm$refsnp_id,'rs')%>%lapply(function(x)x[[2]])%>%unlist%>%as.numeric
# clinvar_inret_bm <- left_join(k,k_bm,by=c('refsnp_id'))

clinvar_sum <- read.table('ref/variant_summary_11.05.18.txt', stringsAsFactors = F, sep = '\t', header = T, comment.char = '', fill = T)%>%
     #dplyr::select(AlleleID,VariationID,Type,GeneSymbol,`RS...dbSNP.`,Origin,Assembly,Chromosome,Start,S)%>%
     filter(complete.cases(.),Assembly=='GRCh38',Type%in%c("single nucleotide variant","deletion","insertion" ,"indel"))

num_chr <- function(seqid){
    ch_seq <- seqid
    ch_seq[ch_seq=='X']='23'
    ch_seq[ch_seq=='Y']='24'
    ch_seq[ch_seq=='MT']='25'
    return(as.numeric(ch_seq))
}
clinvar_sum$Chromosome <- num_chr(clinvar_sum$Chromosome)

load('rmats_stringTie_comparison_all.Rdata')
novel_exon_tab <- final_table
 novel_exon_tab$start <- novel_exon_tab$start-12
 novel_exon_tab$end <- novel_exon_tab$end+12
clinvar_sum$ljid <- 1:nrow(clinvar_sum)
novel_exon_range <- GRanges(seqnames =novel_exon_tab$chr, ranges = IRanges(novel_exon_tab$start, novel_exon_tab$end))
clinvar_range <- GRanges(seqnames = clinvar_sum$Chromosome,ranges =  IRanges(clinvar_sum$Start, clinvar_sum$Stop))
com <- findOverlaps (novel_exon_range,clinvar_range)%>%as.data.frame()
clinvar_sum$nx_id <- NA
clinvar_sum$nx_id[com$subjectHits] <- com$queryHits
retnet <- read.csv('ref/RetNet_genes.csv', header = F, stringsAsFactors = F)[,1]
com_in_ret_net <- filter(com, queryHits%in%which(novel_exon_tab$gene%in%retnet))
clinvar_in_rn <- clinvar_sum[com_in_ret_net$subjectHits,]%>%filter(RS...dbSNP.!=-1)
mart <- useMart('ENSEMBL_MART_SNP',dataset = 'hsapiens_snp')
at <- c('refsnp_id',  'validated' ,'clinical_significance' ,'pmid', 'consequence_type_tv')
bmres <- getBM(attributes = at,filters = 'snp_filter',values = list(paste0('rs',clinvar_in_rn$RS...dbSNP.)), mart = mart)
clinvar_in_rn <- rename(clinvar_in_rn,RS...dbSNP.='refsnp_id')%>%dplyr::select(nx_id,Type,Name,GeneSymbol,refsnp_id,Chromosome,Start,Stop)
bmres$refsnp_id <-  strsplit(bmres$refsnp_id,'rs')%>%lapply(function(x)x[[2]])%>%unlist%>%as.numeric
cvr_irn_bm <- left_join(clinvar_in_rn,bmres,by=c('refsnp_id'))#hs
hs_IntVar <- filter(cvr_irn_bm,!grepl('p.',Name),consequence_type_tv=='intron_variant')
hs_noIntVar <- filter(cvr_irn_bm, !grepl('p.',Name),consequence_type_tv!='intron_variant')
hs_allvar <- rbind(hs_IntVar,hs_noIntVar,stringsAsFactors=F)
hs_allvar <- hs_allvar[!duplicated(hs_allvar$refsnp_id),]
##hand scanned all variants 
sel_vars <- hs_allvar
novel_exon_tab$ljid <- 1:nrow(novel_exon_tab)
#sel_vars$nx_id <- novel_exon_tab$ [sel_vars$nx_id]
sel_vars <- dplyr::rename(sel_vars, ljid=nx_id)
sel_var_exp <- left_join(sel_vars,novel_exon_tab[,c(1:8,12,16:21)],by=c('ljid'))
load('diff_exp_by_tissue.Rdata')
ttf <- topTableF(efit,number = 30000000,adjust.method = 'fdr')
ttf$transcript_id <- rownames(ttf)
sel_var_exp_de <- left_join(sel_var_exp,ttf, by=c('transcript_id'))%>%filter(Retina_vs_body>=1)


#variants_clinvar_nochange <- sel_var_exp_de
#variants_clinvar_plus3 <- sel_var_exp_de
#variants_clinvar_plus6 <- sel_var_exp_de
variants_clinvar_plus12 <- sel_var_exp_de
save(variants_clinvar_nochange,variants_clinvar_plus12,variants_clinvar_plus3,variants_clinvar_plus6,file = 'variants_clinvar_exon Change.Rdata')
##################################################################################################################
# starting with data from david
load('remapped_clinvar_all.Rdata')
clinvar_unSig <- updated_clinvar
clinvar_unSig$chrom[clinvar_unSig$chrom=='X']='23'
clinvar_unSig$chrom <- as.numeric(clinvar_unSig$chrom)
load('rmats_stringTie_comparison_all.Rdata')
load('ref/refgtf.Rdata')


novel_exon_tab <- final_table
# expand to splice sites
novel_exon_tab$start <- novel_exon_tab$start-12
novel_exon_tab$end <- novel_exon_tab$end+12

novel_exon_range <- GRanges( seqnames =novel_exon_tab$chr, ranges = IRanges(novel_exon_tab$start, novel_exon_tab$end))
cv_unSig_range <-   GRanges(seqnames = clinvar_unSig$chrom, ranges =  IRanges(clinvar_unSig$start,clinvar_unSig$stop))
overlap <- GenomicRanges::findOverlaps(novel_exon_range,cv_unSig_range)%>%as.data.frame()
clinvar_unSig$subjectHits <- 1:nrow(clinvar_unSig)
novel_exon_tab$queryHits <- 1:nrow(novel_exon_tab)
longboi <- left_join(overlap,clinvar_unSig,by=c('subjectHits'))%>%left_join(novel_exon_tab,by=c('queryHits'))
clinvar_david_VUS <- longboi[,c(3:12,14:15,19,31,42:49,53,57:61)]
load('diff_exp_by_tissue.Rdata')
ttf <- topTableF(efit,number = 30000000,adjust.method = 'fdr')
ttf$transcript_id <- rownames(ttf)
cv_david_VUS_de <- left_join(clinvar_david_VUS,ttf,by=c('transcript_id'))
# #2 exons match; one is 2K long, but matches to 3 samples, exon is labeled as F3 in my data but ABCA4 is davids?
# f3_vc <- clinvar_unSig[30:32,]
# f3_gtf <- filter(gtf_final,gene_name=='F3')
# f3_cts <- txi.counts_by_tissue[rownames(txi.counts_by_tissue)%in%f3_gtf$transcript_id,]
# 
# #from overlap, this one is the simplest to look at
# rax_vc <- clinvar_unSig[532,]
# rax_nx <- novel_exon_tab[18607,]
# rax_stgtf <- filter(gtf_final, gene_name=='RAX2')
# rax_cts <- txi.counts_by_tissue[rownames(txi.counts_by_tissue)%in%rax_stgtf$transcript_id,]
# #####
# novel_exon_tab$ljid <- 1:nrow(novel_exon_tab)
# vars <- clinvar_unSig[c(30:32,532),c(1,3:6,12:13,17)]
# vars$ljid <- overlap$queryHits
# vars_exp <- left_join(vars,novel_exon_tab[,c(1:8,12,16:21)], by=c('ljid'))
# load('diff_exp_by_tissue.Rdata')
# ttf <- toptablef(efit,)
# vars_exp_de <-  left_join(vars_exp,ttf, by=c('transcript_id'))
# save(vars_exp_de,file='VUS_david_expDE.Rdata')
#######
novel_exon_range <- GRanges( seqnames =novel_exon_tab$chr, ranges = IRanges(novel_exon_tab$start, novel_exon_tab$end))
clinvar_path <- updated_clinvar_path
clinvar_path$chrom[clinvar_path$chrom=='X']='23'
clinvar_path$chrom <- as.numeric(clinvar_path$chrom)
cv_path_range <- GRanges(seqnames = clinvar_path$chrom, ranges = IRanges(start = clinvar_path$start, end = clinvar_path$stop))
overlap_path <- GenomicRanges::findOverlaps(novel_exon_range, cv_path_range)%>%as.data.frame()
clinvar_path$subjectHits <- 1:nrow(clinvar_path)
novel_exon_tab$queryHits <- 1:nrow(novel_exon_tab)
longboi <- left_join(overlap_path,clinvar_path,by=c('subjectHits'))%>%left_join(novel_exon_tab,by=c('queryHits'))
cv_david_solved <- longboi[,c(3:11,14:15,19,31,42:49,53,57:61)]
load('diff_exp_by_tissue.Rdata')
ttf <- topTableF(efit,number = 30000000,adjust.method = 'fdr')
ttf$transcript_id <- rownames(ttf)
cv_david_solved <- left_join(cv_david_solved,ttf,by=c('transcript_id'))

#########
uk10k <-  read_csv('ref/uk10k_intronic.csv')
k <- strsplit(uk10k$Variant_genomic,':|>|T|A|C|G')%>%lapply(as.numeric)%>%lapply(function(x) c(x[1:2],x[[2]]+1))
k[[4]][[1]] <- 23
k[[11]][[1]] <- 23
uk_rang <- GRanges(seqnames = updated_clinvar_uk10k$chrom , ranges = IRanges(start = updated_clinvar_uk10k$start,
                                                                                   end = updated_clinvar_uk10k$stop))
findOverlaps(novel_exon_range,uk_rang)# nada pt 2


##################################
var_sig <- cv_david_VUS_de[c(2,15,13,3),]%>%rbind(cv_david_solved)
load('variants_on_novel_exons.Rdata') 
retnet_in_clinvar <- retnet_in_clinvar[4:6,]
save(retnet_in_clinvar,var_sig, file = 'vars_on_nx_sig.Rdata')
