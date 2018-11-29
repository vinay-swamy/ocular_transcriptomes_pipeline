setwd('~/NIH/eyeintegration_splicing/')
library(dplyr)
library(rtracklayer)
library(tximport)
# there was a bug in the clean gtf script that was only using protein coding genes to determine novelty, but we used
# the full pc-gtf for guiding assembly, so gotta keep it in for consistensy its removed in the script, but currently using old
# data so once its run again you can remove the filtering when reading files



load('ref/gtf_unformatted.Rdata')
load('rdata/plotting_header.Rdata')
ref_gtf <- readGFF('ref/gencodeAno_bsc.gtf')%>%filter(gene_type == "protein_coding")
gtf_final <- filter(gtf_final,gene_name%in%ref_gtf$gene_name)
rpe_at <- read.table("rmats_analysis/RPE_Adult.Tissue_PE_VS_body_PE/SE.MATS.JC.txt" ,header = T, sep = '\t',stringsAsFactors = F)%>%filter(geneSymbol%in%ref_gtf$gene_name)
ret_at <- read.table('rmats_analysis/Retina_Adult.Tissue_VS_body/SE.MATS.JC.txt',header = T, sep = '\t',stringsAsFactors = F)%>%filter(geneSymbol%in%ref_gtf$gene_name)
corn_at <- read.table('rmats_analysis/Cornea_Adult.Tissue_SE_VS_body_SE/SE.MATS.JC.txt',header = T, sep = '\t',stringsAsFactors = F)%>%filter(geneSymbol%in%ref_gtf$gene_name)
colnames(rpe_at) <- colnames(ret_at) <- colnames(corn_at) <- plotting_header
rpe_at$rpe <- paste0('rpe-',1:nrow(rpe_at))
ret_at$ret <- paste0('ret-',1:nrow(ret_at))
corn_at$corn <- paste0('corn-', 1:nrow(corn_at))
gtf_final$st <-paste0('st-',1:nrow(gtf_final))# create a distinct value for easy joining
num_chr <- function(seqid){
    ch_seq <- strsplit(seqid,'chr')%>%lapply(function(x) x[[2]])%>%unlist
    ch_seq[ch_seq=='X']='23'
    ch_seq[ch_seq=='Y']='24'
    ch_seq[ch_seq=='M']='25'
    return(as.numeric(ch_seq))
}
gtf_final$seqid <- num_chr(gtf_final$seqid)
rpe_at$chr <- num_chr(rpe_at$chr)
ret_at$chr <- num_chr(ret_at$chr)
corn_at$chr <- num_chr(corn_at$chr)
novel_exons <- filter(gtf_final,type=='exon', is_novel_exon=='y')
#novel_exons=data.frame(novel_exons,ret_loc=F,rpe_loc=F,corn_loc=F)
#we need to figure out how many of the so called novel tx's are actually novel.
#rmats uses 0base convention
novel_exons$start <- novel_exons$start-1
ref_gtf$start <- ref_gtf$start-1
gtf_final$start <- gtf_final$start-1
# need to check each gene, and figure out some stats about it.
genes_to_examine <- c(corn_at$geneSymbol,ret_at$geneSymbol,rpe_at$geneSymbolnovel_exons$gene_name)%>%unique
# data frame gene start end found in stringtie found in ret found in rpe found in cornea
# filter novel events for a gene, and rmats events for a gene;

pull_rmats_exons <- function(rmats, gtf, gene){
    rmats <- filter(rmats, geneSymbol==gene)
    if(nrow(rmats)==0)return(NULL)
    rmats_exons <- split(rmats[,c(3,5:6)],1:nrow(rmats))%>%lapply(as.numeric)
    names(rmats_exons) <- rmats[,29]
    gtf_exons <- split(gtf[,c('seqid','start','end')],1:nrow(gtf))%>%lapply(as.numeric)
    keep <- lapply(rmats_exons,function(x) any( lapply ( gtf_exons, function(y) all(x==y))))%>%unlist
    rm_novel_exons <- rmats_exons[!keep]
    rm_novel_exons <- rm_novel_exons[!duplicated(rm_novel_exons)]
    
    return(rm_novel_exons)
}
final_table <- list(data.frame())
l=1
#this works beautifully
for (gene in genes_to_examine){
    ref <- filter(ref_gtf,type=='exon',gene_name==gene)
    stringTie_novel <- filter(novel_exons,gene_name==gene)
    st_exons <- split(stringTie_novel[,c('seqid','start','end')],1:nrow(stringTie_novel))%>%lapply(as.numeric)
    names(st_exons) <- stringTie_novel[,'st']
    st_exons <- st_exons[!duplicated(st_exons)]
    if(nrow(stringTie_novel)==0)st_exons <- NULL
    all_exons <- list(pull_rmats_exons(ret_at,ref,gene),pull_rmats_exons(rpe_at,ref,gene),pull_rmats_exons(corn_at,ref,gene),st_exons)
    names(all_exons) <- c('ret','rpe','corn','st')
#   all_exons <- all_exons[!(lapply(all_exons,is.null)%>%unlist)]
    #pick the shortest list; for each exon in the other lists, check whether the exon is present in the others;
    #if it is, remove it. mark accordingly where each exon is or isnt; continue until list is empty.
    exon_table <- matrix(data = NA,nrow = sum(lapply(all_exons,length)%>%unlist),ncol=11)%>%as.data.frame()
    colnames(exon_table) <- c('gene','chr' ,'start','end','ret','rpe','corn','st','nearest_exon_start','nearest_exon_end','n_exon_mp_dist')
    e <- 1
    all_exons <- all_exons[which(sapply(all_exons,length)>0)]
    while(length(all_exons)>0){
        # find tissue w/ lowest number of exons per gene
        m <- sapply(all_exons,length)%>%which.min()
        #stack.peep
        curr_exon <- all_exons[[m]][[1]]
        res <- data.frame(gene=gene,chr=curr_exon[[1]],start=curr_exon[[2]],end=curr_exon[[3]],ret=NA,rpe=NA,corn=NA,st=NA, 
                          nearest_exon_start=NA,nearest_exon_end=NA,n_exon_mp_dist=NA,stringsAsFactors = F)
        #identify id exons are in other tissues
        search <- lapply(all_exons,function(y)lapply(y, function(x) all(x==curr_exon))%>%unlist%>%which)
        search <- search[(lapply(search,length)%>%unlist)>0]
        t_names <- names(search)
        
        for(name in t_names) {
            #print(names(search[[name]]))
            res[,name] <- names(search[[name]])
            #stack.pop
            all_exons[[name]] <- all_exons[[name]][-search[[name]]]
        }
        #update table
        
        #find nearest exon
        gtf_exons <- split(ref[,c('start','end')],1:nrow(ref))%>%lapply(as.numeric)
        curr_exon_mp <- mean(curr_exon)
        gtf_exons_mp <- lapply(gtf_exons,mean)
        exon_diffs <- lapply(gtf_exons_mp,function(x) abs(x - curr_exon_mp ))
        closest_exon <- which.min(exon_diffs)
        res[,c('nearest_exon_start','nearest_exon_end')] <-gtf_exons[[closest_exon]]
        res[,'n_exon_mp_dist'] <- exon_diffs[[closest_exon]]
        exon_table[e,] <- res
        e <- e+1
        # remove empty lists
        all_exons <- all_exons[which(sapply(all_exons,length)>0)]
    }
    exon_table <- exon_table[!is.na(exon_table$gene),]

    final_table[[l+1]] <- exon_table
    l <- l+1
}
final_table <- do.call(rbind,final_table)
save(final_table,file = 'rmats_StringTieExonComp.Rdata')

#TODO add  add information about nearest exon, and count info for each sample
header <- colnames(final_table)
final_table <- left_join(final_table,gtf_final[,c('st','transcript_id')],by=c('st'))
tst <-left_join(final_table,ret_at[,c('ret','IJC_SAMPLE_1_med')],by=c('ret'))%>%
    left_join(rpe_at[,c('rpe', 'IJC_SAMPLE_1_med')],by=c('rpe'))%>%
    left_join(corn_at[,c('corn','IJC_SAMPLE_1_med')],by=c('corn'))

tissues <- c('Retina','RPE','Cornea','ESC','body')
sample_table <- read.table('samplesrun_1010.tissue',stringsAsFactors = F, header = F,sep = '\t')
colnames(sample_table) <- c('sample','run','paired','tissue','subtissue','origin')
samples_eye <- filter(sample_table,tissue%in%c('Retina','RPE','Cornea','ESC'))
body_synth <- read.table('ref/synth_body.txt',stringsAsFactors = F, header = F)
sample_table <- rbind(samples_eye,filter(sample_table, sample%in%body_synth$V1))
sample_table[!sample_table$tissue%in%c('Retina','RPE','Cornea','ESC'),c('tissue','subtissue')] <- 'body'
files <- paste0('quant_files/',sample_table$sample,'/quant.sf')
txdb <- filter(gtf_final,type=='transcript')%>%select(transcript_id,gene_name,st)
files <- paste0('quant_files/',sample_table$sample,'/quant.sf')
txi <- tximport::tximport(files = files,type = 'salmon',tx2gene = txdb[,1:2],countsFromAbundance = 'lengthScaledTPM',txOut=T)
colnames(txi$counts) <- sample_table$sample
counts <- txi$counts%>%as.data.frame()%>%.[txdb$transcript_id,]
txi.counts_by_tissue <- lapply(tissues, function(x) filter(sample_table,tissue==x)[,'sample']%>%counts[,.])%>%
    lapply(as.data.frame)%>%lapply(rowMeans)%>%do.call(cbind,.)%>%as.data.frame
colnames(txi.counts_by_tissue) <- tissues
txi.counts_by_tissue$transcript_id <- rownames(txi.counts_by_tissue)
final_table <- left_join(tst,txi.counts_by_tissue,by=c('transcript_id'))
colnames(final_table) <- c(header,'transcript_id','med_inc_cts_Retina','med_inc_cts_RPE','med_inc_cts_Cornea','lstpm_Retina','lstpm_RPE','lstpm_Cornea','lstpm_ESC','lstpm_Body')
save(final_table,file = 'rmats_stringTie_comparison_all.Rdata')
##looking at double confidence genes, and genes that contain unshared novel st/rmats exons
boolin <- ((!is.na(final_table$ret)) | (!is.na(final_table$rpe)) | (!is.na(final_table$corn))) & (!is.na(final_table$st))
final_table_double_conf <- final_table[boolin,]
save(final_table_double_conf,file = 'stringTie_double_confidence_genes')

test <- final_table
test <- split(final_table,final_table$gene)
k=test[[1]]
res=list()

for( k in test){
    log <- split(k, 1:nrow(k))%>%lapply(function(x) c((all(is.na(x[,4:6]))&(!is.na(x[,7]))), (any(!is.na(x[,4:6]))&(is.na(x[,7])))))%>%do.call(rbind,.)
    keep <- any(log[,1])&any(log[,2])
    if(keep) res[[length(res)+1]] <- k
}
res <- do.call(rbind, res)
res <- as.data.frame(res)
res_t <- res[!complete.cases(res),]
res <- arrange(res,gene, start)






elon_novel<- c('st-116949','st-116950','st-218195','st-521884')
