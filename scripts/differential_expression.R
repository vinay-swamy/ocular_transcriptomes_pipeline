setwd('~/NIH/eyeintegration_splicing/')
library(dplyr)
load('rdata/plotting_header.Rdata')
load('rmats_stringTie_comparison_all.Rdata')
novel_exon_tab <- final_table
rpe_at <- read.table("rmats_analysis/RPE_Adult.Tissue_PE_VS_body_PE/SE.MATS.JC.txt" ,header = T, sep = '\t',stringsAsFactors = F)
ret_at <- read.table('rmats_analysis/Retina_Adult.Tissue_VS_body/SE.MATS.JC.txt',header = T, sep = '\t',stringsAsFactors = F)
corn_at <- read.table('rmats_analysis/Cornea_Adult.Tissue_SE_VS_body_SE/SE.MATS.JC.txt',header = T, sep = '\t',stringsAsFactors = F)
colnames(rpe_at) <- colnames(ret_at) <- colnames(corn_at) <- plotting_header
rpe_at$rpe <- paste0('rpe-',1:nrow(rpe_at))
ret_at$ret <- paste0('ret-',1:nrow(ret_at))
corn_at$corn <- paste0('corn-', 1:nrow(corn_at))

pull_rmats_exons <- function(rmats, gtf, gene){
    rmats <- filter(rmats, geneSymbol==gene, new_fdr<=.05)
    if(nrow(rmats)==0) return(NULL)
    rmats_exons <- split(rmats[,5:6],1:nrow(rmats))%>%lapply(as.numeric)
    names(rmats_exons) <- rmats[,29]
    rm_novel_exons <- rmats_exons
    rm_novel_exons <- rm_novel_exons[!duplicated(rm_novel_exons)]
    # return start, end
    return(rm_novel_exons)
}
genes_to_examine <- unique(c(ret_at$geneSymbol, rpe_at$geneSymbol,corn_at$geneSymbol))
final_table <- list(data.frame())
l=1
#a little slow, but allows to account for  different novel patterns within a single gene
for (gene in genes_to_examine){
    all_exons <- list(pull_rmats_exons(ret_at,ref,gene),pull_rmats_exons(rpe_at,ref,gene),pull_rmats_exons(corn_at,ref,gene))
    names(all_exons) <- c('in_ret','in_rpe','in_corn')
    #   all_exons <- all_exons[!(lapply(all_exons,is.null)%>%unlist)]
    #pick the shortest list; for each exon in the other lists, check whether the exon is present in the others;
    #if it is, remove it. mark accordingly where each exon is or isnt; continue until list is empty.
    exon_table <- matrix(data = NA,nrow = sum(lapply(all_exons,length)%>%unlist),ncol=9)%>%as.data.frame()
    colnames(exon_table) <- c('gene','start','end','in_ret','nov_ret','in_rpe','nov_rpe','in_corn','nov_corn')
    e <- 1
    all_exons <- all_exons[which(sapply(all_exons,length)>0)]
    while(length(all_exons)>0){
        # find tissue w/ lowest number of exons per gene
        m <- sapply(all_exons,length)%>%which.min()
        #stack.peep
        curr_exon <- all_exons[[m]][[1]]
        res <- data.frame(gene=gene,start=curr_exon[[1]],end=curr_exon[[2]],in_ret=NA,nov_ret='n',in_rpe=NA,nov_rpe='n',in_corn=NA,
                          nov_corn='n',stringsAsFactors = F)
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
        tmp <- res
        for( tis in c('in_ret', 'in_rpe','in_corn')){
            i=res[,tis]
            if(is.na(res[,tis]))  i <- NULL
            t <- strsplit(tis,'_')[[1]][2]
            n_t <- paste0('nov_',t)
            if(any(novel_exon_tab[,t]%in%i)) res[,n_t] <- 'y'
            
        }
        exon_table[e,] <- res
        e <- e+1
        # remove empty lists
        all_exons <- all_exons[which(sapply(all_exons,length)>0)]
    }
    exon_table <- exon_table[!is.na(exon_table$gene),]
    
    final_table[[l+1]] <- exon_table
    l <- l+1
}
rmats_sig_table <- do.call(rbind,final_table)
save(rmats_sig_table,file = 'rmats_sig_nov_by_tiss.Rdata')
#########################################################################################
#limma 
library(limma)
library(edgeR)
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
counts <- txi$counts%>%as.data.frame()
keep <- which(rowSums(counts)>2*ncol(counts))
counts <- counts[keep,]
tissue_fac <- as.factor(sample_table$tissue)
design_mat <- model.matrix(~0+tissue_fac)
colnames(design_mat) <- levels(tissue_fac)
counts_norm <- calcNormFactors.DGEList(DGEList(counts))
counts_voom <- voom(counts_norm,design_mat)
conts <- makeContrasts(Retina_vs_body='Retina-body',
              RPE_vs__body='RPE-body',
              Cornea_vs_body='Cornea-body',
              ESC_vs_body='ESC-body',
              levels = design_mat)
lmfit <- lmFit(counts_voom,design_mat)
lmf_cont <- contrasts.fit(lmfit, conts)
efit <- eBayes(lmf_cont)
save(efit, file='diff_exp_by_tissue.Rdata')
deg <- topTableF(efit,number = 300000, adjust.method = 'fdr')
pvals <- efit$p.value 
padj <- apply(pvals,2, function(x)p.adjust(x,method = 'fdr'))
names(padj) <- rownames(pvals)







