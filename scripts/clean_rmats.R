notes <- '
-for a single subtisue, there are different sites present in diffferent comparisons, might want to look into that
combination=i_combination
event <- i_event
files <- i_files
event_header <- i_event_header
- script wont run on its own, need to remove and move some files gotta fix that



> files <- i_files
> event <- i_event
subtissue <- i_subtissue
'
#setwd('~/NIH/autoRNAseq/')
setwd('/data/swamyvs/eyeintegration_splicing')
library(dplyr)
# i_event <- 'MXE.MATS.JC.txt'
#combination=k[[1]]
#somehow this will fail sometimes as a function but will run fine line by line
combine_PE_SE <- function(combination,event,files,event_header){
    target_files <- files[grepl(combination[1],files)]%>%.[grepl(combination[2],.)]
    if(length(target_files)==1){
        countsCol<-c('IJC_SAMPLE_1','SJC_SAMPLE_1','IJC_SAMPLE_2','SJC_SAMPLE_2')
        tmp <- paste('rmats_out',target_files,event, sep = '/')%>%read.table(header = T,sep = '\t',stringsAsFactors = F)
        #tmp[,countsCol] <- apply(tmp[,countsCol],2, function(x) sapply(x,function(y) strsplit(y,',')%>%unlist%>%as.numeric%>%sum) )
        path <- paste0('rmats_out/',combination[1],'_VS_',combination[2])
        dir.create(path = path)
        write.table(tmp,paste(path,event,sep='/'),row.names = F,col.names = T, quote = F,sep = '\t')
        #if(event=='A3SS.MATs.JC.txt') unlink(target_files,recursive = T)# after all events, remove the folders
        #if(event=='A3SS.MATS.JC.txt') unlink(paste0('rmats_out/',target_files),recursive = T)
        return(0)
    }else if(length(target_files)==0){
        print('REEEEEEEEEEEEEE')
        return(1)
    }
    countsCol<-c('IJC_SAMPLE_1','SJC_SAMPLE_1','IJC_SAMPLE_2','SJC_SAMPLE_2')
    names(target_files) <- grepl('_PE',target_files)%>%ifelse('PE','SE')
    samp_PE <- paste('rmats_out',target_files[grep('_PE',target_files)],event,sep = '/')%>%read.table(header = T,sep = '\t',stringsAsFactors = F)
    samp_SE <- paste('rmats_out',target_files[grep('_SE',target_files)],event,sep = '/')%>%read.table(,header = T,sep = '\t',stringsAsFactors = F)
    if(nrow(samp_SE)==0 || nrow(samp_PE)==0){
      #shitty error handling
      print(combination)
      path <- paste0('rmats_out/',combination[1],'_VS_',combination[2],'/',samp_PE,event)
      write.table(samp_PE,ath,row.names = F,col.names = T, quote = F,sep = '\t')
      path <- paste0('rmats_out/',combination[1],'_VS_',combination[2],'/',samp_SE,event)
      write.table(samp_SE,ath,row.names = F,col.names = T, quote = F,sep = '\t')
      return(1)
    }

    # samp_PE[,countsCol] <- apply(samp_PE[,countsCol],2, function(x) sapply(x,function(y) strsplit(y,',')%>%unlist%>%as.numeric%>%sum) )
    # samp_SE[,countsCol] <- apply(samp_SE[,countsCol],2, function(x) sapply(x,function(y) strsplit(y,',')%>%unlist%>%as.numeric%>%sum) )

    #the first tissue is the first one in combination, the second tissue is the second

    st1_se <- paste0(c('IJC_SAMPLE_','SJC_SAMPLE_'),  grep(combination[1],strsplit(target_files['SE'],'VS')%>%unlist))
    st2_se <- paste0(c('IJC_SAMPLE_','SJC_SAMPLE_'),  grep(combination[2],strsplit(target_files['SE'],'VS')%>%unlist))
    st1_pe <- paste0(c('IJC_SAMPLE_','SJC_SAMPLE_'),  grep(combination[1],strsplit(target_files['PE'],'VS')%>%unlist))
    st2_pe <- paste0(c('IJC_SAMPLE_','SJC_SAMPLE_'),  grep(combination[2],strsplit(target_files['PE'],'VS')%>%unlist))
    #test <- full_join(samp_SE,samp_PE, by=c("chr","strand","exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE"))
    good_cols <- c("GeneID","geneSymbol",event_header[[event]],'IJC_SAMPLE_1','SJC_SAMPLE_1','IJC_SAMPLE_2','SJC_SAMPLE_2',"PValue","FDR")
    samp_PE <- samp_PE[,c("GeneID","geneSymbol",event_header[[event]],st1_pe,st2_pe,"PValue","FDR")]
    colnames(samp_PE) <- good_cols
    samp_SE <- samp_SE[,c("GeneID","geneSymbol",event_header[[event]],st1_se,st2_se,"PValue","FDR")]
    colnames(samp_SE) <- good_cols

    event_header[event]
    #z_merge <- full_join(samp_PE,samp_SE,by=c("chr","strand","exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE"))
    z_merge <- full_join(samp_PE,samp_SE,by=event_header[[event]])
    mergeCols.x <- c("GeneID.x","geneSymbol.x",'IJC_SAMPLE_1.x','SJC_SAMPLE_1.x','IJC_SAMPLE_2.x','SJC_SAMPLE_2.x')
    mergeCols.y <- c("GeneID.y","geneSymbol.y",'IJC_SAMPLE_1.y','SJC_SAMPLE_1.y','IJC_SAMPLE_2.y','SJC_SAMPLE_2.y')
    # fill in na values for info
    parse_count_info <- function(df,col){
        t_cols <- df[,grep(col,colnames(df))]
        comb_counts <- apply(t_cols,2,function(x) strsplit(x,','))
        flattened_counts <- sapply(1:nrow(t_cols),function(x)c(comb_counts[[1]][[x]],comb_counts[[2]][[x]]) )%>%lapply(na.omit)%>%
          lapply(function(x) lapply(x,as.numeric)%>%unlist)
        counts <-sapply(flattened_counts, sum)
        t_mean <- sapply(flattened_counts, mean)
        t_median <- sapply(flattened_counts,median)
        t_sd <-  sapply(flattened_counts,sd)
        final <- data.frame(counts,t_mean,t_median,t_sd)
        colnames(final) <- paste(col,c('counts','avg','med','sd'),sep = '_')
        return(final)
    }


    z_merge[is.na(z_merge$IJC_SAMPLE_1.x),mergeCols.x]  <-z_merge[is.na(z_merge$IJC_SAMPLE_1.x),mergeCols.y]
    # fill in na p-values by just replicatng p-value from sample with valid values, so when we average, it will stay the same
    z_merge$PValue.x[is.na(z_merge$PValue.x)] <- z_merge$PValue.y[is.na(z_merge$PValue.x)]
    z_merge$PValue.y[is.na(z_smerge$PValue.y)] <- z_merge$PValue.x[is.na(z_merge$PValue.y)]
    new_pvalue <- rowMeans(z_merge[c('PValue.x','PValue.y')])
    new_fdr <- p.adjust(new_pvalue,method = "BH")
    final <- data.frame(z_merge[,c("GeneID.x","geneSymbol.x",event_header[[event]],'IJC_SAMPLE_1.x','SJC_SAMPLE_1.x','IJC_SAMPLE_2.x','SJC_SAMPLE_2.x')],new_pvalue,new_fdr,stringsAsFactors = F)
    colnames(final) <- good_cols
    path <- paste0('rmats_out/',combination[1],'_VS_',combination[2])
    dir.create(path = path)
    write.table(final,paste(path,event,sep='/'),row.names = F,col.names = T, quote = F,sep = '\t')

}

combine_fromGTF.novel <- function(event,files,first=TRUE){
    for(path in files){
        if (first==TRUE){
            if(nrow(read.table(paste0(path,'/fromGTF.novelEvents.',event,'.txt'),sep = '\t',header = F,stringsAsFactors = F))>1){
                prev <- read.table(paste0(path,'/fromGTF.novelEvents.',event,'.txt'),sep = '\t',header = T,stringsAsFactors = F)
                first <- FALSE
            }
        }else{
            next1 <- read.table(paste0(path,'/fromGTF.novelEvents.',event,'.txt'),sep = '\t',header = T,stringsAsFactors = F)
            prev <- anti_join(next1,prev)%>%rbind(.,prev)# BAAAAAAAAD
        }
    }
    return(prev)
}


i_event_header <- list(SE.MATS.JC.txt=c('chr'	,'strand',	'exonStart_0base',	'exonEnd',	'upstreamES',	'upstreamEE',	'downstreamES',	'downstreamEE'),
                       RI.MATS.JC.txt=c('chr'	,'strand',	'riExonStart_0base',	'riExonEnd'	,'upstreamES'	,'upstreamEE'	,'downstreamES'	,'downstreamEE'),
                       MXE.MATS.JC.txt=c('chr',	'strand',	'X1stExonStart_0base',	'X1stExonEnd',	'X2ndExonStart_0base',	'X2ndExonEnd'	,'upstreamES',	'upstreamEE',	'downstreamES',	'downstreamEE'),
                       A5SS.MATS.JC.txt=c('chr',	'strand',	'longExonStart_0base',	'longExonEnd',	'shortES',	'shortEE',	'flankingES',	'flankingEE'),
                       A3SS.MATS.JC.txt=c('chr',	'strand',	'longExonStart_0base',	'longExonEnd'	,'shortES',	'shortEE'	,'flankingES',	'flankingEE')
                      )


events <- names(i_event_header)
i_files <- dir('rmats_out')
subtissues_PE <- c("Retina_Adult.Tissue", "RPE_Cell.Line", "ESC_Stem.Cell.Line" , "RPE_Adult.Tissue","body" )# add body back in  at some point
k <- combn(subtissues_PE,2,simplify = F)
for (i in 1:length(k)){
  i_combination <- k[[i]]
  for(j in 1:length(events)){
    i_event <- events[j]
    combine_PE_SE(combination = i_combination,event = i_event,files = i_files,event_header = i_event_header)
  }
}# add PESE


# generate tables with all novel events
# t <- c('SE','RI','MXE','A5SS','A3SS')
# files <- dir('~/NIH/autoRNAseq/old_rmats_out',full.names = T)
# for(i in t){
#     all_ev <-  combine_fromGTF.novel('SE',files)
#     write.table(all_ev,paste0('all.',i,'.novelevents.txt'),quote = F,col.names = T,row.names = F,sep = '\t')
# }



# for (combination in k){
#   target_files <- i_files[grepl(combination[1],i_files)]%>%.[grepl(combination[2],.)]%>%paste0('rmats_out/',.)
#   print(target_files)
#   unlink(target_files,recursive = T)
# } # remove stuff
#system2('mv rmats_out/* rmats_comb/')
#couldnt get that^ to work, might have to just run it separately
#nowcombine everything together
#hcekc to see if body files are alive

#strat: combine all novel events per type in one master file, then select specific ones
