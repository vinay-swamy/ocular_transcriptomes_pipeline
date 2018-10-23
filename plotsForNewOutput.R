library(rtracklayer)
library(dplyr)
library(ggplot2)
library(data.table)
#use  these for solving equations for edges
test_fun <- function(x,d,fstart,fend,h){
  mp=mean(c(fstart, fend))
  a=h/((mp-fstart)*(mp-fend))
  a*(x-fstart)*(x-fend)
}
quad <- function(x) return(c(x^2,x,1))
solver <- function(p1,p2,p3){
  mat <- rbind(quad(p1[1]),quad(p2[1]),quad(p3[1]))
  y <- matrix(c(p1[2],p2[2],p3[2]),nrow = 3,ncol = 1)
  return(solve(mat,y))
}
eq <- function(x,cons){
  a=cons[1]
  b=cons[2]
  c=cons[3]
  a*x^2 + b*x + c
}

#setwd('~/NIH/autoRNAseq/spleen_vs_skin')
#load('plotting_curvs.Rdata')
#assuming '+' strand
#gdata::keep(gtf,sure=T)

setwd("~/NIH/eyeintegration_splicing/")
gtf <- readGFF('~/NIH/eyeintegration_splicing/ref/combined_final.gtf')
gtf <- load('ref/gtf_unformatted.Rdata')
all_events=read.table('SE.MATS.JC.txt',stringsAsFactors = F,header = T)
#all_events[,c('IJC_SAMPLE_1','SJC_SAMPLE_1')] <- apply(all_events[,c('IJC_SAMPLE_1','SJC_SAMPLE_1')],2, function(x) sapply(x,function(y) strsplit(y,',')%>%unlist%>%as.numeric%>%sum) )
#only protein coding tx's
#gtf_exons<- gtf%>%filter(type=='exon')

splice_graph<- function(gene,gtf,all_events){
  genes <- filter(all_events, geneSymbol==gene)
  gtf.gene <- filter(gtf,gene_name==gene,type=='exon')%>%arrange(start)#%>%select(1:7)
  gtf.gene$start=gtf.gene$start-1
  #exons ends are ok, exon starts are off by 1
  tx_ids <- unique(gtf.gene$transcript_id)
  i=tx_ids[3]
  gtf.tocomp <- 
  #plot every tx in a gene
  for(i in tx_ids){
    tx_name=i
    gtf.tx <- filter(gtf.gene,transcript_id==tx_name,type=='exon')%>%arrange(start)#%>%select(1:7)
    gtf.tocomp <- filter(gtf.gene,transcript_id==tx_ids[1],type=='exon')%>%arrange(start)
    #pick rmats events in which all three exons match thsoe in GTF
    exon_coords <- as.matrix(gtf.tx[,c('start','end')]) 
    exon_col_names <- c("exonStart_0base","exonEnd" ,"upstreamES", "upstreamEE" ,"downstreamES" ,"downstreamEE")
    rmats_exon_coords <- as.matrix(genes[,exon_col_names])
    #there has to be a better way to do this
    
    to_plot=logical(nrow(rmats_exon_coords))
    # for(j in 1:nrow(rmats_exon_coords)){
    #     #b's are exon start and ends in rmats data, b is true if rmats exon start & end match gtf 
    #     #now that we are using rmats to find new stuff, might want to consider keeping all exons
    #     b1 <- (apply(exon_coords,1 ,function (x) all(x==rmats_exon_coords[j,1:2])))
    #     b2 <- any(apply(exon_coords,1 ,function (x) all(x==rmats_exon_coords[j,3:4])))
    #     b3 <- any(apply(exon_coords,1 ,function (x) all(x==rmats_exon_coords[j,5:6])))
    #     to_plot[j] <- b1 && b2 && b3
    # }
    
    for(j in 1:nrow(rmats_exon_coords)){
      #b's are exon start and ends in rmats data, b is true if rmats exon start & end match gtf 
      #now that we are using rmats to find new stuff, might want to consider keeping all exons
        b1 <- any(apply(exon_coords,1 ,function (x) isTRUE(all.equal(x,rmats_exon_coords[j,1:2],tolerance = 15))))
        b2 <- any(apply(exon_coords,1 ,function (x) isTRUE(all.equal(x,rmats_exon_coords[j,1:2],tolerance = 15))))
        b3 <- any(apply(exon_coords,1 ,function (x) isTRUE(all.equal(x,rmats_exon_coords[j,1:2],tolerance = 15))))
      to_plot[j] <- b1 && b2 && b3
    }
    rmats_to_plot <- genes
    rmats_to_plot$in_gtf <- to_plot
    novel_exons <- genes[!to_plot,]
    
    #prep_plots
    if(nrow(rmats_to_plot)<2) {
      print('plop')  
      return(0)
    }
    start <- min(gtf.tx$start[1]/100, min(rmats_to_plot$upstreamES)/100,min(novel_exons$upstreamES)/100)
    #scale genome coordinates
    gtf.tx[,c('start','end')] <- gtf.tx[,c('start','end')]/100 -start
    gtf.tocomp[c('start','end')] <- gtf.tocomp[,c('start','end')]/100 -start
    rmats_to_plot[,exon_col_names] <- rmats_to_plot[,exon_col_names]/100-start
    novel_exons[,c('exonStart_0base','exonEnd')] <- novel_exons[,c('exonStart_0base','exonEnd')]/100-start
    #current idea: each event contains three edges: a single long edge in the skipped isoform and 2 short edges in the included
    #isoform. RMATS gives counts for the the exclusion event, and inclusion event. So the weight of the long edge 
    #is just the relative counts for the exlcusion event. for the short edges the inclusion counts are distributed over 2 edges,
    #so split the weight over the edges. remember to account for split merging edges.
    # also incusion edges counts are going to be duplicated when merged below, so cut their counts by 2 again(so divide by 4)
    rmats_to_plot[,"IJC_SAMPLE_1_avg"] <-  rmats_to_plot[,"IJC_SAMPLE_1_avg"]/4
    
    
    prep_for_plots <- function(rmats_to_plot,in_gtf){
    edges <- rbind(as.matrix(rmats_to_plot[,c("upstreamEE","downstreamES", "SJC_SAMPLE_1_avg")]),# skipped junction
                   as.matrix(rmats_to_plot[,c("upstreamEE","exonStart_0base","IJC_SAMPLE_1_avg")]),#included junction
                   as.matrix(rmats_to_plot[,c("exonEnd","downstreamES","IJC_SAMPLE_1_avg")]))#included junction 
    edgeList <- split(edges,seq(nrow(edges)))
    
    edges <- as.data.frame(edges)
    rmats_to_plot <- aggregate(edges,list(edges$upstreamEE,edges$downstreamES),sum)
    rmats_to_plot <-  dplyr::select(rmats_to_plot,-c('upstreamEE','downstreamES'))
    colnames(rmats_to_plot) <- c('start','end','counts')
    rmats_to_plot$in_gtf <- in_gtf
    return(rmats_to_plot)
    }
    rmats_to_plot <- rbind(prep_for_plots(rmats_to_plot = filter(rmats_to_plot,in_gtf==TRUE),in_gtf = T),
                           prep_for_plots(rmats_to_plot = filter(rmats_to_plot,in_gtf==F),in_gtf = F),stringsAsFactors=F)
    #duplicated(rmats_to_plot[])
    
    total_counts <- sum(rmats_to_plot[,'counts'])
    rmats_to_plot$total_weight=round(rmats_to_plot[,'counts']/total_counts,3)
    rmats_to_plot <- rmats_to_plot[rmats_to_plot$total_weight>0,]
    #draw_plots
    if(nrow(rmats_to_plot)%%2 == 0){
      rmats_to_plot$flip <- rep(c(-1,1),nrow(rmats_to_plot)/2)
    } else {
      rmats_to_plot$flip <- c(rep(c(-1,1),(nrow(rmats_to_plot)-1)/2),-1)}
    
    
    rmats_to_plot$size <- rmats_to_plot[,'end']- rmats_to_plot[,'start']
    rmats_to_plot <- arrange(rmats_to_plot,size)
    rmats_to_plot$height <- NA
    n=5
    rmats_to_plot$height[rmats_to_plot$flip==-1] <- seq(0,n*length(rmats_to_plot$height[rmats_to_plot$flip==-1]),n)[-1]
    rmats_to_plot$height[rmats_to_plot$flip==1] <- seq(0,-n*length(rmats_to_plot$height[rmats_to_plot$flip==1]),-n)[-1]
    rmats_to_plot$mp <- rowMeans(rmats_to_plot[,c('start','end')])
    ###plotting parameters
    pheight=max(rmats_to_plot$height)*2+10# plot height
    #pwidth=400# plotwidth, gonna have to play arounf with this a lil
    #pwidth <- max(gtf.tx$end)+20
    pwidth=2*pheight
    #if (pwidth<400) pwidth <- 400
    rheight=10/170*pheight/4 
    
    draw_plots <- function(rmats_to_plot,plot,col,pheight,pwidth,rheight){
        for(i in 1:nrow(rmats_to_plot)){
            weight=rmats_to_plot[i,'total_weight']*15
            start=rmats_to_plot[i,'start']
            end=rmats_to_plot[i,'end']
            flip=rmats_to_plot[i,'flip']
            height=rmats_to_plot[i,'height']
            k <- rheight/2
            if (height<0) k <- -rheight/2
            mp=mean(c(start,end))
            plot <- plot+stat_function(size=weight,color=col ,xlim = c(start,end),fun = eq,args = list(cons=solver(p1=c(mp,height),p2=c(start,k),p3=c(end,k))))
        }
        return(plot)  
    }
    
    base_plot <- ggplot(data = rmats_to_plot, aes(x=start)) + xlim(min(gtf.tx$start), max(gtf.tx$end))+ ylim(-pheight/2,pheight/2)
    plot <- draw_plots(filter(rmats_to_plot,in_gtf==F), base_plot,'red',pheight,pwidth,rheight)%>%
        draw_plots(filter(rmats_to_plot,in_gtf==T), . ,'black',pheight,pwidth,rheight)
    
        #rm(start,end,height,mp)
    plot <- plot + 
      geom_label(data = rmats_to_plot,aes(x=mp,y=height,label=round(counts,3)))+
      geom_rect(xmin=min(gtf.tx$start),xmax=max(gtf.tx$end),ymin=-rheight/2,ymax=rheight/2,color='white',fill='white')+
      geom_hline(yintercept=0,alpha=.9, size=1)+ggtitle(paste0(gene,'-',tx_name,' Splicing Graph'))+
    geom_rect(data = novel_exons,aes(xmin=exonStart_0base,xmax=exonEnd,ymin=-rheight/2,ymax=rheight/2),color='red',fill='red')+
        geom_rect(data = gtf.tocomp,aes(xmin=start,xmax=end,ymin=-rheight/2,ymax=rheight/2),color='green',fill='green')+
         geom_rect(data = gtf.tx,aes(xmin=start,xmax=end,ymin=-rheight/2,ymax=rheight/2),color='black',fill='black')+theme_void()+
     labs(subtitle='Retina.Adult.Tissue')+theme(plot.title=element_text(size=pwidth/10),plot.subtitle = element_text(size = pwidth/10 -10))
    
    print(plot)

    ggsave(paste0( gene,'_',tx_name,'.png'),height = pheight, width = pwidth,units = 'mm',limitsize = F )
  }
}

splice_graph('SLC37A4',gtf_exons,all_events)
splice_graph('ABCA4',gtf_exons,all_events)
test_genes <- c('SLC37A4','ATG4D','ANXA2','ATE1','CFAP36','CLN3','GCN1','IFT52','RCBTB2','SHANK3')
for(l in test_genes) splice_graph(l,gtf_exons,all_events)


# for(i in 1:nrow(rmats_to_plot)){
#   weight=rmats_to_plot[i,5]*15
#   start=rmats_to_plot[i,1]
#   end=rmats_to_plot[i,2]
#   flip=rmats_to_plot[i,6]
#   height=rmats_to_plot[i,8]
#   plot <- plot+stat_function(size=weight,xlim = c(start,end),fun = test_fun,args = list(fstart=start,fend=end, d=flip,h=height))
# 
# }