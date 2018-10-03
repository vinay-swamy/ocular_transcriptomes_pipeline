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

setwd('~/NIH/autoRNAseq/spleen_vs_skin')
#load('plotting_curvs.Rdata')
#assuming '+' strand
#gdata::keep(gtf,sure=T)


gtf <- readGFF('../ref/gencodeAno_bsc.gtf')
all_events=read.table('SE.MATS.JC.txt',stringsAsFactors = F,header = T)
#only protein coding tx's
gtf_exons<- filter(gtf, transcript_type=='protein_coding')%>%filter(type=='exon')

splice_graph<- function(gene,gtf,all_events){
  genes <- filter(all_events, geneSymbol==gene)
  gtf.gene <- filter(gtf_exons,gene_name==gene,type=='exon')%>%arrange(start)#%>%select(1:7)
  gtf.gene$start=gtf.gene$start-1
  #exons ends are ok, exon starts are off by 1
  tx_ids <- unique(gtf.gene$transcript_id)
  #i=tx_ids[1]
  #plot every tx in a gene
  for(i in tx_ids){
    tx_name=i
    gtf.tx <- filter(gtf.gene,transcript_id==i,type=='exon')%>%arrange(start)#%>%select(1:7)
    #pick rmats events in which all three exons match thsoe in GTF
    exon_coords <- as.matrix(gtf.tx[,4:5]) 
    rmats_exon_coords <- as.matrix(genes[,6:11])
    #there has to be a better way to do this
    
    to_plot=logical(nrow(rmats_exon_coords))
    for(j in 1:nrow(rmats_exon_coords)){
      #b's are exon start and ends in rmats data, b is true if rmats exon start & end match gtf 
      b1 <- any(apply(exon_coords,1 ,function (x) all(x==rmats_exon_coords[j,1:2])))
      b2 <- any(apply(exon_coords,1 ,function (x) all(x==rmats_exon_coords[j,3:4])))
      b3 <- any(apply(exon_coords,1 ,function (x) all(x==rmats_exon_coords[j,5:6])))
      to_plot[j] <- b1 && b2 && b3
    }
    rmats_to_plot <- genes[to_plot,]
    if(nrow(rmats_to_plot)<2) {
      print('plop')  
      return(0)
    }
    start <- gtf.tx$start[1]/100
    #scale genome coordinates
    gtf.tx[,c('start','end')] <- gtf.tx[,c('start','end')]/100 -start
    rmats_to_plot[,6:11] <- rmats_to_plot[,6:11]/100-start
    #current idea: each event contains three edges: a single long edge in the skipped isoform and 2 short edges in the included
    #isoform. RMATS gives counts for the the exclusion event, and inclusion event. So the weight of the long edge 
    #is just the relative counts for the exlcusion event. for the short edges the inclusion counts are distributed over 2 edges,
    #so split the weight over the edges. remember to account for split merging edges.
    # also incusion edges counts are going to be duplicated when merged below, so cut their counts by 2 again(so divide by 4)
    rmats_to_plot[,c("IJC_SAMPLE_1", "IJC_SAMPLE_2")] <-  rmats_to_plot[,c("IJC_SAMPLE_1", "IJC_SAMPLE_2")]/4
    
    edges <- rbind(as.matrix(rmats_to_plot[,c("upstreamEE","downstreamES", "SJC_SAMPLE_1","SJC_SAMPLE_2")]),# skipped junction
                   as.matrix(rmats_to_plot[,c("upstreamEE","exonStart_0base","IJC_SAMPLE_1", "IJC_SAMPLE_2")]),#included junction
                   as.matrix(rmats_to_plot[,c("exonEnd","downstreamES","IJC_SAMPLE_1", "IJC_SAMPLE_2")]))#included junction 
    
    edges <- as.data.frame(edges)
    rmats_to_plot <- aggregate(edges,list(edges$upstreamEE,edges$downstreamES),sum)
    total_counts <- sum(rowSums(rmats_to_plot[,5:6]))
    rmats_to_plot$total_weight= rowSums(rmats_to_plot[,5:6])/total_counts
    if(nrow(rmats_to_plot)%%2 == 0){
      rmats_to_plot$flip <- rep(c(-1,1),nrow(rmats_to_plot)/2)
    } else {
      rmats_to_plot$flip <- c(rep(c(-1,1),(nrow(rmats_to_plot)-1)/2),-1)}
    
    rmats_to_plot <-  rmats_to_plot[,c(1:2,5:8)]
    rmats_to_plot$size <- rmats_to_plot[,2]- rmats_to_plot[,1]
    rmats_to_plot <- arrange(rmats_to_plot,size)
    rmats_to_plot$height <- NA
    n=40
    rmats_to_plot$height[rmats_to_plot$flip==-1] <- seq(0,n*length(rmats_to_plot$height[rmats_to_plot$flip==-1]),n)[-1]
    rmats_to_plot$height[rmats_to_plot$flip==1] <- seq(0,-n*length(rmats_to_plot$height[rmats_to_plot$flip==1]),-n)[-1]
    colnames(rmats_to_plot) <- c('start','end','sample_1','sample_2','total_weight','flip','size','height')
    rmats_to_plot$mp <- rowMeans(rmats_to_plot[,1:2])
    ###plotting parameters
    pheight=max(rmats_to_plot$height)*2+10# plot height
    #pwidth=400# plotwidth, gonna have to play arounf with this a lil
    #pwidth <- max(gtf.tx$end)+20
    pwidth=2*pheight
    #if (pwidth<400) pwidth <- 400
    rheight=10/170*pheight  
    
    plot <- ggplot(data = rmats_to_plot, aes(x=start)) + xlim(min(gtf.tx$start), max(gtf.tx$end))+ ylim(-pheight/2,pheight/2)
    for(i in 1:nrow(rmats_to_plot)){
      weight=rmats_to_plot[i,5]*15
      start=rmats_to_plot[i,1]
      end=rmats_to_plot[i,2]
      flip=rmats_to_plot[i,6]
      height=rmats_to_plot[i,8]
      k <- rheight/2
      if (height<0) k <- -rheight/2
      mp=mean(c(start,end))
      plot <- plot+stat_function(size=weight,xlim = c(start,end),fun = eq,args = list(cons=solver(p1=c(mp,height),p2=c(start,k),p3=c(end,k))))
    }
    rm(start,end,height,mp)
    plot <- plot +  geom_hline(yintercept = 0,size=7,color='white')+ 
      geom_label(data = rmats_to_plot,aes(x=mp,y=height,label=round(total_weight,3)))+
      geom_rect(xmin=min(gtf.tx$start),xmax=max(gtf.tx$end),ymin=-rheight/2,ymax=rheight/2,color='white',fill='white')+
      geom_hline(yintercept=0,alpha=.9, size=1)+ggtitle(paste0(gene,'-',tx_name,' Splicing Graph'))+
      geom_rect(data = gtf.tx,aes(xmin=start,xmax=end,ymin=-rheight/2,ymax=rheight/2),color='black',fill='black')+theme_void()
      #geom_segment(data = gtf.tx,aes(x=start,xend=end,y=0, yend=0),size=7)+theme_minimal()
    #print(plot)
    ggsave(paste0('plots/', gene,'_',i,'.png'),height = pheight, width = pwidth,units = 'mm',limitsize = F )
  }
}


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