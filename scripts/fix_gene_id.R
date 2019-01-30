library(tidyverse)
args=commandArgs(trailingOnly = T)
infile=args[1]
outfile=args[2]
write_gtf3 <- function(df, file){
    to_write <- lapply(colnames(df[9:ncol(df)]), function(x) pull(df,x) %>% paste0(' "', . ,'"') %>% paste(x,., sep = '')) %>%
        do.call(cbind,.) %>% split(.,1:nrow(df)) %>%lapply( function(x) grep('"NA"',x,invert = T) %>% x[.]) %>%
        sapply( function(x) paste0(x,collapse = ';'))
    res <- cbind(df[,1:8],to_write,stringsAsFactors=F)
    res[is.na(res)] <- '.'
    write('##gtf-version-3',file = file)
    write.table(res,file, append = T, quote=F , sep = '\t',row.names = F , col.names = F, qmethod = 'escape' )
}
in_gtf <- rtracklayer::readGFF(infile) %>% mutate(gffc_loc=gene_id, gene_id=transcript_id)
write_gtf3(in_gtf,outfile)



