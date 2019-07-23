library(tidyverse)
library(parallel)
cores=15
notes <- '
"I think tab-delimited files are a minor evil in the world." - a dumbass
-all txs only have 1 pc form so that
'


domain_file='results/hmmer/domain_hits.tsv'
gtf_file <- 'results/all_tissues.combined.gtf'

cn=c('target_name', 'target_accession' ,'tlen' ,'query_name', 'accession' ,'qlen' ,'evalue', 'total_score', 'total_bias', 'dom_num' , 
     'total_doms' , 'c_evalue' , 'i_evalue' , 'domain_score' , 'domain_bias' , 'start_hmm_coord' ,'end_hmmcoord','start_seq', 
     'end_seq' , 'from_env' ,'to_env',  'post_prob', 'description')
domain_hits <- read_delim(domain_file,skip = 3, delim = ' ' , col_names = cn ) %>%
    mutate(query_name=query_name %>% strsplit('\\.') %>% sapply(function(x) x[1]))
num_cols <- cn[-c(1,2,4,5,23)]
for(col in num_cols) domain_hits[,col] <- pull(domain_hits, col) %>% as.numeric
f_domain_hits <- filter(domain_hits, evalue < 1e-5 )

gtf <- rtracklayer::readGFF(gtf_file)
tx_tab <- filter(gtf, type=='transcript', grepl('MSTRG',oId)) %>% select(transcript_id, oId, cmp_ref)
ref <- pull(tx_tab, cmp_ref) %>% unique %>% {filter(gtf, oId%in%., type=='transcript')} %>% 
    select(cmp_ref, transcript_id) %>% rename(ref_tx_id=transcript_id)
tx_tab <- left_join(tx_tab, ref, by='cmp_ref') %>% filter(!is.na(cmp_ref)) %>% filter(transcript_id%in%f_domain_hits$query_name)



n_tx <- tx_tab$transcript_id[1]
r_tx <- tx_tab$ref_tx_id[1]
novel<- filter(f_domain_hits, query_name==n_tx)
ref <- filter(f_domain_hits, query_name==r_tx)
check_diff_domains <- function(pair, dom){
    subj=pair[1,1]
    qu=pair[1,2]
    n_s <- filter(dom, query_name==subj)
    n_q <- filter(dom, query_name==qu)
    if(nrow(n_s)==0) return(0)# novel tx not found in dom file
    if(nrow(n_s)==nrow(n_q)){ if(n_s$description == n_q$description) return(2); return(0)}
    # if both tx's have the same number of domain but different domains, we want ot know that
    return(1)
}
tx_tab_list <- split(tx_tab%>%select(transcript_id, ref_tx_id,), 1:nrow(tx_tab))
diff_dom_f <- mclapply(tx_tab_list,function(x) check_diff_domains(x, f_domain_hits),mc.cores = cores)
sum(grepl('MSTRG', gtf$oId))
c_domain_hits <- filter(f_domain_hits, i_evalue<=.05)
diff_dom_c <- mclapply(tx_tab_list,function(x) check_diff_domains(x, c_domain_hits),mc.cores = cores)
diff_dom_c <- diff_dom_c %>% unlist
sum(diff_dom_c==1)
sum(diff_dom_c==0)
sum(diff_dom_f==2)

tx_tab_diff_dom <- tx_tab[!diff_dom,]
diff_dom_hits <- filter(f_domain_hits, query_name%in%tx_tab_diff_dom$transcript_id)
diff_dom_hits <- filter(diff_dom_hits, i_evalue<=.05)



