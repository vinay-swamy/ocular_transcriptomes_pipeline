library(tidyverse)
library(data.table)
library(matrixStats)
library(RBedtools)
library(DBI)
library(parallel)
library(argparse)
library(yaml)

parser <- ArgumentParser()
parser$add_argument('--workingDir', action = 'store', dest = 'working_dir')
parser$add_argument('--ddStem', action = 'store', dest = 'dd_stem')
parser$add_argument('--filesYaml', action = 'store', dest = 'files_yaml')

#######
# working_dir <- '/data/swamyvs/ocular_transcriptomes_pipeline/'
#files_yaml <- '/data/swamyvs/ocular_transcriptomes_pipeline/files.yaml'
# dd_stem <- 'data/gtfs/final_tissue_gtfs/REPLACE.detdf'
#######
list2env(parser$parse_args(), .GlobalEnv)
source('~/scripts/write_gtf.R')

files <- read_yaml(files_yaml)

setwd(working_dir)

process_det_files <- function(det_file, tissue){
    df <- read_tsv(det_file) %>% select(transcript_id, everything())
    num_det <- df %>% 
      mutate(!!tissue:= rowSums(.[,-(1:2)])) %>% 
      select(transcript_id, !!tissue)
      
    return(num_det)
}

calc_isoform_percentage <- function(t_tissue){
    df <- counts_by_tissue %>% select(transcript_id, gene_name, new_tx_id, !!t_tissue)
    tt_sym <- as.symbol(t_tissue)
    df_gene_sums <- df %>% 
        select(-transcript_id) %>% 
        group_by(gene_name) %>%  
        summarise(gene_sums:= sum(!!tt_sym)) %>% 
        left_join(df, .) %>% 
        mutate(piu = .[[t_tissue]] / .[['gene_sums']], !!t_tissue :=NULL ) %>% 
        select(transcript_id, gene_name, new_tx_id, !!t_tissue:=piu)
    return(df_gene_sums)
    
}

replace_nan <- function(df) {
    df_fixed <- lapply(colnames(df),function(col) pull(df, col) %>%  
                           {replace(., is.nan(.), 0)}) %>% bind_cols %>% as_tibble
    colnames(df_fixed) <- colnames(df)
    return(df_fixed)
    
}

message('tidying data')
gtf <- rtracklayer::readGFF(files$anno_all_tissue_gtf)
convtab <- fread(files$all2tissue_convtab) %>% as_tibble
dntx2enst <- convtab %>% filter(class_code == '=') %>% select(transcript_id, new_tx_id = refid)
gtf <- gtf %>% left_join(dntx2enst) %>% mutate(new_tx_id= replace(new_tx_id, is.na(new_tx_id), transcript_id[is.na(new_tx_id)]))

### despite my best efforts, its possiblt for multiple DNTX transcripts to match a single reference tx
double_tx <- gtf %>% 
    filter(type == 'exon') %>% 
    group_by(new_tx_id) %>% 
    summarise(bad=sum(exon_number == 1)) %>% 
    filter(bad > 1) %>% pull(new_tx_id) %>% 
    {filter(gtf, new_tx_id %in% ., type == 'transcript')} %>% 
    select(seqid, strand, start, end, transcript_id, new_tx_id)
ref_gtf <- rtracklayer::readGFF(files$ref_GTF) 
ref_double_tx <- ref_gtf %>%  filter(transcript_id %in% double_tx$new_tx_id, type == 'transcript') %>% 
    select(seqid, strand, start, end, new_tx_id=transcript_id)
matching_tx <- inner_join(double_tx, ref_double_tx) 

# replace transcript_id that are either 
## replace non matching doubles in double that have one transcript match the reference, 
## replace both transcripts in doubles that have neither match the reference 
rep_bool <- ((gtf$new_tx_id %in% matching_tx$new_tx_id) & (!gtf$transcript_id %in% matching_tx$transcript_id) ) |
    ((gtf$transcript_id %in% double_tx$transcript_id ) & (!gtf$new_tx_id %in% matching_tx$new_id) )

gtf$new_tx_id[rep_bool] <- gtf$transcript_id[rep_bool]


tx2strand <- gtf %>% select(transcript_id, strand) %>% distinct

sample_table <- read_tsv(files$sample_table_file) %>% filter(subtissue != 'synth')
subtissues <- unique(sample_table$subtissue)

t2g <- gtf %>% filter(type == 'transcript') %>% select(transcript_id, gene_name,new_tx_id) %>% distinct

unique_exons <- gtf %>% filter(type == 'exon') %>%
    select(seqid, strand, start, end) %>% distinct %>% 
    mutate(exon_id = paste0('exon_', 1:nrow(.)))

uexon_bed <- unique_exons %>% mutate(score = 999) %>% select(seqid,start, end, exon_id, score, strand) %>% 
    from_data_frame %>% RBedtools('sort', i=.)
test_bed <- unique_exons %>% sample_n(1000) %>% mutate(score = 999) %>% select(seqid,start, end, exon_id, score, strand)
write_tsv(test_bed, 'testing/gtf_bed_for_testing.bed', col_names = F)

exon_pp_bed <- RBedtools('intersect', options = '-wa -wb -sorted', output = 'stdout', a=uexon_bed, b=files$phylop_scores) %>%
        RBedtools('groupby', options ='-g 1,2,3,4,5,6 -c 11 -o mean',i= .) %>%
        to_data_frame %>%
        select(exon_id=X4, mean_phylop_score=X7)


# intersect with strand specific
exon_snp_bed <- RBedtools('intersect', options = ' -s -wa -wb -sorted ', output = 'stdout', a=uexon_bed, b=files$snps_locs) %>%
    RBedtools('groupby', options ='-g 1,2,3,4,5,6 -c 10 -o collapse ',i= .) %>%
    to_data_frame %>%
    select(exon_id=X4, snps=X7)
save(exon_pp_bed, exon_snp_bed, file = 'testing/exons_snps_phylop.Rdata')
#load('testing/exons_snps_phylop.Rdata')
all_det <- lapply(subtissues, function(tis) 
                    gsub(pattern = 'REPLACE', replacement = tis, x = dd_stem) %>% {process_det_files(.,tis)}) %>% 
    reduce(full_join) %>% inner_join(t2g, .)
all_det[is.na(all_det)] <- 0    

message('calculating PIU')
samples_per_tissue <-  sample_table %>% group_by(subtissue) %>% summarise(count=n()) %>% filter(subtissue %in% subtissues)# REMOVE ME
samples_per_tissue_list <- 1/samples_per_tissue$count
names(samples_per_tissue_list) <- samples_per_tissue$subtissue

frac_samp_det <- as.matrix(all_det[,names(samples_per_tissue_list)]) %*% diag(samples_per_tissue_list) 
colnames(frac_samp_det) <- names(samples_per_tissue_list)
frac_samp_det <- frac_samp_det %>% as_tibble() %>% bind_cols(all_det[,c('transcript_id','gene_name', 'new_tx_id')], .)



load(files$all_tissue_quant)
all_quant[is.na(all_quant)] <- 0
message('The following samples are missing from the quant file:')
print(filter(sample_table, !sample %in% colnames(all_quant)))
sample_table <- filter(sample_table, sample %in% colnames(all_quant))
counts_by_tissue <- lapply(subtissues,
                               function(tis) filter(sample_table, subtissue == tis) %>% pull(sample) %>%
                                   {all_quant[,c('transcript_id', .)]} %>%
                                   mutate(!!tis := rowMeans(.[,-1] %>% as.matrix)) %>%
                                   select(transcript_id, !!tis)
                           ) %>%
    reduce(left_join) %>% inner_join(t2g, .)

piu_raw <- lapply(colnames(counts_by_tissue)[-(1:3)], calc_isoform_percentage) %>% reduce(left_join)
piu <-replace_nan(piu_raw)



# fix transcript id column
gtf <- gtf %>% select(-transcript_id) %>% rename(transcript_id = new_tx_id) %>% left_join(unique_exons) %>% 
    mutate(seqid=as.character(seqid))
frac_samp_det <- frac_samp_det %>% select(-transcript_id) %>% rename(transcript_id = new_tx_id) 
all_det <-  all_det %>% select(-transcript_id) %>% rename(transcript_id = new_tx_id) 
piu <- piu %>% select(-transcript_id) %>% rename(transcript_id = new_tx_id) 
convtab <- convtab %>% inner_join(t2g, .) %>% select(-transcript_id) %>% rename(transcript_id = new_tx_id) 
tissue_det <- convtab %>% select(subtissues) %>% apply(2, function(x) !is.na(x)) %>% as_tibble %>% bind_cols(convtab[,1:3],.) 
cds_df <- gtf %>%  
    mutate(seqid=as.character(seqid)) %>% filter(type == 'CDS') %>% select(-is.last,-is.singleExon, -is.first, -contains('novel'))
exon_info_df <- left_join(exon_pp_bed, exon_snp_bed) 
## fix gtf gene_name column

gtf_t2g <- gtf %>% filter(type == 'transcript') %>% select(transcript_id, gene_name) %>% distinct %>% 
  mutate(gene_name = replace(gene_name , is.na(gene_name), transcript_id[is.na(gene_name)]))
gtf <- gtf %>% select(-gene_name) %>% inner_join(gtf_t2g)



#---- precalculating data for shiny app
# I wanted to add thick lines to indicate where the CDS of transcripts is, like commonly seen in genome browsers

## First for each transcript I identify which exon the CDS start is on, and which exon the CDS end is on 
cds_se <- cds_df %>% filter(type == 'CDS') %>% 
    group_by(transcript_id) %>% 
    summarise(seqid=first(seqid), strand=first(strand), cds_start=min(start), cds_end=max(end), gene_name=first(gene_name))

cds_start_bed <- cds_se %>% mutate(end=cds_start + 1, score= 888) %>% 
    select(seqid, cds_start, end, transcript_id, score, strand) %>% 
    from_data_frame %>% 
    RBedtools('sort', i=.)
cds_end_bed <- cds_se %>% mutate(end=cds_end + 1, score= 999) %>% 
    select(seqid, cds_end, end, transcript_id, score, strand) %>% 
    from_data_frame %>% 
    RBedtools('sort', i=.)
exon_df <- gtf %>% filter(type == 'exon') %>% mutate(score= 111) %>%
    select(seqid, start, end, transcript_id, score, strand)
exon_bed <- exon_df %>% from_data_frame %>% 
    RBedtools('sort', i=.)

## an cds start/end can either overlap and exon, or land exactly on start/end of an exon 

which_exon_overlaps_start <- RBedtools('intersect', options = '-s -wa -wb', a=exon_bed, b=cds_start_bed) %>% 
    to_data_frame %>% 
    filter(X4 == X10) %>% 
    select(seqid=X1, start=X2, end=X3, transcript_id=X4, strand=X6) %>% 
    mutate(is_CDS_start=T)
which_exon_overlaps_start <- cds_se %>% 
    select(seqid, strand, start=cds_start, transcript_id) %>% 
    inner_join(exon_df) %>% 
    select(seqid, start, end, transcript_id, strand) %>% 
    mutate(is_CDS_start=T) %>% 
    bind_rows(which_exon_overlaps_start)
which_exon_overlaps_start <- cds_se %>% 
    select(seqid, strand, end=cds_start, transcript_id) %>% 
    inner_join(exon_df) %>% 
    select(seqid, start, end, transcript_id, strand) %>% 
    mutate(is_CDS_start=T) %>% 
    bind_rows(which_exon_overlaps_start) %>% distinct 

which_exon_overlaps_end <- RBedtools('intersect', options = '-s -wa -wb', a=exon_bed, b=cds_end_bed) %>% 
    to_data_frame %>% 
    filter(X4 == X10) %>% 
    select(seqid=X1, start=X2, end=X3, transcript_id=X4, strand=X6) %>% 
    mutate(is_CDS_end=T)
which_exon_overlaps_end <- cds_se %>% 
    select(seqid, strand, end=cds_end, transcript_id) %>% 
    inner_join(exon_df) %>% 
    select(seqid, start, end, transcript_id, strand) %>% 
    mutate(is_CDS_end=T) %>% 
    bind_rows(which_exon_overlaps_end)
which_exon_overlaps_end <- cds_se %>% 
    select(seqid, strand, start=cds_end, transcript_id) %>% 
    inner_join(exon_df) %>% 
    select(seqid, start, end, transcript_id, strand) %>% 
    mutate(is_CDS_end=T) %>% 
    bind_rows(which_exon_overlaps_end) %>% 
    distinct



gtf <- gtf %>% left_join(which_exon_overlaps_start) %>% left_join(which_exon_overlaps_end) %>% 
    mutate(is_CDS_start=replace_na(is_CDS_start, F), is_CDS_end=replace_na(is_CDS_end, F))
## Next, this function scales the genomic lengthsload for each gene, and adds in thick starts for protein coding transcripts 
save.image('testing/pre_plotting_gtf_image.Rdata')
message('running app data precalculation')
# NOTE:
#   The some cases of plotting gtf will produce extra exons, specifcallt when the CDS starts/ends in the middle of an exon. 
#   This is the right behavior, and exon numbers are preserved
#   Reworked 08/04/2020 - add CDS, then take use squarte distance. Intron scaling is now handled by the app
make_plotting_gtf_by_gene <- function(gtf, cds_df, cds_se,which_exon_overlaps_end, which_exon_overlaps_start,  gene){
  
  gtf_gene <- filter(gtf, gene_name == gene)
  cds_se_gene <- filter(cds_se, gene_name == gene)
  cds_df_gene <- filter(cds_df, gene_name == gene)
  if(nrow(cds_df) == 0) {
    return()
  }
  ### this adds the thick start for each transcript. For most genes its easy, but there are many corner cases
  merge_CDS_scaled_exons <- function(gtf_gene, which_exon_overlaps_start, which_exon_overlaps_end,
                                     cds_df_gene, cds_se_gene, t_tx){
    #print(t_tx) 
    # in order to create thick lines for CDS, we have to first get the CDS locations, and when the start and end 
    # are in the middle of an exon, split that exon into 2 features.
    scaled_exons_tx <- filter(gtf_gene, type == 'exon', transcript_id == t_tx) %>% 
      mutate(Xmax=0, Xmin=0, Ymin=-.5, Ymax=.5)
    if(any(duplicated(scaled_exons_tx$exon_number)) ){ #in a single transcript no transcript should have the same start or end 
      
      write(t_tx, file = 'data/shiny_data/debug/multi_tx_under_same_id.txt', append = T, sep = '\n')# multiple transcripts under same id
      return(tibble())
    }
    cds_df_tx <- filter(cds_df_gene, transcript_id ==t_tx ) %>% select(seqid, strand, start, end)
    cds_se_tx <- filter(cds_se_gene, transcript_id  == t_tx)
    which_exon_CDS_start <- which_exon_overlaps_start %>% filter(transcript_id == t_tx)
    which_exon_CDS_end <-   which_exon_overlaps_end %>% filter(transcript_id == t_tx)
    # this indicates the transcript id was not found in the CDS df, and so must not be protein coding
    if(nrow(cds_se_tx) == 0){
      #write(t_tx, file = 'data/shiny_data/debug/no_cds_bu_marked_as_pc.txt', append = T, sep = '\n')
      return(scaled_exons_tx)
    }
    # this means we failed to identify which exon the CDS start /end is on int the previous step
    if(nrow(which_exon_CDS_end) == 0 | nrow(which_exon_CDS_start) == 0){
      write(t_tx, file = 'data/shiny_data/debug/failed_to_find_CDS_start_or_end.txt', append = T, sep = '\n')
      return(tibble())
    }
    
    # if the cds starts and ends on the same exon, which is always the case for single exon transcripts 
    if(nrow(anti_join(which_exon_CDS_start, which_exon_CDS_end)) == 0 ){# if the cds starts and ends on the same exon
      special_exon <- scaled_exons_tx %>% inner_join(which_exon_CDS_start)
      if(nrow(cds_se_tx) == 0){
        write(t_tx, file = 'data/shiny_data/debug/special_is_funky.txt', append = T, sep = '\n')
        return(tibble())
      }
      # the CDS starts/ ends on the exact same location as the exon its on
      if(special_exon$start == cds_se_tx$cds_start | special_exon$end == cds_se_tx$cds_end){
        # single exon and CDS are the exact same 
        if(special_exon$start == cds_se_tx$cds_start & special_exon$end == cds_se_tx$cds_end){
          complete_df <- special_exon %>% mutate(Ymin=-1, Ymax =1)  
        } else if(special_exon$start == cds_se_tx$cds_start){ # both of these are of either the CDS start /end are on the location
          #print('here')
          cds_mid <- special_exon %>% mutate(start=cds_df_tx$start, end=cds_df_tx$end, 
                                             Xmax= Xmin + sqrt(end-start), 
                                             Ymax=1, Ymin=-1)
          
          split_end_nc_start <- special_exon  %>% mutate(start = cds_se_tx$cds_end+ 1, 
                                                         Xmin=cds_mid$Xmax, 
                                                         Xmax=Xmin + sqrt(end-start))
          complete_df <- bind_rows( cds_mid, split_end_nc_start)
          
        }else{#special_exon$end == cds_se_tx$cds_end
          offset <- (cds_se_tx$cds_start - special_exon$start - 1 )
          split_start_nc_end <- special_exon %>% mutate(end= start + offset, Xmax= Xmin + sqrt(end-start))
          cds_mid <- special_exon %>% mutate(start=cds_df_tx$start, end=cds_df_tx$end, 
                                             Xmin=split_start_nc_end$Xmax, Xmax= Xmin + sqrt(end-start), 
                                             Ymax=1, Ymin=-1)
          complete_df <- bind_rows(split_start_nc_end, cds_mid)
        }
      }else{
        offset <- (cds_se_tx$cds_start - special_exon$start - 1 )
        split_start_nc_end <- special_exon %>% mutate(end= start + offset, Xmax= Xmin + sqrt(offset))
        cds_mid <- special_exon %>% mutate(start=cds_df_tx$start, end=cds_df_tx$end, 
                                           Xmin=split_start_nc_end$Xmax, Xmax= Xmin + sqrt(end-start), 
                                           Ymax=1, Ymin=-1)
        offset <- special_exon$end - cds_se_tx$cds_end
        split_end_nc_start <- special_exon  %>% mutate(start = cds_se_tx$cds_end+ 1, 
                                                       Xmin=cds_mid$Xmax, 
                                                       Xmax=Xmin + sqrt(offset))
        
        complete_df <- bind_rows(split_start_nc_end, cds_mid, split_end_nc_start) %>% arrange(start)
      }
      return(complete_df)
    }
    
    both <- bind_rows(which_exon_CDS_start, which_exon_CDS_end)
    both[is.na(both)] <- F
    exact_match <-  inner_join(scaled_exons_tx, cds_df_tx) %>% mutate(Ymin=-1, Ymax=1)
    non_match <- anti_join(scaled_exons_tx, cds_df_tx) %>%  inner_join(both) %>% arrange(start)
    not_cds <- anti_join(scaled_exons_tx, cds_df_tx) %>%  anti_join(both)
    
    # there are three cases for the spliting start and end problem for multi exon CDSs
    if(nrow(non_match) == 2){
      # both the start and the end land in the middle of an exon, so we have to adjust both
      offset <- (cds_se_tx$cds_start - non_match$start[1] -1)
      split_start_nc_end <- non_match[1,] %>% mutate(end= start + offset, Xmax= Xmin + sqrt(offset))
      split_start_cds_start <- non_match[1,] %>% mutate(start= start +offset + 1, Xmin = Xmin + sqrt(offset), Ymin=-1, Ymax=1)
      merged_starts <- bind_rows(split_start_nc_end, split_start_cds_start)
      
      
      offset <- non_match$end[2] - cds_se_tx$cds_end
      split_end_cds_end <- non_match[2,] %>% mutate(end = start + offset, Xmax = Xmin + sqrt(offset), Ymin=-1, Ymax=1)
      split_end_nc_start <- non_match[2,] %>% mutate(start = start + offset + 1, Xmin=Xmin + sqrt(offset))
      merged_ends <- bind_rows(split_end_cds_end, split_end_nc_start)
      
      complete_df <- bind_rows(not_cds, exact_match, merged_starts, merged_ends) %>% arrange(start)
    } else if(nrow(non_match) == 1){
      #either the start or the end matches perfectly with an exon, and the other is in the middle 
      start_in_middle <- non_match %>% inner_join(which_exon_CDS_start)
      end_in_middle <- non_match %>% inner_join(which_exon_CDS_end)
      if(nrow(start_in_middle)!= 0){
        offset <- (cds_se_tx$cds_start - non_match$start[1] -1)
        split_start_nc_end <- non_match[1,] %>% mutate(end= start + offset, Xmax= Xmin + sqrt(offset))
        split_start_cds_start <- non_match[1,] %>% mutate(start= start +offset + 1, Xmin = Xmin + sqrt(offset), Ymin=-1, Ymax=1)
        merged <- bind_rows(split_start_nc_end, split_start_cds_start)
        
      }else{
        offset <- non_match$end[1] - cds_se_tx$cds_end
        split_end_cds_end <- non_match[1,] %>% mutate(end = start + offset, Xmax = Xmin + sqrt(offset), Ymin=-1, Ymax=1)
        split_end_nc_start <- non_match[1,] %>% mutate(start = start + offset + 1, Xmin=Xmin + sqrt(offset))
        merged <- bind_rows(split_end_cds_end, split_end_nc_start)
        
      }
      complete_df <- bind_rows(not_cds, exact_match, merged)
      
    } else if(nrow(non_match) == 0){
      #both the start and the end exact match the exons, so nothing to do but bind them all together 
      complete_df <- bind_rows(not_cds, exact_match)
    } else {# nothing should ever get here 
      write(gene,file = 'data/shiny_data/debug/bad_genes.txt', sep = '\n', append = T)
      return(tibble())
    }
    
    return(complete_df)
  }
  PC_tx <- unique(gtf_gene$transcript_id)
  res <- lapply(PC_tx, function(tx) merge_CDS_scaled_exons(gtf_gene = gtf_gene, 
                                                           #scaled_exons = scaled_exons, 
                                                           which_exon_overlaps_start = which_exon_overlaps_start, 
                                                           which_exon_overlaps_end = which_exon_overlaps_end,
                                                           cds_df_gene = cds_df_gene, 
                                                           cds_se_gene = cds_se_gene, 
                                                           t_tx = tx)) %>% bind_rows %>% select(-Xmax, -Xmin)
  
  t_gtf_exons <- filter(res, type == 'exon') %>% select(seqid,strand, start, end) %>% 
    mutate(length = sqrt(end-start))
  gap=mean(t_gtf_exons$length)/10
  scale_exon_lengths <- function(df, scale_factor, ymin, ymax){
    df <- t_gtf_exons
    scale_factor = gap 
    n_gtf_exons <- df %>% 
      select(seqid, strand, start, end) %>% 
      distinct %>% 
      arrange(start)
    n_gtf_exons_scaled <- n_gtf_exons %>% mutate(Xmin=start-min(start), Xmax=end-min(start))
    nge_scale_sqrt <- n_gtf_exons_scaled %>% mutate(Xmin=sqrt(Xmin), Xmax=sqrt(Xmax), length = Xmax-Xmin)
    
    
    return(nge_scale_sqrt)
  }
  scaled_exons <- scale_exon_lengths(t_gtf_exons, gap,-.5,.5 )
  final_gtf <- inner_join(res, scaled_exons)
  
  return(final_gtf)
}


set.seed(10021)
all_genes <- gtf %>%  pull(gene_name) %>% unique()
plotting_gtf <- mclapply(all_genes, function(gene) 
      make_plotting_gtf_by_gene(gtf = gtf, cds_df = cds_df, cds_se = cds_se, 
                                 which_exon_overlaps_end = which_exon_overlaps_end,
                                 which_exon_overlaps_start = which_exon_overlaps_start,
                                 gene = gene
                                 ),
      mc.cores = 12) #%>% bind_rows 
save.image('testing/plotting_gtf_list.Rdata')
plotting_gtf <- bind_rows(plotting_gtf)

#add transcript0_ids back in, sort and add exon tooltips
plotting_gtf <- plotting_gtf %>% 
    bind_rows(gtf %>% filter(type == 'transcript', transcript_id %in% plotting_gtf$transcript_id )) %>% 
    arrange(seqid, start, transcript_id, type) %>% 
    left_join(exon_info_df) %>%  
    mutate(ttip=paste0('genomic start: ', start, '\n', 'genomic end: ', end,'\n',
                         'mean phylop score: ', mean_phylop_score, '\n', 'snps: ', snps), 
           ttip=gsub('snps: NA', '', ttip), 
           ttip=gsub('mean phylop score: NA\n', '', ttip))
 
message('plotting gtf successfully made. Generating DB')
#save(plotting_gtf, file = 'data/shiny_data/plotting_GTF.Rdata')
#----
save.image('testing/all_prepped_shiny_data.Rdata')
con <- dbConnect(RSQLite::SQLite(), files$shiny_db)
dbWriteTable(con, 'gtf', gtf)
dbSendQuery(con, 'CREATE INDEX gtf_ind ON gtf (gene_name, transcript_id)')
dbWriteTable(con, 'frac_samp_det', frac_samp_det)
dbSendQuery(con, 'CREATE INDEX frac_samp_det_ind ON frac_samp_det (gene_name, transcript_id)')
dbWriteTable(con, 'all_det', all_det)
dbSendQuery(con, 'CREATE INDEX all_det_ind ON all_det (gene_name, transcript_id)')
dbWriteTable(con, 'piu', piu)
dbSendQuery(con, 'CREATE INDEX piu_ind ON piu (gene_name, transcript_id)')
dbWriteTable(con, 'tc2m', convtab)
dbSendQuery(con, 'CREATE INDEX tc2m_ind ON tc2m (gene_name, transcript_id)')
dbWriteTable(con, 'tissue_det', tissue_det)
dbSendQuery(con, 'CREATE INDEX tissue_det_ind ON tissue_det (gene_name, transcript_id)')
dbWriteTable(con, 'exon_info_df', exon_info_df)
dbSendQuery(con, 'CREATE INDEX exon_info_df_ind ON exon_info_df (exon_id)')
dbWriteTable(con, 'cds_df', cds_df)
dbSendQuery(con, 'CREATE INDEX cds_df_ind ON cds_df (gene_name, transcript_id)')
dbWriteTable(con, 'plotting_gtf', plotting_gtf)
dbSendQuery(con, 'CREATE INDEX plotting_gtf_ind ON plotting_gtf (gene_name, transcript_id)')
dbDisconnect(con)
all_gene_names <- unique(plotting_gtf$gene_name)
all_tissues <- subtissues 
#save(gtf, all_det, frac_samp_det, piu, tc2m,tissue_det, file=outfile)

save(all_gene_names, all_tissues, file = files$shiny_rdata)

#now make data for downloading from shiny app
#----
##panbody
write_gtf3(gtf, paste0(files$shiny_dl_dir, '/panbody.gtf'))
##paneye
eye_subtissues <- c('Retina_Adult.Tissue', 'Retina_Fetal.Tissue', 
                 'RPE_Adult.Tissue', 'RPE_Fetal.Tissue',
                 'Cornea_Adult.Tissue', 'Cornea_Fetal.Tissue')
eye_tissues <- c('Retina', 'Cornea', 'RPE')
get_tx_in_tissues <- function(tissues){
    convtab %>% 
        select(transcript_id, !!tissues) %>%
        mutate(det=rowSums(.[,-1] %>% apply(2, function(x) !is.na(x)))) %>%
        filter(det >0) %>%
        pull(transcript_id)
}

get_tx_in_one_tissue <- function(tissues){
    convtab %>% 
        select(transcript_id, !!tissues) %>%
        filter(!is.na(.[,2])) %>%
        pull(transcript_id)
}

paneye_gtf <- filter(gtf, transcript_id %in% get_tx_in_tissues(eye_subtissues))
write_gtf3(paneye_gtf, paste0(files$shiny_dl_dir, '/paneye.gtf'))

lapply(eye_subtissues, function(tis) 
        {t_gtf <- filter(gtf, transcript_id %in% get_tx_in_one_tissue(tis)) 
        dim(t_gtf) 
        write_gtf3(t_gtf, paste0(files$shiny_dl_dir,'/', tis, '.gtf' ))
        }
    )


lapply(eye_tissues, function(tis)
        {
    cols <- convtab %>% select(contains(tis)) %>% colnames
    ad_col <- which(grepl('Adult', cols))
    fet_col <- which(grepl('Fetal', cols))
    t_gtf <- filter(gtf,transcript_id %in% get_tx_in_tissues(cols) )
        write_gtf3(t_gtf, paste0(files$shiny_dl_dir,'/',cols[ad_col], '-', cols[fet_col], '.gtf' ))
                
        }
    )


