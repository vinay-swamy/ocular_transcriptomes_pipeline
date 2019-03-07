import os
import sys
import pandas as pd
from  dfply import * # this is required for it work properly
from Bio import SeqIO, pairwise2
args=['','/Users/vinayswamy/NIH/eyeintegration_splicing', 'results_38/all_tissues.combined.gtf', 'ref/gencodeGFF3.gff.tsv',\
      'results_b38/best_orfs.transdecoder.pep', 'ref/gencodeProtSeq.fa']
#args=sys.argv
working_dir=args[1]
st_gtf_file=args[2]
ref_gff_file=args[3]
novel_prot_fa_file=args[4]
ref_prot_fa_file=args[5]


os.chdir(working_dir)

## stole these from an old assignment
def file_len(fname):
    with open(fname) as f:
        os=0
        for i, l in enumerate(f):
            if l[0]=='#':
               os+=1
    return ((i + 1 -os , os))

def parse_attributes(st):
    att_list=st.replace('"','').replace('=', ' ').split(';')[:-1]
    out_dict=dict()
    for att in att_list:
        att=att.lstrip(' ').split(' ')
        out_dict[att[0]]=att[1]
    return(out_dict)
    

def read_GTF(file):
    file_length, os=file_len(file)
    outdf=[None]*file_length
    with open(file) as gtf:
        gtf_header=['seqid','source','type','start','stop','score','strand','frame']
        for k in range(os):
            gtf.readline().strip('\n')
        for i in range(file_length):
            line= gtf.readline().strip('\n') 
            line_list=line.split('\t')
            line_dict=dict()
            [line_dict.update({value:line_list[j]}) for j,value in enumerate(gtf_header)]
            line_dict.update(parse_attributes(line_list[-1]))
            outdf[i]=line_dict
    df=[k for k in outdf if k is not None]
    df=pd.DataFrame(outdf)
    col=df.columns.tolist()
    col=[cn for cn in col if cn not in gtf_header]
    df=df[gtf_header+col]
    return(df)


def overlap(el,loc):
   for i, tup in enumerate(el):
       if  loc >= tup[0]   and loc <= tup[1]:
           return(i) 
           
   return(-1)
           
        


def scale_genomic_length(start, stop):
    f_list=list(zip(start, stop))
    last=-1
    res=[None]* len(f_list)
    for i, tup in enumerate(f_list):
        start=last+1
        end=start+(tup[1]-tup[0])
        last=end
        res[i]=(start,end)
    return(pd.DataFrame(res, columns=['scaled_start','scaled_stop']))



def cds_tot(df):
    df=df >> mutate(length=X.stop - X.start)
    return(df.length)        

def make_gene_model(ref, novel, exon_df, cds_df, n=10):
    '''
    ref: seqrecord obj of ref pc tx
    novel seqrecord obj of novel pc tx
    n: number of continous matches to determine start
    -at some point, should throw away tx that don't match a certain threshold
    -right now the best way to match seems like
    three cases : seqs line up perfectly; ref match after novel start; ref match before
    '''
    #part onee do the matching
    match='|'* n 
    aln=pairwise2.align.globalms(ref.seq,novel.seq, 2, -1, -1, -1, penalize_end_gaps=False, one_alignment_only=True)   
    aln=pairwise2.format_alignment(*aln[0])
    r_aln_score=aln.count('|')/len(ref)
    n_aln_score=aln.count('|')/len(novel)
    score=(r_aln_score+ n_aln_score)/2
    #print('{}:{}'.format(novel.id, str(score)))
    r_seq, m_seq, n_seq= aln.split('\n')[0:3]    
    #part 2 figure out where on the genome
    start_codon=cds_df.iloc[[0]]
    if m_seq.find(match)==-1:
        #print('No good match found')
        return('NA', -3)
    # find location of the relatuve  
    ref_start=start_codon['start'].iloc[0]
    gen_loc_exon_list=list(zip(exon_df.start, exon_df.stop))
    ref_start_idx=overlap(gen_loc_exon_list, ref_start)# locartion
    ref_cds_start= exon_df.scaled_start.iloc[ref_start_idx]+ (ref_start - exon_df.start.iloc[ref_start_idx])
    if r_seq[0] != '-' and n_seq[0]!='-' :
        # if the starts are the same, then the end is the length*3 nt away    
        nvl_cds_sclst=ref_cds_start
    elif r_seq[0]=='-':
        # if the start of the novel occurs before the start of the ref, find the ref start relative to the novel
        # and subtract the distance from the genomic ref start ot get the novel ref start
        nvl_cds_sclst=ref_cds_start - m_seq.find(match) * 3    
    else:
        # if the start of the novel occurs after the start of the novel ref, find the distance to the match point
        # and then add it to the genomic start to determine the novel CDS start
        nvl_cds_sclst=ref_cds_start + m_seq.find(match) * 3

    # finding the end of the cds is just adding the length of the prot seq, regard less of
    nvl_cds_sclend= nvl_cds_sclst + len(novel)*3 + 2 # the plus 2 accounts for the stop codon;
    scl_loc_exon_list=list(zip(exon_df.scaled_start, exon_df.scaled_stop))
    
    nv_start_idx= overlap(scl_loc_exon_list, nvl_cds_sclst)
    nv_end_idx=overlap(scl_loc_exon_list, nvl_cds_sclend)
    if nv_end_idx == -1:
        print('Unable to find end of CDS within transcript (cds end > exon end). There may be a problem with annotation')
        return(pd.DataFrame(), -4)
    if nv_start_idx == -1:
        print('Unable to find start of CDS within transcript (cds start < exon  start). There may be a problem with annotation')
        return(pd.DataFrame(), -5)
    ##now remap scaled location to genomic location
    gen_cds_start = exon_df.start.iloc[nv_start_idx] + (nvl_cds_sclst - exon_df.scaled_start.iloc[nv_start_idx])
    gen_cds_end = exon_df.stop.iloc[nv_end_idx] + (nvl_cds_sclend - exon_df.scaled_stop.iloc[nv_end_idx])
    #HOLY SHIT THERE ARE PIPES IN PYTHON  
    if abs(nv_end_idx - nv_start_idx) > 0 :      
        first_row= pd.DataFrame({'type':'CDS', 'start':gen_cds_start, 'stop':exon_df.stop.iloc[nv_start_idx]}, index=[0])
        if abs((nv_start_idx+1)-nv_end_idx) > 0:
            middle_row= exon_df >> mutate(type='CDS') >> select(X.type, X.start, X.stop) 
            middle_row=middle_row.iloc[(nv_start_idx+1):nv_end_idx].reset_index(drop=True)
        else :
            middle_row=pd.DataFrame(columns=['type', 'start','stop'])
        last_row=pd.DataFrame({'type':'CDS','start':exon_df.start.iloc[nv_end_idx], 'stop':gen_cds_end}, index=[nv_end_idx])
        cds_complete= first_row >> bind_rows(middle_row, join='inner') >> bind_rows(last_row, join='inner')
    else:# sneaky little single exon transcripts 
        cds_complete=pd.DataFrame({'type':'CDS', 'start':gen_cds_start, 'stop':gen_cds_end}, index=[0])
    
    cds_complete= cds_complete.reset_index()
    ### build the remaing annotations
    ##build five prime uttr
    
    fputr_end_idx=overlap(scl_loc_exon_list, nvl_cds_sclst - 1)
    top_rows=pd.DataFrame(columns=['type', 'start','stop'])
    if cds_complete.start.iloc[0] == exon_df.start.iloc[0]:# if the cds starts exactly at the first exon, don't build the fputr
        fputr=pd.DataFrame(columns=['type', 'start','stop'])
    else:
        if fputr_end_idx!=0:
            top_rows=exon_df.iloc[0:fputr_end_idx] >> mutate(type='five_prime_UTR') >> select(X.type, X.start, X.stop)
        fputr_end_coord= exon_df.start.iloc[fputr_end_idx]+ ((nvl_cds_sclst - 1) - exon_df.scaled_start.iloc[fputr_end_idx])
        
        bottom_rows=pd.DataFrame({'type':'five_prime_UTR', 'start':exon_df.start.iloc[fputr_end_idx], 'stop':fputr_end_coord}, index=[0])
        
        fputr=top_rows >> bind_rows(bottom_rows, join='inner')    
    #god bless reverse indexing

    if cds_complete.stop.iloc[-1] == exon_df.stop.iloc[-1]: # cds ends on end of last exon, ie no tputr
        tputr=pd.DataFrame(columns=['type', 'start','stop'])
    else:
    ##build three prime utr
        tputr_start_idx = overlap(scl_loc_exon_list, nvl_cds_sclend+1) # t
        bottom_rows=pd.DataFrame(columns=['type', 'start','stop'])
        l_row=len(exon_df.index)
        if tputr_start_idx < l_row-1:
            bottom_rows = exon_df.iloc[(tputr_start_idx+1):l_row] >> mutate(type='three_prime_UTR') >> select(X.type, X.start, X.stop)
        tputr_start_coord = exon_df.start.iloc[tputr_start_idx] + ((nvl_cds_sclend+1) - exon_df.scaled_start.iloc[tputr_start_idx])
        top_rows= pd.DataFrame({'type':'three_prime_UTR', 'start':tputr_start_coord, 'stop':exon_df.stop.iloc[tputr_start_idx]}, index=[0])
        
        tputr= top_rows >> bind_rows(bottom_rows, join='inner')  
    
    full_feat_df= fputr >> bind_rows(cds_complete, join='inner') >> bind_rows(tputr,join='inner') 
    return(full_feat_df.reset_index(drop=True), score)
 

def make_cds_annotation(ref_gff, st_gtf,novel_fasta, ref_fasta, tx ):
    '''
    function to make gene model of novel transcript using information from the reference transcript
    ref_gff: df of gencode reference gff
    dn_gtf: stringtie gtf subset to novel exon transcripts only, with annotation on which exon is novel
    dn_fast: fasta records of transdecoder protein seqs
    ref_fasta: gencode protein coding transcripts
    tx: tx to examine
    return codes:
    0 - no refernce tx found for current tx
    -1 - reference tx is not protien coding
    -2 - novel tx is not protein coding
    -3 - no good match between novel and reference
    '''

    gtf_tx= st_gtf.query('transcript_id == @tx')
    try:
        ref_tx=gtf_tx[gtf_tx['cmp_ref'].str.contains('ENST')]['cmp_ref'].iloc[0]# i have to be doing this w
    except:
       # print('except')
        return(pd.DataFrame, 0)
    if ref_tx not in ref_fasta.keys() :
        # 0 signifies not found, -1 is no match 
        return(pd.DataFrame, -1)
    elif tx not in novel_fasta.keys():
        return(pd.DataFrame(),-2)
        #return('')
    
    
    ref_seq=ref_fasta[ref_tx]
    nov_seq=novel_fasta[tx]
    gtf_tx=gtf_tx.assign(start=pd.to_numeric(gtf_tx.start),stop=pd.to_numeric(gtf_tx.stop)).sort_values(['start']).reset_index(drop=True)
    gtf_exon=gtf_tx.query('type == "exon"')[['start','stop']].reset_index(drop=True)
#    gtf_exon=gtf_exon.assign(start=pd.to_numeric(gtf_exon.start), stop=pd.to_numeric(gtf_exon.stop)).sort_values(['start'])
    ref_cds=ref_gff.query('transcript_id==@ref_tx & type=="CDS"').sort_values(['start'])[['start','stop']]
    scl_df=scale_genomic_length(gtf_exon.start, gtf_exon.stop)
    gtf_exon=pd.concat([gtf_exon, scl_df], axis=1)
    
    
    gene_model, score = make_gene_model(ref_seq, nov_seq, gtf_exon, ref_cds)
    if score <= -3:
        return(pd.DataFrame, score)
    tx_line= gtf_tx.reset_index(drop=True).query('type == "transcript"') >> select(X.seqid, X.source,X.type, X.start, X.stop, X.score,X.strand, X.frame, X.transcript_id )
    gene_model_full = gene_model >> mutate(seqid=tx_line.seqid.iloc[0], source= tx_line.source.iloc[0], score= tx_line.score.iloc[0], strand=tx_line.strand.iloc[0], frame=tx_line.frame.iloc[0], transcript_id=tx_line.transcript_id.iloc[0]) >> select(X.seqid, X.source,X.type, X.start, X.stop, X.score,X.strand, X.frame, X.transcript_id )
    gtf_exon=gtf_tx.query('type=="exon"') >> select(X.seqid, X.source,X.type, X.start, X.stop, X.score,X.strand, X.frame, X.transcript_id )
    all_features_df= tx_line >> bind_rows(gene_model_full, join='inner') >> bind_rows(gtf_exon, join='inner') 
    
    return(all_features_df.reset_index(drop=True), score)    

#
#st_gtf=read_GTF('results_b38/all_tissues.combined.gtf')
#st_gtf=st_gtf.where(pd.notnull(st_gtf),'.')
#novel_exon_gtf=pd.read_csv('novel_exon_gtf.tsv', sep='\t')
#ref_gff=pd.read_csv('ref/gencodeGFF3.gff.tsv', sep='\t')
#ref_gff=ref_gff.rename(columns={'end':'stop'})
#novel_fasta=SeqIO.to_dict(SeqIO.parse(novel_prot_fa_file,'fasta'))
#ref_fasta=SeqIO.to_dict(SeqIO.parse(ref_prot_fa_file, 'fasta'))
#
#
#######
##ipy scratch 
#
#tx=nv_id
#ref_tx=ref_id
#exon_df=gtf_exon
#ref=ref_seq
#novel=nov_seq
#cds_df=ref_cds
#n=10
#########
#
##################
###write tests here
#pc_ref_tx= ref_gff.query('type=="transcript" & transcript_type=="protein_coding" & strand == "+" ') >> pull('transcript_id')
#new_id_to_oId = st_gtf[st_gtf.oId.isin(pc_ref_tx)] >> select(X.transcript_id, X.oId)
#n=100
#test_state_1 = 4224# 87/100 all errors with trd
#test_state_2 = 1824# 88/100 1 UTR error not on my end, rest are trd
#test_state_3 = 7439
#t2g_samp=new_id_to_oId.sample(n,random_state=test_state_3)
#
#res=[0] * n
#for i in range(n):
#    feats=['CDS', 'five_prime_UTR','three_prime_UTR']
#    ref_id=t2g_samp.oId.iloc[i]
#    nv_id=t2g_samp.transcript_id.iloc[i]
#    ref_model=ref_gff.query('transcript_id == @ref_id & type!="start_codon" & type !="stop_codon"')
#    new_model, score=make_cds_annotation(ref_gff, st_gtf,novel_fasta, ref_fasta, nv_id)
#    if score != 1:
#        res[i]=score
#        continue
#    k=list()
#    for feature in feats:
#        rm=ref_model.query('type == @feature').reset_index(drop=True) >> select(X.start, X.stop)
#        nm=new_model.query('type == @feature').reset_index(drop=True) >> select(X.start, X.stop)
#        j=rm.equals(nm.assign(start=pd.to_numeric(nm.start), stop=pd.to_numeric(nm.stop)))
#        k.append(j)
#    #print(k)
#    res[i]=sum(k)
## scores:
## 3 = perfect
## 2 = 1 missing feature
## 1 = 2 missing features 
## 0 = no reference for st tx
## -1 = ref is not protein coding
## -2 = novel is not protein coding
## -3 no good match for between ref and novel 
#res_S=pd.Series(res)
#res_S.value_counts()
#
#my_bad_2=list(res_S[res_S ==2].index)
#
#
#
#j=ref_gff.query('transcript_id == "ENST00000641515.2" & type== "CDS"') >> select( X.start , X.stop)
#n=cds_tot(j)






















