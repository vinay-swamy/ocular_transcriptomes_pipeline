import os
import sys
import pandas as pd
from  dfply import *
from Bio import SeqIO, pairwise2
args=['','/Users/vinayswamy/NIH/eyeintegration_splicing', 'results_38/all_tissues.combined.gtf', 'ref/gencodeGFF3.gff.tsv',\
      'results_b38/best_orfs.transdecoder.pep', 'ref/gencodeProtSeq.fa']
args=sys.argv
working_dir=args[1]
st_gtf=args[2]
ref_gff=args[3]
novel_prot_fa=args[4]
ref_prot_fa=args[5]


os.chdir('/Users/vinayswamy/NIH/eyeintegration_splicing')

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
    print('{}:{}'.format(novel.id, str(score)))
    r_seq, m_seq, n_seq= aln.split('\n')[0:3]    
    #part 2 figure out where on the genome
    start_codon=cds_df.iloc[[0]]
    if m_seq.find(match)==-1:
        print('No good match found')
        #return('')
    # find location of the relatuve  
    ref_start=start_codon['start'].iloc[0]
    gen_loc_exon_list=list(zip(exon_df.start, exon_df.stop))
    ref_start_idx=overlap(gen_loc_exon_list, start)# locartion
    ref_cds_start= exon_df.scaled_start.iloc[ref_start_idx]+ (ref_start - exon_df.start.iloc[ref_start_idx])
    if r_seq[0] != '-' and n_seq[0]!='-' :
        # if the starts are the same, then the end is the length*3 nt away    
        nvl_cds_sclst=cds_start
    elif r_seq[0]=='-':
        # if the start of the novel occurs before the start of the ref, find the ref start relative to the novel
        # and subtract the distance from the genomic ref start ot get the novel ref start
        nvl_cds_sclst=ref_cds_start - m_seq.find(match) * 3    
    else:
        # if the start of the novel occurs after the start of the novel ref, find the distance to the match point
        # and then add it to the genomic start to determine the novel CDS start
        nvl_cds_sclst=cds_start + m_seq.find(match) * 3

    # finding the end of the cds is just adding the length of the prot seq, regard less of
    nvl_cds_sclend=cds_start+len(novel)*3
    scl_loc_exon_list=list(zip(exon_df.scaled_start, exon_df.scaled_stop))
    
    nv_start_idx= overlap(scl_loc_exon_list, nvl_cds_sclst)
    nv_end_idx=overlap(scl_loc_exon_list, nvl_cds_sclend)
    ##now remap scaled location to genomic location
    gen_cds_start = exon_df.start.iloc[nv_start_idx] + (nvl_cds_sclst - exon_df.scaled_start.iloc[nv_start_idx])
    gen_cds_end = exon_df.stop.iloc[nv_end_idx] + (nvl_cds_sclend - exon_df.scaled_stop.iloc[nv_end_idx])
    #HOLY SHIT THERE ARE PIPES IN PYTHON
    first_row= pd.DataFrame({'type':'CDS', 'start':gen_cds_start, 'stop':exon_df.stop.iloc[nv_start_idx]}, index=[0])
    middle_row= exon_df >> mutate(type='CDS') >> select(X.type, X.start, X.stop) 
    middle_row=middle_row.iloc[(nv_start_idx+1):nv_end_idx].reindex(range(1, nv_end_idx))
    
    last_row=pd.DataFrame({'type':'CDS','start':exon_df.start.iloc[nv_end_idx], 'stop':gen_cds_end}, index=[nv_end_idx])
    
    cds_complete= first_row >> bind_rows(middle_row, join='inner') >> bind_rows(last_row, join='inner')
    ### build the remaing annotations
    ##build five prime uttr
    top_rows=pd.DataFrame(columns=['type', 'start','stop'])
    if nv_start_idx!=0:
        top_rows=exon_df.iloc[0:nv_end_idx] >> mutate(type='five_prime_utr') >> select(X.type, X.start, X.stop)
    bottom_rows=pd.DataFrame({'type':'five_prime_utr', 'start':exon_df.start.iloc[nv_start_idx], 'stop':gen_cds_start -1}, index=[0])
    
    fputr=top_rows >> bind_rows(bottom_rows, join='inner')    
    ##build three prime utr
    bottom_rows=pd.DataFrame(columns=['type', 'start','stop'])
    l_row=len(exon_df.index)
    if nv_end_idx != l_row-1:
        bottom_rows = exon_df.iloc[nv_end_idx:l_row] >> mutate(type='three_prime_utr') >> select(X.type, X.start, X.stop)
    top_rows= pd.DataFrame({'type':'three_prime_utr', 'start':gen_cds_end+1, 'stop':exon_df.stop.iloc[end_idx]}, index=[0])
    
    tputr= top_rows >> bind_rows(bottom_rows, join='inner')  
    full_feat_df= fputr >> bind_rows(cds_complete, join='inner') >> bind_rows(tputr,join='inner') 
    return(full_feat_df.reset_index(drop=True))
 

def make_cds_annotation(ref_gff, st_gtf,novel_fasta, ref_fasta, tx ):
    '''
    function to make gene model of novel transcript using information from the reference transcript
    ref_gff: df of gencode reference gff
    dn_gtf: stringtie gtf subset to novel exon transcripts only, with annotation on which exon is novel
    dn_fast: fasta records of transdecoder protein seqs
    ref_fasta: gencode protein coding transcripts
    tx: tx to examine
    '''

    gtf_tx= st_gtf.query('transcript_id == @tx')
    try:
        ref_tx=gtf_tx[gtf_tx['cmp_ref'].str.contains('ENST')]['cmp_ref'].iloc[0]# i have to be doing this w
    except:
        print('No ref tx in strigntie gtf')
        #return('')
    if ref_tx not in ref_fasta.keys():
        print('No Protein Coding Reference Transcript')
        #return('')
    ref_seq=ref_fasta[ref_tx]
    nov_seq=novel_fasta[tx]

    gtf_exon=gtf_tx.query('type == "exon"').reset_index(drop=True).sort_values(['start'])[['start','stop']]
    gtf_exon=gtf_exon.assign(start=pd.to_numeric(gtf_exon.start), stop=pd.to_numeric(gtf_exon.stop))
    ref_cds=ref_gff.query('transcript_id==@ref_tx & type=="CDS"').sort_values(['start'])[['start','stop']]
    scl_df=scale_genomic_length(gtf_exon.start, gtf_exon.stop)
    gtf_exon=pd.concat([gtf_exon, scl_df], axis=1)
    
    gene_model = make_gene_model(ref_seq, nov_seq, gtf_exon, ref_cds)
    tx_line= gtf_tx.reset_index(drop=True).query('type == "transcript"') >> select(X.seqid, X.source,X.type, X.start, X.stop, X.score,X.strand, X.frame, X.transcript_id )
        
        
    # now we know both the ref     
    
    return(0)    


######
#tx='TCONS_00010321'
tx='TCONS_00011242'    
    
ref_tx=st_gtf[st_gtf['cmp_ref'].str.contains('ENST')].query('transcript_id == @tx')['cmp_ref'].iloc[0]# i have to be doing this w

ref=ref_fasta[ref_tx]
novel=novel_fasta[tx]
i=k.split('\n')
[print(j[0:10]) for j in i[0:3]]
loc_tup=match_tx_find_gloc(ref_seq, nov_seq, start_codon )

start=pd.to_numeric(gtf_exon['start'])
stop=pd.to_numeric(gtf_exon['stop'])


######3



st_gtf=read_GTF('results_b38/all_tissues.combined.gtf')
st_gtf=st_gtf.where(pd.notnull(st_gtf),'.')
novel_exon_gtf=pd.read_csv('novel_exon_gtf.tsv', sep='\t')
ref_gff=pd.read_csv('ref/gencodeGFF3.gff.tsv', sep='\t')
ref_gff=ref_gff.rename(columns={'end':'stop'})
novel_fasta=SeqIO.to_dict(SeqIO.parse(novel_prot_fa,'fasta'))
ref_fasta=SeqIO.to_dict(SeqIO.parse(ref_prot_fa, 'fasta'))


































