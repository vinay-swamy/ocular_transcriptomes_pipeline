import pandas as pd
import os
import sys
from itertools import chain

args=sys.argv
wd=args[1]
infile=args[2]
outfile=args[3]

def spread_coverage(exon):
    tdf=cov_df[cov_df['window'] == exon ].reset_index(drop=True)
    min_loc=tdf['start'].min()
    max_loc=tdf['end'].max()
    wstart=tdf['wstart'][0]
    wend=tdf['wend'][0]
    lol=[ [tdf['cov'][i]] *tdf['length'][i] for i in range(tdf.shape[0]) ]
    flat=list(chain(*lol))
    start_offset=wstart-min_loc
    end_offset=len(flat) -(max_loc-wend)
    comp_line=list(tdf.iloc[0,4:8])  + flat[start_offset:end_offset]
    return(comp_line)

os.chdir(wd)
global cov_df 
cov_df=pd.read_csv(infile,sep='\t' , names=['seqid','start', 'end', 'cov','w_seqid','wstart','wend','window']).assign(length= lambda x: x['end']-x['start'])


all_windows=list(cov_df['window'].unique())
res=[spread_coverage(e) for e in all_windows]
fin=pd.DataFrame(res)

fin.to_csv(outfile,sep='\t',header=False, index=False)
