#!/usr/bin/env python3
import sys
args=sys.argv
infasta= args[1]
outfasta=args[2]
len_cor_tab=args[3]
def fix_seq(seq):
    aa=seq[0]
    os=0
    while aa !='M':
        os+=1
        try:
            aa=seq[os]
        except:
            return('no Met found\n', 'NA')
    if len(seq[os:]) <=15:
        return ('too short\n','NA')
    return(seq[os:]+'\n', len(seq[os:]))

with open(infasta) as fasta, open(outfasta, 'w+') as of, open(len_cor_tab, 'w+') as lct :
    od=dict
    seqs=['M']
    header='toss'
    count=0
    for line in fasta:
        if '>' == line[0]:
            header_fixed=header[1:].strip('\n').split('.')[0]
            full_seq=''.join(seqs).replace('\n', '').replace('*','')
            full_seq_clean,seq_len=fix_seq(full_seq)
            of.write('>'+header_fixed+'\n')
            of.write(full_seq_clean)
            lct.write('{}\t{}\n'.format(header_fixed,seq_len))
            header=line
            seqs=[]
        else:
            seqs.append(line)
