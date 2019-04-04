#!/usr/bin/env python3
import sys
from Bio import SeqIO, pairwise2
try:
    args=sys.argv
    ref=args[1]
    novel=args[2]
    ref_name=ref.split('/')[1].split('.')[0]
    novel_name=novel.split('/')[1].split('.')[0]

    seq1=[k for k in SeqIO.parse(ref, "fasta")][0].seq
    seq2=[k for k in SeqIO.parse(novel, "fasta")][0].seq

    aln=pairwise2.align.globalms(seq1,seq2, 1, -1, -2, -2, penalize_end_gaps=False)

    match=pairwise2.format_alignment(*aln[0]).count('|')
    score=match/len(seq2)
    print('{ref}:{novel}\t{rlen}\t{nlen}\t{score}\n'.format(ref=ref_name,novel=novel_name,rlen=len(seq1),nlen=len(seq2), score=score))
except:
    exit()
