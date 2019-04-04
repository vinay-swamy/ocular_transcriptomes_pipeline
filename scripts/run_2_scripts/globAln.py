#!/usr/bin/env python3
import sys
from Bio import SeqIO, pairwise2
args=sys.argv
ref_fasta_fl=args[1]
novel_fasta_fl=args[2]
ref_tx=args[3]
novel_tx=args[4]

ref_fa=SeqIO.to_dict(SeqIO.parse(ref_fasta_fl, 'fasta'))
novel_fa=SeqIO.to_dict(SeqIO.parse(novel_fasta_fl, 'fasta'))
seq1=ref_fa[ref_tx]
seq2=novel_fa[novel_tx]
aln=pairwise2.align.globalms(seq1.seq,seq2.seq, 2, -1, -1, -2, penalize_end_gaps=True)
print(pairwise2.format_alignment(*aln[0]))
match=pairwise2.format_alignment(*aln[0]).count('|')

print(len(seq2.seq))
print(novel_name+':'+str(match/len(seq2.seq)))
