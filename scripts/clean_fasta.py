import sys
import csv
args=sys.argv
fasta_file=args[1]
outfile=args[2]

with open(fasta_file) as fasta , open(outfile,'w+') as of:
    for line in fasta:
        if line[0]=='>':
            spl=line.split('|')
            line='>' + spl[1] +'\n'
        of.write(line)
