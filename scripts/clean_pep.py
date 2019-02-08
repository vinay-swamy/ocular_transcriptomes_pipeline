import sys
import csv
args=sys.argv
pep_file=args[1]
outfile=args[2]
pep_meta_info=args[3]
with open(pep_file) as pep , open(outfile,'w+') as of, open(pep_meta_info, 'w+', newline='\n') as m_inf:
    tsv=writer = csv.writer(m_inf, delimiter='\t')
    for line in pep:
        if line[0]=='>':
            spl=line.split(' ')
            line=spl[0] +'\n'
            meta_info=spl[1].split('~~')+ spl[3:]
            tsv.writerow(meta_info)
        of.write(line)
