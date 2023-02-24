#!/usr/bin/env python
import sys
import re

if len(sys.argv) != 3:
    print("Usage: python {} ref.bed miniprot.gff3".format(sys.argv[0]))
    print("Usage for Miniprot: miniprot --gff hs.gen.fa.gz mm.pep.fa.gz > miniprot.gff3")
    sys.exit(1)

ref_bed = sys.argv[1]
miniprot_gff3 = sys.argv[2]


infordb = {}
with open(miniprot_gff3, 'r') as infile:
    for line in infile:
        if 'mRNA' in line:
            data = line.strip().split()
            match = re.search(r'Target=(\S+)', line)
            if match:
                gene = match.group(1)
                infordb[gene] = infordb.get(gene, '') + data[0] + '\t'

with open('Allele.ctg.table', 'w') as outfile, open(ref_bed, 'r') as infile:
    for line in infile:
        data = line.strip().split()
        gene = re.sub(r';.*', '', data[3])
        if gene not in infordb:
            continue
        tdb = infordb[gene].split()
        tmpdb = {}
        for t in tdb:
            tmpdb[t] = tmpdb.get(t, 0) + 1
        outfile.write(data[0] + '\t' + data[2] + '\t')
        outfile.write('\t'.join(tmpdb.keys()) + '\n')
