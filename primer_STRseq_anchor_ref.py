# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 14:57:29 2020

@author: snm0205
"""

import re
import os
import time
#remember to change directory
directory = "C:\\Users\\snm0205\\Desktop\\UMI\\Run2_10ng_8samples_300cycles"
os.chdir(directory)

start = time.time()

primer_STRSeq_anchor_ref = []
#input primer file
file_primer = "PrimedAnchors.txt"
with open(file_primer, 'r') as fp:
    for line in fp:
        (val1Locus, val2Chr, keyPos, val3Strand, val4Primer, val5Anchor) = (line.rstrip('\n')).split('\t')
        #print(val3Strand)
        if (val3Strand == "0"):
            #print(val1Locus, val2Chr, keyPos, val3Strand, val4Primer, val5Anchor)
            with open('GCA_000001405.15_GRCh38_no_alt_analysis_set.fastq', 'r') as f:
                for line in f:
                    foo = re.search(r'(%s(.{1,400})%s)' % (val4Primer, val5Anchor), line)
                    if foo is not None:
                        primer_STRSeq_anchor_ref.append((foo.group(1),foo.group(2),val4Primer,val5Anchor))
            #        if re.search(r'CAGTCTCCATAAATATGTGAGTCAATTCCCCAAGTG(.{,200})TCGTCTATCTATCCAGTC', line) is not None:
            #            refseq = re.search(r'CAGTCTCCATAAATATGTGAGTCAATTCCCCAAGTG(.{,100})TCGTCTATCTATCCAGTC', line)
            #            #refseq = re.search(r'AGTGAATTGCCT(.*)CTACCTCCTATTAGTCTGTCTCTGGAGAACATTGAC', line)
            #            #refseq = re.search(r'GACCCTGTCCTAGCCTTCTTATAGCTGCTAT(.*{,600})GTTATAAAAATA', line)
            #            print (refseq.group(1))
        else:
            with open('GCA_000001405.15_GRCh38_no_alt_analysis_set.fastq', 'r') as f:
                for line in f:
                    foo = re.search(r'(%s(.{1,400})%s)' % (val5Anchor, val4Primer), line)
                    if foo is not None:
                        primer_STRSeq_anchor_ref.append((foo.group(1),foo.group(2),val4Primer,val5Anchor))

fho = open("primer_STRSeq_anchor_ref.txt",'w')
for i in primer_STRSeq_anchor_ref:
    fho.write(str('\t'.join(i)) + "\n")
            
end = time.time()
print(end-start)