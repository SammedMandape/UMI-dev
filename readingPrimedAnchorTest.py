# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 10:50:55 2019

@author: snm0205
"""

import os
#import time
import re
import strfuzzy
#import collections
#import sys
#import pandas as pd
#from pytidyverse import *
#from dplython import (DplyFrame, X, diamonds, select, sift, sample_n, sample_frac, head, arrange, mutate, group_by, summarize, DelayFunction)

#remember to change directory
os.chdir("C:\\Users\\snm0205\\Desktop\\UMI\\Run2_10ng_8samples_300cycles")

complement = {'A' : 'T', 'C' : 'G', 'T' : 'A', 'G' : 'C'}

def reverse_complement(seq):
    bases = list(seq)
    bases = ''.join(complement[base] for base in reversed(bases))
    return bases



def dict_for_primer(file_primer, dict_primer_empty):
    #testcount = 0
    if not file_primer:
        raise SystemError("Error: Specify primer file name\n")
    with open(file_primer, 'r') as fh_primer:
        for line in fh_primer:
            #print((line.rstrip('\n')).split('\t'))
            (val1Locus, val2Chr, keyPos, val3Strand, val4Primer, val5Anchor) = (line.rstrip('\n')).split('\t')
            #print (val4Primer, val5Anchor)
            if val3Strand == "1":
                #testcount += 1
                val4Primer = reverse_complement(val4Primer)
                val5Anchor = reverse_complement(val5Anchor)
                #print ("This is the reverse complement: %s and %s" % (val4Primer, val5Anchor))
            else:
                pass
            dict_primer_empty[keyPos] = [val1Locus, val2Chr, val3Strand, val4Primer, val5Anchor]
    #sprint (testcount)
    return dict_primer_empty

dict_primer = {}

file_primer = "PrimedAnchors.txt"
dict_for_primer(file_primer, dict_primer)
#print(dict_primer)

s = "ACGTACGTCCCACACGGCCTGGCAACTTATATGTATTTTTGTATTTCATGTGTACATTCGTATCTATCTATCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATTCCCCACAGTGAAAATAATCTACAGGATAGGTAAATAAATTAAGGCATATTCACGCAATGGGATACGATACAGTGATGAAAATGAACTAATTATAGCTACGTGAAACTATACTCATGAACACAATTTTGTAAAAGAAACAGGACTCCAATTTTCGCTCTTCC"

for items in dict_primer.items():
    #print(items[1][3])
    #myregex = re.escape(items[1][3]) + (r".*") + re.escape(items[1][4])
    if re.search(r'%s(.*)%s' % (items[1][3], items[1][4]), s):
        print(re.search(r'%s(.*)%s' % (items[1][3], items[1][4]), s).group(1))
    else:
        print('no match')