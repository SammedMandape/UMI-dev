# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 15:53:21 2019

@author: snm0205
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor
Author: Sammed Mandape
Purpose: This python code will find UMIs given primers, read1(fastq)
and read2(fastq) as inputs. 
This is a temporary script file.
"""

import os
import time
import re
import collections
import strfuzzy
#import sys
#import pandas as pd
#from pytidyverse import *
#from dplython import (DplyFrame, X, diamonds, select, sift, sample_n, sample_frac, head, arrange, mutate, group_by, summarize, DelayFunction)

#remember to change directory
directory = "C:\\Users\\snm0205\\Desktop\\UMI\\Run2_10ng_8samples_300cycles"
os.chdir(directory)


complement = {'A' : 'T', 'C' : 'G', 'T' : 'A', 'G' : 'C'}

def reverse_complement(seq):
    bases = list(seq)
    bases = ''.join(complement[base] for base in reversed(bases))
    return bases

def getFastqNameR1R2(directory):
#directory = os.chdir("C:\\Users\\snm0205\\Desktop\\UMI\\Run2_10ng_8samples_300cycles")
#directory = os.fsencode('.')
    filedict = {}
    for filename in os.listdir(directory):
        #print(filename)
        #file = os.fsdecode(filename)
        if filename.endswith("001.fastq"):
        #    print(filename)
            if re.match(r'\d+-\d+_.*_R1_\d+\.fastq', filename) is not None:
                #print (filename)
                filebegin = re.match(r'(\d+-\d+_.*)_R1_(\d+\.fastq)', filename)
                #filebegin_r2 = re.match(r'\d+-\d+_.*_R2_.*\.fastq', filename)
                filebeginr1 = filebegin.group(1) + "_R1_001.fastq"
                filebeginr2 = filebegin.group(1) + "_R2_001.fastq"
                filedict[filebeginr1] = filebeginr2
                #print(filebeginr1, filebeginr2)
    return filedict

#directory = "C:\\Users\\snm0205\\Desktop\\UMI\\Run2_10ng_8samples_300cycles"
filedict = getFastqNameR1R2(directory)

#function that takes in a primer file and an empty dictionary and returns 
#dictionary with chr-pos as key and primer as value
def dict_for_primer(file_primer):
    dict_primer_empty = {}
    if not file_primer:
        raise SystemError("Error: Specify primer file name\n")
    with open(file_primer, 'r') as fh_primer:
        for line in fh_primer:
            #(key, val) = line.split()
            #dict_primer_empty[key] = val
            (val1Locus, val2Chr, keyPos, val3Strand, val4Primer, val5Anchor) = (line.rstrip('\n')).split('\t')
            if val3Strand == "1":
                #testcount += 1
                val4Primer = reverse_complement(val4Primer)
                val5Anchor = reverse_complement(val5Anchor)
                #print ("This is the reverse complement: %s and %s" % (val4Primer, val5Anchor))
            else:
                pass
            dict_primer_empty[keyPos] = [val1Locus, val2Chr, val3Strand, val4Primer, val5Anchor]
    return dict_primer_empty

#function that takes in a fastq file and an empty dictionary and returs
#dictionary with seqid as key and seq as values
def dict_for_fastq(file_fastq):
    dict_fastq_empty = {}
    if not file_fastq:
        raise SystemError("Error: Specify fastq file name\n")
    n = 4
    with open(file_fastq, 'r') as fh:
         lines = []
         count = 0
         for line in fh:
             lines.append(line.rstrip())
             if len(lines) == n:
                 count += 1
                 ks = ['name', 'sequence', 'optional', 'quality']
                 record = {k: v for k, v in zip(ks, lines)}
                 #sys.stderr.write("Record: %s\n" % (str(record)))
                 #print(record['name'],record['sequence'])
                 dict_fastq_empty[record['name'].split(' ')[0]] = record['sequence']
                 lines = []
         print(count)        
    return dict_fastq_empty


#input primer file
file_primer = "PrimedAnchors.txt"
dict_primer = dict_for_primer(file_primer)


start = time.time()
counterCS_P = 0
counterCS = 0
counter_noCS_match = 0
#key_count = 0

# CS ATTGGAGTCCT
#UmiSTRLociList = []
#LociList = []
#LociRead2Seq_postCS = []

#input Read1 fastq file
filedict = getFastqNameR1R2(directory)
for key in filedict:
    file_fastq_R1 = key
    file_fastq_R2 = filedict[key]
    dict_fastq_R1 = dict_for_fastq(file_fastq_R1)
    dict_fastq_R2 = dict_for_fastq(file_fastq_R2)
    UmiSTRLociList = []
    #print(file_fastq_R1[:-6])
    #sys.exit(0)
    for key in set(dict_fastq_R1) & set(dict_fastq_R2):
        readR1 = dict_fastq_R1[key]
        readR2 = dict_fastq_R2[key]
        #key_count += 1
        #numMatches=0
        if re.match(r'(.{12})(ATTGGAGTCCT)', readR2) is not None:
            counterCS += 1
            for items in dict_primer.items():
                #re.search(r'%s(.*)', readR1).group(1)
                # do fuzzy matching of anchor
                anchor = items[1][4]
                anchorIndex = strfuzzy.fuzzyFind(readR1, anchor, fuzz=1)
                primer = items[1][3]
                #if re.match(r'%s(.*)%s' % (items[1][3], items[1][4]), readR1):
                if ((readR1.startswith(primer, 0, len(primer))) and (anchorIndex >= 0)):
                    Loci = items[1][0]
                    #R1 = readR1
                    #R2 = readR2
                    #primerIdx = re.match(r'%s.*' % primer, readR1)
                    STRseq =  readR1[len(primer):anchorIndex]
                    #print (STRseq)
                    ##print("lenprimer:%d, anchorIdx:%d, SeqId:%s, primer:%s, anchor:%s" % (len(primer), anchorIndex, key, primer, anchor))
                    #STRseq = re.match(r'%s(.*)%s' % (items[1][3], items[1][4]), readR1).group(1)
                    searchCS = re.match(r'(.{12})(ATTGGAGTCCT)', readR2)
                    UMI = searchCS.group(1)
                    #CommSeq = searchCS.group(2)
                    #readR2Seq = searchCS.group(3)
                    ##print (Loci, UMI, primer, readR1, readR2, STRseq)
                    #print(Loci, readR2Seq)
                    counterCS_P += 1
                    print (counterCS_P)
                    UmiSTRLociList.append((Loci, STRseq, UMI))
                    
                    #LociList.append(Loci)
                    #LociRead2Seq_postCS.append((Loci, STRseq, UMI))
                    ##################################################################
                    #
                    # Do snity check to see if there is more than one match for anchors in 
                    # readR1 using re.findall
                    #
                    #
                    ##################################################################
                    
                    ##if counterCS_P == 10:
                    ##    sys.exit(0)
        else:
            counter_noCS_match += 1

    #print(key_count)    
    UmiSTRLociCount = collections.defaultdict(int)       
    for k in UmiSTRLociList:
        UmiSTRLociCount[k] += 1
        
    with open('UmiSTRLociCount_%s_atStart_LName_fuzzy_1.txt' % (file_fastq_R1[:-6]), 'w') as fh:
        fh.writelines('{}\t{}\n'.format('\t'.join(k),v) for k,v in UmiSTRLociCount.items() )
    

end = time.time()
print(end-start)
#print(UmiLociCount.items())
print("Counter for CS = %d and counter for Primer = %d and counter for no CS match = %d" % (counterCS, counterCS_P, counter_noCS_match))

        