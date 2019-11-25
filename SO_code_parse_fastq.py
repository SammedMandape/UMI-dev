# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 12:45:55 2019

@author: snm0205
"""

#!/usr/bin/env python

import sys
import os

# =============================================================================
# def process(lines=None):
#     ks = ['name', 'sequence', 'optional', 'quality']
#     return {k: v for k, v in zip(ks, lines)}
# =============================================================================
#fn = sys.argv[1]
try:
    fn = sys.argv[1]
except IndexError as ie:
    raise SystemError("Error: Specify file name\n")

if not os.path.exists(fn):
    raise SystemError("Error: File does not exist\n")

# =============================================================================
# n = 4
# #breakpoint
# with open(fn, 'r') as fh:
#     lines = []
#     for line in fh:
#         lines.append(line.rstrip())
#         if len(lines) == n:
#             record = process(lines)
#             #sys.stderr.write("Record: %s\n" % (str(record)))
#             print(record['name'])
#             lines = []
# =============================================================================


#function that takes in a fastq file and an empty dictionary and returs
#dictionary with seqid as key and seq as values
def dict_for_fastq(file_fastq, dict_fastq_empty):
     n = 4
     with open(file_fastq, 'r') as fh:
         lines = []
         #count = 0
         for line in fh:
             lines.append(line.rstrip())
             if len(lines) == n:
                 #count += 1
                 ks = ['name', 'sequence', 'optional', 'quality']
                 record = {k: v for k, v in zip(ks, lines)}
                 #sys.stderr.write("Record: %s\n" % (str(record)))
                 #print(record['name'],record['sequence'])
                 dict_fastq_empty[record['name']] = record['sequence']
                 lines = []
                 #print(count)
         return dict_fastq_empty
                
#define an empty dictionary for Read1 fastq   
dict_fastq_R1 = {}
#file_fastq = fn
dict_for_fastq(fn, dict_fastq_R1)
print(dict_fastq_R1.keys())
