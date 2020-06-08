import os
import time
import re
import collections
import strfuzzy

complement = {'A' : 'T', 'C' : 'G', 'T' : 'A', 'G' : 'C'}

def reverse_complement(seq):
    '''
    This function gives out a reverse complement of sequence.
    @param seq: The input seqeunce to construct a reverse complement of.
    @return: Reverse complement of input sequence.
    '''
    bases = list(seq)
    bases = ''.join(complement[base] for base in reversed(bases))
    return bases

def dict_for_primer(file_primer):
    ''' 
    This function constructs a dictionary of a primer file.
    @param file_primer: The primer file with Locus, Chr, Pos, Strand, Primer, Anchor
    @return: Dictionary of primer file with pos as key and list of values
    '''
    dict_primer_empty = {}
    if not file_primer:
        raise SystemError("Error: Specify primer file name\n")
    with open(file_primer, 'r') as fh_primer:
        for line in fh_primer:
            (val1Locus, val2Chr, keyPos, val3Strand, val4Primer, val5Anchor) = (line.rstrip('\n')).split('\t')
            if val3Strand == "1":
                val4Primer = reverse_complement(val4Primer)
                val5Anchor = reverse_complement(val5Anchor)
            else:
                pass
            dict_primer_empty[keyPos] = [val1Locus, val2Chr, val3Strand, val4Primer, val5Anchor]
    return dict_primer_empty
    
#input primer file
file_primer = "PrimedAnchors.txt"
dict_primer = dict_for_primer(file_primer)

def mainfunc(data, name):
    '''
    This function searches for primer and anchor in reads and pulls out STRseq between them.
    @param data: The tibble/data frame.
    @param name: The output file name.
    @return: Writes a output file with STRseq and other metadata.
    '''
    mypydata = data.set_index('ID').T.to_dict('list')
    UmiSTRLociList = []
    counter_P_A = 0
    counter_P = 0
    for mydatakey, mydataitems in mypydata.items():
        ID = mydatakey
        ReadCount = mydataitems[0]
        readR1 = mydataitems[1]
        for key, items in dict_primer.items():
            Pos = key
            anchor = items[4]
            anchorIndex = strfuzzy.fuzzyFind(readR1, anchor, fuzz=1)
            primer = items[3]
            if readR1.startswith(primer, 0, len(primer)):
                counter_P += 1
            if ((readR1.startswith(primer, 0, len(primer))) and (anchorIndex >= 0)):
                Loci = items[0]
                STRseq =  readR1[len(primer):anchorIndex]
                counter_P_A += 1
                UmiSTRLociList.append((ID, ReadCount, readR1, Loci, STRseq, primer, anchor))
                
    UmiSTRLociCount = collections.defaultdict(int)
    for k in UmiSTRLociList:
      UmiSTRLociCount[k] += 1
    
    outfilename = name + "_noN.tsv"
    with open(outfilename, 'w') as fh:
      fh.writelines("Number of Primer match = %d, Number of Primer and Anchor = %d" % (counter_P, counter_P_A))
      fh.writelines('{}\t{}\n'.format('\t'.join(k),v) for k,v in UmiSTRLociCount.items() )


  
