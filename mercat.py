#!/usr/bin/env python

"""mercat.py: Python code for k-mer counting."""
"""Usage: python mercat.py path-to-input-file kmer-value """
"""Example: To compute all 2-mers -> python mercat.py test.faa 2"""
"""Results are stored in input-file-name.csv and input-file-name_summary.csv
   (test.csv and test_summary.csv in the above example)"""
"""test.csv contains kmer count for kmers in individual sequences
   test_summary.csv contains kmer count for all unique kmers across all sequences in the sample test.fa"""

__author__      = "Ajay Panyala, Richard A. White III"
__copyright__   = "Copyright 2016"

import sys
import re
import os
import glob
import pandas as pd
from collections import OrderedDict

inputfile = sys.argv[1]
kmer = int(sys.argv[2])

bif = os.path.splitext(os.path.basename(inputfile))[0]
#print bif

def get_all_substrings(input_string):
    length = len(input_string)
    return [input_string[i:i + kmer] for i in xrange(length-kmer+1)]

# ipfiles = glob.glob("*.faa")
# ipfiles.extend(glob.glob("*.fa"))

sequences = OrderedDict()
#seq_kmers = dict()

is_fastq = False


with open(inputfile,'r') as f:
    for line in f:
        if line.startswith(">"): break
        elif line.startswith("@"):
            is_fastq = True
            break



with open(inputfile,'r') as f:
    if not is_fastq:
        seq = ""
        sname = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if sname: sequences[sname] = ""
                if seq:
                    sequences[sname] = seq
                    seq = ""
                sname = line[1:]
            else:
                seq += line

        assert sname and seq
        sequences[sname] = seq

    else: #process fastq file
        seq = ""
        sname = ""
        for line in f:
            line = line.strip()
            if line.startswith("@"):
                sname = line[1:].split()[0]
            elif line.startswith("+"):
                if seq:
                    sequences[sname] = seq
                    seq = ""
            else:
                if sname not in sequences: seq = line


#print sequences.keys()[0] + "="+ sequences.values()[0]

print "Number of sequences in " + inputfile + " = "+ str(len(sequences))
kmerlist = dict()
kmerlist_all_seq = dict()



for seq in sequences:
    cseq = sequences[seq] #Get current sequence
    sslist = get_all_substrings(cseq) # get all substrings of current sequence
    kmerlist_all_seq[seq] = dict() #kmer count for each substring of current sequence
    #print sslist
    for ss in sslist:
        if ss not in kmerlist: kmerlist[ss] = 0 #global kmer count for substring ss
        if ss not in kmerlist_all_seq[seq]: kmerlist_all_seq[seq][ss] = 0 #kmer count for substring in current sequence
        count = len(re.findall(r'(?=(%s))' % re.escape(ss), cseq))
        kmerlist[ss] += count #global kmer count for substring ss
        kmerlist_all_seq[seq][ss] = count #kmer count for substring in current sequence

#print kmerlist

dfcol = kmerlist.keys()
dfcol.extend(["length","GC","AT"])

df = pd.DataFrame(0,index=sequences.keys(),columns=dfcol)

for seq in sequences:
    cseq = sequences[seq]
    len_cseq = float(len(cseq))
    df.set_value(seq, "length", int(len_cseq))
    df.set_value(seq, "GC", round(((cseq.count("G")+cseq.count("C")) / len_cseq) * 100.0))
    df.set_value(seq, "AT", round(((cseq.count("A")+cseq.count("T")) / len_cseq) * 100.0))
    for ss in kmerlist_all_seq[seq]:
        df.set_value(seq, ss, kmerlist_all_seq[seq][ss])

df.to_csv(bif+".csv",index_label='Sequence',index=True)

df = pd.DataFrame(0,index=kmerlist.keys(),columns=[bif])
for k,v in kmerlist.items():
    df.set_value(k,bif,v)

df.to_csv(bif+"_summary.csv",index_label='Kmer',index=True)

# import csv
# with open(bif+'.csv', 'wb') as csv_file:
#     writer = csv.writer(csv_file)
#     for key, value in kmerlist.items():
#        writer.writerow([key, value])


