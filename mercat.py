#!/usr/bin/env python

"""mercat.py: Python code for Parallel k-mer counting."""
"""Usage: python mercat.py path-to-input-file kmer-value [no-of-cores] """
"""Example: To compute all 2-mers:
            python mercat.py test.faa 2 (uses all available cores)
            python mercat.py test.faa 2 8 (uses 8 cores)"""
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
import psutil
import timeit
import humanize
import pandas as pd
from collections import OrderedDict
from joblib import Parallel, delayed


inputfile = sys.argv[1]
kmer = int(sys.argv[2])

num_cores = 1
prune_kmer = 10 #remove all kmers whose count is < 10

try:
    num_cores = int(sys.argv[3])
except:
    num_cores = psutil.cpu_count(logical=False)

try:
    prune_kmer = int(sys.argv[4])
except:
    prune_kmer = 10

print "Running mercat using " + str(num_cores) + " cores"

bif = os.path.splitext(os.path.basename(inputfile))[0]

def get_all_substrings(input_string):
    length = len(input_string)
    return [input_string[i:i + kmer] for i in xrange(length-kmer+1)]

# ipfiles = glob.glob("*.faa")
# ipfiles.extend(glob.glob("*.fa"))

sequences = OrderedDict()
is_fastq = False
kmerstring = str(kmer)+"-mers"

start_time = timeit.default_timer()


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

print "Number of sequences in " + inputfile + " = "+ str(humanize.intword(len(sequences)))

def calculateKmerCount(seq):
    kmerlist = dict()
    kmerlist_all_seq = dict()
    cseq = sequences[seq] #Get current sequence
    sslist = get_all_substrings(cseq) # get all substrings of current sequence
    kmerlist_all_seq[seq] = dict() #kmer count for each substring of current sequence
    #print sslist
    for ss in sslist:
        if ss not in kmerlist: kmerlist[ss] = 0 #global kmer count for substring ss
        count = len(re.findall(r'(?=(%s))' % re.escape(ss), cseq))
        kmerlist[ss] += count #global kmer count for substring ss
        if count >= prune_kmer:
            #if ss not in kmerlist_all_seq[seq]: kmerlist_all_seq[seq][ss] = 0  # kmer count for substring in current sequence
            kmerlist_all_seq[seq][ss] = count #kmer count for substring in current sequence

    return [kmerlist,kmerlist_all_seq]

results = Parallel(n_jobs=num_cores)(
    delayed(calculateKmerCount)(seq) for seq in sequences)

kmerlist = dict()
kmerlist_all_seq = dict()

for d in results:
    for k,v in d[0].iteritems():
        if k in kmerlist:
            kmerlist[k] += v
        else: kmerlist[k] = v

for d in results:
    for seq,kdict in d[1].iteritems():
        #assert seq not in kmerlist_all_seq
        kmerlist_all_seq[seq] = kdict.copy()

print "Time to compute " + kmerstring +  ": " + str(round(timeit.default_timer() - start_time,2)) + " secs"

significant_kmers = []
for k in kmerlist:
    if kmerlist[k] >= prune_kmer: significant_kmers.append(k)

print "Total number of " + kmerstring +  " found: " + str(humanize.intword(len(kmerlist)))
print kmerstring +  " with count >= " + str(prune_kmer) + ": " + str(humanize.intword(len(significant_kmers)))

#df = df.ix[df[bif] >= prune_kmer]

df = pd.DataFrame(0,index=significant_kmers,columns=[bif])
for k in significant_kmers:
    df.set_value(k,bif,kmerlist[k])
df.to_csv(bif+"_summary.csv",index_label='Kmer',index=True)

dfcol = significant_kmers
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

df = df.loc[:,df.max() >= prune_kmer]
df.to_csv(bif+".csv",index_label='Sequence',index=True)

print "Total time: " + str(round(timeit.default_timer() - start_time,2)) + " secs"


#Debug
# sname = '515620.EUBELI_01521'
# print df.loc[sname,"length"]
# print df.loc[sname,"GC"]
# print df.loc[sname,"AT"]

