#!/usr/bin/env python

"""mercat.py: Python code for Parallel k-mer counting."""
"""Usage: python mercat.py -i path-to-input-file -k kmer-value [-n no-of-cores] [-c kmer-min-count] """
"""Example: To compute all 3-mers:
            python mercat.py -i test.fa -k 3 -n 8 -c 10"""
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
import subprocess
import pandas as pd
from collections import OrderedDict
from joblib import Parallel, delayed

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter


from metrics import *

def check_command(cmd):
    cmd1 = cmd
    if cmd == 'trimmomatic': cmd1 = 'trimmomatic -version'
    with open(os.devnull, 'w') as FNULL:
        try:
            subprocess.check_call(cmd1, stdout=FNULL, stderr=FNULL, shell=True)
        except subprocess.CalledProcessError as e:
            # print e.output -- null since we suppressed output in check_call
            print "Mercat Error: %s not found, please setup prodigal using: conda install %s" %(cmd,cmd)




def parseargs(argv=None):

    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    num_cores = psutil.cpu_count(logical=False)
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('-i', type=str, required = True, help='path-to-input-file') #default=nucleotide
    parser.add_argument('-k', type=int, required = True, help='kmer length')
    parser.add_argument('-n', type=int, default=num_cores, help='no of cores [default = all]')  # no of cores to use
    parser.add_argument('-c', type=int, default=10, help='minimum kmer count [default = 10]')  # minimum kmer count to report
    parser.add_argument('-pro', action='store_true', help='protein input file')
    #parser.add_argument('-nuc', action='store_true', help='nucleotide input file')#default
    parser.add_argument('-q', action='store_true', help='fastQ input file')
    parser.add_argument('-p', action='store_true', help='run prodigal on fasta file')
    parser.add_argument('-t',type=int,nargs='?',const=30,required=False,help='Trimmomatic options')

    # Process arguments
    args = parser.parse_args()
    if os.path.exists(args.i) < 1:
        path_error = "file " + args.i + " does not exist.\n"
        parser.error(path_error)

    if args.p: check_command('prodigal')
    if args.t: check_command('trimmomatic')

    return args


def get_all_substrings(input_string):
    length = len(input_string)
    return [input_string[i:i + kmer] for i in xrange(length-kmer+1)]

def calculateKmerCount(seq,cseq, prune_kmer):
    kmerlist = dict()
    kmerlist_all_seq = dict()
    #cseq = sequences[seq] #Get current sequence
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


if __name__ == "__main__":
    __args__ = parseargs()

    kmer = __args__.k
    num_cores = __args__.n
    inputfile = __args__.i
    prune_kmer = __args__.c
    mflag_fastq = __args__.q
    mflag_prodigal = __args__.p
    mflag_trimmomatic = __args__.t
    mflag_protein = __args__.pro
    #mflag_nucleotide = __args__.nuc

    #print __args__

    bif = os.path.splitext(os.path.basename(inputfile))[0]

    '''trimmomatic SE -phred33 test.fq Out.fastq ILLUMINACLIP:TruSeq2-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:50'''
    if mflag_trimmomatic:
        swq = mflag_trimmomatic
        trimmed_file = bif+"_trimmed.fq"
        prod_cmd = "trimmomatic SE -phred33 %s %s ILLUMINACLIP:TruSeq2-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:%s MINLEN:50" %(inputfile,trimmed_file,swq)
        with open(os.devnull, 'w') as FNULL:
            subprocess.call(prod_cmd, stdout=FNULL, stderr=FNULL, shell=True)
        inputfile = trimmed_file


    "Run prodigal if specified"
    '''prodigal -i test_amino-acid.fa -o output.gff -a output.orf_pro.faa  -f gff -p meta -d output.orf_nuc'''
    if mflag_prodigal:
        mflag_protein = True
        gen_protein_file = bif+"_pro.faa"
        prod_cmd = "prodigal -i %s -o %s -a %s -f gff -p meta -d %s" %(inputfile,bif+".gff",gen_protein_file,bif+"_nuc.ffn")
        with open(os.devnull, 'w') as FNULL:
            subprocess.call(prod_cmd, stdout=FNULL, stderr=FNULL, shell=True)
        inputfile = gen_protein_file

    print "Running mercat using " + str(num_cores) + " cores"
    print "input file: " + inputfile


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


    results = Parallel(n_jobs=num_cores)(
        delayed(calculateKmerCount)(seq, sequences[seq], prune_kmer) for seq in sequences)

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


    if mflag_protein:
        df = pd.DataFrame(0.0, index=significant_kmers, columns=[bif,"PI","MW","Hydro"])
        for k in significant_kmers:
            df.set_value(k, bif, kmerlist[k])
            df.set_value(k,"PI", predict_isoelectric_point_ProMoST(k))
            df.set_value(k, "MW", calculate_MW(k))
            df.set_value(k, "Hydro", calculate_hydro(k))
    else:
        df = pd.DataFrame(0, index=significant_kmers, columns=[bif,"GC_Percent","AT_Percent"])
        for k in significant_kmers:
            c_kmer = k
            df.set_value(k, bif, kmerlist[k])
            len_cseq = float(len(c_kmer))
            df.set_value(k, "GC_Percent", round(((c_kmer.count("G")+c_kmer.count("C")) / len_cseq) * 100.0))
            df.set_value(k, "AT_Percent", round(((c_kmer.count("A")+c_kmer.count("T")) / len_cseq) * 100.0))

    df_summ_sort = df.sort_values(bif, ascending=False)
    df_summ_sort.to_csv(bif + "_summary.csv", index_label=kmerstring, index=True)

    dfcol = significant_kmers


    if not mflag_protein:
        dfcol.extend(["length","GC_Percent","AT_Percent"])

        df = pd.DataFrame(0,index=sequences.keys(),columns=dfcol)

        for seq in sequences:
            cseq = sequences[seq]
            len_cseq = float(len(cseq))
            df.set_value(seq, "length", int(len_cseq))
            df.set_value(seq, "GC_Percent", round(((cseq.count("G")+cseq.count("C")) / len_cseq) * 100.0))
            df.set_value(seq, "AT_Percent", round(((cseq.count("A")+cseq.count("T")) / len_cseq) * 100.0))
            for ss in kmerlist_all_seq[seq]:
                df.set_value(seq, ss, kmerlist_all_seq[seq][ss])


    else:

        dfcol.extend(["length", "PI", "MW","Hydro"])

        df = pd.DataFrame(0, index=sequences.keys(), columns=dfcol)

        for seq in sequences:
            cseq = sequences[seq]
            cseq=cseq.translate(None, '*')
            len_cseq = float(len(cseq))
            df.set_value(seq, "length", int(len_cseq))
            df.set_value(seq, "PI", predict_isoelectric_point_ProMoST(cseq))
            df.set_value(seq, "MW", calculate_MW(cseq))
            df.set_value(seq, "Hydro", calculate_hydro(cseq))
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

    if mflag_protein:
        mercat_scatter_plots(bif,'PI',df_summ_sort,kmerstring)
        mercat_scatter_plots(bif,'MW', df_summ_sort, kmerstring)
        mercat_scatter_plots(bif,'Hydro', df_summ_sort, kmerstring)
    else:
        mercat_scatter_plots(bif,'GC_Percent',df_summ_sort,kmerstring)
        mercat_scatter_plots(bif,'AT_Percent', df_summ_sort, kmerstring)