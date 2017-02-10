#!/usr/bin/env python

"""mercat.py: Python code for Parallel k-mer counting."""

__author__      = "Ajay Panyala, Richard A. White III"
__copyright__   = "Copyright 2016"

import sys
import re
import os
import glob
import psutil
import shutil
import timeit
import humanize
import subprocess
import pandas as pd
from collections import OrderedDict
from joblib import Parallel, delayed

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

import dask.dataframe as dd

from metrics import *
from Chunker import mercat_chunker


def check_command(cmd):
    cmd1 = cmd
    if cmd == 'trimmomatic': cmd1 = 'trimmomatic -version'
    with open(os.devnull, 'w') as FNULL:
        try:
            subprocess.check_call(cmd1, stdout=FNULL, stderr=FNULL, shell=True)
        except subprocess.CalledProcessError as e:
            # print e.output -- null since we suppressed output in check_call
            print("Mercat Error: %s not found, please setup %s using: conda install %s" %(cmd,cmd,cmd))
            sys.exit(1)



protein_file_ext = ['.fa','.fna','.ffn','.fasta']

def parseargs(argv=None):

    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    num_cores = psutil.cpu_count(logical=False)
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('-i', type=str, required=False, help='path-to-input-file') #default=nucleotide
    parser.add_argument('-f', type=str, required=False, help='path-to-folder-containing-input-files')
    parser.add_argument('-k', type=int, required = True, help='kmer length')
    parser.add_argument('-n', type=int, default=num_cores, help='no of cores [default = all]')  # no of cores to use
    parser.add_argument('-c', type=int, default=10, help='minimum kmer count [default = 10]')  # minimum kmer count to report
    parser.add_argument('-pro', action='store_true', help='protein input file')
    #parser.add_argument('-nuc', action='store_true', help='nucleotide input file')#default. We do not need this
    parser.add_argument('-q', action='store_true', help='fastQ input file')
    parser.add_argument('-p', action='store_true', help='run prodigal on fasta file')
    parser.add_argument('-t',type=int,nargs='?',const=30,required=False,help='Trimmomatic options')

    # Process arguments
    args = parser.parse_args()

    if args.i and args.f:
        parser.error("Can only specify either an input file (-i) or path to folder containing input files (-f) at a time")

    if args.i:
        if os.path.exists(args.i) < 1:
            path_error = "file " + args.i + " does not exist.\n"
            parser.error(path_error)
    elif args.f:
        if os.path.exists(args.f) < 1:
            path_error = "folder " + args.f + " does not exist.\n"
            parser.error(path_error)
    else:
        parser.error("Please provide either an input file (-i) or an input folder (-f)")


    if args.p:
        if args.pro: parser.error("Can only provide one of -p or -pro option at a time")
        check_command('prodigal')

    if args.t:
        if not args.q:
            parser.error("-t option has to be provided along with -q")

    if args.q:
        if args.t: check_command('trimmomatic')
        if args.pro:
            parser.error("Can only specify one of -q and -pro options at a time")
        if args.p:
            if not args.t:
                parser.error("Should specify -t option as well when both -q and -p options are provided")
        # if args.nuc:
        #     if not args.t:
        #         parser.error("Should specify -t option as well when both -q and -nuc options are provided")

    return [args,parser]


def get_all_substrings(input_string,kmer):
    length = len(input_string)
    return [input_string[i:i + kmer] for i in range(length-kmer+1)]

def calculateKmerCount(seq,cseq, prune_kmer,kmer):
    kmerlist = dict()
    kmerlist_all_seq = dict()
    #cseq = sequences[seq] #Get current sequence
    sslist = get_all_substrings(cseq,kmer) # get all substrings of current sequence
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


def check_args(ipfile,args,def_option,m_parser):
    given_ext = (os.path.splitext(ipfile)[1]).strip()
    if def_option:
        if given_ext not in protein_file_ext:
            m_parser.error("Input file provided should be one of the following formats: " + str(protein_file_ext))

    if args.p:
        if not args.q and given_ext not in protein_file_ext:
            m_parser.error("Input file provided should be one of the following formats: " + str(protein_file_ext))

    if args.q:
        if given_ext != ".fq":
            m_parser.error("Input file provided should be in .fq format")

    if args.pro:
        if given_ext != ".faa":
            m_parser.error("Input file provided should be in .faa format")


def mercat_main():
    __args__, m_parser = parseargs()

    kmer = __args__.k
    num_cores = __args__.n
    m_inputfile = __args__.i
    m_inputfolder = __args__.f
    prune_kmer = __args__.c
    mflag_fastq = __args__.q
    mflag_prodigal = __args__.p
    mflag_trimmomatic = __args__.t
    mflag_protein = __args__.pro

    np_string = "nucleotide"
    if mflag_protein or mflag_prodigal: np_string = "protein"
    def_option =  not __args__.p and not __args__.q and not __args__.pro

    all_ipfiles = []
    if m_inputfolder:
        m_inputfolder = os.path.abspath(m_inputfolder)
        os.chdir(m_inputfolder)
        #Assume all have same ext
        for fname in os.listdir(m_inputfolder):
            mip = os.path.join(m_inputfolder, fname)
            if not os.path.isdir(mip):
                # skip directories
                all_ipfiles.append(mip)

    else:
        m_inputfolder = os.getcwd()
        all_ipfiles.append(os.path.abspath(m_inputfile))

    for m_inputfile in all_ipfiles:

        os.chdir(m_inputfolder)
        check_args(m_inputfile,__args__,def_option,m_parser)

        m_inputfile = os.path.abspath(m_inputfile)

        basename_ipfile = os.path.splitext(os.path.basename(m_inputfile))[0] + "_" + np_string

        inputfile_size = os.stat(m_inputfile).st_size
        dir_runs = basename_ipfile + "_run"

        if os.path.exists(dir_runs):
            shutil.rmtree(dir_runs)
        os.makedirs(dir_runs)

        all_chunks_ipfile = []
        if inputfile_size >= 104857600: #100MB
            print("Large input file provided: Splitting it into smaller files...\n")
            mercat_chunker(m_inputfile,dir_runs,"100M",">")
            os.chdir(dir_runs)
            all_chunks_ipfile = glob.glob("*")
        else:
            os.chdir(dir_runs)
            all_chunks_ipfile.append(m_inputfile)

        #print all_chunks_ipfile
        #sys.exit(1)

        kmerstring = str(kmer) + "-mers"

        for inputfile in all_chunks_ipfile:

            bif = os.path.splitext(os.path.basename(inputfile))[0] + "_" + np_string

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
                prod_cmd = "prodigal -i %s -o %s -a %s -f gff -p meta -d %s" % (
                inputfile, bif + ".gff", gen_protein_file, bif + "_nuc.ffn")

                if mflag_fastq and mflag_trimmomatic:
                    trimfna = bif + "_trimmed.fna"
                    sequences = OrderedDict()
                    with open(inputfile, 'r') as f:
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

                    fnastring = ""

                    for sname in sequences:
                        fnastring += ">"+sname+"\n"
                        fnastring += sequences[sname]+"\n"

                    with open(trimfna, 'w') as f: f.write(fnastring)

                    gen_protein_file = bif + "_trimmed_pro.faa"
                    prod_cmd = "prodigal -i %s -o %s -a %s -f gff -p meta -d %s" % (
                    trimfna, bif + ".gff", gen_protein_file, bif + "_nuc.ffn")

                print prod_cmd
                with open(os.devnull, 'w') as FNULL:
                    subprocess.call(prod_cmd, stdout=FNULL, stderr=FNULL, shell=True)
                inputfile = gen_protein_file

            print("Running mercat using " + str(num_cores) + " cores")
            print("input file: " + inputfile)

            sequences = OrderedDict()
            is_fastq = False


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
                            sname = sname.split("#",1)[0].strip()
                        else:
                            line = line.replace("*","")
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

            print("Number of sequences in " + inputfile + " = "+ str(humanize.intword(len(sequences))))


            results = Parallel(n_jobs=num_cores)(
                delayed(calculateKmerCount)(seq, sequences[seq], prune_kmer, kmer) for seq in sequences)

            kmerlist = dict()
            kmerlist_all_seq = dict()

            for d in results:
                for k,v in d[0].items():
                    if k in kmerlist:
                        kmerlist[k] += v
                    else: kmerlist[k] = v

            for d in results:
                for seq,kdict in d[1].items():
                    #assert seq not in kmerlist_all_seq
                    kmerlist_all_seq[seq] = kdict#.copy()

            print("Time to compute " + kmerstring +  ": " + str(round(timeit.default_timer() - start_time,2)) + " secs")

            significant_kmers = []
            for k in kmerlist:
                if kmerlist[k] >= prune_kmer: significant_kmers.append(k)

            print("Total number of " + kmerstring +  " found: " + str(humanize.intword(len(kmerlist))))
            print(kmerstring +  " with count >= " + str(prune_kmer) + ": " + str(humanize.intword(len(significant_kmers))))

            #df = df.ix[df[bif] >= prune_kmer]


            if mflag_protein:
                df = pd.DataFrame(0.0, index=significant_kmers, columns=['Count',"PI","MW","Hydro"])
                for k in significant_kmers:
                    df.set_value(k, 'Count', kmerlist[k])
                    df.set_value(k,"PI", predict_isoelectric_point_ProMoST(k))
                    df.set_value(k, "MW", calculate_MW(k))
                    df.set_value(k, "Hydro", calculate_hydro(k))
            else:
                df = pd.DataFrame(0, index=significant_kmers, columns=['Count',"GC_Percent","AT_Percent"])
                for k in significant_kmers:
                    c_kmer = k
                    df.set_value(k, 'Count', kmerlist[k])
                    len_cseq = float(len(c_kmer))
                    df.set_value(k, "GC_Percent", round(((c_kmer.count("G")+c_kmer.count("C")) / len_cseq) * 100.0))
                    df.set_value(k, "AT_Percent", round(((c_kmer.count("A")+c_kmer.count("T")) / len_cseq) * 100.0))

            # df_summ_sort = df.sort_values('Count', ascending=False)
            # df_summ_sort.to_csv(bif + "_summary.csv", index_label=kmerstring, index=True)
            df.to_csv(bif + "_summary.csv", index_label=kmerstring, index=True)

            dfcol = significant_kmers


            if not mflag_protein:
                dfcol.extend(["length","GC_Percent","AT_Percent"])

                df = pd.DataFrame(0,index=list(sequences.keys()),columns=dfcol)

                for seq in sequences:
                    cseq = sequences[seq]
                    len_cseq = float(len(cseq))
                    df.set_value(seq, "length", int(len_cseq))
                    df.set_value(seq, "GC_Percent", round(((cseq.count("G")+cseq.count("C")) / len_cseq) * 100.0))
                    df.set_value(seq, "AT_Percent", round(((cseq.count("A")+cseq.count("T")) / len_cseq) * 100.0))
                    for ss in kmerlist_all_seq[seq]:
                        df.set_value(seq, ss, kmerlist_all_seq[seq][ss])

                    #df = df.loc[:, df.max() >= prune_kmer]
                    df1 = df.ix[:,['length','GC_Percent','AT_Percent']]
                    del df['length']
                    del df['GC_Percent']
                    del df['AT_Percent']
                    df = df.loc[:, df.max() >= prune_kmer]
                    df.loc[:, 'length'] = df1.ix[:,'length']
                    df.loc[:, 'GC_Percent'] = df1.ix[:,'GC_Percent']
                    df.loc[:, 'AT_Percent'] = df1.ix[:,'AT_Percent']


            else:

                dfcol.extend(["length", "PI", "MW","Hydro"])

                df = pd.DataFrame(0, index=list(sequences.keys()), columns=dfcol)

                for seq in sequences:
                    cseq = sequences[seq]
                    cseq=cseq.replace('*','')
                    len_cseq = float(len(cseq))
                    df.set_value(seq, "length", int(len_cseq))
                    df.set_value(seq, "PI", predict_isoelectric_point_ProMoST(cseq))
                    df.set_value(seq, "MW", calculate_MW(cseq))
                    df.set_value(seq, "Hydro", calculate_hydro(cseq))
                    for ss in kmerlist_all_seq[seq]:
                        df.set_value(seq, ss, kmerlist_all_seq[seq][ss])

                    #df = df.loc[:,df.max() >= prune_kmer]
                    df1 = df.ix[:,['length','PI','MW','Hydro']]
                    del df['length']
                    del df['PI']
                    del df['MW']
                    del df['Hydro']
                    df = df.loc[:, df.max() >= prune_kmer]
                    df.loc[:, 'length'] = df1.ix[:,'length']
                    df.loc[:, 'PI'] = df1.ix[:,'PI']
                    df.loc[:, 'MW'] = df1.ix[:,'MW']
                    df.loc[:, 'Hydro'] = df1.ix[:, 'Hydro']

            df.to_csv(bif+".csv",index_label='Sequence',index=True)

            print("Total time: " + str(round(timeit.default_timer() - start_time,2)) + " secs")

            #Debug
            # sname = '515620.EUBELI_01521'
            # print df.loc[sname,"length"]
            # print df.loc[sname,"GC"]
            # print df.loc[sname,"AT"]

        num_chunks = len(all_chunks_ipfile)
        df = dd.read_csv(basename_ipfile+"*_summary.csv")
        dfgb = df.groupby(kmerstring).sum()
        df10 = dfgb.nlargest(10,'Count').compute()
        dfsum = dfgb.sum(0).compute()



        if mflag_protein:
            df10[['PI', 'MW', 'Hydro']] = df10[['PI', 'MW', 'Hydro']]/num_chunks
            mercat_scatter_plots(basename_ipfile,'PI',df10,kmerstring)
            mercat_scatter_plots(basename_ipfile,'MW', df10, kmerstring)
            mercat_scatter_plots(basename_ipfile,'Hydro', df10, kmerstring)
            mercat_stackedbar_plots(basename_ipfile, 'Count', df10, kmerstring, dfsum.Count)
        else:
            df10[['GC_Percent', 'AT_Percent']] = df10[['GC_Percent','AT_Percent']]/num_chunks
            mercat_scatter_plots(basename_ipfile,'GC_Percent',df10,kmerstring)
            mercat_scatter_plots(basename_ipfile,'AT_Percent', df10, kmerstring)
            mercat_stackedbar_plots(basename_ipfile, 'Count', df10, kmerstring, dfsum.Count)


        all_counts = dfgb.Count.values.compute().astype(int)
        mercat_compute_alpha_beta_diversity(all_counts,basename_ipfile)


if __name__ == "__main__":
    mercat_main()
