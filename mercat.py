import sys
import re
import os
import glob
import csv

inputfile = sys.argv[1]
kmer = int(sys.argv[2])


def get_all_substrings(input_string):
    length = len(input_string)
    return [input_string[i:i + kmer] for i in xrange(length-kmer)]

# ipfiles = glob.glob("*.faa")
# ipfiles.extend(glob.glob("*.fa"))


sequences = []

with open(inputfile,'r') as f:
    seq = ""
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if len(seq)>0:
                sequences.append(seq)
                seq = ""
        else:
            seq += line

print len(sequences)
kmerlist = dict()

for seq in sequences:
    sslist = get_all_substrings(seq)
    for ss in sslist:
        if ss not in kmerlist: kmerlist[ss] = 0
        kmerlist[ss] += len(re.findall(r'(?=(%s))' % re.escape(ss), seq))


#print kmerlist

# bif = os.path.basename(inputfile).split[0]
# print bif


with open('output.csv', 'wb') as csv_file:
    writer = csv.writer(csv_file)
    for key, value in kmerlist.items():
       writer.writerow([key, value])

