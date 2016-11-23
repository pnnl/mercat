Mercat: Python code for Parallel k-mer counting
================================================
  
Usage: python mercat.py -i path-to-input-file -k kmer-length [-n no-of-cores] [-c kmer-min-count]
> Example: To compute all 3-mers:
            python mercat.py -i test.fa -k 3 -n 8 -c 10
            
- Supported Formats: `faa`, `fasta` and `fastQ`
- If core count is not specified, mercat will use all available physical cores on the machine.
- Results are stored in input-file-name.csv and input-file-name_summary.csv  
   (test.csv and test_summary.csv in the above example)  
- test.csv contains kmer count for kmers in individual sequences  
- test_summary.csv contains kmer count for all unique kmers across all sequences in the sample test.fa

New Items:
1. take a multiple file input then take output csv do MDS in matlibplot. 
 - MDS with each sample equaling one point
 - MDS each kmer from each sample as many points with different colors per sample 
2. calculate these metrics: from a .faa file
- Add pI calculation
http://isoelectric.ovh.org/

-Hydro, MW and pH
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3283891/



Allows user to clean up fastq file prior to using mercat
java -jar trimmomatic-0.33.jar SE -phred33 In.fastq Out.fastq ILLUMINACLIP:adapters/TruSeq2-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:50

To format nucleotide.fa/.fasta to .faa
add prodigal program for converting .fa, .fasta, fastq nucleotides to .faa amino acid file
shell command 
"prodigal -i {input} -o {output.gff} -a {output.orf_pro} -f gff -p meta -d {output.orf_nuc} 2> {log}"
input = nucleotide.fa or nucleotide.fasta
outputs = .faa file amino acid fasta, a gff, and nucleotide orf called .ffn
We will use the .faa for amino acid kmer counting if the person requests it!



Dependencies
------------
Mercat needs the following python libraries:   
> pip install pandas joblib
>pip3 install humanize --user
Add humanize 


Current Limitations
--------------------
- Mercat only works with 1 input file at a time.
- multiple file summary

