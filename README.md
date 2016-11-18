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
 
 
Dependencies
------------
Mercat needs the following python libraries:   
> pip install pandas joblib
  

Current Limitations
--------------------
- Mercat only works with 1 input file at a time.
- multiple file summary

