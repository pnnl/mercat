Mercat: Python code for Parallel k-mer counting
================================================
  
Usage: python mercat.py path-to-input-file kmer-value [no-of-cores]
> Example: To compute all 2-mers use:  
 python mercat.py test.faa 2 (uses all available cores)  
 python mercat.py test.faa 2 8 (uses 8 cores)

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