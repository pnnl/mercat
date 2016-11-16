Mercat: Python code for k-mer counting.
=======================================
  
>Usage: python mercat.py path-to-input-file kmer-value   
Example: To compute all 2-mers -> python mercat.py test.faa 2  

- Results are stored in input-file-name.csv and input-file-name_summary.csv  
   (test.csv and test_summary.csv in the above example)  
- test.csv contains kmer count for kmers in individual sequences  
- test_summary.csv contains kmer count for all unique kmers across all sequences in the sample test.fa
 
 
Dependencies
------------
Mercat needs pandas library: pip install pandas

TODO
----
- Mercat currently takes only 1 file at a time and supports only `faa` and `fasta` formats. FastQ format is not yet supported.