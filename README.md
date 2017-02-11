Mercat: Python code for Parallel k-mer counting
================================================

![GitHub Logo](mercat.jpg){:height="36px" width="36px"}

  
Installing MerCat: 
 - Available in BioConda: Enable BioConda repo and run "conda install mercat"  
 
Usage:

optional arguments:
  -h, --help  show this help message and exit
  -i I        path-to-input-file
  -f F        path-to-folder-containing-input-files
  -k K        kmer length
  -n N        no of cores [default = all]
  -c C        minimum kmer count [default = 10]
  -pro        protein input file
  -q          fastQ input file
  -p          run prodigal on fasta file one of ['.fa', '.fna', '.ffn', '.fasta']
  -t [T]      Trimmomatic options


> Example: To compute all 3-mers:
            python mercat.py -i test.fa -k 3 -n 8 -c 10
            
- Supported Formats: `faa`, `fasta` and `fastQ`
- If core count is not specified, mercat will use all available physical cores on the machine.
- Results are stored in input-file-name.csv and input-file-name_summary.csv  
   (test.csv and test_summary.csv in the above example)  
- test.csv contains kmer count for kmers in individual sequences  
- test_summary.csv contains kmer count for all unique kmers across all sequences in the sample test.fa


- Usage: python mercat.py -i path-to-input-file -k kmer-value [-n no-of-cores] [-c kmer-min-count]  
- Example: To compute all 3-mers:  
            python mercat.py -i test.fa -k 3 -n 8 -c 10
- Results are stored in input-file-name.csv and input-file-name_summary.csv  
   (test.csv and test_summary.csv in the above example)
- test.csv contains kmer count for kmers in individual sequences  
   test_summary.csv contains kmer count for all unique kmers across all sequences in the sample test.fa 
   
   
