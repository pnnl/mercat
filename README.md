MerCat: python code for versatile k-mer counter and diversity estimation for database independent property analysis (DIPA)  obtained from metagenomic and/or metatranscriptomic sequencing data
================================================

![GitHub Logo](mercat.jpg)

  
Installing MerCat: 
 - Will be Available soon via Anaconda: Enable BioConda repo and run `conda install mercat`  
 - We do not have a pip installer available as of now. If you would like to use pip, please install the 
   modules listed in `dependencies.txt` via pip and run `python setup.py install` for setting up mercat.
 
Usage:
-----
 * -i I        path-to-input-file
 * -f F        path-to-folder-containing-input-files
 * -k K        kmer length
 * -n N        no of cores [default = all]
 * -c C        minimum kmer count [default = 10]
 * -pro        run mercat on protein input file specified as .faa 
 * -q          tell mercat that input file provided are raw nucleotide reads as [.fq, .fastq]
 * -p          run prodigal on nucleotide assembled contigs. Must be one of ['.fa', '.fna', '.ffn', '.fasta']
 * -t [T]      Trimmomatic options
 * -h, --help  show this help message


By default mercat assumes that inputs provided is one of ['.fa', '.fna', '.ffn', '.fasta']

> Example: To compute all 3-mers, run `mercat -i test.fa -k 3 -n 8 -c 10 -p`          
 
 The above command:
* Runs prodigal on `test.fa`, then runs mercat on the resulting protein file.            
* Results are generally stored in input-file-name_{protein|nucleotide}.csv and input-file-name_{protein|nucleotide}_summary.csv  
   * `test_protein.csv` and `test_protein_summary.csv` in this example  
* `test_protein.csv` contains kmer frequency count, pI, Molecular Weight, and Hydrophobicity metrics for individual sequences.  
* `test_protein_summary.csv` contains kmer frequency count, pI, Molecular Weight, and Hydrophobicity metrics for all unique kmers across all sequences in `test.fa`
* `test_protein_diversity_metrics.txt` containing the alpha diversity metrics.
  
Other usage examples:
---------------------

* `mercat -i test.fq -k 3 -n 8 -c 10 -q`  
   Runs mercat on raw nucleotide read (.fq or .fastq) 
   
*  `mercat -i test.fq -k 3 -n 8 -c 10 -q -t`  
   Runs trimmomatic on raw nucleotide reads (.fq or .fastq), then runs mercat on the trimmed nucleotides
    
*  `mercat -i test.fq -k 3 -n 8 -c 10 -q -t 20`  
   Same as above but can provide the quality option to trimmomatic
   
*  `mercat -i test.fq -k 3 -n 8 -c 10 -q -t 20 -p`
   Run trimmomatic on raw nucleotide reads, then run prodigal on the trimmed read to produce a protein file which is then processed by mercat
      
*  `mercat -i test.fna -k 3 -n 8 -c 10`  
   Run mercat on nucleotide input - one of ['.fa', '.fna', '.ffn', '.fasta']
    
*   `mercat -i test.fna -k 3 -n 8 -c 10 -p`  
    Run prodigal on nucleotide input, generate a .faa protein file and run mercat on it
    
*   `mercat -i test.faa -k 3 -n 8 -c 10 -pro`  
    Run mercat on a protein input (.faa)

* All the above examples can also be used with  `-f input-folder` instead of `-i input-file` option
  -  Example:  `mercat  -f /path/to/input-folder -k 3 -n 8 -c 10` --- Runs mercat on all inputs in the folder
  
  
Citing Mercat
-------------
If you are publishing results obtained using MerCat, please cite:

Richard Allen White III, Ajay Panyala, Kevin Glass, Sean Colby, Kurt R Glaesemann, Christer Jansson, Janet K Jansson. (2017).
MerCat: a versatile k-mer counter and diversity estimator for database-independent property analysis obtained from metagenomic and/or metatranscriptomic sequencing data. PeerJ Pre-Print.

Richard Allen White III, Ajay Panyala, Kevin Glass, Sean Colby, Kurt R Glaesemann, Christer Jansson, Janet K Jansson. (2017).
MerCat: a versatile k-mer counter and diversity estimator for database-independent property analysis obtained from metagenomic and/or metatranscriptomic sequencing data. Bioinfomatics. Submitted


CONTACT
-------

Please send all queries to Ajay Panyala <[ajay.panyala@gmail.com](ajay.panyala@gmail.com)> and Richard White <[richard.white@pnnl.gov](richard.white@pnnl.gov)> 
 or <[raw937@gmail.com](raw937@gmail.com)> 
