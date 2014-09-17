Version
===================
This software was last updated 7/15/2014. Questions, comments, suggestions, etc?
Send emails to nick.zhangyuan@gmail.com

Installation
===================
1. Untar the source file called "SALT.tar.gz".
2. g++ compiler is required in your Unix system. To install SALT, run the Makeme file using the following command:
make
3. HMMER (http://hmmer.janelia.org/) should be installed. The bin file hmmsearch should be in the path.
4. Python version 2.5 or later is needed.
5. Python libraries of NetworkX (http://networkx.lanl.gov/) and Biopython (http://biopython.org/wiki/Main_Page) should be installed.

Run SATL
===================
To run SALT, use the following command:

./SALT.sh -m <HMMER3 HMM file> -f <fasta file> [options]
  Options:
    -h:  show this message
    -g:  gamma (position-specific score threshold, in the range of [0,1], default: 0.3)
    -k:  overlap threshold, default: (average read length) / 2
    -K:  number of contigs generated, default: the number of sinks
    -E:  E-value threshold for contigs, default: 1e-6
    -o:  output file name, default: standard error

The hmm file can contain multiple hmm models and should be in HMMER3.0's hmm file format. These files can be downloaded from Pfam ftp. The nucleotide sequence file should in be fasta format.
 
Output
===================
The output are the reads that are classified into each input family. The family names are the accessions of the Pfam families in the HMM files. Each family will have one block in the output file, with a header line that begins with a ">" symbol, which is followed by the family name. All the lines after the header line are names of reads that are classified into the family in the header line. If a family does not have any classified read, it will only have a header line. In the following example, PF00006 has 3 reads classified and PF00009 has 5 reads classified:

>PF00006
read1
read3
read5
>PF00009
read4
read6
read8
read9
read10

Referencing SALT
=======================
SALT can be referenced as:

Yuan Zhang, Yanni Sun, and James Cole, A Sensitive and Accurate protein domain cLassification Tool (SALT) for short reads, Bioinformatics, 29(17): 2103-2111, 2013

License
===================
Copyright (C) 2013 Yuan Zhang, Yanni Sun.
