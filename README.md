VAPOR: influenza Virus samPle clAssification frOm WGS Reads

usage: vapor.py [-h] [-q] [-k K] [-s S] [-fa FA]
                              [-fq FQ [FQ ...]]

optional arguments:
    -h, --help          Show this help message and exit
    -q, --quiet         Suppresses output to stderr
    --return_seqs       Returns a fasta of sequences, instead of hits       

    -k K                Kmer length [21]
    -s S                Pre-Filtering Score threshold [0.7]
    -fa FA              Fasta file
    -fq FQ [FQ ...]     Fastq file/files

Example:
    vapor.py -k 21 -s 0.7 -fa HA_sequences.fa -fq reads_1.fq

Author: Joel Southgate
Email (for inquiries): southgateJA@cardiff.ac.uk

