DESCRIPTION

VAPOR is a tool for classification of Influenza samples from raw short read sequence data for downstream bioinformatics analysis. VAPOR is provided with a fasta file of full-length sequences (> 20,000) for a given segment, a set of reads, and attempts to retrieve a reference that is closest to the sample strain.

REQUIREMENTS

Python 3.x
NumPy >= 1.5.1

INSTALLATION

The VAPOR module executable script can be installed with pip:

    sudo pip3 install git+https://github.com/connor-lab/vapor

Alternatively, clone the repo, and add to your PATH:

    git clone https://github.com/connor-lab/vapor

TESTING

A test dataset is provided in the repository/tests. To test, clone the repo as above, and run:

    vapor.py -fq tests/test_reads.fq -fa tests/HA_sample.fa

which should yield:

    0.9782480893592005  191190.0    1701    112.39858906525573  1000    >cds:ADO12563 A/Chile/3935/2009 2009/07/07 HA H1N1 Human

Where the tab-delimited fields correspond to: approximate fraction of query bases found in reads; total score; query length; mean score; number of reads surviving culling; query description

USAGE:

    usage: vapor.py [-h] [-q] [-k K] [-s S] [-fa FA]
                              [-fq FQ [FQ ...]]

    optional arguments:
        -h, --help          Show this help message and exit
        -q, --quiet         Suppresses output to stderr
        --return_seqs       Returns a fasta of sequences, instead of hits       
        --return_best_n     Returns the highest scoring n queries

        -o                  Combined output to files with prefix O, none by default
        -k K                Kmer length [21]
        -t T                Pre-Filtering score threshold [0.2]
        -s S                Number of reads to sub-sample
        -c, --min_kmer_cov  Minimum kmer coverage for culling [5]
        -m, --min_kmer_prop
                            Minimum proportion of kmers required for query [0.1]
        -fa FA              Fasta file
        -fq FQ [FQ ...]     Fastq file/files, can be gzipped
        -f, --top_seed_frac
                            Fraction of best seeds to extend [0.2]
        -v, --version       Show version

Example:

    vapor.py -fa HA_sequences.fa -fq reads_1.fq.gz reads_2.fq.gz

Author: Joel Southgate

Email (for inquiries): southgateJA@cardiff.ac.uk

COPYRIGHT

Copyright 2018 Joel Southgate

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.


