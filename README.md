DESCRIPTION

VAPOR is provided with a large set of reference fastas (> 20,000) for a given segment, a set of reads, and retrieves one that is closest to the sample strain.

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
                            Minimum proportion of kmers required for query [0.25]
        -fa FA              Fasta file
        -fq FQ [FQ ...]     Fastq file/files, can be gzipped

Example:

    vapor.py -fa HA_sequences.fa -fq reads_1.fq.gz

Author: Joel Southgate

Email (for inquiries): southgateJA@cardiff.ac.uk

REQUIREMENTS

Python 3.x
NumPy >= 1.51.x

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

TESTING

To test, run:

    vapor.py -fq tests/11.fq -fa tests/HA_sample.fa

which should yield:

    0.9782480893592005  191190.0    1701    112.39858906525573  1000    >cds:ADO12563 A/Chile/3935/2009 2009/07/07 HA H1N1 Human

Where the tab-delimited fields correspond to: approximate fraction of query bases found in reads; total score; query length; mean score; number of reads surviving culling; query description
