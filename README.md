DESCRIPTION

VAPOR is provided with a large set of reference fastas (> 20,000) for a given segment, a set of reads, and retrieves one that is closest to the sample strain.

USAGE:

usage: vapor.py [-h] [-q] [-k K] [-s S] [-fa FA]
                              [-fq FQ [FQ ...]]

optional arguments:
    -h, --help          Show this help message and exit
    -q, --quiet         Suppresses output to stderr
    --return_seqs       Returns a fasta of sequences, instead of hits       

    -o                  Combined output to files with prefix O, none by default
    -k K                Kmer length [21]
    -s S                Pre-Filtering Score threshold [0.0]
    -fa FA              Fasta file
    -fq FQ [FQ ...]     Fastq file/files, can be gzipped

Example:
    vapor.py -k 21 -s 0.7 -fa HA_sequences.fa -fq reads_1.fq.gz

Author: Joel Southgate
Email (for inquiries): southgateJA@cardiff.ac.uk

REQUIREMENTS

Python 3.x
NumPy 1.51.x

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

