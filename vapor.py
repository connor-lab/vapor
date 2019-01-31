#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" 
usage: vapor.py [-h] [-q] [-k K] [-s S] [-fa FA]
                              [-fq FQ [FQ ...]]

required arguments:
    -fa FA              Input fasta file
    -fq FQ [FQ ...]     Input fastq file/files

optional arguments:
    -h, --help          Show this help message and exit
    -q, --quiet         Suppresses output to stderr
    --return_seqs       Returns a fasta of sequences, instead of hits       
    --return_best_n     returns the best n hits [1]
    -o                  Combined output to files with prefix O, none by default
    -k                  Kmer length [21]
    -t, --threshold     Pre-Filtering Score threshold [0.0]
    -s, --subsample     Number of reads to subsample, no subsampling by default
    -m, --min_kmer_prop
                        Minimum proportion of kmers required [0.8]

Example:
    vapor.py -fa HA_sequences.fa -fq reads_1.fq

Author: Joel Southgate
Email (for inquiries): southgateJA@cardiff.ac.uk

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

"""

import sys
import argparse
import os
import vapor as vp

def blockErr():
    """ Block std err for quiet mode """
    sys.stderr = open(os.devnull, 'w')

#def main(quiet, k, score_threshold, subsample_amount, return_seqs, fasta, fastqs, min_path_weight):
def main(args):

    # If quiet, don't output anything to stderr 
    if args.quiet:
        blockErr()

    sys.stderr.write("Loading database sequences\n")
    seqsh, seqs = vp.parse_fasta_uniq(args.fa)
    sys.stderr.write("Got %d unique sequences\n" % len(seqs))

    # Get database kmers for filtering
    sys.stderr.write("Getting database kmers\n")
    dbkmers = vp.get_kmers(seqs, args.k)
    dbkmersset = set()
    for dbkmer in dbkmers:
        dbkmersset.update(set(dbkmer))
    sys.stderr.write("Got %d database kmers\n" % len(dbkmersset))

    # Parse and pre-filter reads
    sys.stderr.write("Filtering reads\n")
    reads, nrawreads = vp.parse_and_prefilter(args.fq, dbkmersset, args.threshold, args.k)
    sys.stderr.write("%d of %d reads survived\n" % (len(reads),nrawreads))

    # Subsample reads
    sys.stderr.write("Subsampling reads\n")
    if args.subsample != None:
        reads = vp.subsample(reads, args.subsample)

    # Check there are still sequences remaining
    if len(reads) == 0:
        sys.stderr.write("Exiting. Is there any virus in your sequences? Try a lower filtering threshold.\n")
        sys.exit(1)

    # Build the wDBG from reads
    wdbg = vp.wDBG(reads, args.k, dbkmersset)
    sys.stderr.write("Got %d wdbg kmers\n" % len(wdbg.edges))
    if len(wdbg.edges) == 0:
        sys.stderr.write("Zero kmers remaining! None of the kmers in your reads were found in the database. More reads or a lower -k could help. \n")
        sys.exit(1)

    # Ask the wdbg to classify
    sys.stderr.write("Classifying\n")
    path_results = wdbg.classify(dbkmers, args.min_kmer_prop)
    results = path_results[-args.return_best_n:]
    results = results[::-1]

    # Output results
    if args.return_seqs == True:
        for c, score in results:
            if score != -1:
                print(seqsh[c])
                print(seqs[c])
    elif args.output_prefix != None:
        scores_outf = open(args.output_prefix + ".out", "w")
        for c, score in results:
            if score != -1:
                slen = len(seqs[c])
                prop = str(score/slen)
                scores_outf.write(str(score)+"\t"+str(slen)+"\t" +str(prop) + "\t"+ str(len(reads)) + "\t"+seqsh[c] + "\n")
        scores_outf.close()
        seqs_outf = open(args.output_prefix + ".fa", "w")
        for c, score in results:
            if score != -1:
                seqs_outf.write(seqsh[c]+"\n")
                seqs_outf.write(seqs[c]+"\n")
        seqs_outf.close()
    else:
         for c, score in results:
            if score != -1:
                slen = len(seqs[c])
                prop = str(score/slen)
                print(str(score)+"\t"+str(slen)+"\t" +str(prop) + "\t"+ str(len(reads)) + "\t"+seqsh[c])

if __name__ == '__main__':
    # CLI
    parser = argparse.ArgumentParser(description="Do some viral classification!")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-q", "--quiet", action="store_true")
    group2 = parser.add_mutually_exclusive_group()
    group2.add_argument("--return_seqs", action="store_true")
    group2.add_argument("-o", "--output_prefix", type=str, help="Prefix to write full output to, stout by default", nargs='?', default=None)

    parser.add_argument("--return_best_n", type=int, default=1)
    parser.add_argument("-m", "--min_kmer_prop", type=float, help="Minimum proportion of mismatched kmers allowed [default=0.3]", nargs='?', default=0.9)
    parser.add_argument("-k", type=int, help="Kmer Length [15 > int > 30, default=21]", nargs='?', default=21)
    parser.add_argument("-t", "--threshold", type=float, help="Kmer filtering threshold [0 > float > 1, default=0.0]", nargs='?', default=0.0)
    parser.add_argument("-fa", type=str, help="Fasta file")
    parser.add_argument("-fq", nargs='+', type=str, help="Fastq file/files")
    parser.add_argument("-s", "--subsample", type=int, help="Number of reads to subsample [default=all reads]", nargs='?', default=None)

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    # Set some thresholds for user input
    max_kmer = 30
    min_kmer = 5
    max_thres = 1
    min_thres = 0

    if args.k <= max_kmer and args.k >= min_kmer:
        if args.threshold <= max_thres and args.threshold >= min_thres:
            main(args)
        else:
            sys.stderr.write("\nPlease input valid score threshold for prefiltering ({} to {}) \n \n".format(min_thres, max_thres))
            parser.print_help(sys.stderr)
            sys.exit(1)
    else:
        sys.stderr.write("\nPlease input valid kmer length ({} to {}) \n \n".format(min_kmer, max_kmer))
        parser.print_help(sys.stderr)
        sys.exit(1)
