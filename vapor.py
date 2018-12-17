#!/usr/bin/python3
# -*- coding: utf-8 -*-

""" 
usage: vapor.py [-h] [-q] [-k K] [-s S] [-fa FA]
                              [-fq FQ [FQ ...]]

optional arguments:
    -h, --help          Show this help message and exit
    -q, --quiet         Suppresses output to stderr
    --return_seqs       Returns a fasta of sequences, instead of hits       

    -k                  Kmer length [21]
    -t, --threshold     Pre-Filtering Score threshold [0.7]
    -s, --subsample     Number of reads to subsample, no subsampling by default
    -w, --weight        Min path weight to consider [20]
    -fa FA              Fasta file
    -fq FQ [FQ ...]     Fastq file/files

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
import numpy as np
import gzip
import random
import argparse
import os

class Path():
    """ Path object for traversing DBG """
    def __init__(self, start_edge, start_score):
        # Initialize from a start position in a DBG
        # Recount, for each edge, the next edge
        # For resolving cycles
        self.edges = {start_edge:[]}
        self.path = [start_edge]
        self.score = start_score
        self.cyclic = False
    def add_edge(self, edge, score):
        if edge in self.edges:
            print("Cycling")
            self.cyclic = True
        if edge not in self.edges:
            self.edges[edge] = []
        # Also add the last base to the previous edge
        # To avoid cycles
        last_kmer = self.path[-1]
        self.edges[last_kmer].append(edge[-1])
        self.path.append(edge)
        self.score += score
    def get_string(self):
        s = self.path[0]
        for edge in self.path[1:]:
            s += edge[-1]
        return s

def remove_overlaps(paths):
    """ Greedily resolving conflicting paths """
    # Removes the lowest scoring of a pair with nonzero kmer set intersection """
    flagarr = [0 for p in paths]
    for i in range(len(paths)-1):
        for j in range(i+1, len(paths)):
            if flagarr[i] == 0 and flagarr[j] == 0:
                pi = paths[i]
                pj = paths[j]
                if len(set(pi.edges.keys()) & set(pj.edges.keys())) > 0:
                    if pi.score > pj.score:
                        flagarr[j] = 1
                    else:
                        flagarr[i] = 1
    for i in range(len(flagarr)):
        if flagarr[i] == 0:
            yield paths[i]
 
class wDBG():
    """ Basic DBG with associated edge weights """
    def __init__(self, strings, k):
        # Only explicitly store edges
        self.edges = {}
        self.k = k
        self.start_positions = set()
        self._build(strings)

    def _build(self, strings):
        sys.stderr.write("Building wDBG\n")
        newkmerc = 0
        for si, string in enumerate(strings):
            sys.stderr.write(str(si) + "      \r")
            kmers = [string[i:i+self.k] for i in range(len(string)-self.k+1)]
            newkmerc = 0
            for kmer in kmers:
                if kmer in self.edges:
                    self.edges[kmer] += 1
                else:
                    newkmerc += 1
                    self.edges[kmer] = 1

    def get_statistics(self):
        """ Gets statistics from the graph, such as: """
        """ Degree distribution, branch ratios, weight distribution"""
        """ Number of kmers, weight distribution """
        """ The should be more thinly spread out for coinfection """
        n_kmers = len(self.edges)
        weights = np.array([w for kmer, w in self.edges.items()])
        total_weight = sum(weights)
        mean_weight = np.mean(weights)
        std_weight = np.std(weights)
        degrees = []
        branch_ratios = []
        for kmer, weight in self.edges.items():
            d = 0
            indegree_weights = []
            outdegree_weights = []
            for b in "ATCG":
                new_outkmer = kmer[1:] + b
                new_inkmer = b + kmer[:-1]
                if new_inkmer in self.edges:
                    d += 1
                    indegree_weights.append(self.edges[new_inkmer])
                if new_outkmer in self.edges:
                    d += 1
                    outdegree_weights.append(self.edges[new_outkmer])
            if len(outdegree_weights) > 1:
                branch_ratio = max(outdegree_weights)/sum(outdegree_weights)
                branch_ratios.append(branch_ratio)
            if len(indegree_weights) > 1:
                branch_ratio = max(indegree_weights)/sum(indegree_weights)
                branch_ratios.append(branch_ratio)
            degrees.append(d)
        mean_degree = np.mean(degrees)
        std_degree= np.std(degrees)
        mean_bratio = np.mean(branch_ratios)
        std_bratio = np.std(branch_ratios)
        return n_kmers, total_weight, mean_weight, std_weight, mean_degree, std_degree, mean_bratio, std_bratio

    def get_start_positions(self):
        """ Gets positions in the graph that have no inward edge """
        for kmer, weight in self.edges.items():
            edge = True
            for b in "ATCG":
                if b+kmer[:-1] in self.edges:
                    edge = False
            if edge == True:
                self.start_positions.add(kmer)

    def cull(self, kmers):
        """ Removes any kmer in kmers from edges """
        todel = []
        for kmer in self.edges:
            if kmer not in kmers:
                todel.append(kmer)
        for kmer in todel:
             del self.edges[kmer]

    def get_paths(self):
        """ Traverses the wDBG heuristically, getting a list of paths """
        paths = [Path(p, self.edges[p]) for p in self.start_positions]
        pstrings = [p.get_string() for p in paths]
        assert len(pstrings) == len(set(pstrings))
        final_paths = []
        sys.stderr.write("Building paths, from %d start points\n" % len(paths))
        tmp_paths = []
        while paths != []:
            switch = True
            assert len(paths) + len(final_paths) == len(self.start_positions)
            tmp_paths = []
            for path in paths:
                assert paths.count(path) == 1
                localswitch = False
                tmp_additions = {}
                for b in "ATCG":
                    # Check to see if the path has taken the route before
                    last_kmer = path.path[-1]
                    if b not in path.edges[last_kmer]:
                        tmpkmer = path.path[-1][1:] + b
                        if tmpkmer in self.edges:
                            tmp_additions[b] = self.edges[tmpkmer]
                            localswitch = True
                if localswitch == True:
                    maxtmp = max(tmp_additions.items(), key = lambda x:x[1])[1]
                    maxtmps = [t[0] for t in tmp_additions.items() if t[1] == maxtmp]
                    assert maxtmp != 0
                    best_tmp = maxtmps[0]
                    path.add_edge(path.path[-1][1:]+best_tmp, maxtmp)
                    tmp_paths.append(path)
                else:
                    final_paths.append(path)                
            paths = tmp_paths
            assert len(tmp_paths) <= len(paths)
            assert len(paths) <= len(self.start_positions)
            paths = tmp_paths
            sys.stderr.write("%d             \r" % (len(final_paths)))
        sys.stderr.write("\n")
        return final_paths          

class cDBG():
    """ Basic DBG with additional color information """
    def __init__(self, k):
        self.edges = {}
        self.k = k
        self.n = None
        self.color_classes = set()
    
    @classmethod
    def from_strings(cls, strings, k):
        """ Builds the cDBG from strings """
        c = cls(k)
        c.n = len(strings)
        c._build(strings)
        return c

    @classmethod
    def from_strings_and_subgraph(cls, strings, k, subgraph):
        """ As before, but only for parts intersecting with subgraph """
        c = cls(k)
        c.n = len(strings)
        c._build(strings, subgraph)
        c.build_color_classes()
        return c

    def _build(self, strings, subgraph=False):
        sys.stderr.write("Adding strings:\n")
        for si, string in enumerate(strings):
            sys.stderr.write(str(si)+"        \r")
            kmers = (string[i:i+self.k] for i in range(len(string)-self.k+1))
            for kmer in kmers:
                if subgraph==False or kmer in subgraph.edges:
                    if kmer in self.edges:
                        self.edges[kmer].append(si)
                    else:
                        self.edges[kmer] = [si]
        sys.stderr.write("\n")

    def build_color_classes(self, subkmers=None):
        """ Builds the color classes using binary strings """
        sys.stderr.write("Building color classes:\n")
        z = 0
        if subkmers == None:
            kmers = self.edges
        else:
            kmers = subkmers
        for kmer, ccl in kmers.items():
            sys.stderr.write(str(z) + "          \r")
            # Create empty bit string with leading 1
            cc = (1 << (self.n))
            for c in ccl:
                cc = cc | (1 << c)
            self.color_classes.add(cc)
            self.edges[kmer] = cc
            z += 1
        sys.stderr.write("\n")

    def classify(self, wdbg, seqs, min_path_weight=20.):
        """ Classifies a wdbg """

        # First get paths
        paths = [p for p in wdbg.get_paths()]
        paths = [p for p in remove_overlaps(paths) if p.score > min_path_weight]
        sys.stderr.write("Got %d fragments\n" % len(paths))
        if len(paths) == 0:
            sys.stderr.write("No paths with a greater weight than %d found. Please try a lower threshold (-w) than %d \n" % (min_path_weight, min_path_weight))
            sys.exit(1)
        
        mpl = max([len(p.get_string()) for p in paths])
        sys.stderr.write("%d paths found\n" % len(paths))
        sys.stderr.write("Maximum path length: %d\n" % mpl)

        total_colors = []
        # Next get arrays of color classes for paths
        for pi, path in enumerate(paths):
            colors = []   
            assert path.score > 0
            seq = path.get_string()
            kmers = (seq[i:i+self.k] for i in range(len(seq)-self.k+1))
            for kmer in kmers:
                assert kmer in wdbg.edges or kmer in self.edges
                if kmer in self.edges:
                    colors.append(self.edges[kmer])
            assert len(colors) > 0
            total_colors.append(colors)

        # Lastly contiguize the colors for scoring
        # Need to do it individually for path, and then aggregate them
        path_scores = []
        for colors in total_colors:
            sumo = np.zeros(len(seqs))
            color_flag = np.zeros(self.n)
            for c in colors:
                arr = np.fromstring(np.binary_repr(c), dtype='S1').astype(int)[1:]
                sumo += (self.k * (1-color_flag) + color_flag) * arr
            path_scores.append(sumo)

        # Now aggregate the scores for each path and rank
        aggsumo = np.zeros(len(seqs))
        for sumo in path_scores:
            aggsumo += sumo
        maxs = max(aggsumo)
        # Take the highest indices, noting they are backwards wrt 
        maxi = [i for i in range(len(aggsumo)) if aggsumo[i] == maxs]            
        # Finally, for the purposes of coinfection detection, we are going to 
        # Examine the worst rank that this score gets, in each path
        ranks = []
        for pi,sumo in enumerate(path_scores):
            maxi_path_score = sumo[maxi[0]]
            sortind = np.argsort(sumo)[::-1]
            sorted_scores = sumo[sortind]
            path_weight = paths[pi].score
            for i in range(len(sorted_scores)):
                if sorted_scores[i] <= maxi_path_score:
                    ranks.append([i,path_weight])
                    break

        # Order of the classes is different 
        maxi_cls = [len(aggsumo)-i-1 for i in range(len(aggsumo)) if aggsumo[i] == maxs]       
        yield len(path.path), path.score, maxi_cls, maxs, ranks

def get_kmers(strings,k):
    """ Takes strings and returns a set of kmers """
    kmers = set()
    for string in strings:
        for i in range(len(string)-k+1):
            kmers.add(string[i:i+k])
    return kmers

def rev_comp(read):
    read = read.replace("T", "a")
    read = read.replace("A", "t")
    read = read.replace("C", "g")
    read = read.replace("G", "c")
    return read.upper()[::-1]

def parse_and_prefilter(fqs, dbkmers, threshold, k):
    """ Parses fastq files fqs, and filters them """
    c = 0
    M = float(len(dbkmers))
    seen = set()
    for fq in fqs:
        if fq.endswith(".gz"):
            f = gzip.open(fq,'rt')
        else:
            f = open(fq,'r')
        for line in f:
            if c == 1:
                tmpseq = line.strip()
                kcount = 0
                # Don't allow Ns in read
                # Don't allow reads < k
                if "N" not in tmpseq and len(tmpseq) >= k and tmpseq not in seen:
                    seen.add(tmpseq)
                    for i in range(0, len(tmpseq), k):
                        if tmpseq[i:i+k] in dbkmers:
                            kcount += 1
                    # Only allow reads with a given number of words in dbkmers
                    if k*kcount/M < threshold:
                        yield tmpseq
            c += 1                  
            if c == 4:
                c = 0
        f.close()

def parse_fasta_uniq(fasta, filter_Ns=True):
    """ Gets unique sequences from a fasta, with filtering of Ns"""
    tmph = ""
    tmps = ""
    hs = []
    ss = []
    sseen = set()
    with open(fasta) as f:
        for line in f:
            l = line.strip()
            if l[0] == ">":
                if tmps not in sseen:
                    if ((filter_Ns == True) and "N" not in tmps) or filter_Ns == False:
                        hs.append(tmph)
                        ss.append(tmps) 
                        sseen.add(tmps)
                tmph = l
                tmps = ""
            else:
                tmps += l
    hs.append(tmph)
    ss.append(tmps)
    return hs, ss 

def subsample(reads, n):
    """ Takes a sample of n from reads """
    if n >= len(reads):
        return reads
    else:
#        random.shuffle(reads)
#        return reads[:n]
        return random.sample(reads, n)

def blockErr():
    """ Block std err for quiet mode """
    sys.stderr = open(os.devnull, 'w')

#def main(quiet, k, score_threshold, subsample_amount, return_seqs, fasta, fastqs, min_path_weight):
def main(args):

    # If quiet, don't output anything to stderr 
    if args.quiet:
        blockErr()

    sys.stderr.write("Loading database sequences\n")
    seqsh, seqs = parse_fasta_uniq(args.fa)
    sys.stderr.write("Got %d unique sequences\n" % len(seqs))

    # Get database kmers for filtering
    sys.stderr.write("Getting database kmers\n")
    dbkmers = get_kmers(seqs, args.k)

    sys.stderr.write("Filtering reads\n")
    # Parse and pre-filter reads
    reads = [r for r in parse_and_prefilter(args.fq, dbkmers, args.threshold, args.k)]
    sys.stderr.write("Subsampling reads\n")
    if args.subsample != None:
        reads = subsample(reads, args.subsample)
    sys.stderr.write(str(len(reads)) + " reads survived\n")

    # Check there are still sequences remaining
    if len(reads) == 0:
        sys.stderr.write("Exiting. Is there any virus in your sequences? Try a lower filtering threshold.\n")
        sys.exit(1)

    # Build the wDBG from reads
    wdbg = wDBG(reads, args.k)

    # Cull any kmers that are not present in the reference kmers
    sys.stderr.write("Culling kmers, beginning with %s\n" % len(wdbg.edges))
    wdbg.cull(dbkmers)
    sys.stderr.write("%d kmers remaining\n" % len(wdbg.edges))
    if len(wdbg.edges) == 0:
        sys.stderr.write("Zero kmers remaining! None of the kmers in your reads were found in the database. More reads or a lower -k could help. \n")
        sys.exit(1)
        
    sys.stderr.write("Getting start positions\n")
    # Get start positions for paths
    wdbg.get_start_positions()
    sys.stderr.write("%d start positions found\n" % len(wdbg.start_positions))
    mew = max(wdbg.edges.items(), key = lambda x : x[1])[1]
    sys.stderr.write("Largest edge weight: %d \n" % mew)

    # Build cDBG; don't build nodes that are not used, though
    cdbg = cDBG.from_strings_and_subgraph(seqs, args.k, wdbg)

    # Finally, classify
    path_results = cdbg.classify(wdbg, seqs, args.weight)
    if args.return_seqs == True:
        for length, weight, cls, score, ranks in path_results:
            for c in cls:
                    print(seqsh[c])
                    print(seqs[c])
    else:
        for length, weight, cls, score, ranks in path_results:
            print(str(length)+"\t"+str(weight)+"\t"+str(score)+"\t"+str(len(reads))+"\t"+",".join([seqsh[c] for c in cls]))
            for rank in ranks:
                print(rank)

    sys.stderr.write("\nClassification Complete\n")

if __name__ == '__main__':
    # CLI
    parser = argparse.ArgumentParser(description="Do some sweet viral classification")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-q", "--quiet", action="store_true")
    group2 = parser.add_mutually_exclusive_group()
    group2.add_argument("--return_seqs", action="store_true")

    parser.add_argument("-w", "--weight", type=int, help="Minimum Path Weight [default=20]", nargs='?', default=20)
    parser.add_argument("-k", type=int, help="Kmer Length [15 > int > 30, default=21]", nargs='?', default=21)
    parser.add_argument("-t", "--threshold", type=float, help="Kmer filtering threshold [0 > float > 1, default=0.7]", nargs='?', default=0.7)
    parser.add_argument("-fa", type=str, help="Fasta file")
    parser.add_argument("-fq", nargs='+', type=str, help="Fastq file/files")
    parser.add_argument("-s", "--subsample", type=int, help="Number of reads to subsample [default=all reads]", nargs='?', default=None)

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    # Set some thresholds for user input
    max_kmer = 30
    min_kmer = 15
    max_thres = 1
    min_thres = 0

    if args.k < max_kmer and args.k > min_kmer:
        if args.threshold < max_thres and args.threshold > min_thres:
            ###########Run main
#            main(args.quiet, args.k, args.s, args.r, args.return_seqs, args.fa, args.fq, args.w)
            main(args)
        else:
            sys.stderr.write("\nPlease input valid score threshold for prefiltering ({} to {}) \n \n".format(min_thres, max_thres))
            parser.print_help(sys.stderr)
            sys.exit(1)
    else:
        sys.stderr.write("\nPlease input valid kmer length ({} to {}) \n \n".format(min_kmer, max_kmer))
        parser.print_help(sys.stderr)
        sys.exit(1)
