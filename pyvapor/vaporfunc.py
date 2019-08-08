""" This file contains general functions used by vapor """

import numpy as np
import gzip
import random

def kmers2str(kmers):
    """ Takes a set off kmers and extracts their string """
    s = kmers[0]    
    for k in kmers[1:]:
        s += k[-1]
    return s

def get_kmers(strings, k):
    """ Takes strings and returns a set of kmers """
    kmers = []
    for string in strings:
        kmers.append([string[i:i+k] for i in range(len(string)-k+1)])
    return kmers

def get_kmers_set(strings, k):
    """ Takes strings and returns a set of kmers """
    kmers = set()
    for string in strings:
        for kmer in [string[i:i+k] for i in range(len(string)-k+1)]:
            kmers.add(kmer)
    return kmers

def rev_comp(read):
    """ Basic (slow) reverse complement """
    read = read.replace("T", "a")
    read = read.replace("A", "t")
    read = read.replace("C", "g")
    read = read.replace("G", "c")
    return read.upper()[::-1]

def parse_and_prefilter(fqs, dbkmers, threshold, k):
    """ Parses fastq files fqs, and filters them """
    nraw = 0
    reads = []
    c = 0
    for fq in fqs:
        if fq.endswith(".gz"):
            f = gzip.open(fq,'rt')
        else:
            f = open(fq,'r')
        for line in f:
            if c == 1:
                nraw += 1
                stripped = line.strip()
                ktotal = int(len(stripped)/k)
                # Don't allow Ns in read
                # Don't allow reads < k
                rev = rev_comp(stripped)
                for tmpseq in [stripped, rev]: 
                    kcount = 0
                    if "N" not in tmpseq and len(tmpseq) >= k:
                        for i in range(0, len(tmpseq)-k+1, k):
                            if tmpseq[i:i+k] in dbkmers:
                                kcount += 1
                                if kcount/ktotal > threshold:
                                    reads.append(tmpseq)
                                    # As soon as our threshold is exceeded, break
                                    break
            c += 1                  
            if c == 4:
                c = 0
        f.close()
    return reads, nraw

def parse_fasta_uniq(fasta, filter_Ns=True):
    """ Gets unique sequences from a fasta, with filtering of Ns"""
    tmph = ""
    tmps = ""
    hs = []
    ss = []
    sseen = set()
    with open(fasta) as f:
        for li, line in enumerate(f):
            l = line.strip()
            if not l.startswith(">"):
                l = l.upper()
            if len(l) == 0:
                continue
            elif l[0] == ">":
                if tmps not in sseen and li > 0:
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
        return random.sample(reads, n)


