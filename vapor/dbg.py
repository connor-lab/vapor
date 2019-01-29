""" DBG classes and related objects """
import sys
import numpy as np
from collections import deque

class wDBG():
    """ Basic DBG with associated edge weights """
    def __init__(self, strings, k, ref_kmers):
        """ Initialized with strings, k, and reference kmers """
        # Only explicitly store edges
        self.edges = {}
        self.k = k
        self.start_positions = set()
        self._build(strings, ref_kmers)

    def _build(self, strings, ref_kmers):
        # Builds by taking a set of strings (reads), reference kmers
        # Any kmer not present in references is discarded
        sys.stderr.write("Building wDBG\n")
        newkmerc = 0
        for si, string in enumerate(strings):
            sys.stderr.write(str(si) + "      \r")
            kmers = [string[i:i+self.k] for i in range(len(string)-self.k+1)]
            newkmerc = 0
            for kmer in kmers:
                if kmer in ref_kmers:
                    if kmer in self.edges:
                        self.edges[kmer] += 1
                    else:
                        newkmerc += 1
                        self.edges[kmer] = 1

    def brute_score_bases(self, kmers):
        # Takes ordered kmers
        scores = np.zeros(len(kmers))
        for ki, kmer in enumerate(kmers):
            maxw = 0
#            print(max(ki-self.k, 0), ki, self.edges[kmer])
            c = 0
            for kj in range(max(ki-self.k+1, 0), ki+1):
                c += 1
                kmerj = kmers[kj]
                if kmerj in self.edges:
                    maxw = max([maxw, self.edges[kmerj]])
    
            assert c <= self.k
            scores[ki] = maxw
#            assert maxw >= self.edges[kmer]
        return scores

    def deque_score_bases(self, kmers):
        """ For each contiguous stretch of kmers that are present
            uses a deque to find local maxima
            to approximate coverage of a base """
        local_minima = []
        array = []
        for kmer in kmers:
            if kmer in self.edges:
                array.append(self.edges[kmer])
            else:
                array.append(0)
        # local window alg using deque
        # window size needs to be a function of k
        deq = deque()
        deq.append(0) 
        for ki in range(1, len(array)):
            # append the prev
            local_minima.append(array[deq[0]])
            while deq and deq[0] <= ki-self.k:
                deq.popleft()
            while deq and array[ki] >= array[deq[-1]]:
                deq.pop()
            deq.append(ki)
        local_minima.append(array[deq[0]])
        return local_minima            

    def query(self, seq):
        kmers = [seq[ki:ki+self.k] for ki in range(len(seq)-self.k+1)]
        # get brute score for debugging
        deq_scores = self.deque_score_bases(kmers)            
        score = sum(deq_scores)
        return score, deq_scores

    def classify(self, seqs, seqsh, min_missing_kmers=100.):
        # ALLOW FOR PARTIAL MATCHES, WHEN WALKING
        # IF A KMER IS NOT FOUND, WE CAN CHECK HAMMING DISTANCES OF KMERS?
        # ATTEMPT TO `BRIDGE THE GAP' (FOR LATER)
        scores = []
        for si, seq in enumerate(seqs):
            print(si)
            scores.append(self.query(seq))
    
        maxs = max(scores, key = lambda x:x[0])[0]
        maxcls = [si for si in range(len(seqs)) if scores[si][0] == maxs]
        inds = np.argsort([s[0] for s in scores])[::-1]
        return maxs, maxcls


