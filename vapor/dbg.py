""" DBG classes and related objects """
import sys
from collections import deque
import numpy as np
from Bio import pairwise2

class wDBG():
    """ Basic DBG with associated edge weights """
    def __init__(self, strings, k, ref_kmers):
        """ Initialized with strings, k, and reference kmers """
        # Only explicitly store edges
        self.edges = {}
        self.k = k
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

    def search_gap(self, kmers, gapl, gapr):
        # heuristically build an optimal path in the wdbg
        if gapl > 0:
            stringo = kmers[gapl-1]
        else:
           stringo = kmers[gapr]
        assert stringo in self.edges
        for i in range(gapr-gapl):
            if gapl > 0:
                poss_edges = [stringo[-self.k+1:] + b for b in "ATCG"]
            else:
                poss_edges = [b + stringo[:self.k-1] for b in "ATCG"]
            assert len(poss_edges[0]) == self.k
            max_score = 0
            max_base = None
            for pe in poss_edges:
                if pe in self.edges:
                    tmpscore = self.edges[pe]
                    if tmpscore > max_score:
                        if gapl == 0:
                            max_base = pe[0]
                        else:
                            max_base = pe[-1]
            if max_base != None:
                if gapl > 0:
                    stringo += max_base
                else:
                    stringo = max_base + stringo
            else:
                break
    
        extra_scores = []
        stringk = kmers[gapl] + "".join([kmer[-1] for kmer in kmers[gapl+1:gapr]])
        aln = pairwise2.align.globalxx(stringk, stringo, one_alignment_only=True)[0]
        for i in range(len(aln[0])-self.k+1):
            bi = aln[0][i]
            bj = aln[1][i]
            if bi == bj:
                kmerj = ""
                ki = i
                while len(kmerj) < self.k and ki < len(aln[1]):
                    if aln[1][ki] != "-":
                        kmerj += aln[1][ki]
                    ki += 1 
                if len(kmerj) == self.k:
                    extra_scores.append(self.edges[kmerj])
        # Now pad the extra scores (this is in case the bridge is not as long as the query gap)
        if gapl == 0:
            extra_scores = [0]*(gapr-gapl-len(extra_scores)) + extra_scores
        else:
             extra_scores += [0]*(gapr-gapl-len(extra_scores)) 
        return extra_scores
    
    def get_weight_array(self, kmers, min_kmer_prop):
        array = []
        in_gap = False
        gapl = -1
        gapr = -1
        gaps = []
        for ki, kmer in enumerate(kmers):
            if kmer in self.edges:
                if in_gap == True:
                    gapr = ki
                    gaps.append((gapl, gapr))
  
                in_gap = False
                array.append(self.edges[kmer])
            else:
                if in_gap == False:
                    gapl = ki
                in_gap = True
                array.append(0)
        if in_gap == True:
            gapr = len(kmers)
            gaps.append((gapl, gapr))
        # Only fill in gaps if it is considered a close match
        gapsum = 0
        for gapl, gapr in gaps:
            gapsum += gapr-gapl
        if gapsum/len(kmers) < 1-min_kmer_prop:
            for gapl, gapr in gaps:
                search = self.search_gap(kmers, gapl, gapr)
                array[gapl:gapr] = search
            return array        
        else:
            return False

    def deque_score_bases(self, array):
        """ For each contiguous stretch of kmers that are present
            uses a deque to find local maxima
            to approximate coverage of a base """
        local_minima = []
        deq = deque()
        deq.append(0) 
        for ki in range(1, len(array)):
            # append the prev
            if array[ki] > 0:
                local_minima.append(array[deq[0]])
            else:
                local_minima.append(0)
            while deq and deq[0] <= ki-self.k:
                deq.popleft()
            while deq and array[ki] >= array[deq[-1]]:
                deq.pop()
            deq.append(ki)
        local_minima.append(array[deq[0]])
        return local_minima            

    def query(self, kmers, min_kmer_prop=0.9):
        weight_array = self.get_weight_array(kmers, min_kmer_prop)
        if weight_array != False:
            deq_scores = self.deque_score_bases(weight_array)            
            score = sum(deq_scores)
            return score
        else:
            return -1

    def classify(self, kmersets, min_kmer_prop):
        scores = []
        for si, kmers in enumerate(kmersets):
            query = self.query(kmers, min_kmer_prop)
            scores.append(query)
        # Sort and return results
        results = []
        inds = np.argsort(scores)
        for ind in inds:
            results.append((ind, scores[ind]))
        return results


