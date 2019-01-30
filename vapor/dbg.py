""" DBG classes and related objects """
import sys
import numpy as np
from collections import deque
from Bio import pairwise2

class wDBG():
    """ Basic DBG with associated edge weights """
    def __init__(self, strings, k, ref_kmers):
        """ Initialized with strings, k, and reference kmers """
        # Only explicitly store edges
        self.edges = {}
        # Path cache to prevent re-computing
        # Of paths between nodes
        self.path_cache = {}
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

    def search_gap(self, kmers, gapl, gapr):
        # heuristically build an optimal path in the wdbg
        if gapl > 0:
            stringo = kmers[gapl-1]
        else:
           stringo = kmers[gapr]
        assert stringo in self.edges
        for i in range(gapr-gapl):
#            print("stringo", stringo)
            if gapl > 0:
                poss_edges = [stringo[-self.k+1:] + b for b in "ATCG"]
            else:
                poss_edges = [b + stringo[:self.k-1] for b in "ATCG"]
            assert len(poss_edges[0]) == self.k
            max_score = 0
            max_base = None
            for pe in poss_edges:
                if pe in self.edges:
#                    print("pe", pe)
                    tmpscore = self.edges[pe]
                    if tmpscore > max_score:
                        if gapl == 0:
                            max_base = pe[0]
                        else:
                            max_base = pe[-1]
#            print(stringo, max_base, pe)
            if max_base != None:
                if gapl > 0:
                    stringo += max_base
                else:
                    stringo = max_base + stringo
            else:
                break
    
        extra_scores = []
        # compare that to the kmer string
        stringk = kmers[gapl] + "".join([kmer[-1] for kmer in kmers[gapl+1:gapr]])
        aln = pairwise2.align.globalxx(stringk, stringo, one_alignment_only=True)[0]
#        print(aln, len(aln[1].replace("-","")))
#        print(aln[0])
#        print(aln[1])
        # finally need to pad
#        if 
        for i in range(len(aln[0])-self.k+1):
            bi = aln[0][i]
            bj = aln[1][i]
            if bi == bj:
                kmerj = aln[1][i]
                ki = i+1
                while len(kmerj) < self.k:
                    if aln[1][ki] != "-":
                        kmerj += aln[1][ki]
                    ki += 1                
                extra_scores.append(self.edges[kmerj])
        # Now pad the extra scores (this is in case the bridge is not as long as the query gap)
        if gapl == 0:
            extra_scores = [0]*(gapr-gapl-len(extra_scores)) + extra_scores
        else:
             extra_scores += [0]*(gapr-gapl-len(extra_scores)) 
        return extra_scores
    
    def get_weight_array(self, kmers):
        array = []
        array2 = []
        in_gap = False
        gapl = -1
        gapr = -1
        for ki, kmer in enumerate(kmers):
            if kmer in self.edges:
                if in_gap == True:
                    gapr = ki
                    search = self.search_gap(kmers, gapl, gapr)
#                    print(array)
                    array[gapl:gapr] = search
#                    print(search)
#                    print(array)
#                    print(len(search), gapr-gapl)
#                    print()
                    assert len(search) == gapr-gapl
#                    assert gapr-gapl != self.k
   
                in_gap = False
                array.append(self.edges[kmer])
            else:
                if in_gap == False:
                    gapl = ki
                in_gap = True
                array.append(0)
        return array        

    def deque_score_bases(self, array):
        """ For each contiguous stretch of kmers that are present
            uses a deque to find local maxima
            to approximate coverage of a base """
        local_minima = []
               # local window alg using deque
        # window size needs to be a function of k
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

    def query(self, kmers, name):
#        kmers = [seq[ki:ki+self.k] for ki in range(len(seq)-self.k+1)]
        # get brute score for debugging
        weight_array = self.get_weight_array(kmers)
#        if "AIM56458" in name or "Chile" in name:
#            print(name, weight_array, sum(weight_array), len(weight_array))
        deq_scores = self.deque_score_bases(weight_array)            
        score = sum(deq_scores)
        return score

    def classify(self, kmersets, seqsh, n_threads):
        # FILTER FIRST

        # ALLOW FOR PARTIAL MATCHES, WHEN WALKING
        # IF A KMER IS NOT FOUND, WE CAN CHECK HAMMING DISTANCES OF KMERS?
        # ATTEMPT TO `BRIDGE THE GAP' (FOR LATER)
        scores = []
        for si, kmers in enumerate(kmersets):
            scores.append(self.query(kmers, seqsh[si]))

        inds = np.argsort(scores)
        for i in inds:
            print(scores[i], seqsh[i])
        maxs = max(scores, key = lambda x:x)
        maxcls = [si for si in range(len(seqsh)) if scores[si] == maxs]

        # resolve ties by completion
        return maxs, maxcls


