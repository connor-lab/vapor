""" DBG classes and related objects """
import sys
import numpy as np
from collections import deque
from Bio import pairwise2
from vapor.vaporfunc import *

class SearchResult():
    def __init__(self):
        self.i = -1
        self.prop_kmers = -1
        self.gap_positions = []
        self.raw_array = []
        self.filled_array = []
        self.filled_deque_array = []
        self.raw_score = -1
        self.score = -1
    def compare(self, sr):
        for i in range(len(self.filled_deque_array)):
            if self.filled_deque_array[i] != sr.filled_deque_array[i]:
                print(i, self.filled_deque_array[i], sr.filled_deque_array[i], "*")
            else:
                print(i, self.filled_deque_array[i])

class wDBG():
    """ Basic DBG with associated edge weights """
    def __init__(self, strings, k):
        """ Initialized with strings, k, and reference kmers """
        # Only explicitly store edges
        self.edges = {}
        self.k = k
        self._build(strings)

    def _build(self, strings):
        # Builds by taking a set of strings (reads), reference kmers
        # Any kmer not present in references is discarded
        for si, string in enumerate(strings):
            kmers = [string[i:i+self.k] for i in range(len(string)-self.k+1)]
            for kmer in kmers:
                if kmer in self.edges:
                    self.edges[kmer] += 1
                else:
                    self.edges[kmer] = 1

    def cull_low(self, perc=5):
        # Provide a percentile p;
        # Cull any kmers below p
        vals = np.array(list(self.edges.values()))
        percentile = np.percentile(vals, perc)
        keyvals = [(k,v) for k,v in self.edges.items()]
        for key, val in keyvals:
            if val <= percentile:
                del self.edges[key]            

    def score_against_bridge(self, query, bridge, bridge_scores):
        # hamming distance, for now
        scores = []
        for i in range(len(query)):
            if query[i] == bridge[i]:
                scores.append(bridge_scores[i])
            else:
                scores.append(-1)
        return scores

    def extend_bridge(self, kmer, n, direction=1):
        string = kmer
        scorearr = np.empty(n)
        scorearr.fill(-1)
        if direction == 1:
            si = 0
        else:
            si = -1
        while len(string) < n+self.k:
            if direction == 1:
                poss_edges = [string[-self.k+1:] + b for b in "ATCG"]
            elif direction == -1:
                poss_edges = [b + string[:self.k-1] for b in "ATCG"]
            max_score = 0
            max_base = None
            for pe in poss_edges:
                if pe in self.edges:
                    tmpscore = self.edges[pe]
                    if tmpscore > max_score:
                        max_score = tmpscore
                        if direction == 1:
                            max_base = pe[-1]
                        elif direction == -1:
                            max_base = pe[0]
            if max_base != None:
                scorearr[si] = max_score
                if direction == 1:
                    si += 1
                    string += max_base
                elif direction == -1:
                    si -= 1
                    string = max_base + string
            else:
                break
        if direction == 1:
            string = string[len(kmer):]
            string += "X"*(n-len(string))
        elif direction == -1:
            string = string[:-len(kmer)]
            string = "X"*(n-len(string)) + string
        return string, scorearr

    def get_raw_weight_array(self, kmers):
        array = np.zeros(len(kmers))
        for ki, kmer in enumerate(kmers):
            if kmer in self.edges:
                array[ki] = self.edges[kmer]
            else:
                array[ki] = 0
        return array
    
    def get_weight_array_gaps(self, array):
        in_gap = False
        gapl = -1
        gapr = -1
        gaps = []
        for ki, val in enumerate(array):
            if val != 0:
                if in_gap == True:
                    gapr = ki
                    gaps.append((gapl, gapr))
 
                in_gap = False
            else:
                if in_gap == False:
                    gapl = ki
                in_gap = True
        if in_gap == True:
            gapr = len(array)
            gaps.append((gapl, gapr))
        return gaps    

    def deque_score_bases(self, array):
        """ For each contiguous stretch of kmers that are present
            uses a deque to find local maxima
            to approximate coverage of a base """
        local_maxima = np.zeros(len(array))
        deq = deque(maxlen=self.k)
        deq.append(0)
        for ki in range(1, len(array)):
            if array[ki-1] == -1:
                local_maxima[ki-1] = 0
            else:
                local_maxima[ki-1] = array[deq[0]]
            while deq and deq[0] <= ki-self.k:
                deq.popleft()
            while deq and array[ki] >= array[deq[-1]]:
                deq.pop()
            deq.append(ki)
        if array[-1] == -1:
            local_maxima[-1] = 0
        else:
            local_maxima[-1] = array[deq[0]]
        return local_maxima            

    def query(self, kmers, seqsh, min_kmer_prop):
        sr = SearchResult()
        raw_weight_array = self.get_raw_weight_array(kmers)
        kmer_cov = np.count_nonzero(raw_weight_array)/len(raw_weight_array)
        if kmer_cov > min_kmer_prop:
            filled_weight_array = raw_weight_array
            gaps = self.get_weight_array_gaps(raw_weight_array)
            for gapl, gapr in gaps:
                if gapl != 0 and gapr != len(kmers):
                    gapstring = kmers2str(kmers[gapl:gapr])[self.k-1:]
                    bridge, bridge_scores = self.extend_bridge(kmers[gapl-1], gapr-gapl)
                    bridge_rev, bridge_scores_rev = self.extend_bridge(kmers[gapr], gapr-gapl, -1)
                    gapstring_rev = kmers2str(kmers[gapl:gapr])[:-self.k+1]
                    if sum(bridge_scores_rev) > sum(bridge_scores):
                        extra_scores = self.score_against_bridge(gapstring_rev, bridge_rev, bridge_scores_rev)
                        filled_weight_array[gapl:gapr] = extra_scores
                    else:
                        extra_scores = self.score_against_bridge(gapstring, bridge, bridge_scores)
                        filled_weight_array[gapl:gapr] = extra_scores

                elif gapr != len(kmers) and gapl == 0:
                    gapstring = kmers2str(kmers[gapl:gapr])[self.k-1:]
                    bridge, bridge_scores = self.extend_bridge(kmers[gapr], gapr-gapl, -1)
                    extra_scores = self.score_against_bridge(gapstring, bridge, bridge_scores)
                    filled_weight_array[gapl:gapr] = extra_scores
                    # only one direction
                elif gapl > 0 and gapr == len(kmers):
                    # only one direction
                    for kmer in kmers[gapl:gapr]:
                        assert kmer not in self.edges
                    gapstring = kmers2str(kmers[gapl:gapr])[self.k-1:]
                    bridge, bridge_scores = self.extend_bridge(kmers[gapl-1], gapr-gapl)
                    extra_scores = self.score_against_bridge(gapstring, bridge, bridge_scores)
                    filled_weight_array[gapl:gapr] = extra_scores
            filled_deque_array = self.deque_score_bases(filled_weight_array)
        else:
            filled_deque_array = raw_weight_array
            sr.est_pid = -1
            sr.score = -1
            return sr
        sr.filled_deque_array = filled_deque_array
        np.set_printoptions(threshold=np.nan)
        nonzeros = [i for i in filled_deque_array if i > 0]
        est_pid = len(nonzeros)/len(filled_deque_array)
        sr.est_pid = est_pid
        score = sum(nonzeros)
        sr.score = score
        return sr

    def classify(self, kmersets, seqsh, min_kmer_prop):
        scores = []
        debug_srs = []
        for si, kmers in enumerate(kmersets):
            sr = self.query(kmers, seqsh[si], min_kmer_prop)
            debug_srs.append(sr)
            scores.append((sr.est_pid, sr.score))
#        debug_srs[0].compare(debug_srs[1])

        # Sort and return results
        results = []
        argsort = lambda x : sorted(range(len(x)), key=x.__getitem__)

        inds = argsort(scores)
        for ind in inds:
            results.append((ind, scores[ind][0], scores[ind][1]))
        return results


