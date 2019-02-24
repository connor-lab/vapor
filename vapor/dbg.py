""" DBG classes and related objects """
import sys
import numpy as np
np.set_printoptions(threshold=np.nan)
from collections import deque
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
        for i in range(min(len(self.raw_array), len(sr.raw_array))):
            if self.filled_deque_array[i] != sr.filled_deque_array[i]:
                print(i, self.raw_array[i], sr.raw_array[i], self.filled_deque_array[i], sr.filled_deque_array[i], "*")
            else:
                print(i, self.raw_array[i], sr.raw_array[i], self.filled_deque_array[i], sr.filled_deque_array[i])

class wDBG():
    """ Basic DBG with associated edge weights """
    def __init__(self, strings, k):
        """ Initialized with strings, k, and reference kmers """
        # Only explicitly store edges
        self.edges = {}
        self.k = k
        self._build(strings)
        self.caching = True
        self.path_cache = {}
        self.max_trim_size = self.k+1

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

    def cull_low(self, perc=5):
        # Provide a percentile p;
        # Cull any kmers below p
        vals = np.array(list(self.edges.values()))
        percentile = np.percentile(vals, perc)
        keyvals = [(k,v) for k,v in self.edges.items()]
        for key, val in keyvals:
            if val <= percentile:
                del self.edges[key]            

    def mask_against_bridge(self, query, bridge, gapl):
        mask = []
        for i in range(len(query)):
            if query[i] != bridge[i]:
                mask.append(gapl+i)
        return mask

    def extend_bridge(self, kmer, n, direction=1):
        # First check the cache
        if self.caching == True:
            if (kmer, n) in self.path_cache:
                return self.path_cache[(kmer, n)]
        string = kmer
        scorearr = np.zeros(n)
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
        if self.caching == True:
            self.path_cache[(kmer, n)] = (string, scorearr)
        return string, scorearr

    def get_raw_weight_array(self, kmers):
        array = np.zeros(len(kmers))
        for ki, kmer in enumerate(kmers):
            if kmer in self.edges:
                array[ki] = self.edges[kmer]
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
                    gaps.append([gapl, gapr])
 
                in_gap = False
            else:
                if in_gap == False:
                    gapl = ki
                in_gap = True
        if in_gap == True:
            gapr = len(array)
            gaps.append([gapl, gapr])
        return gaps    

    def deque_score_bases(self, array):
        """ For each contiguous stretch of kmers that are present
            uses a deque to find local maxima
            to approximate coverage of a base """
        local_maxima = np.zeros(len(array))
        deq = deque(maxlen=self.k)
        deq.append(0)
        for ki in range(1, len(array)):
            local_maxima[ki-1] = array[deq[0]]
            while deq and deq[0] <= ki-self.k:
                deq.popleft()
            while deq and array[ki] >= array[deq[-1]]:
                deq.pop()
            deq.append(ki)
        local_maxima[-1] = array[deq[0]]
        return local_maxima            

    def get_suboptimal_branches(self, kmers):
        suboptimal_branches = set()
        bases = list("ATCG")
        for ki, kmer in enumerate(kmers):
            if kmer in self.edges:
                score = self.edges[kmer]
                alts = [kmer[:-1] + b for b in bases]
                altscore = max([self.edges[amer] for amer in alts if amer in self.edges])
                if score < altscore:
                    suboptimal_branches.add(ki)
        return suboptimal_branches 

    def expand_gaps(self, gaps, suboptimal_branches, max_gapr):        
        """ 
        Acts in place to expand gaps to suboptimal branch positions
        takes the must distant sub branch within self.max_trim_size
         """
        for gapi, gap in enumerate(gaps):
            gapl, gapr = gap
            for li in range(max(0, gapl-self.max_trim_size), gapl):
                if li in suboptimal_branches:
                    gaps[gapi][0] = li
                    break
            for ri in range(min(gapr + self.max_trim_size, max_gapr), gapr, -1):
                if ri in suboptimal_branches:
                    gaps[gapi][1] = ri
                    break

    def query(self, kmers, seqsh, min_kmer_prop, debug=False):
        sr = SearchResult()
        # First obtain the raw weight array for kmers of a sequence
        raw_weight_array = self.get_raw_weight_array(kmers)
        kmer_cov = np.count_nonzero(raw_weight_array)/len(raw_weight_array)
        # Next trim the raw weight array
        if debug==True:
            sr.raw_array = raw_weight_array
        if kmer_cov > min_kmer_prop:
            # Get the gaps
            gaps = self.get_weight_array_gaps(raw_weight_array)
            # Get the suboptimal branches
            sub = self.get_suboptimal_branches(kmers)
            self.expand_gaps(gaps, sub, len(kmers))
            # Copy the raw weight array to modify
            filled_weight_array = [r for r in raw_weight_array]
#            print(gaps)
            all_masks = []
            for gapl, gapr in gaps:
                if gapl != 0 and gapr != len(kmers):
                    gapstring = kmers2str(kmers[gapl:gapr])[self.k-1:]
                    bridge, bridge_scores = self.extend_bridge(kmers[gapl-1], gapr-gapl)
                    bridge_rev, bridge_scores_rev = self.extend_bridge(kmers[gapr], gapr-gapl, -1)
                    gapstring_rev = kmers2str(kmers[gapl:gapr])[:-self.k+1]
                    if sum(bridge_scores_rev) > sum(bridge_scores):
                        mask = self.mask_against_bridge(gapstring_rev, bridge_rev, gapl)
                        filled_weight_array[gapl:gapr] = bridge_scores
                    else:
                        mask = self.mask_against_bridge(gapstring, bridge, gapl)
                        filled_weight_array[gapl:gapr] = bridge_scores

                elif gapr != len(kmers) and gapl == 0:
                    gapstring = kmers2str(kmers[gapl:gapr])[self.k-1:]
                    bridge, bridge_scores = self.extend_bridge(kmers[gapr], gapr-gapl, -1)
                    mask = self.mask_against_bridge(gapstring, bridge, gapl)
                    filled_weight_array[gapl:gapr] = bridge_scores

                elif gapl > 0 and gapr == len(kmers):
                    gapstring = kmers2str(kmers[gapl:gapr])[self.k-1:]
                    bridge, bridge_scores = self.extend_bridge(kmers[gapl-1], gapr-gapl)
                    mask = self.mask_against_bridge(gapstring, bridge, gapl)
                    filled_weight_array[gapl:gapr] = bridge_scores
                for mi in mask:
                    assert mi in range(gapl, gapr)
                all_masks += mask
            # Add the last k - 1 bases
            filled_weight_array = np.concatenate((filled_weight_array, np.zeros(self.k-1)))
            # Deque score
            filled_deque_array = self.deque_score_bases(filled_weight_array)
            for maski in all_masks:
#                print(maski, raw_weight_array[maski])
#                assert raw_weight_array[maski] == 0
                filled_deque_array[maski] = 0
        else:
            sr.filled_deque_array = raw_weight_array
            sr.est_pid = -1
            sr.score = -1
            return sr
        sr.filled_deque_array = filled_deque_array
        nonzeros = [i for i in filled_deque_array if i > 0]
        est_pid = len(nonzeros)/len(filled_deque_array)
        sr.est_pid = est_pid
        score = sum(nonzeros)
        sr.score = score
        return sr

    def classify(self, seqs, seqsh, min_kmer_prop, debug_query=None):
        scores = []
        for si, seq in enumerate(seqs):
            kmers = [seq[i:i+self.k] for i in range(len(seq)-self.k+1)]
            sr = self.query(kmers, seqsh[si], min_kmer_prop, (debug_query != False))
            sr.header = seqsh[si]
            sr.index = si
            scores.append(sr)

        # Sort the results
        results = sorted(scores, key = lambda x: (x.score, x.est_pid), reverse=True)

        if debug_query != None:
            for hi, h in enumerate(seqsh):
                if h == debug_query:
                    seq = seqs[hi]
                    kmers = [seq[i:i+self.k] for i in range(len(seq)-self.k+1)]
                    sr = self.query(kmers, seqsh[hi], min_kmer_prop, True)
                    sr.compare(results[0])
                    print(sr.score)

        return results
