""" DBG classes and related objects """
import sys
from collections import deque
import numpy as np
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
    def print(self):
        print("i", self.i)
        print("kmer prop", self.prop_kmers)
        print("gap_positions", self.gap_positions)
        print("rawarray", self.raw_array)
        print("filledarray", self.filled_array)
        print("score", self.score)
        print("raw_score", self.raw_score)
    def compare(self, sr):
        gap_posls = [i[0] for i in self.gap_positions] + [i[0] for i in sr.gap_positions]
        gap_posrs = [i[1] for i in self.gap_positions] + [i[1] for i in sr.gap_positions]
        for i in range(len(self.raw_array)):
            if i in gap_posls:
                print("/====/")
            if i in gap_posrs:
                print("/====")
#            if (self.filled_deque_array[i] < sr.filled_deque_array[i]):
            print(i, self.raw_array[i], self.filled_array[i], self.filled_deque_array[i], sr.raw_array[i], sr.filled_array[i], sr.filled_deque_array[i],"*")
#            else:
#                print(i, self.raw_array[i], self.filled_array[i], self.filled_deque_array[i], sr.raw_array[i], sr.filled_array[i], sr.filled_deque_array[i])
#                pass


        print(self.score, sr.score)

class wDBG():
    """ Basic DBG with associated edge weights """
    def __init__(self, strings, k):
        """ Initialized with strings, k, and reference kmers """
        # Only explicitly store edges
        self.edges = {}
        self.k = k
        self._build(strings)
        self._cull_low()

    def _build(self, strings):
        # Builds by taking a set of strings (reads), reference kmers
        # Any kmer not present in references is discarded
        sys.stderr.write("Building wDBG\n")
        newkmerc = 0
        for si, string in enumerate(strings):
            sys.stderr.write(str(si) + "      \r")
            kmers = [string[i:i+self.k] for i in range(len(string)-self.k+1)]
            newkmerc = 0
            for kmer in kmers:
#                if kmer in ref_kmers:
                if kmer in self.edges:
                    self.edges[kmer] += 1
                else:
                    newkmerc += 1
                    self.edges[kmer] = 1

    def _cull_low(self, threshold=5):
        keyvals = [(i,q) for i,q in self.edges.items()]
        for key, val in keyvals:
            if val <= threshold:
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
        scores = []
        while len(string) < n+self.k:
            if direction == 1:
                poss_edges = [string[-self.k+1:] + b for b in "ATCG"]
            elif direction == -1:
                poss_edges = [b + string[:self.k-1] for b in "ATCG"]
            max_score = 0
            max_base = None
            for pe in poss_edges:
                assert len(pe) == self.k
                if pe in self.edges:
                    tmpscore = self.edges[pe]
                    if tmpscore > max_score:
                        max_score = tmpscore
                        if direction == 1:
                            max_base = pe[-1]
                        elif direction == -1:
                            max_base = pe[0]
            if max_base != None:
                scores.append(max_score)
                if direction == 1:
                    string += max_base
                elif direction == -1:
                    string = max_base + string
            else:
                break
        if direction == 1:
            string = string[len(kmer):]
            scores += [-1 for i in range(n-len(string))]
            string += "X"*(n-len(string))
        elif direction == -1:
            string = string[:-len(kmer)]
            scores = [-1 for i in range(n-len(string))] + scores
            string = "X"*(n-len(string)) + string

        # return early
        assert len(string) == n
        assert len(scores) == n
        return string, scores


    def get_raw_weight_array(self, kmers):
        array = []
        for ki, kmer in enumerate(kmers):
            if kmer in self.edges:
                array.append(self.edges[kmer])
            else:
                array.append(0)    
        # for each gap starting on the left, should change to 0, not -1?
        # 0
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
                array.append(-1)
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
#                print("search",search)
                if search != False:
#                    print(len(search), gapr-gapl, gapl, gapr)
                    assert len(search) == gapr-gapl
                    array[gapl:gapr] = search
            return array
        else:
            return False

    def deque_score_bases(self, array):
        """ For each contiguous stretch of kmers that are present
            uses a deque to find local maxima
            to approximate coverage of a base """
        local_maxima = []
        deq = deque()
        deq.append(0)
        # MASK; special flags for array:
            # -1 gives no information, also not counted
            # (-1,x) gives information, but is not counted
            # we note down these positions
            # and reset them
 #       mask_inds = []
#        print(array)
#        for ai in range(len(array)):
#            w = array[ai]
#            if type(w) == int:
 #               if w == -1:
 #                   mask_inds.append(ai)
 #           else:
 #               mask_inds.append(ai)
 #               array[ai] = w[1]

#        print(mask_inds)
#        print(array)
        for ki in range(1, len(array)):
            # append the prev
            if array[ki-1] == -1:
                local_maxima.append(-1)
            else:
                local_maxima.append(array[deq[0]])
            while deq and deq[0] <= ki-self.k:
                deq.popleft()
            while deq and array[ki] >= array[deq[-1]]:
                deq.pop()
            deq.append(ki)
        if array[-1] == -1:
            local_maxima.append(-1)
        else:
            local_maxima.append(array[deq[0]])

        # Finally use the mask indices
#        for ai in mask_inds:
#            local_maxima[ai] = 0
        return local_maxima            

    def query(self, kmers, seqsh, min_kmer_prop=0.8, max_gap_size=500):
        sr = SearchResult()
        raw_weight_array = self.get_raw_weight_array(kmers)
        gaps = self.get_weight_array_gaps(raw_weight_array)
        filled_weight_array = [r for r in raw_weight_array]
        kmer_cov = sum([min(1,i) for i in range(len(raw_weight_array))])
        if kmer_cov > min_kmer_prop:
            for gapl, gapr in gaps:
                if gapl != 0 and gapr != len(kmers):
    #                print("gap_positions", gapl, gapr, gapr-gapl)
                    for kmer in kmers[gapl:gapr]:
                        assert kmer not in self.edges
                    gapstring = kmers2str(kmers[gapl:gapr])[self.k-1:]
                    assert kmers[gapl-1] in self.edges
                    assert kmers[gapr] in self.edges
                    bridge, bridge_scores = self.extend_bridge(kmers[gapl-1], gapr-gapl)
                    assert gapr-gapl == len(gapstring)
                    extra_scores = self.score_against_bridge(gapstring, bridge, bridge_scores)
                    assert len(bridge) == gapr-gapl
                    assert len(bridge_scores) == gapr-gapl
    #                print(gapl, gapr, extra_scores)
                    filled_weight_array[gapl:gapr] = extra_scores
                elif gapr != len(kmers) and gapl == 0:
                    for kmer in kmers[gapl:gapr]:
                        assert kmer not in self.edges
                    gapstring = kmers2str(kmers[gapl:gapr])[self.k-1:]
                    assert kmers[gapr] in self.edges
                    bridge, bridge_scores = self.extend_bridge(kmers[gapr], gapr-gapl, -1)
                    assert gapr-gapl == len(gapstring)
                    extra_scores = self.score_against_bridge(gapstring, bridge, bridge_scores)
                    assert len(bridge) == gapr-gapl
                    assert len(bridge_scores) == gapr-gapl
    #                print(gapl, gapr, extra_scores)
                    filled_weight_array[gapl:gapr] = extra_scores
                    # only one direction
                elif gapl > 0 and gapr == len(kmers):
                    # only one direction
                    for kmer in kmers[gapl:gapr]:
                        assert kmer not in self.edges
                    gapstring = kmers2str(kmers[gapl:gapr])[self.k-1:]
                    assert kmers[gapl-1] in self.edges
                    bridge, bridge_scores = self.extend_bridge(kmers[gapl-1], gapr-gapl)
                    assert gapr-gapl == len(gapstring)
                    extra_scores = self.score_against_bridge(gapstring, bridge, bridge_scores)
                    assert len(bridge) == gapr-gapl
                    assert len(bridge_scores) == gapr-gapl
    #                print(gapl, gapr, extra_scores)
                    filled_weight_array[gapl:gapr] = extra_scores

        raw_deque_array = self.deque_score_bases(raw_weight_array)
        filled_deque_array = self.deque_score_bases(filled_weight_array)

#        print(seqsh)
#        print(raw_weight_array)
#        print(filled_weight_array)
#        print(filled_deque_array)
#        print(raw_deque_array)
#        print([z for z in zip(raw_weight_array, filled_weight_array, filled_deque_array)])
        assert sum([i for i in filled_weight_array if i > 0]) >= sum(raw_weight_array)
        assert sum(filled_deque_array) >= sum(filled_weight_array)
        
        sr.gap_positions = gaps
        sr.filled_array = filled_weight_array
        sr.raw_array = raw_weight_array
        sr.filled_deque_array = filled_deque_array
#            sr.filled_deque_array = filled_deque_array
#            sr.raw_deque_array = raw_deque_array
        sr.prop_kmers = len([i for i in raw_deque_array if i > 0])/len(kmers)            
#        print(raw_weight_array)
#        print(sr.prop_kmers, len([i for i in raw_weight_array if i > 0])/len(kmers))

#        print(filled_deque_array)
        est_pid = sum([1 for i in filled_deque_array if i > 0])/len(filled_deque_array)
        sr.est_pid = est_pid
        score = sum([i for i in filled_deque_array if i > 0])
        sr.score = score
#            print(seqsh)
#            sr.print()
        return sr

    def classify(self, kmersets, seqsh, min_kmer_prop):
        scores = []
        debug_srs = []
        for si, kmers in enumerate(kmersets):
#            if "ADK87309" in seqsh[si] or "ACT34407" in seqsh[si]:
#            print(seqsh[si])
            sr = self.query(kmers, seqsh[si], min_kmer_prop)
            scores.append((sr.est_pid, sr.score))
#            debug_srs.append(sr)
#            else:
#                scores.append((-1,-1))
#        debug_srs[0].compare(debug_srs[1])
        # Sort and return results
        results = []
#        print(scores)
        argsort = lambda x : sorted(range(len(x)), key=x.__getitem__)

        inds = argsort(scores)
        for ind in inds:
            results.append((ind, scores[ind][0], scores[ind][1]))
        return results


