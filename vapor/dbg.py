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
            poss_edges = [string[-self.k+1:] + b for b in "ATCG"]
            max_score = 0
            max_base = None
            for pe in poss_edges:
                assert len(pe) == self.k
                if pe in self.edges:
                    tmpscore = self.edges[pe]
                    if tmpscore > max_score:
                        max_score = tmpscore
                        max_base = pe[-1]
            if max_base != None:
                scores.append(max_score)
                string += max_base
            else:
                print("break", len(string)-self.k)
                break
        print("finished", len(string))
        string = string[len(kmer):]
        scores += [-1 for i in range(n-len(string))]
        string += "X"*(n-len(string))
        # return early
        print(kmer, n)
        print(string)
        print(scores)
        print(len(string), len(scores), n)
        assert len(string) == n
        assert len(scores) == n
        return string, scores

    def search_gap(self, kmers, gapl, gapr):
        if gapl == 0 and gapr == len(kmers):
            # no kmers matched at all
            return False
        if gapl > 0:
            # BRIDGE EXTENSION
            early_break = False
            # BOUNDARY ASSERTIONS
            assert kmers[gapl-1] in self.edges
            assert kmers[gapl] not in self.edges
            if gapr != len(kmers):
                assert kmers[gapr] in self.edges
            # INITIALIZE STRINGO, STRINGK
            # E.g. ATG CTAT, ATG CCAT
            stringo = kmers[gapl-1]
            to_join = [kmer[-1] for kmer in kmers[gapl:gapr]]
            joined = "".join(to_join)
            stringk = kmers[gapl-1] + joined 
#            print(len(to_join), gapr-gapl, gapr-gapl+self.k, gapr, len(kmers))
            assert len(to_join) == gapr-gapl
            assert gapr <= len(kmers)
            assert len(stringk) == gapr-gapl+self.k
#            print("INIT", gapr, gapl)
#            if gapr == len(kmers):
#                print(kmers[gapl])
#            else:
#                print(kmers[gapl], kmers[gapr])
#            print("stringo", stringo)
#            print("stringk", stringk)
            # EXTEND NOW; GAPR-GAPL bases
            for i in range(gapr-gapl):
                poss_edges = [stringo[-self.k+1:] + b for b in "ATCG"]
                max_score = 0
                max_base = None
                for pe in poss_edges:
                    if pe in self.edges:
                        tmpscore = self.edges[pe]
                        if tmpscore > max_score:
                            max_base = pe[-1]
                if max_base != None:
                    stringo += max_base
                else:
                    if i == 0:
                        # no extension possible, return False
                        return False
                    early_break = True
#                    print("early break")
                    break            
#            print(gapl, gapr)
#            print("BRIDGE")
#            print("stringo", stringo)
#            print("stringk", stringk)
#            print("TRIM1")
            # TRIM 1
            stringo = stringo[1:]
            assert len(stringk) == gapr-gapl+self.k
            stringk = stringk[1:]
#            print(len(stringk), gapr-gapl-1+self.k, gapr, gapl, len(kmers))
            assert len(stringk) == gapr-gapl-1+self.k
            assert len(stringo) <= gapr-gapl-1+self.k
            stringk = stringk[:len(stringo)]
            assert len(stringo) == len(stringk)
#            print("stringo", stringo)
#            print("stringk", stringk)
#            print("SCORE")
            # SCORE
            bridge_scores = [-1 for i in range(gapr-gapl)]
            for i in range(len(stringo)-self.k+1):
                bridge_scores[i] = self.edges[stringo[i:i+self.k]]
#            bridge_scores = self.deque_score_bases(bridge_scores)
#            print(bridge_scores)
#            print("TRIM2")
            # TRIM 2
            stringo = stringo[self.k-1:gapr-gapl]
            stringk = stringk[self.k-1:gapr-gapl]
            if len(stringk) == 0:
                return False
#            print("stringo", stringo)
#            print("stringk", stringk)
            assert len(stringk) < gapr-gapl+self.k-1
            assert len(stringo) < gapr-gapl+self.k-1
            # ADJUST SCORES
            aln = pairwise2.align.globalxx(stringk, stringo, one_alignment_only=True)[0]
#            print("ALIGNMENT")
#            print(aln[0])
#            print(aln[1])
            bi = self.k-1
#            print(len(bridge_scores), len(stringo), len(stringk), gapl, gapr)
        
            for i in range(len(aln[0])):
                ci = aln[0][i]
                cj = aln[1][i]
                if cj != "-":
                    if ci != cj:
                        bridge_scores[bi] = (-1, bridge_scores[bi])
                    bi += 1
            # Now pad the extra scores (this is in case the bridge is not as long as the query gap)
#            print(len(bridge_scores), gapr-gapl, gapr, gapl, len(kmers))
            assert len(bridge_scores) == gapr-gapl
            return bridge_scores

        else:
            # INITIALIZE
            stringo = kmers[gapr]
            stringk = "".join([kmer[0] for kmer in kmers[gapl:gapr]]) + kmers[gapr]
#            print("INIT", gapr, gapl)
#            print(kmers[gapl], kmers[gapr])
#            print("stringo", stringo)
#            print("stringk", stringk)
            early_break = False
            for i in range(gapr-gapl):
                poss_edges = [b + stringo[:self.k-1] for b in "ATCG"]
                max_score = 0
                max_base = None
#                print(poss_edges)
                for pe in poss_edges:
                    if pe in self.edges:
                        tmp_score = self.edges[pe]
                        if tmp_score > max_score:
                            max_base = pe[0]
                if max_base != None:
                    stringo = max_base + stringo
                else:
                    early_break = True
#                    print("early break")
                    stringk = stringk[-len(stringo):]
                    if i == 0:
                        # we were not able to extend at all; return False
                        return False
                    break

#            print(len(stringk), gapr-gapl)
#            print("stringo", stringo)
#            print("stringk", stringk)
#            print("TRIM 1")
            stringo = stringo[:-1]
            stringk = stringk[:-1]
            stringk = stringk[-len(stringo):]
#            print("stringo", stringo)
#            print("stringk", stringk)
            # calculate scores for the whole bridge
            assert len(stringo) <= gapr-gapl-1+self.k
            assert len(stringo) <= gapr-gapl-1+self.k
#            print("SCORE")
            bridge_scores = []
            for i in range(len(stringo)-self.k+1):
                bridge_scores.append(self.edges[stringo[i:i+self.k]])
#            bridge_scores = self.deque_score_bases(bridge_scores)
#            print(bridge_scores)
#            print("TRIM 2")
            stringk = stringk[:-self.k+1]
            stringo = stringo[:-self.k+1]
            if len(stringk) == 0:
                return False
#            print("stringo", stringo)
#            print("stringk", stringk)
#            print("ADJUST SCORES")
            aln = pairwise2.align.globalxx(stringk, stringo, one_alignment_only=True)[0]
            bi = 0
            for i in range(len(aln[0])):
                ci = aln[0][i]
                cj = aln[1][i]
                if cj != "-":
                    if ci != cj:
                        bridge_scores[bi] = (-1, bridge_scores[bi])
                    bi += 1
            # Now pad the extra scores (this is in case the bridge is not as long as the query gap)
#            print(bridge_scores)
#            print("FINAL PAD")
            if early_break == True:
                bridge_scores = [-1 for i in range(gapr-gapl-len(bridge_scores))] + bridge_scores
#            print(bridge_scores)
            return bridge_scores 

        if len(stringo) == self.k:
            # abandon, no bridge exists
            return [0] * len(gapr-gapl)

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
        mask_inds = []
#        print(array)
        for ai in range(len(array)):
            w = array[ai]
            if type(w) == int:
                if w == -1:
                    mask_inds.append(ai)
            else:
                mask_inds.append(ai)
                array[ai] = w[1]

#        print(mask_inds)
#        print(array)
        for ki in range(1, len(array)):
            # append the prev
#            if array[ki-1] == -1:
#                local_maxima.append(-1)
#            else:
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
        for ai in mask_inds:
            local_maxima[ai] = 0
#        print(local_maxima)
        return local_maxima            

    def query(self, kmers, seqsh, min_kmer_prop=0.9):
        sr = SearchResult()
        raw_weight_array = self.get_raw_weight_array(kmers)
        gaps = self.get_weight_array_gaps(raw_weight_array)
        filled_weight_array = raw_weight_array
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
#                print("bounds", gapr, gapl)
#                print("bridge", bridge)
#                print("gapstr", gapstring)
#                print(bridge_scores)
#                print(extra_scores)
                assert len(bridge) == gapr-gapl
                assert len(bridge_scores) == gapr-gapl
                filled_weight_array[gapl:gapr] = extra_scores
        raw_deque_array = self.deque_score_bases(raw_weight_array)
        filled_deque_array = self.deque_score_bases(filled_weight_array)

#        print(raw_weight_array)
#        print(raw_deque_array)
#        print([z for z in zip(raw_weight_array, filled_weight_array)])
        assert sum(filled_weight_array) >= sum(raw_weight_array)
        assert sum(filled_deque_array) >= sum(raw_deque_array)
        
        sr.gap_positions = gaps
        sr.filled_array = filled_weight_array
        sr.raw_array = raw_weight_array
#            sr.filled_deque_array = filled_deque_array
#            sr.raw_deque_array = raw_deque_array
        sr.prop_kmers = len([i for i in raw_deque_array if i > 0])/len(kmers)            
#        print(raw_weight_array)
#        print(sr.prop_kmers, len([i for i in raw_weight_array if i > 0])/len(kmers))

        score = sum([i for i in raw_deque_array if i != -1])
        sr.score = score
#            print(seqsh)
#            sr.print()
        return score

    def classify(self, kmersets, seqsh, min_kmer_prop):
        scores = []
        for si, kmers in enumerate(kmersets):
            query = self.query(kmers, seqsh[si], min_kmer_prop)
            scores.append(query)
        # Sort and return results
        results = []
        inds = np.argsort(scores)
        for ind in inds:
            results.append((ind, scores[ind]))
        return results


