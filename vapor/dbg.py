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
            print(len(to_join), gapr-gapl, gapr-gapl+self.k, gapr, len(kmers))
            assert len(to_join) == gapr-gapl
            assert gapr <= len(kmers)
            assert len(stringk) == gapr-gapl+self.k
            print("INIT", gapr, gapl)
            if gapr == len(kmers):
                print(kmers[gapl])
            else:
                print(kmers[gapl], kmers[gapr])
            print("stringo", stringo)
            print("stringk", stringk)
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
                    print("early break")
                    break            
            print(gapl, gapr)
            print("BRIDGE")
            print("stringo", stringo)
            print("stringk", stringk)
            print("TRIM1")
            # TRIM 1
            stringo = stringo[1:]
            assert len(stringk) == gapr-gapl+self.k
            stringk = stringk[1:]
            print(len(stringk), gapr-gapl-1+self.k, gapr, gapl, len(kmers))
            assert len(stringk) == gapr-gapl-1+self.k
            assert len(stringo) <= gapr-gapl-1+self.k
            stringk = stringk[:len(stringo)]
            assert len(stringo) == len(stringk)
            print("stringo", stringo)
            print("stringk", stringk)
            print("SCORE")
            # SCORE
            bridge_scores = [0 for i in range(gapr-gapl)]
            for i in range(len(stringo)-self.k+1):
                bridge_scores[i] = self.edges[stringo[i:i+self.k]]
#            bridge_scores = self.deque_score_bases(bridge_scores)
            print(bridge_scores)
            print("TRIM2")
            # TRIM 2
            stringo = stringo[self.k-1:gapr-gapl]
            stringk = stringk[self.k-1:gapr-gapl]
            if len(stringk) == 0:
                return False
            print("stringo", stringo)
            print("stringk", stringk)
            assert len(stringk) < gapr-gapl+self.k-1
            assert len(stringo) < gapr-gapl+self.k-1
            # ADJUST SCORES
            aln = pairwise2.align.globalxx(stringk, stringo, one_alignment_only=True)[0]
            print("ALIGNMENT")
            print(aln[0])
            print(aln[1])
            bi = self.k-1
            print(len(bridge_scores), len(stringo), len(stringk), gapl, gapr)
        
            for i in range(len(aln[0])):
                ci = aln[0][i]
                cj = aln[1][i]
                if cj != "-":
                    if ci != cj:
                        bridge_scores[bi] = -1
                    bi += 1
            # Now pad the extra scores (this is in case the bridge is not as long as the query gap)
            print(len(bridge_scores), gapr-gapl, gapr, gapl, len(kmers))
            assert len(bridge_scores) == gapr-gapl
            return bridge_scores

        else:
            # INITIALIZE
            stringo = kmers[gapr]
            stringk = "".join([kmer[0] for kmer in kmers[gapl:gapr]]) + kmers[gapr]
            print("INIT", gapr, gapl)
            print(kmers[gapl], kmers[gapr])
            print("stringo", stringo)
            print("stringk", stringk)
            early_break = False
            for i in range(gapr-gapl):
                poss_edges = [b + stringo[:self.k-1] for b in "ATCG"]
                max_score = 0
                max_base = None
                print(poss_edges)
                for pe in poss_edges:
                    if pe in self.edges:
                        tmp_score = self.edges[pe]
                        if tmp_score > max_score:
                            max_base = pe[0]
                if max_base != None:
                    stringo = max_base + stringo
                else:
                    early_break = True
                    print("early break")
                    stringk = stringk[-len(stringo):]
                    if i == 0:
                        # we were not able to extend at all; return False
                        return False
                    break

            print(len(stringk), gapr-gapl)
            print("stringo", stringo)
            print("stringk", stringk)
            print("TRIM 1")
            stringo = stringo[:-1]
            stringk = stringk[:-1]
            stringk = stringk[-len(stringo):]
            print("stringo", stringo)
            print("stringk", stringk)
            # calculate scores for the whole bridge
            assert len(stringo) <= gapr-gapl-1+self.k
            assert len(stringo) <= gapr-gapl-1+self.k
            print("SCORE")
            bridge_scores = []
            for i in range(len(stringo)-self.k+1):
                bridge_scores.append(self.edges[stringo[i:i+self.k]])
#            bridge_scores = self.deque_score_bases(bridge_scores)
            print(bridge_scores)
            print("TRIM 2")
            stringk = stringk[:-self.k+1]
            stringo = stringo[:-self.k+1]
            if len(stringk) == 0:
                return False
            print("stringo", stringo)
            print("stringk", stringk)
            print("ADJUST SCORES")
            aln = pairwise2.align.globalxx(stringk, stringo, one_alignment_only=True)[0]
            bi = 0
            for i in range(len(aln[0])):
                ci = aln[0][i]
                cj = aln[1][i]
                if cj != "-":
                    if ci != cj:
                        bridge_scores[bi] = -1
                    bi += 1
            # Now pad the extra scores (this is in case the bridge is not as long as the query gap)
            print(bridge_scores)
            print("FINAL PAD")
            if early_break == True:
                bridge_scores = [-1 for i in range(gapr-gapl-len(bridge_scores))] + bridge_scores
            print(bridge_scores)
            return bridge_scores 

        if len(stringo) == self.k:
            # abandon, no bridge exists
            return [0] * len(gapr-gapl)
    
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
                print("search",search)
                if search != False:
                    print(len(search), gapr-gapl, gapl, gapr)
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
        return local_maxima            

    def query(self, kmers, seqsh, min_kmer_prop=0.9):
        weight_array = self.get_weight_array(kmers, min_kmer_prop)
        if weight_array != False:
            deq_scores = self.deque_score_bases(weight_array)            
            if "Chile" in seqsh or "A/Cambodia/NHRCC00010/2009" in seqsh:
                print(seqsh)
                print([z for z in zip(range(len(weight_array)), weight_array)])
                print([z for z in zip(range(len(deq_scores)), deq_scores)]) 
                print(len(deq_scores))
                print(deq_scores.count(-1))
            score = sum([i for i in deq_scores if i != -1])
            return score
        else:
            return -1

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


