""" DBG classes and related objects """
import sys
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
from collections import deque
from pyvapor.vaporfunc import *

class SearchResult():
    """
    Class to hold search results and methods
    Used for debugging
    """
    def __init__(self):
        # Most variables are unused unless debugging is specified
        self.i = -1
        self.kmers = []
        self.prop_kmers = -1
        self.gap_positions = []
        self.bridges = {}
        self.suboptimal_branches = None
        self.raw_array = []
        self.filled_array = []
        self.filled_deque_array = []
        self.score = -1

    def compare(self, sr, wdbg):
        """
        For debugging; function to compare a search result with
        a specific query; outputs per-base query seeding,
        trimming, bridging, and scoring information
        """
        gapset1 = set()
        for gapl, gapr in self.gap_positions:
            for gi in range(gapl, gapr):
                gapset1.add(gi)
        gapset2 = set()
        for gapl, gapr in sr.gap_positions:
            for gi in range(gapl, gapr):
                gapset2.add(gi)

        for i in range(min(len(self.raw_array), len(sr.raw_array))):
            tag1 = ""
            tag2 = ""
            if self.kmers[i] != sr.kmers[i]:
                tag1 += "*"
                tag2 += "*"
            if i in gapset1:
                tag1 += "g"
            if i in gapset2:
                tag2 += "g"
            if i in self.suboptimal_branches:
                tag1 += "b"
            if i in sr.suboptimal_branches:
                tag2 += "b"
            if i in sr.bridges:
                tag2 += sr.bridges[i]
            if i in self.bridges:
                tag1 += self.bridges[i]
            print(i, self.kmers[i], self.kmers[i] in wdbg.nodes, self.raw_array[i], self.filled_array[i], self.filled_deque_array[i], tag1, "\t", sr.kmers[i], sr.kmers[i] in wdbg.nodes, sr.raw_array[i], sr.filled_array[i], sr.filled_deque_array[i], tag2)
        self_cycles = [kmer for kmer in self.kmers if self.kmers.count(kmer) > 1]
        print("self cycles:", len(self_cycles))
        sr_cycles = [kmer for kmer in self.kmers if sr.kmers.count(kmer) > 1]
        print("sr cycles:", len(sr_cycles))

class wDBG():
    """ Basic DBG with associated edge weights """
    def __init__(self, strings, k):
        """ Initialized with strings, k, and reference kmers """
        # Only explicitly store nodes
        self.nodes = {}
        self.k = k
        self._build(strings)
        self.caching = True
        self.path_cache = {}
        self.max_trim_size = self.k+1

    def _build(self, strings):
        for si, string in enumerate(strings):
            kmers = [string[i:i+self.k] for i in range(len(string)-self.k+1)]
            for kmer in kmers:
                if kmer in self.nodes:
                    self.nodes[kmer] += 1
                else:
                    self.nodes[kmer] = 1

    def cull_low(self, min_cov=5):
        """
        Culls kmers with a coverage less than min_cov
        """
        keyvals = [(k, v) for k, v in self.nodes.items()]
        for key, val in keyvals:
            if val <= min_cov:
                del self.nodes[key]            

    def mask_against_bridge(self, query, bridge, gapl):
        """
        Takes two strings; returns indices where they mismatch
        offset by gapl
        """
        mask = []
        for i in range(len(query)):
            if query[i] != bridge[i]:
                mask.append(gapl+i)
        return mask

    def extend_bridge(self, kmer, n, direction, debug=True):
        """
        Walks along the wDBG n positions
        making heuristic locally optimal decisions
        at branches. Returns string, score array
        """
        # First check the cache
        if self.caching == True:
            if (kmer, n, direction) in self.path_cache:
                return self.path_cache[(kmer, n, direction)]
        string = kmer
        scorearr = np.zeros(n)
        if direction == 1:
            si = 0
        else:
            si = -1
        while len(string) < n+self.k:
            if direction == 1:
                poss_nodes = [string[-self.k+1:] + b for b in "ATCG"]
            elif direction == -1:
                poss_nodes = [b + string[:self.k-1] for b in "ATCG"]
            max_score = 0
            max_base = None
            for pe in poss_nodes:
                if pe in self.nodes:
                    tmpscore = self.nodes[pe]
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
            self.path_cache[(kmer, n, direction)] = (string, scorearr)
        return string, scorearr

    def get_raw_weight_array(self, kmers):
        """
        Takes a sequence of kmers, and returns
        an array of weights, zero if a kmer is
        not present in the graph
        """
        array = np.zeros(len(kmers))
        for ki, kmer in enumerate(kmers):
            if kmer in self.nodes:
                array[ki] = self.nodes[kmer]
        return array
    
    def get_weight_array_gaps(self, array):
        """
        Obtains positions of gaps (sequences of zeroes)
        in an array
        """
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
        """ 
        For a weight array corresponding to a 
        contiguous stretch of kmers,
        uses a deque to find max weight that contains a base,
        returning another array
        """
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
        """
        Walks backward along an array of kmers until a
        step is detected that is suboptimal
        returns positions of these branches
        """
        suboptimal_branches = set()
        bases = list("ATCG")
        for ki, kmer in enumerate(kmers):
            if kmer in self.nodes:
                score = self.nodes[kmer]
                alts = [kmer[:-1] + b for b in bases]
                altscores = [self.nodes[amer] for amer in alts if amer in self.nodes]
                if altscores != []:
                    if score < max(altscores):
                        suboptimal_branches.add(ki)
        return suboptimal_branches 

    def expand_gaps(self, gaps, suboptimal_branches, max_gapr):        
        """ 
        Acts in place to expand gaps to suboptimal branch positions
        takes the most distant sub branch within self.max_trim_size
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

    def seed(self, kmers):
        """
        Takes a query sequence of kmers,
        seeds it in the wDBG
        return SearchResult object
        """
        sr = SearchResult()
        # First obtain the raw weight array for kmers of a sequence
        raw_weight_array = self.get_raw_weight_array(kmers)
        kmer_cov = np.count_nonzero(raw_weight_array)/len(raw_weight_array)
        sr.raw_weight_array = raw_weight_array
        sr.kmer_cov = kmer_cov
        return sr

    def complete_query(self, sr, kmers=None, string=None, debug=False):
        if kmers == None:
            kmers = get_kmers([string], self.k)[0]
        kmer_cov = sr.kmer_cov
        raw_weight_array = sr.raw_weight_array
        # Get the gaps
        gaps = self.get_weight_array_gaps(raw_weight_array)
        # Get the suboptimal branches
        sub = self.get_suboptimal_branches(kmers)
        self.expand_gaps(gaps, sub, len(kmers))
        if debug == True:
            # Copy the raw weight as a record
            sr.suboptimal_branches = sub
            sr.gap_positions = gaps
            filled_weight_array = [r for r in raw_weight_array]
        else:
            # Don't copy, not debugging
            filled_weight_array = raw_weight_array
        all_masks = []
        for gapl, gapr in gaps:
            if gapl != 0 and gapr != len(kmers):
                gapstring = kmers2str(kmers[gapl:gapr])[self.k-1:]
                bridge, bridge_scores = self.extend_bridge(kmers[gapl-1], gapr-gapl, 1, debug)
                bridge_rev, bridge_scores_rev = self.extend_bridge(kmers[gapr], gapr-gapl, -1, debug)
                gapstring_rev = kmers2str(kmers[gapl:gapr])[:-self.k+1]
                if sum(bridge_scores_rev) > sum(bridge_scores):
                    mask = self.mask_against_bridge(gapstring_rev, bridge_rev, gapl)
                    filled_weight_array[gapl:gapr] = bridge_scores_rev
                else:
                    mask = self.mask_against_bridge(gapstring, bridge, gapl)
                    filled_weight_array[gapl:gapr] = bridge_scores

            elif gapr != len(kmers) and gapl == 0:
                gapstring = kmers2str(kmers[gapl:gapr])[self.k-1:]
                bridge, bridge_scores = self.extend_bridge(kmers[gapr], gapr-gapl, -1, debug)
                mask = self.mask_against_bridge(gapstring, bridge, gapl)
                filled_weight_array[gapl:gapr] = bridge_scores

            elif gapl > 0 and gapr == len(kmers):
                gapstring = kmers2str(kmers[gapl:gapr])[self.k-1:]
                bridge, bridge_scores = self.extend_bridge(kmers[gapl-1], gapr-gapl, 1, debug)
                mask = self.mask_against_bridge(gapstring, bridge, gapl)
                filled_weight_array[gapl:gapr] = bridge_scores
            all_masks += mask
            if debug == True:
                for i in range(gapl, gapr):
                    sr.bridges[i] = bridge[i-gapl]
        if debug == True:
            sr.filled_array = filled_weight_array
        # Deque score
        filled_weight_array = np.concatenate((filled_weight_array, np.zeros(self.k-1)))
        filled_deque_array = self.deque_score_bases(filled_weight_array)
        for maski in all_masks:
            filled_deque_array[maski] = 0
        if debug == True:
            sr.filled_deque_array = filled_deque_array
        # Sum, also get estimated pid
        nonzeros = [i for i in filled_deque_array if i > 0]
        est_pid = len(nonzeros)/len(filled_deque_array)
        sr.est_pid = est_pid
        score = sum(nonzeros)
        sr.score = score
        return sr

    def query(self, kmers, min_kmer_prop, debug=False):
        # DEPRECIATED
        """ 
        Takes a query set of kmers,
        param min_kmer_prop and debug flag 
        for recording information in search result.
        Returns a SearchResult object
        """
        sr = SearchResult()
        # First obtain the raw weight array for kmers of a sequence
        raw_weight_array = self.get_raw_weight_array(kmers)
        kmer_cov = np.count_nonzero(raw_weight_array)/len(raw_weight_array)
        sr.kmer_cov = kmer_cov
        # Next trim the raw weight array
        if debug==True:
            sr.raw_array = raw_weight_array
            sr.kmers = kmers
        if kmer_cov > min_kmer_prop:
            # Get the gaps
            gaps = self.get_weight_array_gaps(raw_weight_array)
            # Get the suboptimal branches
            sub = self.get_suboptimal_branches(kmers)
            self.expand_gaps(gaps, sub, len(kmers))
            if debug == True:
                # Copy the raw weight as a record
                sr.suboptimal_branches = sub
                sr.gap_positions = gaps
                filled_weight_array = [r for r in raw_weight_array]
            else:
                # Don't copy, not debugging
                filled_weight_array = raw_weight_array
            all_masks = []
            for gapl, gapr in gaps:
                if gapl != 0 and gapr != len(kmers):
                    gapstring = kmers2str(kmers[gapl:gapr])[self.k-1:]
                    bridge, bridge_scores = self.extend_bridge(kmers[gapl-1], gapr-gapl, 1, debug)
                    bridge_rev, bridge_scores_rev = self.extend_bridge(kmers[gapr], gapr-gapl, -1, debug)
                    gapstring_rev = kmers2str(kmers[gapl:gapr])[:-self.k+1]
                    if sum(bridge_scores_rev) > sum(bridge_scores):
                        mask = self.mask_against_bridge(gapstring_rev, bridge_rev, gapl)
                        filled_weight_array[gapl:gapr] = bridge_scores_rev
                    else:
                        mask = self.mask_against_bridge(gapstring, bridge, gapl)
                        filled_weight_array[gapl:gapr] = bridge_scores

                elif gapr != len(kmers) and gapl == 0:
                    gapstring = kmers2str(kmers[gapl:gapr])[self.k-1:]
                    bridge, bridge_scores = self.extend_bridge(kmers[gapr], gapr-gapl, -1, debug)
                    mask = self.mask_against_bridge(gapstring, bridge, gapl)
                    filled_weight_array[gapl:gapr] = bridge_scores

                elif gapl > 0 and gapr == len(kmers):
                    gapstring = kmers2str(kmers[gapl:gapr])[self.k-1:]
                    bridge, bridge_scores = self.extend_bridge(kmers[gapl-1], gapr-gapl, 1, debug)
                    mask = self.mask_against_bridge(gapstring, bridge, gapl)
                    filled_weight_array[gapl:gapr] = bridge_scores
                all_masks += mask
                if debug == True:
                    for i in range(gapl, gapr):
                        sr.bridges[i] = bridge[i-gapl]
            if debug == True:
                sr.filled_array = filled_weight_array
            # Deque score
            filled_weight_array = np.concatenate((filled_weight_array, np.zeros(self.k-1)))
            filled_deque_array = self.deque_score_bases(filled_weight_array)
            for maski in all_masks:
                filled_deque_array[maski] = 0
        else:
            if debug == True:
                sr.filled_deque_array = raw_weight_array
            sr.est_pid = -1
            sr.score = -1
            return sr
        if debug == True:
            sr.filled_deque_array = filled_deque_array
        # Sum, also get estimated pid
        nonzeros = [i for i in filled_deque_array if i > 0]
        est_pid = len(nonzeros)/len(filled_deque_array)
        sr.est_pid = est_pid
        score = sum(nonzeros)
        sr.score = score
        return sr

    def classify(self, seqs, seqsh, min_kmer_prop, top_seed_frac, debug_query=None, low_mem=False):
        """
        Queries a set of sequences seqs, with headers seqsh,
        parameter min_kmer_prop
        and if debugging, a debug query
        """
        seeds = []
        for si, seq in enumerate(seqs):
            kmers = [seq[i:i+self.k] for i in range(len(seq)-self.k+1)]
            seed = self.seed(kmers)
            seed.index = si
            if seed.kmer_cov > min_kmer_prop:
                if low_mem == False:
                    seed.kmers = kmers
                seeds.append(seed)

        topseeds = sorted(seeds, key=lambda x:x.kmer_cov, reverse=True)[:int(np.ceil(top_seed_frac*len(seqs)))]

        scores = []
        for seed in topseeds:
            if low_mem == False:
                sr = self.complete_query(seed, seed.kmers, debug_query)
            else:
                sr = self.complete_query(seed, None, seqs[seed.index], debug_query)
            scores.append(sr)

        # Sort the results
        results = sorted(scores, key = lambda x: x.score * x.est_pid, reverse=True)

        if debug_query != None:
            # For debugging
            for hi, h in enumerate(seqsh):
                if h == debug_query:
                    seq = seqs[hi]
                    kmers = [seq[i:i+self.k] for i in range(len(seq)-self.k+1)]
                    sr = self.query(kmers, min_kmer_prop, True)
                    sr.compare(results[0], self)
                    print(sr.score)

        return results
