""" DBG classes and related objects """
import sys

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

    def query(self, seq):
        score = 0
        for ki in range(len(seq)-self.k+1):
            kmer = seq[ki:ki+self.k]
            if kmer in self.edges:
                score += self.edges[kmer]
        return score

    def classify(self, seqs, min_missing_kmers=100.):
        # ALLOW FOR PARTIAL MATCHES, WHEN WALKING
        # IF A KMER IS NOT FOUND, WE CAN CHECK HAMMING DISTANCES OF KMERS?
        # ATTEMPT TO `BRIDGE THE GAP' (FOR LATER)
        scores = [self.query(seq) for seq in seqs]
        maxs = max(scores)
        maxcls = [si for si in range(len(seqs)) if scores[si] == maxs]
        return maxs, maxcls


