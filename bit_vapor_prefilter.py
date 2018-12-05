import sys
from Bio import SeqIO
import numpy as np
import math
import random
import argparse

class wDBG():
    # DBG with additional weights on each edge
    def __init__(self, strings, k):
        self.edges = {}
        self.k = k
        self.start_positions = set()
        self._build(strings)

    def _build(self, strings):
        sys.stderr.write("Building wDBG\n")
        for si, string in enumerate(strings):
            sys.stderr.write(str(si) + "             \r")
            kmers = [string[i:i+self.k] for i in range(len(string)-self.k+1)]
            for kmer in kmers:
                if kmer in self.edges:
                    self.edges[kmer] += 1
                else:
                    self.edges[kmer] = 1
            # Add the first kmer to start positions, since they can represent the start of biological sequences
#            start = string[:self.k]
#            self.start_positions.add(start)            
        # Now get the start positions, this is inefficient
        for kmer, weight in self.edges.items():
            edge = True
            for b in "ATCG":
                if b+kmer[:-1] in self.edges:
                    edge = False
            if edge == True:
                self.start_positions.add(kmer)
        sys.stderr.write("\n")

    def get_n_branches(self):
        nodes = {}
        for e in self.edges:
            n1 = e[1:]
            n2 = e[:-1]
            nodes[n1] = 0
            nodes[n2] = 0
        for node in nodes:
            for b in "ATCG":
                if node+b in self.edges:
                    nodes[node] += 1
        n = 0
        for node in nodes:
            if nodes[node] > 1:
                 n += 1
        return n

    def get_paths(self):
        paths = [[p, self.edges[p]] for p in self.start_positions]
        pstrings = [p[0] for p in paths]
        assert len(pstrings) == len(set(pstrings))
        final_paths = []
        sys.stderr.write("Building paths, from %d start points\n" % len(paths))
        while paths != []:
            switch = True
            tmp_paths = []
            for path in paths:
                assert paths.count(path) == 1
                localswitch = False
                best_tmps = []
                for b in "ATCG":
                    tmpkmer = path[0][-self.k+1:] + b
                    if tmpkmer in self.edges:
                        tmp_path = [path[0], path[1]]
                        tmp_path[0] = tmp_path[0] + b
                        tmp_path[1] += self.edges[tmpkmer]
#                        tmp_paths.append(tmp_path)
                        best_tmps.append(tmp_path)
                        localswitch = True
                if localswitch == True:
                    # Never a tie at a branch! Assert it
                    maxtmp = max(best_tmps, key = lambda x:x[1])[1]
                    maxtmps = [t for t in best_tmps if t[1] == maxtmp]
    #                assert len(maxtmps) == 1
                    best_tmp = maxtmps[0]
    #                if best_tmp not in tmp_paths:
                    tmp_paths.append(best_tmp)
                else:
                    final_paths.append(path)
            # For our heuristic walk, we cant haave more paths than starts
            assert len(tmp_paths) <= len(paths)
            assert len(paths) <= len(self.start_positions)
            paths = tmp_paths
            sys.stderr.write("%d %d             \r" % (len(final_paths), len(tmp_paths)))
        sys.stderr.write("\n")
        return final_paths           

class cDBG():
    def __init__(self, k):
        self.edges = {}
        self.k = k
        self.n = None
        self.color_classes = set()
    
    @classmethod
    def from_strings(cls, strings, k):
        c = cls(k)
        c.n = len(strings)
        c._build(strings)
        return c

    @classmethod
    def from_strings_and_subgraph(cls, strings, k, subgraph):
        c = cls(k)
        c.n = len(strings)
        c._build(strings, subgraph)
        c.build_color_classes()
        return c

    def _build(self, strings, subgraph=False):
        # Build the color classes
        sys.stderr.write("Adding strings:\n")
        for si, string in enumerate(strings):
            sys.stderr.write(str(si)+"        \r")
            kmers = (string[i:i+self.k] for i in range(len(string)-self.k+1))
            for kmer in kmers:
                if subgraph==False or kmer in subgraph.edges:
                    if kmer in self.edges:
                        self.edges[kmer].append(si)
                    else:
                        self.edges[kmer] = [si]
        sys.stderr.write("\n")

    def score(self, string):
        score = 0
        for kmer in (string[i:i+self.k] for i in range(len(string)-self.k+1)):
            if kmer in self.edges:
                score += 1
        return score            

    def build_color_classes(self, subkmers=None):
        sys.stderr.write("Building color classes:\n")
        z = 0
        if subkmers == None:
            kmers = self.edges
        else:
            kmers = subkmers
        for kmer, ccl in kmers.items():
            sys.stderr.write(str(z) + "          \r")
            # Create empty bit string with leading 1
            cc = (1 << (self.n))
            for c in ccl:
                cc = cc | (1 << c)
            self.color_classes.add(cc)
            self.edges[kmer] = cc
            z += 1
        sys.stderr.write("\n")

    def classify(self, wdbg, seqs):

        ## Note - Previously seqs arg was not specified but was being used below to define sumo variable from main
        # not sure if you realised this or not

        # First build a dbg, second find shortest path in it, thirdly parse the color information
#        wdbg = wDBG(reads, self.k)
        paths = wdbg.get_paths()
        lens = [len(p) for p in paths]
        scores = [p[1] for p in paths]
        pathseqs = [p[0] for p in paths]
        maxpath = max(paths, key = lambda x:x[1])
        colors = []
        for seq in [maxpath[0]]:
            kmers = (seq[i:i+self.k] for i in range(len(seq)))
            for kmer in kmers:
                if kmer in self.edges:
                    colors.append(self.edges[kmer])

        sumo = np.zeros(len(seqs))
        sys.stderr.write("Summing colors")
        for c in colors:
            arr = np.fromstring(np.binary_repr(c), dtype='S1').astype(int)
            sumo += arr[1:]

        maxi = max(sumo)
        maxs = [len(sumo)-i-1 for i in range(len(sumo)) if sumo[i] == maxi]
        return maxs

def get_kmers(strings,k):
    kmers = set()
    for string in strings:
        for i in range(len(string)-k+1):
            kmers.add(string[i:i+k])
    return kmers

def rev_comp(read):
    read = read.replace("T", "a")
    read = read.replace("A", "t")
    read = read.replace("C", "g")
    read = read.replace("G", "c")
    return read.upper()[::-1]

def kmer_prefilter(reads, filterkmers, threshold, k=21):
    for raw_read in reads:
        for r in [raw_read, rev_comp(raw_read)]:
            kmers = set([r[i:i+k] for i in range(len(r))])
            prop = float(len(kmers & filterkmers)) / len(r)
            if prop > threshold:
                yield r

def choose_paths(paths):
    maxi = max(paths, key=lambda x: x[1])[1]
    for path in paths:
        if path[1] == maxi:
            yield path

def main():

    random.seed(100)

    sys.stderr.write("Loading sequences\n")

    seqsr = [r for r in SeqIO.parse(sys.argv[3], "fasta")]
    seqs = [str(r.seq) for r in seqsr]
    seqsh = [r.description for r in seqsr]

    K = int(sys.argv[1])
    score_threshold = float(sys.argv[2])

    reads = []
    for f in sys.argv[4:]:
        reads += [str(r.seq) for r in SeqIO.parse(f, "fastq")]

    sys.stderr.write("Prefiltering reads and taking revComp where necessary\n")

    dbkmers = get_kmers(seqs, K)
    reads = [r for r in kmer_prefilter(reads, dbkmers, score_threshold, K)]
    sys.stderr.write(str(len(reads)) + " survived\n")

    if len(reads) == 0:
        sys.stderr.write("Exiting. Is there any virus in your sequences? Try a lower filtering threshold.\n")
        sys.exit()

    # prefilter
    wdbg = wDBG(reads, K)
    cdbg = cDBG.from_strings_and_subgraph(seqs, K, wdbg)

    ### not sure if you want seqs below?
    cls = cdbg.classify(wdbg, seqs)

    for c in cls:
        print(len(reads), c, seqsh[c])


parser = argparse.ArgumentParser(description="Do some sweet viral classification")
group = parser.add_mutually_exclusive_group()
group.add_argument("-v", "--verbose", action="store_true")
group.add_argument("-q", "--quiet", action="store_true")
parser.add_argument("Kmer", type=int, help="Kmer Length")
parser.add_argument("Score_thres", type=float, help="Score thresehold")
parser.add_argument("Fasta", type=str, help="Fasta file")
parser.add_argument("Fastq", type=str, help="Fastq file")

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

if args.quiet:
    # Do something
    print (str("WTF"))
    main()
elif args.verbose:
    # Do something else
    answer = "YOLO FOOL"
    print ("{} is the name of the game and the answer is {} ({})".format(args.Fasta, args.Fastq, answer))
    main()
else:
    # Still do something
    answer = "YOLO FOOL"
    print ("{} is way better than {} but make sure you {}".format(args.Fasta, args.Fastq, answer))
    main()








#    for i in range(200):
#        roll = random.randint(0,len(seqs)-1)
#        ref = seqs[roll]
#        quasispecies = generate_quasispecies(ref, 10, 0.01)
#        reads = readize(quasispecies)
#        wdbg = wDBG(reads, K)
#        print("nreads", len(reads), nmut)
#        print("nbranches", wdbg.get_n_branches())
#        paths = wdbg.get_paths()
#        best_paths = choose_paths(paths)
#        path_kmers = set()
#        print()
#        for path in best_paths:
#            print(path)
#            pkmers = [path[0][i:i+K-1] for i in range(len(path[0])-K+2)]
#            for kmer in pkmers:
#                path_kmers.add(kmer)
#        print(len(reads), len(wdbg.edges))
#        dbgplotter.plot_dbg(wdbg, path_kmers)
#        cls = cdbg.classify(reads)
#        print(roll, cls, roll in cls)
