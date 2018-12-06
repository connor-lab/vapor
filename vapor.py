import sys
import numpy as np
import random
import argparse
import os

class Path():
    def __init__(self, start_edge, start_score):
        self.edges = {start_edge}
        self.path = [start_edge]
        self.score = start_score
        self.cyclic = False
    def add_edge(self, edge, score):
        if edge in self.edges:
            self.cyclic = True
        self.edges.add(edge)
        self.path.append(edge)
        self.score += score
    def get_string(self):
        s = self.path[0]
        for edge in self.path[1:]:
            s += edge[-1]
        return s

def remove_overlaps(paths):
    # Greedy approach to resolving overlaps
    flagarr = [0 for p in paths]
    for i in range(len(paths)-1):
        for j in range(i+1, len(paths)):
            pi = paths[i]
            pj = paths[j]
            if len(pi.edges | pj.edges) > 0:
                if pi.score > pj.score:
                # The paths overlap at at least one kmer
                    flagarr[j] = 1
                else:
                    flagarr[i] = 1
    for i in range(len(flagarr)):
        if flagarr[i] == 0:
            yield paths[i]
 
class wDBG():
    # DBG with additional weights on each edge
    def __init__(self, strings, k):
        self.edges = {}
        self.k = k
        self.start_positions = set()
        self._build(strings)

    def _build(self, strings):

        sys.stderr.write("Building wDBG\n")
        newkmerc = 0
        for si, string in enumerate(strings):
            sys.stderr.write(str(si) + "      \r")
            kmers = [string[i:i+self.k] for i in range(len(string)-self.k+1)]
            newkmerc = 0
            for kmer in kmers:
                if kmer in self.edges:
                    self.edges[kmer] += 1
                else:
                    newkmerc += 1
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

    def cull(self, kmers):
        todel = []
        for kmer in self.edges:
            if kmer not in kmers:
                todel.append(kmer)
        for kmer in todel:
            if kmer not in self.start_positions:
                del self.edges[kmer]

    def get_paths(self):
        paths = [Path(p, self.edges[p]) for p in self.start_positions]
        pstrings = [p.get_string() for p in paths]
        assert len(pstrings) == len(set(pstrings))
        final_paths = []
        sys.stderr.write("Building paths, from %d start points\n" % len(paths))
        tmp_paths = []
        while paths != []:
            switch = True
            assert len(paths) + len(final_paths) == len(self.start_positions)
            tmp_paths = []
            for path in paths:
                if path.cyclic == True:
                    final_paths.append(path)
                    print("cycling")
                else:
                    assert paths.count(path) == 1
                    localswitch = False
                    tmp_additions = {}
                    for b in "ATCG":
                        tmpkmer = path.path[-1][1:] + b
                        if tmpkmer in self.edges:
                            tmp_additions[b] = self.edges[tmpkmer]
#                            tmp_path.add_edge(tmpkmer, self.edges[tmpkmer])
    #                        tmp_path[0] = tmp_path[0] + b
    #                        tmp_path[1] += self.edges[tmpkmer]
    #                        tmp_paths.append(tmp_path)
#                            best_tmps.append(tmp_path)
                            localswitch = True
                    if localswitch == True:
                        # Never a tie at a branch! Assert it
                        maxtmp = max(tmp_additions.items(), key = lambda x:x[1])[1]
                        maxtmps = [t[0] for t in tmp_additions.items() if t[1] == maxtmp]
                        best_tmp = maxtmps[0]
                        path.add_edge(path.path[-1][1:]+best_tmp, maxtmp)
                        tmp_paths.append(path)
                    else:
                        final_paths.append(path)                
            paths = tmp_paths

            # For our heuristic walk, we cant haave more paths than starts
            assert len(tmp_paths) <= len(paths)
            assert len(paths) <= len(self.start_positions)
            paths = tmp_paths
            sys.stderr.write("%d             \r" % (len(final_paths)))
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
        paths = remove_overlaps(paths)
        colors = []
        for path in paths:
            seq = path.get_string()
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
        return (maxs, maxi)

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

def parse_and_prefilter(fqs, dbkmers, threshold, k):
    c = 0
    M = float(len(dbkmers))
    for fq in fqs:
        with open(fq) as f:
            for line in f:
                if c == 1:
                    tmpseq = line.strip()
                    kcount = 0
                    for i in range(0, len(tmpseq), k):
                        if tmpseq[i:i+k] in dbkmers:
                            kcount += 1
                    if k*kcount/M < threshold:
                        yield tmpseq
                c += 1                  
                if c == 4:
                    c = 0

def parse_fasta_uniq(fasta, filter_Ns=True):
    tmph = ""
    tmps = ""
    hs = []
    ss = []
    sseen = set()
    with open(fasta) as f:
        for line in f:
            l = line.strip()
            if l[0] == ">":
                if tmps not in sseen:
                    if ((filter_Ns == True) and "N" not in tmps) or filter_Ns == False:
                        hs.append(tmph)
                        ss.append(tmps) 
                        sseen.add(tmps)
                tmph = l
                tmps = ""
            else:
                tmps += l
    hs.append(tmph)
    ss.append(tmps)
    return hs, ss 

def choose_paths(paths):
    maxi = max(paths, key=lambda x: x[1])[1]
    for path in paths:
        if path[1] == maxi:
            yield path

def subsample(reads, n):
    if n >= len(reads):
        return reads
    else:
        return random.sample(reads, n)
    

def blockErr():
    sys.stderr = open(os.devnull, 'w')

# Restore
def enablePrint():
    sys.stdout = sys.__stdout__

def main(quiet, K, score_threshold, fasta, fastqs):

    random.seed(100)

    ## maybe something like this if you can be bothered.
    if quiet:
        blockErr()

    sys.stderr.write("Loading database sequences\n")
    seqsh, seqs = parse_fasta_uniq(fasta)
    sys.stderr.write("Got %d unique sequences\n" % len(seqs))

    sys.stderr.write("Getting database kmers\n")
    dbkmers = get_kmers(seqs, K)

    sys.stderr.write("Filtering reads\n")
    reads = [r for r in parse_and_prefilter(fastqs, dbkmers, score_threshold, K)]
    reads = subsample(reads, 5000)
    sys.stderr.write(str(len(reads)) + " reads survived\n")

    if len(reads) == 0:
        enablePrint()
        sys.stderr.write("Exiting. Is there any virus in your sequences? Try a lower filtering threshold.\n")
        sys.exit(1)

    # prefilter
    wdbg = wDBG(reads, K)
    # cull any kmers that are not present in the reference kmers; these do not give additional information
    sys.stderr.write("Culling kmers\n")
    wdbg.cull(dbkmers)

    cdbg = cDBG.from_strings_and_subgraph(seqs, K, wdbg)

    ### not sure if you want seqs below?
    cls,score = cdbg.classify(wdbg, seqs)
    for c in cls:
        print(seqsh[c]+","+str(score))

    sys.stderr.write("\nClassification Complete\n")


##### Command line Interface

parser = argparse.ArgumentParser(description="Do some sweet viral classification")
group = parser.add_mutually_exclusive_group()
group.add_argument("-q", "--quiet", action="store_true")

parser.add_argument("-k", type=int, help="Kmer Length")
parser.add_argument("-s", type=float, help="Score threshold")
parser.add_argument("-fa", type=str, help="Fasta file")
parser.add_argument("-fq", nargs='+', type=str, help="Fastq file/files")

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

## set thresholds for user input
max_kmer = 30
min_kmer = 2
max_thres = 1
min_thres = 0

if args.k < max_kmer and args.k > min_kmer:
    if args.s < max_thres and args.s > min_thres:
        ###########Run main
        main(args.quiet, args.k, args.s, args.fa, args.fq)

    else:
        sys.stderr.write("\nPlease input correct score threshold ({} to {}) \n \n".format(min_thres, max_thres))
        parser.print_help(sys.stderr)
        sys.exit(1)
else:
    sys.stderr.write("\nPlease input correct kmer length ({} to {}) \n \n".format(min_kmer, max_kmer))
    parser.print_help(sys.stderr)
    sys.exit(1)
