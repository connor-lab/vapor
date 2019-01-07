""" DBG classes and related objects """

import numpy as np
import sys

class Path():
    """ Path object for traversing DBG """
    def __init__(self, start_edge, start_score):
        # Initialize from a start position in a DBG
        # Recount, for each edge, the next edge
        # For resolving cycles
        self.edges = {start_edge:[]}
        self.path = [start_edge]
        self.score = start_score
        self.cyclic = False
    def add_edge(self, edge, score):
        if edge in self.edges:
            sys.stderr.write("Cycle detected\n")
            self.cyclic = True
        if edge not in self.edges:
            self.edges[edge] = []
        # Also add the last base to the previous edge
        # To avoid cycles
        last_kmer = self.path[-1]
        self.edges[last_kmer].append(edge[-1])
        self.path.append(edge)
        self.score += score
    def get_string(self):
        s = self.path[0]
        for edge in self.path[1:]:
            s += edge[-1]
        return s

def remove_overlaps(paths):
    """ Greedily resolving conflicting paths """
    # Removes the lowest scoring of a pair with nonzero kmer set intersection """
    flagarr = [0 for p in paths]
    for i in range(len(paths)-1):
        for j in range(i+1, len(paths)):
            if flagarr[i] == 0 and flagarr[j] == 0:
                pi = paths[i]
                pj = paths[j]
                if len(set(pi.edges.keys()) & set(pj.edges.keys())) > 0:
                    if pi.score < pj.score:
                        flagarr[j] = 1
                    else:
                        flagarr[i] = 1

    for i in range(len(flagarr)):
        if flagarr[i] == 0:
            yield paths[i]
 
class wDBG():
    """ Basic DBG with associated edge weights """
    def __init__(self, strings, k):
        # Only explicitly store edges
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

    def get_branch_statistics(self):
        """ Gets weights at each branch for modelling purposes """
        branch_weights = []
        for kmer, weight in self.edges.items():
            indegree_weights = []
            outdegree_weights = []
            for b in "ATCG":
                new_outkmer = kmer[1:] + b
                new_inkmer = b + kmer[:-1]
                if new_inkmer in self.edges:
                    indegree_weights.append(self.edges[new_inkmer])
                if new_outkmer in self.edges:
                    outdegree_weights.append(self.edges[new_outkmer])
            if len(outdegree_weights) > 1:
                branch_weights.append(outdegree_weights)
            if len(indegree_weights) > 1:
                branch_weights.append(indegree_weights)
        return sorted(branch_weights, key = lambda x:sum(x))

    def get_statistics(self):
        """ Gets statistics from the graph, such as: """
        """ Degree distribution, branch ratios, weight distribution"""
        """ Number of kmers, weight distribution """
        """ The should be more thinly spread out for coinfection """
        n_kmers = len(self.edges)
        weights = np.array([w for kmer, w in self.edges.items()])
        total_weight = sum(weights)
        mean_weight = np.mean(weights)
        std_weight = np.std(weights)
        degrees = []
        branch_ratios = []
        branch_losses = []
        for kmer, weight in self.edges.items():
            d = 0
            indegree_weights = []
            outdegree_weights = []
            for b in "ATCG":
                new_outkmer = kmer[1:] + b
                new_inkmer = b + kmer[:-1]
                if new_inkmer in self.edges:
                    d += 1
                    indegree_weights.append(self.edges[new_inkmer])
                if new_outkmer in self.edges:
                    d += 1
                    outdegree_weights.append(self.edges[new_outkmer])
            if len(outdegree_weights) > 1:
                branch_ratio = max(outdegree_weights)/sum(outdegree_weights)
                branch_ratios.append(branch_ratio)
                branch_losses.append(sum(outdegree_weights)-max(outdegree_weights))
            if len(indegree_weights) > 1:
                branch_ratio = max(indegree_weights)/sum(indegree_weights)
                branch_ratios.append(branch_ratio)
                branch_losses.append(sum(indegree_weights)-max(indegree_weights))
            degrees.append(d)
        mean_branch_loss = np.mean(branch_losses)
        std_branch_loss = np.std(branch_losses)
        max_branch_loss = max(branch_losses)
        mean_degree = np.mean(degrees)
        std_degree= np.std(degrees)
        mean_bratio = np.mean(branch_ratios)
        std_bratio = np.std(branch_ratios)
        return n_kmers, total_weight, mean_weight, std_weight, mean_degree, std_degree, mean_bratio, std_bratio, mean_branch_loss, std_branch_loss, max_branch_loss

    def get_start_positions(self):
        """ Gets positions in the graph that have no inward edge """
        for kmer, weight in self.edges.items():
            edge = True
            for b in "ATCG":
                if b+kmer[:-1] in self.edges:
                    edge = False
            if edge == True:
                self.start_positions.add(kmer)

    def cull(self, kmers):
        """ Removes any kmer in kmers from edges """
        todel = []
        for kmer in self.edges:
            if kmer not in kmers:
                todel.append(kmer)
        for kmer in todel:
             del self.edges[kmer]

    def get_paths(self):
        """ Traverses the wDBG heuristically, getting a list of paths """
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
                assert paths.count(path) == 1
                localswitch = False
                tmp_additions = {}
                for b in "ATCG":
                    # Check to see if the path has taken the route before
                    last_kmer = path.path[-1]
                    if b not in path.edges[last_kmer]:
                        tmpkmer = path.path[-1][1:] + b
                        if tmpkmer in self.edges:
                            tmp_additions[b] = self.edges[tmpkmer]
                            localswitch = True
                if localswitch == True:
                    maxtmp = max(tmp_additions.items(), key = lambda x:x[1])[1]
                    maxtmps = [t[0] for t in tmp_additions.items() if t[1] == maxtmp]
                    assert maxtmp != 0
                    best_tmp = maxtmps[0]
                    path.add_edge(path.path[-1][1:]+best_tmp, maxtmp)
                    tmp_paths.append(path)
                else:
                    final_paths.append(path)                
            paths = tmp_paths
            assert len(tmp_paths) <= len(paths)
            assert len(paths) <= len(self.start_positions)
            paths = tmp_paths
            sys.stderr.write("%d             \r" % (len(final_paths)))
        sys.stderr.write("\n")
        return final_paths          

    def get_kmer_frac(self, seq):
        c = 0
        for i in range(len(seq)-self.k+1):
            kmer = seq[i:i+self.k]
            if kmer in self.edges:
                c += 1
        return c/float(len(self.edges))

    def get_weighted_kmer_frac(self, seq):
        c = 0
        for i in range(len(seq)-self.k+1):
            kmer = seq[i:i+self.k]
            if kmer in self.edges:
                c += self.edges[kmer]
        total_weight = sum([w for e,w in self.edges.items()])
        return c/float(total_weight)

class cDBG():
    """ Basic DBG with additional color information """
    def __init__(self, k):
        self.edges = {}
        self.k = k
        self.n = None
        self.color_classes = set()
    
    @classmethod
    def from_strings(cls, strings, k):
        """ Builds the cDBG from strings """
        c = cls(k)
        c.n = len(strings)
        c._build(strings)
        return c

    @classmethod
    def from_strings_and_subgraph(cls, strings, k, subgraph):
        """ As before, but only for parts intersecting with subgraph """
        c = cls(k)
        c.n = len(strings)
        c._build(strings, subgraph)
        c.build_color_classes()
        return c

    def _build(self, strings, subgraph=False):
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

    def build_color_classes(self, subkmers=None):
        """ Builds the color classes using binary strings """
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

    def classify(self, wdbg, seqs, min_path_weight=100., cull_overlapping=True):
        """ Classifies a wdbg """

        # First get paths
        paths = [p for p in wdbg.get_paths() if p.score > min_path_weight]
        n_paths_initial = len(paths)
        if cull_overlapping == True:
            paths = [p for p in remove_overlaps(paths)]
        n_paths_surviving = len(paths)
        sys.stderr.write("Got %d paths from %d\n" % (n_paths_surviving, n_paths_initial))
        if len(paths) == 0:
            sys.stderr.write("No paths with a greater weight than %d found. Please try a lower threshold (-w) than %d \n" % (min_path_weight, min_path_weight))
            sys.exit(1)
        
        mpl = max([len(p.get_string()) for p in paths])
        sys.stderr.write("%d paths found\n" % len(paths))
        sys.stderr.write("Maximum path length: %d\n" % mpl)
        sys.stderr.write("Coloring paths...\n")
        total_colors = []
        # Next get arrays of color classes for paths
        for pi, path in enumerate(paths):
            sys.stderr.write(str(pi)+"        \r")
            colors = []   
            seq = path.get_string()
            kmers = (seq[i:i+self.k] for i in range(len(seq)-self.k+1))
            for kmer in kmers:
                if kmer in self.edges:
                    colors.append(self.edges[kmer])
            assert len(colors) > 0
            total_colors.append(colors)

        sys.stderr.write("Scoring colors...\n")
        # Lastly contiguize the colors for scoring
        # Need to do it individually for path, and then aggregate them
        path_scores = []
        for ci, colors in enumerate(total_colors):
            sys.stderr.write(str(ci)+"        \r")

            sumo = np.zeros(len(seqs))
            color_flag = np.zeros(self.n)
            for c in colors:
                arr = np.fromstring(np.binary_repr(c), dtype='S1').astype(int)[1:]
                sumo += (self.k * (1-color_flag) + color_flag) * arr
            path_scores.append(sumo)

        # Now aggregate the scores for each path and rank
        aggsumo = np.zeros(len(seqs))
        for sumo in path_scores:
            aggsumo += sumo
        maxs = max(aggsumo)
        # Take the highest indices, noting they are backwards wrt 
        maxi = [i for i in range(len(aggsumo)) if aggsumo[i] == maxs]            

        # Order of the classes is different 
        maxi_cls = [len(aggsumo)-i-1 for i in range(len(aggsumo)) if aggsumo[i] == maxs]       
        return maxs, maxi_cls


