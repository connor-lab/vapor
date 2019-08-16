DESCRIPTION

VAPOR is a tool for classification of Influenza samples from raw short read sequence data for downstream bioinformatics analysis. VAPOR is provided with a fasta file of full-length sequences (> 20,000) for a given segment, a set of reads, and attempts to retrieve a reference that is closest to the sample strain.

REQUIREMENTS

Python 3.x
NumPy >= 1.5.1

INSTALLATION

The VAPOR module executable script can be installed with pip:

    sudo pip3 install git+https://github.com/connor-lab/vapor

Alternatively, clone the repo, and add to your PATH:

    git clone https://github.com/connor-lab/vapor

TESTING

A test dataset is provided in the repository/tests. To test, clone the repo as above, and run:

    vapor.py -fq tests/test_reads.fq -fa tests/HA_sample.fa

which should yield:

    0.9782480893592005  186719.0    1701    109.77013521457965  1000    >cds:ADO12563 A/Chile/3935/2009 2009/07/07 HA H1N1 Human

Where the tab-delimited fields correspond to: approximate fraction of query bases found in reads; total score; query length; mean score; number of reads surviving culling; query description

USAGE:

    usage: vapor.py [-h] [-q] [-k K] [-s S] [-fa FA]
                              [-fq FQ [FQ ...]]

    optional arguments:
        -h, --help          Show this help message and exit
        -q, --quiet         Suppresses output to stderr
        --return_seqs       Returns a fasta of sequences, instead of hits       
        --return_best_n     Returns the highest scoring n queries

        -o                  Combined output to files with prefix O, none by default
        -k K                Kmer length [21]
        -t T                Pre-Filtering score threshold [0.2]
        -s S                Number of reads to sub-sample
        -c, --min_kmer_cov  Minimum kmer coverage for culling [5]
        -m, --min_kmer_prop
                            Minimum proportion of kmers required for query [0.1]
        -fa FA              Fasta file
        -fq FQ [FQ ...]     Fastq file/files, can be gzipped
        -f, --top_seed_frac
                            Fraction of best seeds to extend [0.2]
        --low_mem           Does not store reference kmer arrays, produces same result, marginally slower but less memory [False]
        -v, --version       Show version

Example:

    vapor.py -fa HA_sequences.fa -fq reads_1.fq.gz reads_2.fq.gz

PARAMETER OPTIMIZATION FOR OTHER VIRUSES:

VAPOR can in principle be used for other viruses (although performance has not yet been comprehensively benchmarked as with influenza, future versions will benchmark generalizability). For example, for ~79,448 HIV env sequences downloaded from https://www.hiv.lanl.gov/components/sequence/HIV/search/search.html, and public read sets with the ENA run accession SRR8389950, derived from HIV BF520.W14M (see https://www.ebi.ac.uk/ena/data/view/SRR8389950). A BF520 reference can be retrieved using the following parameters:

    vapor.py -c 100 --subsample 1000000 --low_mem -m 0.2 -f 0.1 -fq SRR8389950_1.fastq.gz SRR8389950_2.fastq.gz -fa env_db.fasta

Which returns:

    0.8030480656506448  235867226.0 2559    92171.63970300899   11899117    >A1.KE.1994.BF520.W14M.C2.KX168094

Mapping to this reference should result in a mismatch rate of < 2e-03%.

In this case we have customized the parameters:  

- Since the reference space is larger than for example, influenza A HA sequences, we use --low mem to reduce memory (this does not affect the result, but may increase run-time slightly)  

- We also assume, with the larger database, that there are sufficient close sequences to our sample, and use -m 0.2, requiring at least 20% exact matches to improve run-time (this does not affect the result as long as there are enough close references to the sample). If the reference space was very sparse, or our sample very novel, we may need to use -m 0.0  

- Again, due to the larger database, we decrease -f to 0.1 in order to extend fewer sequences with high-scoring exact matches.  

- Because the sample has over 12,000,000 reads, we improve run-time by subsampling 1,000,000 reads (--subsample 1000000). If depth is extremely skewed, sub-sampling may result in zero coverage in some sequence regions. In general, sub-sampling may decrease performance, and should be avoided where possible.

- Since depth is expected to be very high, we can also cull any k-mers with coverage less than 100 (assuming that they are, for example, errors or minor quasispecies variants). In some cases, this can affect perfomance (either increase or decrease), but will reduce memory and improve run-time as well. As with sub-sampling, this may also cull legitimate k-mers, and reduce performance, especially where depth is low.  

- If we had reason to believe our virus includes many k-mers that are also present in non-viral background sequences, we could increase -t, or -k, although the latter may have more implications for performance in general (sensitivity/specificity).

RECOMMENDED CPU AND MEMORY

For .fastq files with up to 10,000,000 reads, and influenza A segment .fasta files with approximately 47,000 references, VAPOR (single thread) requires < 4 Gb of memory and can generally be run within 4 minutes on an Intel(R) Core(TM) i7-6600U CPU @ 2.60GHz. For future datasets, memory requirements may be larger, and parameters may need to be optimized (see the above discussion of HIV). As such, we recommend 8 Gb of RAM for general usage, and any modern CPU.

Author: Joel Southgate

Email (for inquiries): southgateJA@cardiff.ac.uk

COPYRIGHT

Copyright 2018 Joel Southgate

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.


