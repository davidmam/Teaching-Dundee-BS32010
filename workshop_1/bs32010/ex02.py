# ex02.py
#
# Functions and data useful in exercise 2 (k-mer spectra) of
# the BS32010 course at the University of Dundee

import pandas as pd                
from collections import defaultdict
import os
from Bio import SeqIO

bact_datadir = "genome_data/gc_content"
files = {"Mycoplasma genitalium": ("NC_018495.fna",
                                   "NC_018496.fna",
                                   "NC_018497.fna",
                                   "NC_018498.fna"),
         "Mycoplasma pneumoniae": ("NC_000912.fna",
                                   "NC_016807.fna", 
                                   "NC_017504.fna",
                                   "NC_020076.fna"),
         "Nostoc punctiforme": ("NC_010628.fna",),
         "Escherichia coli": ("NC_000913.fna",
                              "NC_002695.fna",
                              "NC_004431.fna",
                              "NC_010468.fna"),
         "Mycobacterium tuberculosis": ("NC_016934.fna",
                                        "NC_017523.fna",
                                        "NC_022350.fna",
                                        "NC_000962.fna")}
bacteria = files.keys()

bact_files = {}
for k, v in files.items():
    bact_files[k] = tuple([os.path.join(bact_datadir, fn) for fn in v])

def add_genome(filename, species):
    ''' when passed a filename and species, checks to see if the file is
    a sequence file and if so, adds it to the list of bacterial files'''
    if os.path.exists(filename):
        try:
            ch = SeqIO.parse(filename, 'fasta')
            
            if species in bact_files:
                bact_files[species] = tuple(set(list(bact_files[species])+[filename]))
            else:
                bact_files[species] = tuple([filename])
            print('file {} added to species {}'.format(filename, species))
        except Exception as e:
            print('Cannot read {} {}'.format(filename,e ))
    else:
        print('cannot find file {}'.format(filename))

def count_genome_kmers(filename, k):
    '''This method counts kmers across all sequences in a fasta file.
    It will check if the file exists and if not search in bact_datadir.
    
    It calls count_str_kmers(seq, k) for each sequence and merges the results.
    These are returned as a Pandas dataframe'''
    seqfile = os.path.join(bact_datadir, filename)
    if os.path.exists(filename):
        seqfile = filename
        print('found file {}'.format(filename))
    try:
        sr = SeqIO.parse(seqfile, 'fasta')
        kdict = defaultdict(int)
        for inseq in sr:
            kdict = count_str_kmers(str(inseq.reverse_complement().seq), k, kdict)
            kdict = count_str_kmers(str(inseq.seq), k, kdict)
        df = pd.DataFrame.from_dict(kdict, orient="index")
        df.columns = ("frequency",)
        return df
    except Exception as e:
        print('Error calculating kmers for {}: {}'.format(filename, e))
    

def count_str_kmers(instr, k, kdict=None):
    """Counts sequences of size k in instr, populating kdict.
    
    Loops over instr with a window of size k, populating the
    dictionary kdict with a count of occurrences of each k-mer.
    Returns the dictionary kdict.
    """
    if kdict is None:
        kdict = defaultdict(int)
    for idx in range(len(instr)-k):
        kdict[instr[idx:idx+k]] += 1
    return kdict

def count_seq_kmers(inseq, k):
    """Counts kmers of size k in the sequence inseq.
    
    Counts kmers in forward and reverse directions, returning
    a Pandas dataframe of k-mer and count.
    """
    kdict = count_str_kmers(str(inseq.seq), k)
    kdict = count_str_kmers(str(inseq.reverse_complement().seq), k, kdict)
    df = pd.DataFrame.from_dict(kdict, orient="index")
    df.columns = ("frequency",)
    return df
