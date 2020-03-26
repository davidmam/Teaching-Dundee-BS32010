# ex01.py
#
# Functions and data useful in exercise 1 (bacterial GC content, etc.) of
# the BS32010 course at the University of Dundee

from Bio import SeqIO  # For working with sequence data
from Bio.Graphics.ColorSpiral import get_color_dict  # For defining colours

import matplotlib.pyplot as plt  # For creating graphics

import pandas as pd  # For working with dataframes

import os  # For working with local files

bact_datadir = "genome_data/gc_content"
bact_files = {"Mycoplasma genitalium": ("NC_018495.fna",
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
bacteria = bact_files.keys()

unknown = pd.DataFrame([dict(species="Unknown", length=4391174, 
                             GC=0.656209, color=(1, 0.2, 0.2)), ])

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

def calc_size_gc(*names):
    """ When passed names corresponding to the bacteria
    listed in bact_files, returns a Pandas dataframe
    representing sequence length and GC content for
    each chromosome.
    DM: adapt to calculate over a genome expressed as contigs
    """
    # Use a Pandas DataFrame to hold data. Dataframes are 
    # useful objects/concepts, and support a number of 
    # operations that we will exploit later.
    df = pd.DataFrame(columns=['species', 'length', 'GC', 'color'])
    # Get one colour for each species, from Biopython's 
    # ColorSpiral module
    colors = get_color_dict(names, a=6, b=0.2)
    # Loop over the passed species names, and collect data
    for name in names:
        try:
            for filename in bact_files[name]:
                    
                seqfile = os.path.join(bact_datadir, filename)
                if os.path.exists(filename):
                    seqfile=filename
                ch = SeqIO.parse(seqfile, 'fasta')
                description = name
                ch_size = 0
                ch_g = 0
                ch_c = 0
                for s in ch:
                    ch_size += len(s.seq)
                    ch_g += s.seq.count('G')
                    ch_c += s.seq.count('C')
                    if s.description:
                        description = s.description
                ch_gc = float(ch_c + ch_g) / ch_size
                df = df.append(pd.DataFrame([dict(species=name, length=ch_size, 
                                                  GC=ch_gc, filename=filename,
                                                  description=description,
                                                  color=colors[name]), ]), 
                               ignore_index=True)
        except KeyError:
            print("Did not recognise species: %s" % name)
            continue
    return df


# Plot chromosome size and GC data
def plot_data(dataframe, filename=None, return_fig=False):
    """ When passed a dataframe corresponding to the output
    of calc_size_gc, renders a scatterplot of chromosome length
    against GC content.
    """
    # One advantage of using a Pandas dataframe is that we can
    # operate on the data by the content of the data. Here we're
    # treating the dataframe as a series of subsets on the basis
    # of named species. This allows us to label our scatterplot
    # by species, too.
    fig = plt.figure(figsize=(8,4))
    ax = fig.add_subplot(111)
    ax.set_position([0.15, 0.15, 0.45, 0.75])
    for k, sub in dataframe.groupby("species"):
        ax.scatter(x=sub.GC, y=sub.length, c=list(sub.color), label=k, s=50)
    ax.set_xlabel("GC content/%")
    ax.set_ylabel("chromosome length/bp")
    ax.set_title("Chr length vs GC%, grouped by species")
    leg = ax.legend(bbox_to_anchor=(1.0, 0.5), loc='center left')
    if filename is not None:
        fig.savefig(filename)
    if return_fig:
        return fig

