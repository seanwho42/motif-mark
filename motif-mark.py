#!/usr/bin/env python

import argparse
import cairo
import re
from bioinfo import oneline_fasta


# ['#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f','#e5c494','#b3b3b3']



# handle argparse
def get_args():
    # TODO: make sure this is clear/descriptive
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-f", "--fasta", help = "Fasta file input with sequences", required = True)
    parser.add_argument("-m", "--motifs", help = "Text file with new motifs each on a new line", required=True)
    parser.add_argument("-o", "--out", help = "Output path for output png", required = True)
    # TODO: add in future support for SVGs as alternative since the assignment isn't chill with PNGs
    return parser.parse_args()

args = get_args()

class FastaRead:
    def __init__(self, seq):
        self.seq = seq
        self.segments = []
        self.motifs = []

    def get_motifs():
        # find motifs in the sequence
        pass

    def get_segments():
        # assign segments to the segments list
        pass

class Motif:
    def __init__(self):
        self.color = None

def main(fasta = args.fasta, motifs_file = args.motifs, out = args.out):
    #TODO: flush out docstring
    '''
    Docstring for main
    
    :param fasta: Description
    :param motifs: Description
    :param out: Description
    '''
    fasta = oneline_fasta(fasta)
    motifs = get_motifs(motifs_file)


def get_motifs(motifs_file):
    '''
    Docstring for get_motifs
    
    :param motifs_file: Description
    '''
    
    color_palette = ['#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f','#e5c494','#b3b3b3']
    motifs = {}

    with open(motifs_file, 'r') as m:
        for n, motif in enumerate(m):
            motif = motif.strip()
            motifs[motif] = color_palette[n]
    
    return motifs