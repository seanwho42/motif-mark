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
    def __init__(self, header, seq):
        self.header = header
        self.seq = seq
        self.segments = []
        self.motifs = None

    def find_motifs(self, motifs):
        # find motifs in the sequence

        # keys are the motif color, values are tuples of their x position span (in base pairs)
        motifs_positions = {}
        
        for motif in motifs.keys():
            # span doesn't work because when finding regex matches using finditer
            # with lookahead, it doesn't capture the string itself, just the start

            # could make this more simple by passing the length to motif dictionary
            # instead of converting back to original length
            motif_len = len(re.sub(r"\[[ACTUG]+\]", "N", motif))

            color = motifs[motif]
            motif_matches = re.finditer(f"(?={motif})", self.seq.upper())

            for match in motif_matches:
                # initialize 
                x_span = (match.start(), match.start() + motif_len)
                if color in motifs_positions.keys():
                    motifs_positions[color].append(x_span)
                else:
                    motifs_positions[color] = [x_span]

        self.motifs = motifs_positions

    def get_segments(self):
        # assign segments to the segments list
        pass

class Motif:
    def __init__(self, color, start_bp, end_bp):
        self.color = color
        self.start_bp = start_bp
        self.end_bp = end_bp
        self.start_x = None
        self.end_x = None

class Segment:
    def __init__(self, start_bp, end_bp):
        self.is_exon = None
        self.start_bp = start_bp
        self.end_bp = end_bp
        self.start_x = None
        self.end_x = None

    def set_x_positions(self, bp):
        '''
        Docstring for set_x_positions
        
        :param self: Description
        :param bp: Description
        '''
        # TODO: make this figure it out within the context of the longest read -- it should be called after all have been read in
        x = bp
        return

class Intron(Segment):
    pass

class Exon(Segment):
    pass


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
    reads = []

    with open(fasta, "r") as f:
        max_length = 0
        for n, line in enumerate(f):
            line = line.strip()
            if n % 2 == 0:
                # header line
                header = line
            else:
                reads.append(FastaRead(header, line))

                # set max read length for purposes of setting the canvas width
                read_length = len(line)

                if read_length > max_length:
                    max_length = read_length
        for read in reads:
            read.find_motifs(motifs)
            print(read.motifs)
        



def get_motifs(motifs_file):
    '''
    Gets motifs from motifs file and replaces characters with appropriate
    regex pattern to match (upper case only) to nucleotides in fasta file.
    
    :param motifs_file: Description
    '''

    color_palette = ['#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f','#e5c494','#b3b3b3']
    motifs = {}
    # list of substitutes for regex patterns to properly match
    iupac_subs = {
        "[TU]" : "[TU]",
        "R" : "[AG]",
        "Y" : "[CTU]",
        "S" : "[GC]",
        "W" : "[ATU]",
        "K" : "[GTU]",
        "M" : "[AC]",
        "B" : "[CGTU]",
        "D" : "[AGTU]",
        "H" : "[ACTU]",
        "V" : "[ACG]",
        "N" : "[ACTUG]"
        }

    with open(motifs_file, 'r') as mf:
        for n, motif in enumerate(mf):
            # make it all upper case -- comparison later will make the FASTA nucleotides lower case to match
            motif = motif.strip().upper()
            for nuc in iupac_subs.keys():
                motif = re.sub(nuc, iupac_subs[nuc], motif)
            motifs[motif] = color_palette[n]
    print(motifs)
    return motifs

main()