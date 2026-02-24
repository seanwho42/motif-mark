#!/usr/bin/env python

import argparse
import cairo
import re
from bioinfo import oneline_fasta


# ['#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f','#e5c494','#b3b3b3']

READ_DRAWING_HEIGHT = 40

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
    def __init__(self, header, seq, motifs_colors):
        self.header = header
        self.seq = seq
        self.motifs_colors = motifs_colors
        self.segments = self.find_segments()
        self.motifs = self.find_motifs(motifs_colors)

    def find_motifs(self, motifs_colors):
        # find motifs in the sequence

        # keys are the motif color, values are tuples of their x position span (in base pairs)
        found_motifs = []
        
        for motif in motifs_colors.keys():
            # span doesn't work because when finding regex matches using finditer
            # with lookahead, it doesn't capture the string itself, just the start

            # could make this more simple by passing the length to motif dictionary
            # instead of converting back to original length
            motif_len = len(re.sub(r"\[[ACTUG]+\]", "N", motif))

            color = motifs_colors[motif]
            motif_matches = re.finditer(f"(?={motif})", self.seq.upper())

            for match in motif_matches:
                # initialize 
                start_bp, end_bp = (match.start(), match.start() + motif_len)
                found_motifs.append(Motif(color, start_bp, end_bp))
        return found_motifs

    def find_segments(self):
        # assign segments to the segments list
        # TODO: test to see if this works properly with multiple introns/exons
        segments = []
        segment_matches = re.findall(r"([actug]+|[ACTUG]+)", self.seq)
        bp = 0
        for segment_str in segment_matches:
            segment_bp = len(segment_str)
            # oneline_fasta() should validate that these are valid fastas? if not, add that. 
            if re.search(r"[actug]", segment_str) is None:
                is_exon = True
            else:
                is_exon = False
            # TODO: check to see how this translates to the drawings
            segment = Segment(bp, bp + segment_bp, is_exon)
            segments.append(segment)
            bp += len(segment_str)
        return segments
    
    def draw_read(self):
        # TODO: add docstring for this
        
        read_width = len(self.seq)
        bounds = cairo.Rectangle(0, 0, read_width, READ_DRAWING_HEIGHT) # type: ignore
        # segments_rs = cairo.RecordingSurface(cairo.Content.COLOR, bounds)
        segments_rs = cairo.RecordingSurface(cairo.Content.COLOR, bounds)
        seg_context = cairo.Context(segments_rs)
        # scale so the line weights and y are relative to the frame of that read as a whole
        seg_context.set_source_rgb(1, 1, 1)
        seg_context.rectangle(0, 0, read_width, READ_DRAWING_HEIGHT)
        seg_context.fill()

        seg_context.scale(1, READ_DRAWING_HEIGHT)

        # drawing motifs first so they don't cover the segments
        seg_context.set_line_width(0.8)
        for motif in self.motifs:
            # have to convert hex codes into normalized r g b
            hex_code = motif.color
            r = int(hex_code[1:3], 16)/255
            g = int(hex_code[3:5], 16)/255
            b = int(hex_code[5:], 16)/255
            # print(f"{r}, {g}, {b}")

            seg_context.move_to(motif.start_bp, 0.5)
            seg_context.set_source_rgb(r, g, b)
            seg_context.line_to(motif.end_bp, 0.5)
            seg_context.stroke()
        
        # drawing segments
        seg_context.set_source_rgb(0, 0, 0)
        for segment in self.segments:
            if segment.is_exon:
                seg_context.set_line_width(0.6)
            else:
                seg_context.set_line_width(0.1)
            seg_context.move_to(segment.start_bp, 0.5)
            seg_context.line_to(segment.end_bp, 0.5)
            seg_context.stroke()
        
        # # drawing staggered 
        # read_rs = cairo.RecordingSurface(cairo.Content.COLOR, bounds)
        # read_context = cairo.Context(segments_rs)
        # # scale so the line weights and y are relative to the frame of that read as a whole
        # read_context.set_source_rgb(1, 1, 1)
        # read_context.rectangle(0, 0, read_width, READ_DRAWING_HEIGHT)
        # read_context.fill()
        
        return segments_rs

class Motif:
    def __init__(self, color, start_bp, end_bp):
        self.color = color
        self.start_bp = start_bp
        self.end_bp = end_bp


class Segment:
    def __init__(self, start_bp, end_bp, is_exon):
        self.is_exon = is_exon
        self.start_bp = start_bp
        self.end_bp = end_bp


def main(fasta = args.fasta, motifs_file = args.motifs, out = args.out):
    #TODO: flush out docstring
    '''
    Docstring for main
    
    :param fasta: Description
    :param motifs: Description
    :param out: Description
    '''
    motifs_colors = get_motifs(motifs_file)

    reads, max_length = read_fasta(fasta, motifs_colors)
        
    fasta_width = max_length + 20
    fasta_height = len(reads) * READ_DRAWING_HEIGHT
        
    with cairo.SVGSurface(f"{out}.svg", fasta_width, fasta_height) as surface:
        context = cairo.Context(surface)
        context.rectangle(0, 0, fasta_width, fasta_height)
        context.set_source_rgb(1, 1, 1)
        context.fill()
        for n, read in enumerate(reads):
            context.set_source_surface(read.draw_read(), 10, n * READ_DRAWING_HEIGHT)
            context.paint()
        surface.write_to_png(out)




def get_motifs(motifs_file):
    # TODO: finish docstring
    '''
    Gets motifs from motifs file and replaces characters with appropriate
    regex pattern to match (upper case only) to nucleotides in fasta file.
    
    :param motifs_file: Description of motifs file parameter
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
    # print(motifs)
    return motifs

def read_fasta(fasta, motifs_colors):
    # TODO: write docstring -- returns reads
    '''
    '''
    reads = []
    fasta = oneline_fasta(fasta)
    with open(fasta, "r") as f:
        max_length = 0
        for n, line in enumerate(f):
            line = line.strip()
            if n % 2 == 0:
                # header line
                header = line
            else:
                reads.append(FastaRead(header, line, motifs_colors))

                # set max read length for purposes of setting the canvas width
                read_length = len(line)

                if read_length > max_length:
                    max_length = read_length
        # for read in reads:
        #     read.find_motifs(motifs)
        #     read.find_segments()
        #     # for segment in read.segments:
        #     #     print(segment.start_bp)
        #     #     print(segment.end_bp)
        #     # print(read.motifs)
    return reads, max_length

main()
