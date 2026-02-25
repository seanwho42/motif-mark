#!/usr/bin/env python

import argparse
import cairo
import re
from bioinfo import oneline_fasta


# ['#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f','#e5c494','#b3b3b3']

READ_DRAWING_HEIGHT = 40
KEY_HEIGHT = 30

# handle argparse
def get_args():
    parser = argparse.ArgumentParser(description="Visualizes motifs found across nucleotides, as well as which regions "\
                                     "are introns or exons. Note that this current iteration does not stagger motifs" \
                                     "to show overlapping of motifs. Provide input fasta file and motifs file.")
    parser.add_argument("-f", "--fasta", help = "Fasta file input with sequences", required = True)
    parser.add_argument("-m", "--motifs", help = "Text file with new motifs each on a new line", required=True)
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

        # print(motifs_colors)
        
        for pattern, color in motifs_colors.values():
            # span doesn't work because when finding regex matches using finditer
            # with lookahead, it doesn't capture the string itself, just the start

            # could make this more simple by passing the length to motif dictionary
            # instead of converting back to original length
            motif_len = len(re.sub(r"\[[ACTUG]+\]", "N", pattern))

            motif_matches = re.finditer(f"(?={pattern})", self.seq.upper())

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
    
    def draw_label(self):
        '''
        Returns RecordingSurface for read label
        
        :returns: returns RecordingSurface for read label
        '''
        cleaner_label = self.header[1:]
        bounds = cairo.Rectangle(0, 0, 1000, 30) # type: ignore
        label_rs = cairo.RecordingSurface(cairo.Content.COLOR, bounds)
        label_context = cairo.Context(label_rs)

        label_context.set_source_rgb(1, 1, 1)
        label_context.paint()
        label_context.set_source_rgb(0, 0, 0)
        
        # font styling
        label_context.set_font_size(15)
        label_context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)

        label_context.move_to(0, 20)
        label_context.show_text(cleaner_label)
        label_context.stroke()
        return label_rs



    def draw_segments(self):
        '''
        Returns RecordingSurface for segments with overlaid motifs
        
        :returns: returns RecordingSurface for read
        '''

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
            r, g, b = get_norm_rgb(motif.color)

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
        
        return segments_rs
    
    # TODO: implement this
    # def draw_motifs(self):
    #     '''
    #     Returns RecordingSurface
        
    #     :returns: returns RecordingSurface for staggered motifs
    #     '''
    #     read_width = len(self.seq)
    #     # TODO: write docstring -- patches together components of an individual read into one recording surface
    #     motifs_rs = cairo.RecordingSurface(cairo.Content.COLOR, None)
    #     motifs_context = cairo.Context(motifs_rs)

    #     motifs_context.set_source_surface(self.draw_segments(), 0, READ_DRAWING_HEIGHT)
    #     motifs_context.paint()
    #     return motifs_rs
    
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



def main(fasta = args.fasta, motifs_file = args.motifs):
    #TODO: flush out docstring
    '''
    Main function to read, visualize, and draw the image
    
    :param fasta: properly formatted fasta file name
    :param motifs: motifs file name (IUPAC ambiguous motifs deliminated by newlines)
    '''
    motifs_colors = get_motifs(motifs_file)

    reads = read_fasta(fasta, motifs_colors)
    
    max_width = 0
    fasta_height = KEY_HEIGHT
    for read in reads:
        x0, y0, width, label_height = read.draw_label().ink_extents()
        fasta_height += label_height
        # TODO: figure out unbound copy things so it isn't fixed at 1000 and we can include this
        # without it looking dumb... right now ludicrously long headers will get cut off.
        # if width > max_width:
        #     max_width = width
        x0, y0, width, segments_height = read.draw_segments().ink_extents()
        fasta_height += segments_height
        if width > max_width:
            max_width = width
    fasta_width = max_width + 20

    
    out = re.sub(r"\.(fa)|(fasta)$", "", fasta)
    with cairo.SVGSurface(f"{out}.svg", fasta_width, fasta_height) as surface:
        context = cairo.Context(surface)
        context.rectangle(0, 0, fasta_width, fasta_height)
        context.set_source_rgb(1, 1, 1)
        context.fill()

        current_height = 0
        for n, read in enumerate(reads):
            label_rc = read.draw_label()
            context.set_source_surface(label_rc, 10, current_height)
            context.paint()
            # add that height to track where to put next canvas
            x0, y0, width, label_height = label_rc.ink_extents()
            current_height += label_height
            # print(label_height)

            context.set_source_surface(read.draw_segments(), 10, current_height)
            context.paint()
            current_height += READ_DRAWING_HEIGHT
        context.set_source_surface(draw_key(motifs_dict = motifs_colors), 0, current_height)
        context.paint()
        surface.write_to_png(f"{out}.png")

def get_motifs(motifs_file):
    '''
    Gets motifs from motifs file and replaces characters with appropriate
    regex pattern to match (upper case only) to nucleotides in fasta file.
    
    :param motifs_file: Motifs file with each motif deliminated by new lines
    :returns: motifs dictionary {motif: (regex_pattern, hex_color)}
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
            motif = motif.strip()
            # make it all upper case -- comparison later will make the FASTA nucleotides lower case to match
            pattern = motif.upper()
            for nuc in iupac_subs.keys():
                pattern = re.sub(nuc, iupac_subs[nuc], pattern)
            motifs[motif] = (pattern, color_palette[n])
    # print(motifs)
    return motifs

def read_fasta(fasta, motifs_colors):
    '''
    Reads in fasta and saves them as FastaRead objects

    :returns: list if FastaRead objects
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
        # for read in reads:
        #     read.find_motifs(motifs)
        #     read.find_segments()
        #     # for segment in read.segments:
        #     #     print(segment.start_bp)
        #     #     print(segment.end_bp)
        #     # print(read.motifs)
    return reads

def get_norm_rgb(hex_code):
    # pycairo requires normalized rgb colors, and it felt like that shouldn't be the case, so here we go
    '''
    Takes hexidecimal color (i.e. #ffffff) and returns normalized r, g, b values (i.e. 1, 1, 1)

    :param hex_code: hexidemical color as a string
    :returns: tuple (r, g, b)
    '''
    r = int(hex_code[1:3], 16)/255
    g = int(hex_code[3:5], 16)/255
    b = int(hex_code[5:], 16)/255
    # print(f"{r}, {g}, {b}")
    return r, g, b

def draw_key(motifs_dict):
    '''
    Takes motifs dictionary, returns RecordingSurface to be added to bottom of the visualization
     
    :param motifs_dict: motifs dictionary

    :returns: returns RecordingSurface for motifs key
    '''
    # can't figure out how to make the RecordingSurfaces work properly when unbounded.. using hard coded size for this for now
    bounds = cairo.Rectangle(0, 0, 1000, KEY_HEIGHT) # type: ignore
    key_rs = cairo.RecordingSurface(cairo.Content.COLOR, bounds)
    key_context = cairo.Context(key_rs)

    key_context.set_source_rgb(1, 1, 1)
    key_context.paint()
    key_context.set_source_rgb(0, 0, 0)

    # make a separator line
    key_context.move_to(0, 0)
    key_context.set_line_width(2)
    key_context.line_to(1000, 0)
    key_context.stroke()
    


    # now draw the colors and corresponding motifs
    key_context.scale(KEY_HEIGHT, KEY_HEIGHT) # makes it easier to have standard square sizes
    key_context.set_line_width(0.5)
    # start at 10 for some padding
    current_x = 10/KEY_HEIGHT

    # font styling
    key_context.set_font_size(0.5)
    key_context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)

    for motif, (pattern, color) in motifs_dict.items():
        r, g, b = get_norm_rgb(color)
        key_context.set_source_rgb(r, g, b)
        key_context.move_to(current_x, 0.5)
        new_x = current_x + 0.5
        key_context.line_to(new_x, 0.5)
        current_x = new_x + 0.05
        key_context.stroke()

        key_context.set_source_rgb(0, 0, 0)
        key_context.move_to(current_x, 0.75)
        key_context.show_text(motif)
        key_context.stroke()
        # now stagger to the next
        (x, y, width, height, dx, dy) = key_context.text_extents(motif)
        current_x += width + 0.5
    return key_rs

main()
