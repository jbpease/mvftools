#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This program takes a DNA MVF alignment and annotates the output into
gene boudaries.
"""

import os
import sys
import argparse
import re
from copy import deepcopy
from mvfbase import MultiVariantFile

_LICENSE = """
MVFtools: Multisample Variant Format Toolkit
James B. Pease and Ben K. Rosenzweig
http://www.github.org/jbpease/mvftools

If you use this software please cite:
Pease JB and BK Rosenzweig. 2016.
"Encoding Data Using Biological Principles: the Multisample Variant Format
for Phylogenomics and Population Genomics"
IEEE/ACM Transactions on Computational Biology and Bioinformatics. In press.
http://www.dx.doi.org/10.1109/tcbb.2015.2509997

This file is part of MVFtools.

MVFtools is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
MVFtools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MVFtools.  If not, see <http://www.gnu.org/licenses/>.
"""

RE_GENEID = re.compile("ID=gene:(.*?);")


def parse_gff(gff_file, contigs, filter_annotation=None):
    """Parses a GFF3 file for exon locations
        Arguments:
            gff_file: path to GFF3 file
            args: passthrough from main args

        Output: triplets for codon locations

    """
    gff_entries = {}
    relabeled_gff_entries = {}
    geneids = {}
    geneid = 0
    with open(gff_file) as gff:
        for line in gff:
            if line[0] == '#':
                continue
            arr = line.rstrip().split()
            if arr[2] != 'gene':
                continue
            if filter_annotation:
                if filter_annotation in arr[8]:
                    continue
            genename = re.findall(RE_GENEID, arr[8])[0]
            contig = arr[0]
            if contig not in gff_entries:
                gff_entries[contig] = {}
            coords = [int(arr[3]), int(arr[4])]
            if genename not in geneids:
                geneids[geneid] = {'label': genename,
                                   'length': max(coords) - min(coords)}
            for j in range(min(coords), max(coords) + 1):
                gff_entries[contig][j] = geneid
            geneid += 1
    for contig in gff_entries:
        matchlabel = False
        for contigid in contigs:
            if contigs[contigid]['label'] == contig:
                relabeled_gff_entries[contigid] = gff_entries[contig].copy()
                matchlabel = True
                break
        if matchlabel is False:
            relabeled_gff_entries[contig] = gff_entries[contig].copy()
    gff_entries = None
    return relabeled_gff_entries, geneids


def generate_argparser():
    parser = argparse.ArgumentParser(
        prog="mvf_annotate.py",
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=_LICENSE)
    parser.add_argument("-i", "--mvf", type=os.path.abspath,
                        help="Input MVF file.")
    parser.add_argument("-g", "--gff", type=os.path.abspath,
                        help="Input gff annotation file.")
    parser.add_argument("-o", "--out", type=os.path.abspath,
                        help="Output annotated MVF file")
    parser.add_argument("-F", "--filter_annotation",
                        help=("Skip entries in the GFF file that "
                              "contain this string in their 'Notes'"))
    parser.add_argument("-N", "--nongenic-mode", action="store_true",
                        help=("Instead of returning annotated genes, "
                              "return the non-genic regions without "
                              "without changing contigs or coordinates"))
    parser.add_argument("-M", "--nongenic-margin", type=int, default=0,
                        help=("for --unnanotated-mode, only retain "
                              "positions that are this number of bp away "
                              "from an annotated region boundary"))
    parser.add_argument("-B", "--linebuffer", type=int, default=100000,
                        help="Number of entries to store in memory at a time.")
    parser.add_argument("--overwrite", action="store_true",
                        help="USE WITH CAUTION: force overwrite of outputs")
    parser.add_argument("--quiet", action="store_true",
                        help="suppress progress meter")
    parser.add_argument("-v", "--version", action="version",
                        version="2017-06-24",
                        help="display version information")
    return parser


def main(arguments=None):
    """Main method"""
    arguments = arguments if arguments is not None else sys.argv[1:]
    parser = generate_argparser()
    args = parser.parse_args(args=arguments)
    mvf = MultiVariantFile(args.mvf, 'read')
    gff, geneids = parse_gff(args.gff, mvf.metadata['contigs'])
    if args.quiet is False:
        print("gff_processed")
    outmvf = MultiVariantFile(args.out, 'write', overwrite=args.overwrite)
    outmvf.metadata = deepcopy(mvf.metadata)
    if args.nongenic_mode is False:
        outmvf.metadata['contigs'] = geneids
    outmvf.write_data(outmvf.get_header())
    entrybuffer = []
    nentry = 0
    for contigid, pos, allelesets in mvf.iterentries(decode=False):
        annotated_pos = False
        if contigid in gff:
            if pos in gff[contigid]:
                annotated_pos = True
            elif args.nongenic_mode is True and args.unmargin > 0:
                for xpos in range(pos - args.unmargin,
                                  pos + args.unmargin + 1):
                    if xpos in gff[contigid]:
                        annotated_pos = True
                        break
        if args.nongenic_mode is False and annotated_pos is True:
            entrybuffer.append((gff[contigid][pos], pos, allelesets))
            nentry += 1
            if nentry == args.linebuffer:
                outmvf.write_entries(entrybuffer)
                entrybuffer = []
                nentry = 0
        elif args.nongenic_mode is True and annotated_pos is False:
            entrybuffer.append((contigid, pos, allelesets))
            nentry += 1
            if nentry == args.linebuffer:
                outmvf.write_entries(entrybuffer)
                entrybuffer = []
                nentry = 0
    if entrybuffer:
        outmvf.write_entries(entrybuffer)
        entrybuffer = []
        nentry = 0
    return ''


if __name__ == "__main__":
    main()
