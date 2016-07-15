#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
MVFtools: Multisample Variant Format Toolkit
http://www.github.org/jbpease/mvftools (Stable Releases)
http://www.github.org/jbpease/mvftools-dev (Latest Testing Updates)

If you use this software please cite:
Pease JB and BK Rosenzweig. 2016.
"Encoding Data Using Biological Principles: the Multisample Variant Format
for Phylogenomics and Population Genomics"
IEEE/ACM Transactions on Computational Biology and Bioinformatics. In press.
http://www.dx.doi.org/10.1109/tcbb.2015.2509997

MVF_translate: Translate MVF dna to protein or codon format
@author: James B. Pease

version: 2015-02-05 - MVF1.2 update
version: 2015-09-04 - minor fixes and style cleanup
version: 2015-12-31 - updates to header and minor fixes
@version: 2016-01-01 - added GFF-less mode for de novo alignments

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

from __future__ import print_function
import sys
import argparse
import re
from copy import deepcopy
from mvfbase import MultiVariantFile
from mvfbiolib import FULL_CODON_TABLE, COMPCODE


PARENTGENE = re.compile("Parent=mRNA:(.*?);")


def crop_to_stop(seq, firststop=""):
    """Crops sequence to first stop codon
        Arguments:
            firststop: none=keep all characters
                       inclusive=all characters up to and including '*'
                       exclusive=all charaters up to '*'
    """
    stop_index = len(seq)
    for i in range(int(len(seq) / 3.0) + 1):
        codon = ''.join(seq[i * 3:(i + 1)*3]).upper().replace('U', 'T')
        if not codon:
            break
        if FULL_CODON_TABLE.get(codon, 'X') == '*':
            stop_index = 3 * (i + int(firststop == "inclusive"))
            break
    return seq[:stop_index]


def translate(seq, firststop=None):
    """Returns translated amino acids from nucleotides

        Arguments:
            firststop: none=keep all characters
                       inclusive=all characters up to and including '*'
                       exclusive=all charaters up to '*'
    """
    aa_seq = []
    for i in range(int(len(seq) / 3.0) + 1):
        codon = ''.join(seq[i * 3:(i + 1)*3]).upper().replace('U', 'T')
        if not codon:
            break
        if codon == '---':
            aa_seq.append('-')
        else:
            aa_seq.append(FULL_CODON_TABLE.get(codon, 'X'))
    if firststop:
        aa_seq = '*' in aa_seq and aa_seq[:aa_seq.index('*') + int(
            firststop == "inclusive")] or aa_seq
    return aa_seq


def parse_gff(gff_file, args):
    """Parses a GFF3 file for exon locations
        Arguments:
            gff_file: path to GFF3 file
            args: passthrough from main args

        Output: triplets for codon locations

    """
    gff_entries = {}
    gff_triplets = {}
    with open(gff_file) as gff:
        for line in gff:
            if line[0] == '#':
                continue
            arr = line.rstrip().split()
            if arr[2] != 'CDS':
                continue
            if args.filter_annotation:
                if args.filter_annotation in arr[8]:
                    continue
            parent = re.findall(PARENTGENE, arr[8])[0]
            if arr[0] not in gff_entries:
                gff_entries[arr[0]] = {}
                gff_triplets[arr[0]] = []
            coords = [int(arr[3]), int(arr[4])]
            strand = arr[6]
            if parent not in gff_entries[arr[0]]:
                gff_entries[arr[0]][parent] = [strand, []]
            gff_entries[arr[0]][parent][1].extend(range(min(coords),
                                                        max(coords) + 1))
    for contigname in gff_entries:
        for gene in gff_entries[contigname]:
            if len(gff_entries[contigname][gene][1]) % 3:
                continue
            strand = gff_entries[contigname][gene][0]
            coords = sorted(gff_entries[contigname][gene][1])
            for j in range(0, len(coords), 3):
                try:
                    gff_triplets[contigname].append((coords[j], coords[j+1],
                                                     coords[j+2], strand))
                except:
                    raise RuntimeError(len(coords), j, strand,
                                       contigname, coords[j])
    gff_entries = None
    return gff_triplets


def iter_codons(inputbuffer, mvf):
    print(inputbuffer)
    for i in range(0, len(inputbuffer), 3):
        alleles = [inputbuffer[i][1][0],
                   inputbuffer[i+1][1][0],
                   inputbuffer[i+2][1][0]]
        print(alleles)
        if all(len(x) == 1 for x in alleles):
            amino_acids = translate(''.join(alleles))[0]
        else:
            decoded_alleles = [mvf.decode(x) for x in alleles]
            print(decoded_alleles)
            amino_acids = [translate(''.join(x))
                           for x in zip(*decoded_alleles)]
            amino_acids = mvf.encode(
                ''.join([x[0] for x in amino_acids]))
        yield inputbuffer[i][0], amino_acids, [
            inputbuffer[i][1][0], inputbuffer[i+1][1][0],
            inputbuffer[i+2][1][0]]


def main(arguments=sys.argv[1:]):
    """Main method for mvf filter"""

    parser = argparse.ArgumentParser(description="""
    Filters and Transforms MVF""")
    parser.add_argument("--mvf", help="Input MAF file", required=True)
    parser.add_argument("--gff", help="""
        Input GFF3 file. If GFF3 not provided, alignments are assumed
        to be in-frame as input.""")
    parser.add_argument("--out", help="Output MAF file", required=True)
    parser.add_argument("--outtype", choices=['protein', 'codon'],
                        default="codon",
                        help="""single data column protein alleles,
                                or four column protein frame1 frame2 frame3""")
    parser.add_argument("--filter_annotation",
                        help="""skip GFF entries with text
                                matching this in their 'Notes' field""")
    parser.add_argument("--linebuffer", type=int, default=100000,
                        help="number of entries to write in a block")
    parser.add_argument("--overwrite", action="store_true",
                        help="USE WITH CAUTION: force overwrite of outputs")
    parser.add_argument("--quiet", action="store_true",
                        help="suppress progress meter")
    parser.add_argument("-v", "--version", action="store_true",
                        help="display version information")
    args = parser.parse_args(args=arguments)
    if args.version:
        print("Version 2015-12-31")
        sys.exit()
    mvf = MultiVariantFile(args.mvf, 'read')
    if not mvf.metadata['flavor'] == 'dna':
        raise RuntimeError("MVF must be flavor=dna to translate")
    if args.gff:
        gff = parse_gff(args.gff, args)
        if not args.quiet:
            print("gff_processed")
    outmvf = MultiVariantFile(args.out, 'write', overwrite=args.overwrite)
    outmvf.metadata = deepcopy(mvf.metadata)
    outmvf.metadata['flavor'] = args.outtype
    outmvf.write_data(outmvf.get_header())
    entrybuffer = []
    nentry = 0
    if not args.gff:
        inputbuffer = []
        current_contig = ''
        for contigid, pos, allelesets in mvf.iterentries(decode=False):
            if current_contig == '':
                current_contig = contigid[:]
            if contigid == current_contig:
                inputbuffer.append((pos, allelesets))
            else:
                for coord, amino_acids, alleles in iter_codons(
                            inputbuffer, mvf):
                    if all([x in '-X' for x in amino_acids]):
                        continue
                    if args.outtype == 'protein':
                        entrybuffer.append(
                            (current_contig, pos, (amino_acids,)))
                    else:
                        entrybuffer.append((
                            current_contig, pos, (
                                amino_acids, alleles[0],
                                alleles[1], alleles[2])))
                    nentry += 1
                    if nentry == args.linebuffer:
                        outmvf.write_entries(entrybuffer)
                        entrybuffer = []
                        nentry = 0
                inputbuffer = [(pos, allelesets)]
                current_contig = contigid[:]
        if inputbuffer:
            for coord, amino_acids, alleles in iter_codons(
                        inputbuffer, mvf):
                if all([x in '-X' for x in amino_acids]):
                    continue
                if args.outtype == 'protein':
                    entrybuffer.append(
                        (current_contig, pos, (amino_acids,)))
                else:
                    entrybuffer.append((
                        current_contig, pos, (
                            amino_acids, alleles[0],
                            alleles[1], alleles[2])))
                nentry += 1
                if nentry == args.linebuffer:
                    outmvf.write_entries(entrybuffer)
                    entrybuffer = []
                    nentry = 0
    else:
        mvf_entries = {}
        for contigid, pos, allelesets in mvf.iterentries(decode=False):
            if contigid not in mvf_entries:
                mvf_entries[contigid] = {}
            mvf_entries[contigid][pos] = allelesets[0]
        for contigname in sorted(gff):
            contigid = mvf.get_contig_id(contigname)
            for coords in sorted(gff[contigname]):
                reverse_strand = False
                if coords[3] == '-':
                    reverse_strand = True
                    alleles = [mvf_entries[contigid].get(x, '-')
                               for x in coords[2::-1]]
                else:
                    alleles = [mvf_entries[contigid].get(x, '-')
                               for x in coords[0:3]]
                if all(len(x) == 1 for x in alleles):
                    if reverse_strand:
                        alleles = [COMPCODE[x] for x in alleles]
                    decoded_alleles = alleles
                    amino_acids = translate(''.join(alleles))[0]
                else:
                    if reverse_strand:
                        decoded_alleles = [[COMPCODE[y] for y in mvf.decode(x)]
                                           for x in alleles]
                        alleles = [mvf.encode(''.join(x))
                                   for x in decoded_alleles]
                    else:
                        decoded_alleles = [mvf.decode(x) for x in alleles]
                    amino_acids = [translate(''.join(x))
                                   for x in zip(*decoded_alleles)]
                    amino_acids = mvf.encode(''.join([x[0]
                                                      for x in amino_acids]))
                if all([x in '-X' for x in amino_acids]):
                    continue
                if args.outtype == 'protein':
                    entrybuffer.append((contigid, coords[0], (amino_acids,)))
                else:
                    entrybuffer.append((
                        contigid, coords[0], (
                            amino_acids, alleles[0], alleles[1], alleles[2])))
                nentry += 1
                if nentry == args.linebuffer:
                    outmvf.write_entries(entrybuffer)
                    entrybuffer = []
                    nentry = 0
    if entrybuffer:
        print(entrybuffer)
        outmvf.write_entries(entrybuffer)
        entrybuffer = []
        nentry = 0
    return ''

if __name__ == "__main__":
    main()
