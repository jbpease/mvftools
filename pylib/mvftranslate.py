# -*- coding: utf-8 -*-
"""
This program takes a DNA MVF alignment and annotates the output into
gene boudaries.

MVFtools: Multisample Variant Format Toolkit
James B. Pease and Ben K. Rosenzweig
http://www.github.org/jbpease/mvftools

If you use this software please cite:
Pease JB and BK Rosenzweig. 2015.
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

import re
from copy import deepcopy
from pylib.mvfbase import MultiVariantFile
from pylib.mvfbiolib import MvfBioLib
MLIB = MvfBioLib()

RE_GENEID = re.compile("ID=gene:(.*?);")
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
        if MLIB.codon_tables['full'].get(codon, 'X') == '*':
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
            aa_seq.append(MLIB.codon_tables['full'].get(codon, 'X'))
    if firststop:
        aa_seq = '*' in aa_seq and aa_seq[:aa_seq.index('*') + int(
            firststop == "inclusive")] or aa_seq
    return aa_seq


def iter_codons(inputbuffer, mvf):
    """Iterate through codons
    """

    for i in range(0, len(inputbuffer), 3):
        alleles = [inputbuffer[i][1][0],
                   inputbuffer[i+1][1][0],
                   inputbuffer[i+2][1][0]]
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


def parse_gff_translate(gff_file, args):
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
                except IndexError:
                    raise RuntimeError(len(coords), j, strand,
                                       contigname, coords[j])
    gff_entries = None
    return gff_triplets


def parse_gff_annotate(gff_file, contigs, filter_annotation=None):
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


def parse_gff_analysis(gffpath):
    """Parse GFF function used for the analsis script"""
    rxpr3 = re.compile(r' \(AHRD.*?\)')
    coordinates = {}
    annotations = {}
    with open(gffpath, 'r') as gff:
        for line in gff:
            arr = line.split("\t")
            if len(arr) < 6 or line[0] == "#":
                continue
            if arr[2] == 'mRNA':
                gcoord = "{!s}:{!s}..{!s}".format(arr[0], arr[3], arr[4])
                notes = arr[8].split(';')
                refid = notes[0][notes[0].find(':') + 1:notes[0].rfind('.')]
                annot = '.'
                for elem in notes:
                    if elem.startswith("Note="):
                        annot = elem[5:]
                annotations[refid] = re.sub(rxpr3, '', annot)
                for (hexchar, asciichar) in (
                        ("%3B", ";"), ("%27", "'"), ("%2C", ","),
                        (" contains Interpro domain(s)  ", " ")):
                    annotations[refid] = annotations[refid].replace(
                        hexchar, asciichar)
                coordinates[refid] = gcoord
    return annotations, coordinates


def annotate_mvf(args):
    """Main method"""
    mvf = MultiVariantFile(args.mvf, 'read')
    gff, geneids = parse_gff_annotate(args.gff, mvf.metadata['contigs'])
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


def translate_mvf(args):
    """Main method"""
    mvf = MultiVariantFile(args.mvf, 'read')
    if mvf.metadata['flavor'] != 'dna':
        raise RuntimeError("MVF must be flavor=dna to translate")
    if args.gff:
        gff = parse_gff_translate(args.gff, args)
        if not args.quiet:
            print("gff_processed")
    outmvf = MultiVariantFile(args.out, 'write', overwrite=args.overwrite)
    outmvf.metadata = deepcopy(mvf.metadata)
    outmvf.metadata['flavor'] = args.output_data
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
                for _, amino_acids, alleles in iter_codons(
                        inputbuffer, mvf):
                    if all([x in '-X' for x in amino_acids]):
                        continue
                    if args.output_data == 'protein':
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
            for _, amino_acids, alleles in iter_codons(
                    inputbuffer, mvf):
                if all([x in '-X' for x in amino_acids]):
                    continue
                if args.output_data == 'protein':
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
                        alleles = [MLIB.complement_bases[x]
                                   for x in alleles]
                    decoded_alleles = alleles
                    amino_acids = translate(''.join(alleles))[0]
                else:
                    if reverse_strand:
                        decoded_alleles = [[MLIB.complement_bases[y]
                                            for y in mvf.decode(x)]
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
                if args.output_data == 'protein':
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
