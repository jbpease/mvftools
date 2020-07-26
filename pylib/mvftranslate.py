# -*- coding: utf-8 -*-
"""
This program takes a DNA MVF alignment and annotates the output into
gene boudaries.

MVFtools: Multisample Variant Format Toolkit
James B. Pease and Ben K. Rosenzweig
http://www.github.org/jbpease/mvftools

If you use this software please cite:
Pease, James B. and Benjamin K. Rosenzweig. 2018.
"Encoding Data Using Biological Principles: the Multisample Variant Format
for Phylogenomics and Population Genomics"
IEEE/ACM Transactions on Computational Biology and Bioinformatics. 15(4) 1231–1238.
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
from pylib.mvfbase import MultiVariantFile
from pylib.mvfbiolib import MvfBioLib
MLIB = MvfBioLib()

# RE_GENEID = re.compile("ID=gene:(.*?);")
# PARENTGENE = re.compile("Parent=mRNA:(.*?);")


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
            if codon not in MLIB.codon_tables['full']:
                print(codon)
            amino = MLIB.codon_tables['full'].get(codon, 'X')
            if amino == 'X':
                print(codon, amino)
            aa_seq.append(amino)
    if firststop:
        aa_seq = aa_seq[:aa_seq.index('*') + int(
            firststop == "inclusive")] if '*' in aa_seq else aa_seq
    return aa_seq


def translate_single_codon(seq):
    """Returns translated amino acid from single codon
        Arguments:
    """
    amino = MLIB.codon_tables['full'].get(
        seq.upper().replace('U', 'T'), 'X')
    if amino == 'X' and all(x in 'ATGCU' for x in seq):
        print(amino, seq)
    return amino


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
            amino_acids = [translate(''.join(x))
                           for x in zip(*decoded_alleles)]
            amino_acids = mvf.encode(
                ''.join([x[0] for x in amino_acids]))
        yield inputbuffer[i][0], amino_acids, [
            inputbuffer[i][1][0], inputbuffer[i+1][1][0],
            inputbuffer[i+2][1][0]]


def parse_gff_translate(gff_file, args, parent_gene_prefix='gene:'):
    """Parses a GFF3 file for exon locations
        Arguments:
            gff_file: path to GFF3 file
            args: passthrough from main args
            parent_gene_prefix: prefix for parent gene (def. "gene:")

        Output: triplets for codon locations

    """
    gff_entries = {}
    gff_triplets = {}
    parent_gene_prefix = (
        '' if parent_gene_prefix == 'none' else parent_gene_prefix)
    regex_parentgene = re.compile(
        "Parent={}([^;$]*)".format(parent_gene_prefix))
    #regex_parentgene = re.compile(
    #    'gene_id "([^"]*)'.format(parent_gene_prefix))
    with open(gff_file) as gff:
        for line in gff:
            if line[0] == '#':
                continue
            arr = line.rstrip().split("\t")
            if arr[2] not in ('CDS'):
                continue
            if args.filter_annotation:
                if args.filter_annotation in arr[8]:
                    continue
            if args.require_annotation:
                if args.require_annotation not in arr[8]:
                    continue
            parent = re.findall(regex_parentgene, arr[8])[0]
            if arr[0] not in gff_entries:
                gff_entries[arr[0]] = {}
                gff_triplets[arr[0]] = []
            coords = [int(arr[3]), int(arr[4])]
            strand = arr[6]
            if parent not in gff_entries[arr[0]]:
                gff_entries[arr[0]][parent] = [strand, []]
            gff_entries[arr[0]][parent][1].extend(range(min(coords),
                                                        max(coords) + 1))
    for contiglabel in gff_entries:
        for gene in gff_entries[contiglabel]:
            if len(gff_entries[contiglabel][gene][1]) % 3:
                continue
            strand = gff_entries[contiglabel][gene][0]
            coords = tuple(sorted(gff_entries[contiglabel][gene][1]))
            for j in range(0, len(coords), 3):
                try:
                    gff_triplets[contiglabel].append((coords[j], coords[j+1],
                                                     coords[j+2], strand))
                except IndexError:
                    raise RuntimeError(len(coords), j, strand,
                                       contiglabel, coords[j])
    gff_entries = None
    return gff_triplets


def parse_gff_annotate(gff_file, contigs, filter_annotation=None,
                       gene_prefix='mRNA:'):
    """Parses a GFF3 file for exon locations
        Arguments:
            gff_file: path to GFF3 file
            args: passthrough from main args
            filter_annotation: ignore entries matching these terms
            gene_prefix: add prefix to gene for standard (def. 'mRNA')

        Output: triplets for codon locations

    """
    gff_entries = {}
    relabeled_gff_entries = {}
    geneids = {}
    geneid = 0
    gene_prefix = '' if gene_prefix == 'none' else gene_prefix
    regex_geneid = re.compile("ID={}([^;$]*)".format(gene_prefix))
    #regex_geneid = re.compile('gene_id "([^"]*)')
    with open(gff_file) as gff:
        for line in gff:
            if line[0] == '#':
                continue
            arr = line.rstrip().split("\t")
            if arr[2] != 'gene':
                continue
            if filter_annotation:
                if filter_annotation in arr[8]:
                    continue
            genename = re.findall(regex_geneid, arr[8])
            if not genename:
                raise RuntimeError("Gene IDs not parsing correctly, "
                                   "your current prefix is '{}'. Check "
                                   "your GFF file.".format(gene_prefix))
            genename = genename[0]
            contig = arr[0]
            if contig not in gff_entries:
                gff_entries[contig] = {}
            coords = [int(arr[3]), int(arr[4])]
            strand = arr[6]
            if genename not in geneids:
                geneids[geneid] = {'id': str(geneid),
                                   'label': genename,
                                   'length': max(coords) - min(coords),
                                   'strand': strand}
            for j in range(min(coords), max(coords) + 1):
                gff_entries[contig][j] = geneid
            geneid += 1
    for contig in gff_entries:
        matchlabel = False
        for contigid in contigs:
            if contigs[contigid]['label'] == contig:
                relabeled_gff_entries[str(contigid)] = (
                    gff_entries[contig].copy())
                matchlabel = True
                break
        if matchlabel is False:
            relabeled_gff_entries[str(contig)] = gff_entries[contig].copy()
    gff_entries = None
    return relabeled_gff_entries, geneids


def parse_gff_analysis(gffpath):
    """Parse GFF function used for the analysis script"""
    rxpr3 = re.compile(r' \(AHRD.*?\)')
    coordinates = {}
    annotations = {}
    with open(gffpath, 'r') as gff:
        for line in gff:
            arr = line.split("\t")
            if len(arr) < 6 or line[0] == "#":
                continue
            if arr[2] == 'mRNA':
                gcoord = "{}:{}..{}".format(arr[0], arr[3], arr[4])
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
    args.qprint("Running AnnotateMVF")
    mvf = MultiVariantFile(args.mvf, 'read')
    args.qprint("Input MVF header processed.")
    args.qprint("MVF flavor: {}".format(mvf.flavor))
    gff, geneids = parse_gff_annotate(args.gff, mvf.contig_data,
                                      gene_prefix=args.gene_prefix)
    args.qprint("GFF processed.")
    outmvf = MultiVariantFile(args.out, 'write', overwrite=args.overwrite,
                              flavor=mvf.flavor)
    outmvf.copy_headers_from(mvf)
    if args.nongenic_mode is False:
        outmvf.contig_data = geneids.copy()
        outmvf.contig_indices = list(range(len(geneids)))
        outmvf.contig_ids = [geneids[x]['id'] for x in
                             outmvf.contig_indices]
        outmvf.contig_labels = [geneids[x]['label'] for x in
                                outmvf.contig_indices]
    outmvf.write_data(outmvf.get_header())
    args.qprint("Output MVF established.")
    entrybuffer = []
    nentry = 0
    args.qprint("Processing MVF entries.")
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
        if annotated_pos and not args.nongenic_mode:
            entrybuffer.append((gff[contigid][pos], pos, allelesets))
        elif args.nongenic_mode and not annotated_pos:
            entrybuffer.append((contigid, pos, allelesets))
        if args.nongenic_mode or annotated_pos:
            nentry += 1
            if nentry == args.line_buffer:
                args.qprint("Writing block of entries.")
                outmvf.write_entries(entrybuffer)
                entrybuffer = []
                nentry = 0
    if entrybuffer:
        outmvf.write_entries(entrybuffer)
        args.qprint("Writing final block of entries.")
        entrybuffer = []
        nentry = 0
    return ''


def translate_mvf(args):
    """Main method"""
    args.qprint("Running TranslateMVF")
    if args.gff:
        args.qprint("Reading and Indexing MVF.")
    else:
        args.qprint("Reading MVF.")
    mvf = MultiVariantFile(args.mvf, 'read', contigindex=bool(args.gff))
    if mvf.flavor != 'dna':
        raise RuntimeError("MVF must be flavor=dna to translate")
    if args.gff:
        args.qprint("Processing MVF Index File.")
        mvf.read_index_file()
        args.qprint("GFF processing start.")
        gff = parse_gff_translate(args.gff, args,
                                  parent_gene_prefix=args.parent_gene_prefix)
        args.qprint("GFF processed.")
    outmvf = MultiVariantFile(args.out, 'write', overwrite=args.overwrite)
    outmvf.copy_headers_from(mvf)
    outmvf.flavor = args.output_data
    outmvf.write_data(outmvf.get_header())
    args.qprint("Output MVF Established.")
    entrybuffer = []
    nentry = 0
    pos = None
    if not args.gff:
        args.qprint("No GFF used, translating sequences as pre-aligned in "
                    "coding frame.")
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
                    if nentry == args.line_buffer:
                        outmvf.write_entries(entrybuffer)
                        entrybuffer = []
                        nentry = 0
                inputbuffer = [(pos, allelesets)]
                current_contig = contigid[:]
        if inputbuffer:
            for _, amino_acids, alleles in iter_codons(
                    inputbuffer, outmvf):
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
                if nentry == args.line_buffer:
                    outmvf.write_entries(entrybuffer)
                    entrybuffer = []
                    nentry = 0
    else:
        args.qprint("Indexing GFF gene names.")
        # mvfid_to_gffname = outmvf.get_contig_reverse_dict()
        for xcontig in outmvf.get_contig_indices():
            mvf_entries = {}
            xcontiglabel = outmvf.get_contig_labels(indices=xcontig)[0]
            xcontigid = outmvf.get_contig_ids(indices=xcontig)[0]
            if xcontiglabel not in gff:
                if args.verbose:
                    print("No entries in GFF, skipping contig: index:{} id:{} label:{}".format(
                        xcontig, xcontigid, xcontiglabel))
                continue
            if not xcontig % 100:
                args.qprint("Processing contig: {} {}".format(
                    xcontigid, xcontiglabel))
            for contigid, pos, allelesets in mvf.itercontigentries(
                    xcontig, decode=False):
                mvf_entries[pos] = allelesets[0]
            for coords in sorted(gff[xcontiglabel]):
                reverse_strand = coords[3] == '-'
                alleles = (tuple(mvf_entries.get(x, '-')
                                 for x in coords[2::-1]) if
                           reverse_strand is True else tuple(mvf_entries.get(x, '-')
                                                     for x in coords[0:3]))
                if xcontig == 0 and 166252 in coords:
                    print(coords)
                    print(alleles)
                if all(len(x) == 1 for x in alleles):
                    if reverse_strand:
                        alleles = tuple(
                            MLIB.complement_bases[x] for x in alleles)
                    decoded_alleles = alleles
                    amino_acids = translate_single_codon(''.join(alleles))
                else:
                    if reverse_strand is True:
                        decoded_alleles = tuple(tuple(MLIB.complement_bases[y]
                                                      for y in mvf.decode(x))
                                                for x in alleles)
                        alleles = tuple(outmvf.encode(''.join(x))
                                        for x in decoded_alleles)
                    else:
                        decoded_alleles = tuple(mvf.decode(x) for x in alleles)
                    amino_acids = tuple(translate_single_codon(''.join(x))
                                        for x in zip(*decoded_alleles))
                    # print("aminx", amino_acids)
                    amino_acids = outmvf.encode(''.join(amino_acids))
                # if all(x in '-X' for x in amino_acids):
                #    continue
                # print("amino", amino_acids)
                # print("translated", amino_acids, alleles)
                if args.output_data == 'protein':
                    entrybuffer.append((xcontig, coords[0], (amino_acids,)))
                else:
                    entrybuffer.append((
                        xcontigid, coords[0], (
                            amino_acids, alleles[0], alleles[1], alleles[2])))
                nentry += 1
                if nentry >= args.line_buffer:
                    args.qprint("Writing a block of {} entries.".format(
                        args.line_buffer))
                    outmvf.write_entries(entrybuffer)
                    entrybuffer = []
                    nentry = 0
    if entrybuffer:
        outmvf.write_entries(entrybuffer)
        entrybuffer = []
        nentry = 0
    return ''
