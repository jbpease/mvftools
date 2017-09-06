#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This program is used to convert a FASTA file into MVF format.
"""

import sys
import re
import os
import argparse
from mvfbase import encode_mvfstring, MultiVariantFile, fasta_iter


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


RE_CONTIG_NAME = re.compile("ID=(.*?),")
RE_CONTIG_LENGTH = re.compile("length=(.*?)>")


def generate_argparser():
    parser = argparse.ArgumentParser(
        prog="fasta2mvf.py",
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=_LICENSE)
    parser.add_argument("-i", "--fasta", nargs='*', required=True,
                        help="input FASTA file(s)")
    parser.add_argument("-o", "--out", help="output MVF file", required=True)
    parser.add_argument("-f", "--flavor", choices=['dna', 'protein'],
                        help="type of file [dna] or protein", default='dna')
    parser.add_argument("-c", "--contig-ids", "--contigids", nargs='*',
                        help=("manually specify one or more contig ids "
                              "as ID:NAME"))
    parser.add_argument("-s", "--sample-replace", "--samplereplace", nargs="*",
                        help=("one or more TAG:NEWLABEL or TAG, items, "
                              "if TAG found in sample label, replace with "
                              "NEW (or TAG if NEW not specified) "
                              "NEW and TAG must each be unique"))
    parser.add_argument("-R", "--ref-label", "--reflabel", default="REF",
                        help="label for reference sample")
    parser.add_argument("-B", "--read-buffer", "--readbuffer",
                        type=int, default=100000,
                        help="number of lines to hold in READ buffer")
    parser.add_argument("-W", "--write-buffer", "--writebuffer",
                        type=int, default=100000,
                        help="number of lines to hold in WRITE buffer")
    parser.add_argument("-F", "--field-sep", "--fieldsep", nargs='*',
                        default=None,
                        choices=['TAB', 'SPACE', 'DBLSPACE',
                                 'COMMA', 'MIXED', 'PIPE', 'AT',
                                 'UNDER', 'DBLUNDER'],
                        help=("FASTA field separator; assumes "
                              "'>database accession locus' "
                              "format"))
    parser.add_argument("--contig-field", "--contigfield", type=int,
                        help=("When headers are split by --field-sep, "
                              "the 0-based index of the contig id."))
    parser.add_argument("--contig-by-file", "--contigbyfile",
                        action="store_true",
                        help=("Contigs are designated by separate files."))
    parser.add_argument("--sample-field", "--samplefield", type=int,
                        help=("when headers are split by --field-sep, "
                              "the 0-based index of the sample id"))
    parser.add_argument("--manual-coord", "--manualcoord", nargs='*',
                        help=("manually specify reference coordinates "
                              "for each file in the format "
                              "CONTIGID:START..STOP, ..."))
    parser.add_argument("--overwrite", action="store_true",
                        help="USE WITH CAUTION: force overwrite of outputs")
    parser.add_argument("--version", action="version",
                        version="2017-09-05",
                        help="display version information")
    return parser


def main(arguments=None):
    """Main method"""
    arguments = sys.argv[1:] if arguments is None else arguments
    parser = generate_argparser()
    args = parser.parse_args(args=arguments)
    sepchars = dict([("PIPE", "\\|"), ("TAB", "\\t"),
                     ("SPACE", "\\s"), ("DBLSPACE", "\\s\\s"),
                     ("COMMA", "\\,"), ("NONE", None),
                     ("AT", "\\@"), ('UNDER', "\\_"), ("DBLUNDER", "\\_\\_")])
    if args.field_sep is None:
        args.field_sep = ''
    else:
        args.field_sep = re.compile("[{}]".format(
            ''.join([sepchars[x] for x in args.field_sep])))
    if args.manual_coord:
        assert len(args.manual_coord) == len(args.fasta)
        args.manual_coord = [
            (x.split(':')[0], int(x.split(":")[1].split('..')[0]),
                int(x.split(':')[1].split('..')[1]))
            for x in args.manual_coord]
    mvf = MultiVariantFile(args.out, 'write', overwrite=args.overwrite)
    fasta = {}
    current_contig = 0
    fsamples = []
    fcontigs = []
    for ifasta, fastapath in enumerate(args.fasta):
        print("Processing {}".format(fastapath))
        for header, seq in fasta_iter(fastapath):
            if args.field_sep is None:
                header = header[:]
            if args.field_sep != '' and args.field_sep is not None:
                header = [str(x) for x in re.split(args.field_sep, header)]
            if args.contig_by_file is True:
                contig = os.path.basename(fastapath[:])
                if args.sample_field is None:
                    sample = header[:]
                else:
                    sample = header[args.sample_field]
            elif (len(header) < max(args.contig_field if
                                    args.contig_field
                                    is not None else 0,
                                    args.sample_field if
                                    args.sample_field is not None else
                                    0) or
                  args.contig_field is None or args.sample_field is None):
                contig = "UNK{}".format(current_contig)
                sample = header[:]
            elif args.manual_coord:
                contig = args.manual_coord[ifasta][0]
            else:
                contig = header[args.contig_field]
                sample = header[args.sample_field]
            if contig not in fcontigs:
                fcontigs.append(contig)
                fasta[contig] = {}
            if sample not in fsamples:
                fsamples.append(sample)
            fasta[contig][sample] = (len(seq), seq)
    reflabel = None
    if args.ref_label:
        for i, samplename in enumerate(fsamples):
            if args.ref_label in samplename:
                reflabel = i
                break
    if reflabel:
        newref = fsamples.pop(i)
        fsamples = [newref] + fsamples
    for i, contig in enumerate(fcontigs):
        mvf.metadata['contigs'][i] = {
            'label': contig,
            'length': max([fasta[contig][x][0] for x in fasta[contig]])}
    mvf.metadata['labels'] = fsamples[:]
    for i, label in enumerate(fsamples[:]):
        mvf.metadata['samples'][i] = {'label': label}
    mvf.metadata['ncol'] = len(mvf.metadata['labels'])
    mvf.metadata['sourceformat'] = 'fasta'
    mvf.metadata['flavor'] = args.flavor
    # WRITE MVF HEADER
    mvf.write_data(mvf.get_header())
    mvfentries = []
    nentry = 0
    mvf_alleles = {}
    for cind, contig in enumerate(fcontigs):
        for pos in range(mvf.metadata['contigs'][cind]['length']):
            mvf_alleles = encode_mvfstring(
                ''.join(samp not in fasta[contig] and '-' or
                        pos >= fasta[contig][samp][0] and '-' or
                        fasta[contig][samp][1][pos]
                        for samp in fsamples))
            if mvf_alleles:
                if args.flavor == 'dna':
                    mvf_alleles = ''.join(["X" if x in 'NXOBDHVnxobdhv' else x
                                           for x in mvf_alleles])
                mvfentries.append(
                    (cind, pos+1, (mvf_alleles,)))
                nentry += 1
                if nentry == args.write_buffer:
                    mvf.write_entries(mvfentries, encoded=True)
                    mvfentries = []
                    nentry = 0
    if mvfentries:
        mvf.write_entries(mvfentries)
        mvfentries = []
    return ''


if __name__ == "__main__":
    main()
