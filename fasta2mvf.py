#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
MVFtools: Multisample Variant Format Toolkit
http://www.github.org/jbpease/mvftools

FASTA2MVF: Variant Call Format (FASTA) to MVF conversion program
@author: James B. Pease
@author: Ben K. Rosenzweig

version 2015-09-05 - First release
@version: 2015-12-16 - Updates to cmd line args additional fixes

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
import re
import os
import argparse
from mvfbase import encode_mvfstring, MultiVariantFile, fasta_iter

RE_CONTIG_NAME = re.compile("ID=(.*?),")
RE_CONTIG_LENGTH = re.compile("length=(.*?)>")


def main(arguments=sys.argv[1:]):
    """Main method for fasta2mvf"""
    parser = argparse.ArgumentParser(description="""
    Converts multisample-FASTA to MVF file with filtering """)
    parser.add_argument("--fasta", nargs='*', required=True,
                        help="input FASTA file(s)")
    parser.add_argument("--out", help="output MVF file", required=True)
    parser.add_argument("--contigids", nargs='*',
                        help=("""manually specify one or more contig ids
                                 as ID:NAME"""))
    parser.add_argument("--samplereplace", nargs="*",
                        help="""one or more TAG:NEWLABEL or TAG, items,
                                if TAG found in sample label, replace with
                                NEW (or TAG if NEW not specified)
                                NEW and TAG must each be unique""")
    parser.add_argument("--reflabel", default="REF",
                        help="label for reference sample (default='REF')")
    parser.add_argument("--allelesfrom", default=None,
                        help="""get additional alignment columns
                                from INFO fields (:-separated)""")
    parser.add_argument("--readbuffer", type=int, default=100000,
                        help="number of lines to hold in READ buffer")
    parser.add_argument("--writebuffer", type=int, default=100000,
                        help="number of lines to hold in WRITE buffer")
    parser.add_argument("--fieldsep", nargs='*', default=None,
                        choices=['TAB', 'SPACE', 'DBLSPACE',
                                 'COMMA', 'MIXED', 'PIPE', 'AT',
                                 'UNDER', 'DBLUNDER'],
                        help="""FASTA field separator; assumes
                                '>database/SEP/accession/SEP/locus'
                                format (default='NONE')""")
    parser.add_argument("--contigfield", type=int,
                        help="""when headers are split by --fieldsep,
                        the 0-based index of the contig id""")
    parser.add_argument("--contigbyfile", action="store_true",
                        help="""Contigs are designated by separate files""")
    parser.add_argument("--samplefield", type=int,
                        help="""when headers are split by --fieldsep,
                        the 0-based index of the sample id""")
    parser.add_argument("--manualcoord", nargs='*',
                        help="""manually specify reference coordinates for each file
                                in the format CONTIGID:START..STOP, ...""")
    parser.add_argument("--overwrite", action="store_true",
                        help="USE WITH CAUTION: force overwrite of outputs")
    parser.add_argument("-v", "--version", action="store_true",
                        help="display version information")
    args = parser.parse_args(args=arguments)
    if args.version:
        print("Version 2015-09-05")
        sys.exit()
    sepchars = dict([("PIPE", "\\|"), ("TAB", "\\t"),
                     ("SPACE", "\\s"), ("DBLSPACE", "\\s\\s"),
                     ("COMMA", "\\,"), ("NONE", None),
                     ("AT", "\\@"), ('UNDER', "\\_"), ("DBLUNDER", "\\_\\_")])
    if args.fieldsep is None:
        args.fieldsep = ''
    else:
        args.fieldsep = re.compile("[{}]".format(
            ''.join([sepchars[x] for x in args.fieldsep])))
    if args.manualcoord:
        assert len(args.manualcoord) == len(args.fasta)
        args.manualcoord = [
            (x.split(':')[0], int(x.split(":")[1].split('..')[0]),
                int(x.split(':')[1].split('..')[1]))
            for x in args.manualcoord]
    mvf = MultiVariantFile(args.out, 'write', overwrite=args.overwrite)
    fasta = {}
    current_contig = 0
    fsamples = []
    fcontigs = []
    for ifasta, fastapath in enumerate(args.fasta):
        print("Processing {}".format(fastapath))
        for header, seq in fasta_iter(fastapath):
            header = [str(x) for x in re.split(args.fieldsep, header)]
            if args.contigbyfile:
                contig = os.path.basename(fastapath[:])
                sample = header[args.samplefield or 0]
            elif (len(header) < max(3, args.contigfield or 0,
                                    args.samplefield or 0) or
                  args.contigfield is None or args.samplefield is None):
                contig = "UNK{}".format(current_contig)
                sample = header[0]
            elif args.manualcoord:
                contig = args.manualcoord[ifasta][0]
            else:
                contig = header[args.contigfield]
                sample = header[args.samplefield]
            if contig not in fcontigs:
                fcontigs.append(contig)
                fasta[contig] = {}
            if sample not in fsamples:
                fsamples.append(sample)
            fasta[contig][sample] = (len(seq), seq)
    reflabel = None
    if args.reflabel:
        for i, samplename in enumerate(fsamples):
            if args.reflabel in samplename:
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
    # WRITE MVF HEADER
    mvf.write_data(mvf.get_header())
    mvfentries = []
    nentry = 0
    mvf_alleles = {}
    for cind, contig in enumerate(fcontigs):
        for pos in range(mvf.metadata['contigs'][cind]['length']):
            mvf_alleles = encode_mvfstring(
                ''.join(samp not in fasta[contig] and '-' or
                        pos > fasta[contig][samp][0] and '-' or
                        fasta[contig][samp][1][pos]
                        for samp in fsamples))
            if mvf_alleles:
                mvfentries.append(
                    (cind, pos+1, (mvf_alleles,)))
                nentry += 1
                if nentry == args.writebuffer:
                    mvf.write_entries(mvfentries, encoded=True)
                    # print(mvfentries[:5])
                    mvfentries = []
                    nentry = 0
    if mvfentries:
        mvf.write_entries(mvfentries)
        mvfentries = []
    return ''


if __name__ == "__main__":
    main()
