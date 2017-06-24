#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This program takes an MVF file and converts the data to a FASTA file
"""

import sys
import argparse
import os
from random import randint
from mvfbase import MultiVariantFile, is_int


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


def parse_regions_arg(regionfilepath, contigs):
    """Parses the regions into coordinates"""
    fmt_regions = []
    region_max_coord = {}
    if regionfilepath is None:
        fmt_regions = [(x, None, None, None) for x in contigs]
        region_max_coord = dict.fromkeys(contigs, None)
    else:
        with open(regionfilepath) as regfile:
            for line in regfile:
                entry = line.rstrip().split(',')
                if len(entry) > 4 or len(entry) < 1 or len(entry[0]) == 0:
                    print("malformed entry ({}), ignoring...".format(entry))
                    continue
                contig = ''
                if entry[0] in contigs:
                    contig = entry[0][:]
                elif is_int(entry[0]):
                    if int(entry[0]) in contigs:
                        contig = int(entry[0])
                if contig == '':
                    for cid in contigs:
                        if contigs[cid]['label'] == entry[0]:
                            contig = cid
                assert contig in contigs
                if len(entry) == 1:
                    fmt_regions.append((entry[0], None, None, '+'))
                elif len(entry) == 2:
                    assert int(entry[1]) > 0
                    fmt_regions.append((entry[0], int(entry[1]), None, '+'))
                else:
                    assert int(entry[1]) > 0
                    assert int(entry[2]) > 0
                    assert int(entry[2]) > int(entry[1])
                    fmt_regions.append((
                        entry[0], int(entry[1]), int(entry[2]),
                        "+" if int(entry[2]) > int(entry[1]) else "-"))
    fmt_regions.sort()
    for contigid, _, maxcoord, _ in fmt_regions:
        if contigid not in region_max_coord:
            region_max_coord[contigid] = maxcoord + 0
        elif maxcoord > region_max_coord[contigid]:
            region_max_coord[contigid] = maxcoord + 0
    regionlabel = ','.join(["{}:{}{}{}{}".format(
        contigs[x[0]]['label'],
        "" if x[1] == -1 else x[1],
        "" if x[2] == -1 else '..',
        "" if x[2] == -1 else x[2],
        "" if x[2] == -1 else "({})".format(x[3])
        ) for x in fmt_regions])
    return fmt_regions, region_max_coord, regionlabel


def generate_argparser():
    parser = argparse.ArgumentParser(
        prog="mvf2fasta.py",
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=_LICENSE)
    parser.add_argument("-i", "--mvf", type=os.path.abspath,
                        help="Input MVF file.", required=True)
    parser.add_argument("-o", "--out", type=os.path.abspath,
                        help="target FASTA file", required=True)
    parser.add_argument("-r", "--regions", nargs='*', required=True,
                        help=("One or more space-separated arguments "
                              "formatted as contigid,start,stop "
                              "coordinates are (inclusive)"))
    parser.add_argument("-l", "--labeltype", choices=('long', 'short'),
                        default='long',
                        help=("Long labels with all metadata or short ids"))
    parser.add_argument("-d", "--outdata", choices=("dna", "rna", "prot"),
                        help="Output dna, rna or prot data.")
    parser.add_argument("-s", "--samples", nargs='*',
                        help="One or more taxon labels, leave blank for all")
    parser.add_argument("-B", "--buffer", type=int, default=10,
                        help="size (Mbp) of write buffer for each sample")
    parser.add_argument("-t", "--tmpdir", default=".",
                        help="directory to write temporary fasta files")
    parser.add_argument("--quiet", action="store_true", default=True,
                        help="suppress screen output")
    parser.add_argument("--version", action="version",
                        version="2017-06-14",
                        help="display version information")
    return parser


def main(arguments=None):
    """Main method"""
    arguments = sys.argv[1:] if arguments is None else arguments
    parser = generate_argparser()
    args = parser.parse_args(args=arguments)
    mvf = MultiVariantFile(args.mvf, 'read')
    flavor = mvf.metadata['flavor']
    if (flavor in ("dna", "rna") and args.outdata == "prot") or (
            flavor == "prot" and args.outdata in ("dna", "rna")):
        raise RuntimeError(
            "--outdata {} incompatiable with '{}' flavor mvf".format(
                args.outdata, flavor))
    regions, max_region_coord, regionlabel = parse_regions_arg(
        args.regions, mvf.metadata['contigs'])
    sample_cols = mvf.get_sample_indices(args.samples or None)
    labels = mvf.get_sample_labels(sample_cols)
    skipcontig = ''
    tmp_files = dict((fn, open("{}-{}.tmp".format(
        fn, randint(1000000, 9999999)), 'w+', args.buffer)) for fn in labels)
    labelwritten = dict.fromkeys(labels, False)
    for contig, pos, allelesets in mvf.iterentries(
            contigs=[x for x in max_region_coord],
            quiet=args.quiet, decode=True):
        if contig == skipcontig:
            continue
        if (contig not in max_region_coord) or (
                max_region_coord[contig] is not None and
                pos > max_region_coord[contig]):
            skipcontig = contig[:]
            continue
        inregion = False
        for rcontig, rstart, rstop, _ in regions:
            if (contig == rcontig):
                if rstart is None or pos >= rstart:
                    if rstop is None or pos <= rstop:
                        inregion = True
                        break
        if inregion is False:
            continue
        for col, label in zip(sample_cols, labels):
            if not labelwritten[label]:
                if args.labeltype == 'long':
                    xlabel = "{} region={}".format(label, regionlabel)
                elif args.labeltype == 'short':
                    xlabel = "{}".format(label)
                tmp_files[label].write(">{}\n".format(xlabel))
                labelwritten[label] = True
            if flavor == 'dna':
                tmp_files[label].write(
                    "N" if allelesets[0][col] == 'X'
                    else allelesets[0][col])
            elif flavor in ('codon', 'prot') and (
                    args.outdata == 'prot'):
                tmp_files[label].write(allelesets[0][col])
            elif flavor == 'codon' and args.outdata == 'dna':
                codon = ["N" if allelesets[x][col] == 'X' else
                         allelesets[x][col] for x in (1, 2, 3)]
                tmp_files[label].write(''.join(codon))
    with open(args.out, 'w') as outfile:
        for filehandler in tmp_files.values():
            filehandler.seek(0, 0)
            buff = filehandler.read(args.buffer)
            while len(buff):
                outfile.write(buff)
                buff = filehandler.read(args.buffer)
            outfile.write("\n")
            filehandler.close()
            os.remove(os.path.join(args.tmpdir, filehandler.name))
    return ''


if __name__ == "__main__":
    main()
