#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This program is used to export MVF data to Phylip format.
"""

import sys
import argparse
import os
from random import randint
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


def parse_regions_arg(regionfilepath, contigs):
    """Parses the regions into coordinates"""
    fmt_regions = {}
    region_max_coord = {}
    with open(regionfilepath) as regfile:
        for line in regfile:
            entry = line.rstrip().split(',')
            if len(entry) > 4 or len(entry) < 1 or len(entry[0]) == 0:
                print("malformed entry ({}), ignoring...".format(entry))
                continue
            contig = ''
            if entry[0] in contigs:
                contig = entry[0][:]
            if contig == '':
                for cid in contigs:
                    if contigs[cid]['label'] == entry[0]:
                        contig = cid
            assert contig in contigs
            if len(entry) == 1:
                regionentry = (None, None, '+')
            elif len(entry) == 2:
                assert int(entry[1]) > 0
                regionentry = (int(entry[1]), None, '+')
            else:
                assert int(entry[1]) > 0
                assert int(entry[2]) > 0
                assert int(entry[2]) > int(entry[1])
                regionentry = (
                    int(entry[1]), int(entry[2]),
                    "+" if int(entry[2]) > int(entry[1]) else "-")
            if contig not in fmt_regions:
                fmt_regions[contig] = []
                region_max_coord[contig] = None
            fmt_regions[contig].append(regionentry)
            if regionentry[1] is not None:
                if region_max_coord[contig] is None:
                    region_max_coord[contig] = regionentry[1] + 0
                elif region_max_coord[contig] < regionentry[1]:
                    region_max_coord[contig] = regionentry[1] + 0
    return fmt_regions, region_max_coord  # , regionlabel


def generate_argparser():
    parser = argparse.ArgumentParser(
        prog="mvf2phy.py",
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=_LICENSE)
    parser.add_argument("-i", "--mvf", type=os.path.abspath,
                        help="Input MVF file.", required=True)
    parser.add_argument("-o", "--out", type=os.path.abspath,
                        help="Output Phylip file.", required=True)
    parser.add_argument("-r", "--region",
                        type=os.path.abspath,
                        help=("Path of a plain text file "
                              "containing one more lines "
                              "with entries 'contigid,stop,start' "
                              "(one per line, inclusive coordinates) "
                              "all data will be returned if left blank."))
    parser.add_argument("-L", "--labeltype", choices=('long', 'short'),
                        default='short',
                        help="Long labels with all metadata or short ids")
    parser.add_argument("-d", "--outdata", choices=("dna", "rna", "prot"),
                        help="Output dna, rna or prot data.")
    parser.add_argument("-s", "--samples", nargs='*',
                        help="One or more taxon labels, leave blank for all.")
    parser.add_argument("-B", "--buffer", type=int, default=100000,
                        help="size (bp) of write buffer for each sample")
    parser.add_argument("-t", "--tmpdir", default=".",
                        help="directory to write temporary fasta files")
    parser.add_argument("-p", "--partition", action="store_true",
                        help=("Output a CSV partitions file with RAxML"
                              "formatting for use in partitioned "
                              "phylogenetic methods."))
    parser.add_argument("--quiet", action="store_true", default=True,
                        help="suppress screen output")
    parser.add_argument("--version", action="version",
                        version="2017-06-24",
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
    regions = None
    max_region_coord = dict.fromkeys(mvf.metadata['contigs'], None)
    if args.region is not None:
        regions, max_region_coord = parse_regions_arg(
            args.region, mvf.metadata['contigs'])
    sample_cols = mvf.get_sample_indices(args.samples or None)
    labels = mvf.get_sample_labels(sample_cols)
    skipcontig = ''
    tmp_files = dict((fn, open("{}-{}.tmp".format(
        fn, randint(1000000, 9999999)), 'w+', args.buffer)) for fn in labels)
    labelwritten = dict.fromkeys(labels, False)
    curcontigname = None
    curcontigstart = 1
    curcontigend = 1
    if args.partition is True:
        partprefix = "PROT" if args.outdata == "prot" else "DNA"
        partitionfile = open("{}.part".format(args.out), 'w')
    for contig, pos, allelesets in mvf.iterentries(
            contigs=(mvf.metadata['contigs'] if args.region is None else
                     [x for x in max_region_coord]),
            quiet=args.quiet, decode=True):
        if contig == skipcontig:
            continue
        if contig not in max_region_coord:
            skipcontig = contig[:]
            continue
        if curcontigname is None:
            curcontigname = contig[:]
        elif contig != curcontigname:
            if args.partition is True:
                if curcontigend > curcontigstart:
                    partitionfile.write("{}, {} = {}-{}\n".format(
                        partprefix, mvf.get_contig_label(curcontigname),
                        curcontigstart, curcontigend - 1))
            curcontigname = contig[:]
            # reset start as one position after end of last
            curcontigstart = curcontigend
            curcontigend = curcontigend + 1
        for col, label in zip(sample_cols, labels):
            if not labelwritten[label]:
                if args.labeltype == 'long':
                    tmp_files[label].write("{}{}".format(
                        label[:100], " "*(100 - len(label[:100]))))
                elif args.labeltype == 'short':
                    tmp_files[label].write("{}{}".format(
                        label[:20], " "*(20 - len(label[:20]))))
                labelwritten[label] = True
            if flavor == 'dna':
                tmp_files[label].write(
                    allelesets[0][col] == 'X' and
                    'N' or allelesets[0][col])
                if label == labels[0]:
                    curcontigend += 1
            elif ((flavor == 'codon' and args.outdata == 'prot') or (
                    flavor == 'prot')):
                tmp_files[label].write(allelesets[0][col])
                if label == labels[0]:
                    curcontigend += 1
            elif flavor == 'codon':
                codon = ["N" if allelesets[x][col] == 'X' else
                         allelesets[x][col] for x in (1, 2, 3)]
                tmp_files[label].write(''.join(codon))
                if label == labels[0]:
                    curcontigend += 3
    first_file = True
    totalseqlen = 0
    with open(args.out, 'w') as outfile:
        for filehandler in tmp_files.values():
            # read first file to establish sequence length for phylip header
            if first_file is True:
                filehandler.seek(0, 0)
                buff = filehandler.read(args.buffer)
                while buff != '':
                    if " " in buff:
                        totalseqlen += len(buff.strip().split(" ")[-1])
                    else:
                        totalseqlen += len(buff.strip())
                    buff = filehandler.read(args.buffer)
                outfile.write("{} {}\n".format(len(labels), totalseqlen))
                first_file = False
            filehandler.seek(0, 0)
            buff = filehandler.read(args.buffer)
            while buff != '':
                if first_file is True:
                    outfile.write("{} {}\n".format(
                        len(labels), len(buff.split()[1])))
                    first_file = False
                outfile.write(buff)
                buff = filehandler.read(args.buffer)
            outfile.write("\n")
            filehandler.close()
            os.remove(os.path.join(args.tmpdir, filehandler.name))
    if args.partition is True:
        if curcontigend > curcontigstart:
            partitionfile.write("{},{},{},{}\n".format(
                partprefix, mvf.get_contig_label(curcontigname),
                curcontigstart, curcontigend - 1))
        partitionfile.close()
    return ''


if __name__ == "__main__":
    main()
