#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MVFtools: Multisample Variant Format Toolkit
http://www.github.org/jbpease/mvftools

If you use this software please cite:
Pease JB and BK Rosenzweig. 2016.
"Encoding Data Using Biological Principles: the Multisample Variant Format
for Phylogenomics and Population Genomics"
IEEE/ACM Transactions on Computational Biology and Bioinformatics. In press.
http://www.dx.doi.org/10.1109/tcbb.2015.2509997

MVF2FASTA: Convert MVF file to a FASTA file
@author: James B. Pease
@author: Ben K. Rosenzweig

version: 2015-02-01 - First Public Release
version: 2015-09-04 - Cleanup
version: 2015-12-31 - New headers and cleanup
version: 2016-03-15 - Major refurbish, changes to args
version: 2016-08-02 - Python3 conversion
@version: 2016-10-25 - Minor fixes to regions lookup
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

import sys
import argparse
import os
from random import randint
from mvfbase import MultiVariantFile, is_int


def parse_regions_arg(regions, contigs):
    """Parses the regions argument into coordinates"""
    fmt_regions = []
    region_max_coord = {}
    for elem in regions:
        row = elem.split(',')
        assert len(row) > 0 and len(row) < 4
        contig = ''
        if row[0] in contigs:
            contig = row[0][:]
        elif is_int(row[0]):
            if int(row[0]) in contigs:
                contig = int(row[0])
        if contig == '':
            for cid in contigs:
                if contigs[cid]['label'] == row[0]:
                    contig = cid
        assert contig in contigs
        if len(row) == 1:
            fmt_regions.append((row[0], -1, -1, '+'))
        elif len(row) == 2:
            assert int(row[1]) > 0
            fmt_regions.append((row[0], int(row[1]), -1, '+'))
        else:
            assert int(row[1]) > 0
            assert int(row[2]) > 0
            if int(row[2]) > int(row[1]):
                fmt_regions.append((row[0], int(row[1]), int(row[2]), '+'))
            else:
                fmt_regions.append((row[0], int(row[1]), int(row[2]), '-'))
    regions.sort()
    for contigid, _, maxcoord, _ in fmt_regions:
        if contigid not in region_max_coord:
            region_max_coord[contigid] = maxcoord + 0
        elif maxcoord > region_max_coord[contigid]:
            region_max_coord[contigid] = maxcoord + 0
    regionlabel = ','.join(["{}:{}{}{}{}".format(
        contigs[x[0]]['label'], 
        x[1] != -1 and x[1] or '', 
        x[2] != -1 and '..' or '',
        x[2] != -1 and x[2] or '', 
        x[2] != -1 and "({})".format(x[3]) or ''
        ) for x in fmt_regions])
    return fmt_regions, region_max_coord, regionlabel


def main(arguments=sys.argv[1:]):
    """Main method for mvf2fasta"""
    parser = argparse.ArgumentParser(description="""
    Process MVF into FASTA alignment""")
    parser.add_argument("--mvf", help="input MVF file", required=True)
    parser.add_argument("--out", help="target FASTA file", required=True)
    parser.add_argument("--regions", nargs='*', required=True,
                        help="""one or more regions contigid,
                                start,stop (inclusive)""")
    parser.add_argument("--labeltype", choices=['long', 'short'],
                        default='long',
                        help="long labels with all metadata or short ids")
    parser.add_argument("--outdata", choices=("dna", "rna", "prot"),
                        help="output dna, rna or prot data")
    parser.add_argument("--samples", nargs='*',
                        help="one or more taxon labels, leave blank for all")
    parser.add_argument("--buffer", type=int, default=10,
                        help="size (Mbp) of write buffer for each sample")
    parser.add_argument("--tmpdir", default=".",
                        help="directory to write temporary fasta files")
    parser.add_argument("--quiet", action="store_true", default=True,
                        help="suppress screen output")
    parser.add_argument("-v", "--version", action="store_true",
                        help="display version information")
    args = parser.parse_args(args=arguments)
    if args.version:
        print("Version 2016-10-25")
        sys.exit()
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
        if contig not in max_region_coord or (
                pos > max_region_coord[contig] and 
                max_region_coord[contig] != -1):
            skipcontig = contig[:]
            continue
        inregion = False
        for (rcontig, rstart, rstop, _) in regions:
            if contig == rcontig and (
                    rstart <= pos <= rstop or (
                        pos >= rstart and rstop == -1)):
                inregion = True
                break
        if not inregion:
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
                    allelesets[0][col] == 'X' and
                    'N' or allelesets[0][col])
            elif flavor in ('codon', 'prot') and (
                    args.outdata == 'prot'):
                tmp_files[label].write(allelesets[0][col])
            elif flavor == 'codon' and args.outdata == 'dna':
                codon = [allelesets[x][col] == 'X' and 'N' or
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
