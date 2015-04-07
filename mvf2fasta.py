#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
MVFtools: Multisample Variant Format Toolkit
http://www.github.org/jbpease/mvftools

MVF2FASTA: Convert MVF file to a FASTA file
@author: James B. Pease
@author: Ben K. Rosenzweig

Version: 2015-02-01 - First Public Release

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
import sys, argparse, os
from mvfbase import  MultiVariantFile

def main(arguments=sys.argv[1:]):
    """Main method for mvf2fasta"""
    parser = argparse.ArgumentParser(description="""
    Process MVF into FASTA alignment""")
    parser.add_argument("--mvf", help="input MVF file", required=True)
    parser.add_argument("--out", help="target FASTA file", required=True)
    parser.add_argument("--labeltype", choices=['long', 'short'],
                        default='long',
                        help="long labels with all metadata or short ids")
    parser.add_argument("--regions", nargs='*',
                        help="one or more regions id,start,stop (inclusive)")
    parser.add_argument("--samples", nargs='*',
                        help="one or more taxon labels, leave blank for all")
    parser.add_argument("--outgroups", nargs="*")
    parser.add_argument("--contigs", nargs='*',
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
        print("Version 2015-02-01: Initial Public Release")
        sys.exit()
    mvf = MultiVariantFile(args.mvf, 'read')
    if args.contigs:
        contigs = dict(mvf.metadata['contigs'][c] for c in args.contigs)
    else:
        contigs = dict(mvf.metadata['contigs'])
    sample_cols = mvf.get_sample_indices(args.samples or None)
    labels = mvf.get_sample_labels(sample_cols)
    current_contig = None
    tmp_files = dict((fn, open(fn+'.tmp', 'w+', args.buffer)) for fn in labels)
    for contig, _, allelesets in mvf.iterentries(
            contigs=args.contigs, subset=sample_cols,
            quiet=args.quiet, decode=True):
        alleles = mvf.decode(allelesets)
        if current_contig != contig:
            current_contig = contig
            for col, label in zip(sample_cols, labels):
                if args.labeltype == 'long':
                    tmp_files[label].write(
                        '\n>{} contig={}  length={}\n{}'.format(
                            label,
                            contigs[current_contig]['label'],
                            contigs[current_contig]['length'],
                            alleles[col]))
                elif args.labeltype == 'short':
                    tmp_files[label].write(
                        '\n>{}_{}\n{}'.format(
                            label, contigs[current_contig]['label'],
                            alleles[col]))
        else:
            for col, label in zip(sample_cols, labels):
                tmp_files[label].write(alleles[col])
    with open(args.out, 'w') as outfile:
        for filehandler in tmp_files.values():
            filehandler.seek(0, 0)
            buff = filehandler.read(args.buffer)
            while  len(buff):
                outfile.write(buff)
                buff = filehandler.read(args.buffer)
            filehandler.close()
            os.remove(os.path.join(args.tmpdir, filehandler.name))
    return ''


if __name__ == "__main__":
    main()
