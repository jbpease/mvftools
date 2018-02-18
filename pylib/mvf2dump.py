#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This program exports the entirety of an MVF to FASTA format,
with many fewer options than mvf2fasta.py.  This is designed
to export large MVF files faster, but with less specific
formatting and region-finding options.
"""

import os
import sys
import argparse
from mvfbase import MultiVariantFile

_LICENSE = """
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


def generate_argparser():
    parser = argparse.ArgumentParser(
        prog="mvf2dump.py",
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=_LICENSE)
    parser.add_argument("-i", "--mvf", type=os.path.abspath,
                        help="Input MVF file.",
                        required=True)
    parser.add_argument("-o", "--outprefix", type=os.path.abspath,
                        help="Target FASTA file", required=True)
    parser.add_argument("-d", "--outdata",
                        choices=("dna", "rna", "prot"),
                        help="output dna, rna or prot data")
    parser.add_argument("-s", "--samples", nargs='*',
                        help="One or more taxon labels, leave blank for all")
    parser.add_argument("-B", "--buffer", type=int, default=10,
                        help="size (Mbp) of write buffer for each sample")
    parser.add_argument("-t", "--tmpdir", default=".",
                        type=os.path.abspath,
                        help="directory to write temporary fasta files")
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
    sample_cols = mvf.get_sample_indices(args.samples or None)
    labels = mvf.get_sample_labels(sample_cols)
    current_contig = ''
    seqs = {}
    for contig, _, allelesets in mvf.iterentries(
            quiet=args.quiet, decode=True):
        if contig != current_contig:
            if seqs:
                with open("{}.{}.fa".format(
                        args.outprefix,
                        mvf.metadata['contigs'][contig]['label']),
                        'wt') as outfile:
                    for seqname in sorted(seqs):
                        outfile.write(">{}\n{}\n".format(
                            seqname, ''.join(seqs[seqname])))
            seqs = None
            seqs = {}
            current_contig = contig[:]
        for col, label in zip(sample_cols, labels):
            if label not in seqs:
                seqs[label] = []
            if flavor in ('dna', 'rna'):
                seqs[label].append(allelesets[0][col] == 'X' and
                                   'N' or allelesets[0][col])
            elif flavor in ('codon', 'prot') and (
                    args.outdata == 'prot'):
                seqs[label].append(allelesets[0][col])
            elif flavor == 'codon' and args.outdata == 'dna':
                seqs[label].extend([allelesets[x][col] == 'X' and 'N' or
                                    allelesets[x][col] for x in (1, 2, 3)])
    if seqs:
        with open("{}.{}.fa".format(
                args.outprefix,
                mvf.metadata['contigs'][contig][
                    'label']), 'wt') as outfile:
            for seqname in sorted(seqs):
                outfile.write(">{}\n{}\n".format(
                    seqname, ''.join(seqs[seqname])))
            seqs = None
            seqs = {}
    return ''


if __name__ == "__main__":
    main()
