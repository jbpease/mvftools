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

MVF2DUMP: Dumps an MVF to multiple FASTA files by contig
@author: James B. Pease
@author: Ben K. Rosenzweig

version: 2016-04-13 - First testing release
@version 2016-08-02 - Python3 release

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
from mvfbase import MultiVariantFile


def main(arguments=sys.argv[1:]):
    """Main method for mvf2fasta"""
    parser = argparse.ArgumentParser(description="""
    Process MVF into FASTA alignment""")
    parser.add_argument("--mvf", help="input MVF file", required=True)
    parser.add_argument("--outprefix", help="target FASTA file", required=True)
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
        print("Version 2016-04-13")
        sys.exit()
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
    for contig, pos, allelesets in mvf.iterentries(
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
