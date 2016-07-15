#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
MVFtools: Multisample Variant Format Toolkit
http://www.github.org/jbpease/mvftools (Stable Releases)
http://www.github.org/jbpease/mvftools-dev (Latest Testing Updates)

If you use this software please cite:
Pease JB and BK Rosenzweig. 2016.
"Encoding Data Using Biological Principles: the Multisample Variant Format
for Phylogenomics and Population Genomics"
IEEE/ACM Transactions on Computational Biology and Bioinformatics. In press.
http://www.dx.doi.org/10.1109/tcbb.2015.2509997
http://www.github.org/jbpease/mvftools

MVF_check: Verify Compliance and Check for Errors in MVF Format
@author: James B. Pease

@version: 2015-12-31 - Updates to headers and cleanup

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
import argparse
from mvfbase import MultiVariantFile


def main(arguments=sys.argv[1:]):
    """Main method for mvf_annotate"""
    parser = argparse.ArgumentParser(description="""
    Reannotates MVF based on genes""")
    parser.add_argument("mvf", help="input MVF file")
    parser.add_argument("--quiet", action="store_true",
                        help="suppress progress meter")
    parser.add_argument("-v", "--version", action="store_true",
                        help="display version information")
    args = parser.parse_args(args=arguments)
    if args.version:
        print("Version 2015-12-31")
        sys.exit()
    mvf = MultiVariantFile(args.mvf, 'read')
    contigs = mvf.metadata['contigs']
    ncol = mvf.metadata['ncol']
    previous_location = (None, None)
    if mvf.metadata['mvftype'] in ('dna', 'protein'):
        if mvf.metadata['mvftype'] == 'protein':
            valid_bases = 'ACDEFGHIKLMNPQRSTVWY'
            valid_characters = 'ACDEFGHIKLMNPQRSTVWYX-'
        else:
            valid_bases = 'ATGCKMRSWY'
            valid_characters = 'ATGCKMRSWYX-'
        for contigid, pos, allelesets in mvf:
            alleles = allelesets[0]
            nonref = False
            if alleles[0] == '@':
                alleles = alleles[1:]
                nonref = True
            errmsg = []
            #  CHECK ALLELES
            if len(alleles) == 1:
                if alleles in 'X-':
                    errmsg.append("no data")
                elif alleles not in valid_bases:
                    errmsg.append("invalid alleles")
            elif len(alleles) == 2:
                if alleles[0] == alleles[1]:
                    errmsg.append("invalid format")
                elif alleles[0] in '-' and not nonref:
                    errmsg.append("empty reference")
                elif (alleles[0] not in valid_characters or
                      alleles[1] not in valid_characters):
                    errmsg.append("invalid alleles")
            elif alleles[1] == '+':
                if alleles[2] == '-':
                    errmsg.append("invalid format")
                elif alleles[0] == '-' and not nonref:
                    errmsg.append("empty reference")
                elif (alleles[0] not in valid_characters or
                      alleles[2] not in valid_characters):
                    errmsg.append("invalid alleles")
                elif int(alleles[3:]) > ncol:
                    errmsg.append("invalid sample number")
            elif alleles[2] == '+':
                if alleles[0] == alleles[1] and alleles[0] == alleles[3]:
                    errmsg.append("invalid format")
                elif any(alleles[x] not in valid_characters
                         for x in (0, 1, 3)):
                    errmsg.append("invalid alleles")
                elif int(alleles[4:]) > ncol:
                    errmsg.append("invalid sample number")
            else:
                if alleles[0] in '-' and not nonref:
                    errmsg.append("empty reference")
                if alleles[0] in '-':
                    errmsg.append("empty reference")
                if any(x not in valid_characters for x in alleles):
                    errmsg.append("invalid alleles")
            #  CHECK POSITION
            if contigid not in contigs:
                errmsg.append("invalid contigid")
            elif pos > contigs[contigid]['length']:
                errmsg.append("invalid position on contig")
            elif contigid != previous_location[0]:
                previous_location = (contigid, pos)
            elif pos <= previous_location[1]:
                errmsg.append("position out of order")
            #  PRINT MESSAGES
            if errmsg:
                print(contigid, pos, allelesets, errmsg)
    elif mvf.metadata['mvftype'] == 'codon':
        print("codon checking coming soon")
    return ''

if __name__ == "__main__":
    main()
