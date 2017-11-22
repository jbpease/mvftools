#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MVFtools: Multisample Variant Format Toolkit
James B. Pease and Ben K. Rosenzweig
http://www.github.org/jbpease/mvftools
"""

import os
import sys
import argparse
from time import time
from pylib.chromoplot import chromoplot, Pallette
from pylib.vcf import vcf2mvf
from pylib.fasta import fasta2mvf, mvf2fasta
from pylib.maf import maf2mvf
from pylib.phylip import mvf2phy

_LICENSE = """
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


class MVFcall(object):

    def __init__(self, arguments=None):
        """Main method for vcf2mvf"""
        self.arguments = arguments if arguments is not None else sys.argv[1:]
        parser = argparse.ArgumentParser(
            prog="mvftools.py",
            usage="""Choose from the following commands:
            CalcCountCharacterWindows
            CalcPairwiseDistances
            CalcPairwiseDistancesWindow
            CalcSampleCoverage
            CalcDAllTrios
            CalcPatternCount
            CalcPatternList
            ConvertFASTAtoMVF
            ConvertMAFtoMVF
            ConvertMVFtoFASTA
            ConvertMVFtoPhylip
            ConvertVCFtoMVF
            PlotChromoplot
            """,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog=_LICENSE)
        parser.add_argument("command", help="MVFtools command to run")
        args = parser.parse_args(self.arguments[:1])
        if not hasattr(self, args.command):
            print('Unrecognized command')
            parser.print_help()
            exit(1)
        getattr(self, args.command)()

    def PlotChromoplot(self):
        def generate_argparser():
            pallette = Pallette()
            parser = argparse.ArgumentParser(
                prog="mvf_chromoplot.py",
                description=__doc__,
                formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                epilog=_LICENSE)
            parser.add_argument("-i", "--mvf", help="Input MVF file.",
                                required=True,
                                type=os.path.abspath)
            parser.add_argument("-o", "--outprefix",
                                help="Output prefix (not required).")
            parser.add_argument("-s", "--samples", nargs='*', required=True,
                                help="3 or more taxa to use for quartets")
            parser.add_argument("-G", "--outgroup", nargs='*', required=True,
                                help="1 or more outgroups to use for quartets")
            parser.add_argument("-w", "--windowsize", type=int, default=100000)
            parser.add_argument("-c", "--contigs", nargs='*',
                                help=("Enter the ids of one or more"
                                      "contigs in the "
                                      "order they will appear "
                                      "in the chromoplot. "
                                      "(defaults to all ids in "
                                      "order present in MVF)"))
            parser.add_argument(
                "-M", "--majority", action="store_true",
                help=("Plot only 100% shading in the majority track "
                      " rather than shaded proportions in all tracks."))
            parser.add_argument(
                "-I", "--infotrack", action="store_true",
                help=("Include an additional coverage information "
                      "track that will show empty, uninformative, "
                      "and informative loci. (Useful for "
                      "ranscriptomes/RAD or other reduced sampling."))
            parser.add_argument(
                "-E", "--emptymask", choices=pallette.colornames,
                default="none",
                help="Mask empty regions with this color.")
            parser.add_argument("-y", "--yscale", default=20, type=int,
                                help=("Height (in number of pixels)"
                                      " for each track"))
            parser.add_argument(
                "-x", "--xscale", default=1, type=int,
                help="Width (in number of pixels) for each window")
            parser.add_argument(
                "-C", "--colors", nargs=3, choices=pallette.colornames,
                help="three colors to use for chromoplot")
            parser.add_argument("-q", "--quiet", action="store_true",
                                help="suppress all output messages")
            parser.add_argument("-P", "--plottype", choices=["graph", "image"],
                                default="image",
                                help=("PNG image (default) "
                                      "or graph via matplotlib (experimental)"))
            parser.add_argument("-v", "--version", action="version",
                                version="2017-08-14",
                                help="display version information")
            return parser
        parser = generate_argparser()
        args = parser.parse_args(self.arguments[1:])
        chromoplot(args)
        return ''

    def ConvertVCFtoMVF(self):
        def generate_argparser():
            parser = argparse.ArgumentParser(description="VCF conversion")
            parser.add_argument("inputfile",
                                help="input file (varies with application)",
                                type=os.path.abspath)
            parser.add_argument("--out", help="output MVF file", required=True)
            parser.add_argument("--outflavor",
                                choices=['dna', 'dnaqual', 'dnaqual-indel',
                                         'dna-indel'], default='dna',
                                help=("choose output MVF flavor to include "
                                      "quality scores and/or indels"))
            parser.add_argument("--maskdepth", type=int, default=1,
                                help="below this read depth mask with N/n")
            parser.add_argument("--lowdepth", type=int, default=3,
                                help=("below this read depth coverage, "
                                      "convert to lower case set to 0 to disable"))
            parser.add_argument("--maskqual", type=int, default=3,
                                help="""low quality cutoff, bases replaced by N/-
                                     set to 0 to disable""")
            parser.add_argument("--lowqual", type=int, default=20,
                                help="""below this quality convert to lower case
                                        set to 0 to disable""")
            parser.add_argument("--contigids", nargs='*',
                                help=("""manually specify one or more contig ids
                                         as ID;VCFLABE;MVFLABEL, note that
                                         VCFLABEL must match EXACTLY the contig string
                                         labels in the VCF file"""))
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
            parser.add_argument("--linebuffer", type=int, default=100000,
                                help="number of lines to hold in read/write buffer")
            parser.add_argument("--no_autoindex", action="store_true",
                                help="do not automatically index contigs from the VCF")
            parser.add_argument("--fieldsep", default="TAB",
                                choices=['TAB', 'SPACE', 'DBLSPACE', 'COMMA', 'MIXED'],
                                help="""VCF field separator (default='TAB')""")
            parser.add_argument("--qual", action="store_true",
                                help="""Include Phred genotype quality (GQ) scores""")
            parser.add_argument("--overwrite", action="store_true",
                                help="USE WITH CAUTION: force overwrite of outputs")
            return parser
        parser = generate_argparser()
        args = parser.parse_args(self.arguments[1:])
        vcf2mvf(args)
        return ''

    def ConvertFASTAtoMVF(self):

        def generate_argparser():
            parser = argparse.ArgumentParser(description="fasta conversion")
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
        parser = generate_argparser()
        args = parser.parse_args(self.arguments[1:])
        fasta2mvf(args)
        return ''
    def ConvertMAFtoMVF(self):
        def generate_argparser():
            parser = argparse.ArgumentParser("MAF to MVF converter")
            parser.add_argument("-i", "--maf", help="input MAF file",
                                required=True, type=os.path.abspath,)
            parser.add_argument("-o", "--out", help="output MVF file",
                                type=os.path.abspath, required=True)
            parser.add_argument("-R", "--ref-tag", "--reftag",
                                help="old reference tag")
            parser.add_argument("-M", "--mvf-ref-label", "--mvfreflabel",
                                default="REF",
                                help=("new label for reference sample "
                                      "(default='REF')"))
            parser.add_argument("-s", "--sample-tags", "--sampletags", nargs="*",
                                help=("one or more TAG:NEWLABEL or TAG, items, "
                                      "if TAG found in sample label, replace with "
                                      "NEW (or TAG if NEW not specified) "
                                      "NEW and TAG must each be unique."),
                                required=True)
            parser.add_argument("-B", "--line-buffer", "--linebuffer",
                                type=int, default=100000,
                                help="number of lines to hold in read/write buffer")
            parser.add_argument("--overwrite", action="store_true")
            parser.add_argument("--version", action="version",
                                    version="2017-06-24",
                                    help="display version information")
            return parser
        parser = generate_argparser()
        args = parser.parse_args(self.arguments[1:])
        maf2mvf(args)
        return ''

    def ConvertMVFtoFASTA(self):
        def generate_argparser():
            parser = argparse.ArgumentParser(description="MVF to FASTA")
            parser.add_argument("-i", "--mvf", type=os.path.abspath,
                                help="Input MVF file.", required=True)
            parser.add_argument("-o", "--out", type=os.path.abspath,
                                help="target FASTA file", required=True)
            parser.add_argument("-r", "--regions", nargs='*', required=True,
                                help=("A file path to a plain-text file with"
                                      "one region per line formatted as"
                                      "formatted as: contigid,start,stop"
                                      "(coordinates are inclusive)"))
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
                                version="2017-09-26",
                                help="display version information")
            return parser
        parser = generate_argparser()
        args = parser.parse_args(self.arguments[1:])
        mvf2fasta(args)
        return ''

    def MVFtoPhylip(self):
        def generate_argparser():
            parser = argparse.ArgumentParser(description="MVFtoPhylip")
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
        parser = generate_argparser()
        args = parser.parse_args(self.arguments[1:])
        mvf2phy(args)
        return ''


    def CalcSampleCoverage():

        def generate_argparser():
            parser = argparse.ArgumentParser()
            parser.add_argument("--mvf", type=os.path.abspath,
                                required=True, help="Input MVF file.")
            parser.add_argument("--out", help="output file", required=True,
                                type=os.path.abspath)
            parser.add_argument("-c", "--contigs", nargs='*',
                                help="limit analyses to these contigs")
            parser.add_argument("-s", "--samples", nargs='*',
                                help="limit analyses to these samples")
            return parser

        parser = generate_argparser()
        args = parser.parse_args(self.arguments[1:])
        calcsamplecoverage(args)
        return ''



    def CalcCountCharWindow():

        def generate_argparser():
            parser = argparse.ArgumentParser()
            parser.add_argument("--mvf", type=os.path.abspath, required=True, help="Input MVF file.")
            parser.add_argument("-o", "--out", help="output file", required=True,
                                type=os.path.abspath)
            parser.add_argument("-c", "--contigs", nargs='*',
                                help="limit analyses to these contigs")
            parser.add_argument("-s", "--samples", nargs='*',
                                help="limit analyses to these samples")
            parser.add_argument("-m", "--mincoverage", type=int,
                                help="mininum sample coverage for site")
            parser.add_argument("-w", "--windowsize", type=int, default=100000,
                                help="""window size, use -1 to use whole contigs""")
            parser.add_argument("--base-match",
                                help=("[BaseCountWindow] string of "
                                      "bases to match (i.e. numerator)."))
            parser.add_argument("--base-total",
                                help=("[BaseCountWindow] string of bases "
                                      "for total (i.e. denominator)."))
            return parser



if __name__ == "__main__":
    time0 = time()
    MVFcall()
    print("Total time: ", time() - time0)
