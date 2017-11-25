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
from pylib.mvffasphy import fasta2mvf, mvf2fasta, mvf2phy
from pylib.mvfmaf import maf2mvf
from pylib import mvfanalysis
from pylib.mvfargo import MvfArgumentParser
from pylib.mvfcheck import check_mvf
from pylib.mvfwindowtree import infer_window_tree
from pylib.mvftranslate import annotate_mvf, translate_mvf
from pylib.mvfjoin import mvf_join
from pylib.mvfgroupallele import calc_group_unique_allele_window
from pylib.mvffilter import filter_mvf

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
            AnnotateMVF
            CalcCountCharacterWindows
            CalcDstatCombinations
            CalcPairwiseDistances
            CalcPairwiseDistancesWindow
            CalcSampleCoverage
            CalcDstatCombinations
            CalcPatternCount
            CheckMVF
            ConvertFASTAtoMVF
            ConvertMAFtoMVF
            ConvertMVFtoFASTA
            ConvertMVFtoPhylip
            ConvertVCFtoMVF
            FilterMVF
            InferWindowTree
            PlotChromoplot
            TranslateMVF
            """,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog=_LICENSE)
        parser.add_argument("command", help="MVFtools command to run")
        parser.add_argument("--version", action="version",
                            version="0.8.0",
                            help="display version information")
        args = parser.parse_args(self.arguments[:1])
        if not hasattr(self, args.command):
            print('Unrecognized command')
            parser.print_help()
            exit(1)
        getattr(self, args.command)()

    def PlotChromoplot(self):
        from pylib.chromoplot import chromoplot, Pallette

        def generate_argparser():
            pallette = Pallette()
            parser = MvfArgumentParser()
            parser.add_argument("--outprefix",
                                help="Output prefix (not required).")
            parser.add_argument("--samples", nargs='*', required=True,
                                help="3 or more taxa to use for quartets")
            parser.add_argument("--outgroup", nargs='*', required=True,
                                help="1 or more outgroups to use for quartets")
            parser.addarg_windowsize()
            parser.add_argument(
                "--contigs", nargs='*',
                help=("Enter the ids of one or more contigs in the "
                      "order they will appear in the chromoplot. "
                      "(defaults to all ids in order present in MVF)"))
            parser.add_argument(
                "--majority", action="store_true",
                help=("Plot only 100% shading in the majority track "
                      " rather than shaded proportions in all tracks."))
            parser.add_argument(
                "--info-track", "--infotrack", action="store_true",
                help=("Include an additional coverage information "
                      "track that will show empty, uninformative, "
                      "and informative loci. (Useful for "
                      "ranscriptomes/RAD or other reduced sampling."))
            parser.add_argument(
                "--empty-mask", "--emptymask", choices=pallette.colornames,
                default="none",
                help="Mask empty regions with this color.")
            parser.add_argument(
                "--yscale", default=20, type=int,
                help=("Height (in number of pixels) for each track"))
            parser.add_argument(
                "--xscale", default=1, type=int,
                help="Width (in number of pixels) for each window")
            parser.add_argument(
                "--colors", nargs=3, choices=pallette.colornames,
                help="three colors to use for chromoplot")
            parser.add_argument(
                "--plot-type", "--plottype",
                choices=["graph", "image"], default="image",
                help=("PNG image (default) or graph via matplotlib "
                      "(experimental)"))
            return parser
        parser = generate_argparser()
        args = parser.parse_args(self.arguments[1:])
        chromoplot(args)
        return ''

    def ConvertVCFtoMVF(self):
        def generate_argparser():
            parser = MvfArgumentParser(description="VCF conversion")
            parser.add_argument(
                "--vcf", type=os.path.abspath,
                help="VCF input file")
            parser.add_argument(
                "--out", required=True,
                help="output MVF file")
            parser.add_argument(
                "--outflavor", default="dna",
                choices=['dna', 'dnaqual', 'dnaqual-indel', 'dna-indel'],
                help=("choose output MVF flavor to include "
                      "quality scores and/or indels"))
            parser.add_argument(
                "--maskdepth", type=int, default=1,
                help="below this read depth mask with N/n")
            parser.add_argument(
                "--lowdepth", type=int, default=3,
                help=("below this read depth coverage, "
                      "convert to lower case set to 0 to disable"))
            parser.add_argument(
                "--maskqual", type=int, default=3,
                help=("low quality cutoff, bases replaced by N/- "
                      "set to 0 to disable"))
            parser.add_argument(
                "--lowqual", type=int, default=20,
                help=("below this quality convert to lower case "
                      "set to 0 to disable"))
            parser.add_argument(
                "--contigids", nargs='*',
                help=("manually specify one or more contig ids "
                      "as ID;VCFLABE;MVFLABEL, note that "
                      "VCFLABEL must match EXACTLY the contig string "
                      "labels in the VCF file"))
            parser.add_argument(
                "--samplereplace", nargs="*",
                help=("one or more TAG:NEWLABEL or TAG, items, "
                      "if TAG found in sample label, replace with "
                      "NEW (or TAG if NEW not specified) "
                      "NEW and TAG must each be unique"))
            parser.add_argument(
                "--reflabel", default="REF",
                help="label for reference sample (default='REF')")
            parser.add_argument(
                "--allelesfrom", default=None,
                help="""get additional alignment columns
                from INFO fields (:-separated)""")
            parser.addarg_linebuffer()
            parser.add_argument(
                "--no_autoindex", action="store_true",
                help="do not automatically index contigs from the VCF")
            parser.add_argument(
                "--fieldsep", default="TAB",
                choices=['TAB', 'SPACE', 'DBLSPACE', 'COMMA', 'MIXED'],
                help="""VCF field separator (default='TAB')""")
            parser.add_argument(
                "--qual", action="store_true",
                help="""Include Phred genotype quality (GQ) scores""")
            parser.addarg_overwrite()
            return parser
        parser = generate_argparser()
        args = parser.parse_args(self.arguments[1:])
        from pylib.vcf import vcf2mvf
        vcf2mvf(args)
        return ''

    def ConvertFASTAtoMVF(self):

        def generate_argparser():
            parser = MvfArgumentParser(description="fasta conversion")
            parser.add_argument(
                "--fasta", nargs='*', required=True,
                help="input FASTA file(s)")
            parser.add_argument(
                "--out", required=True,
                help="output MVF file")
            parser.add_argument(
                "--flavor", choices=['dna', 'protein'],
                help="type of file [dna] or protein", default='dna')
            parser.add_argument(
                "--contig-ids", "--contigids", nargs='*',
                help=("manually specify one or more contig ids "
                      "as ID:NAME"))
            parser.add_argument(
                "--sample-replace", "--samplereplace", nargs="*",
                help=("one or more TAG:NEWLABEL or TAG, items, "
                      "if TAG found in sample label, replace with "
                      "NEW (or TAG if NEW not specified) "
                      "NEW and TAG must each be unique"))
            parser.add_argument(
                "--ref-label", "--reflabel", default="REF",
                help="label for reference sample")
            parser.add_argument(
                "--read-buffer", "--readbuffer",
                type=int, default=100000,
                help="number of lines to hold in READ buffer")
            parser.add_argument(
                "--write-buffer", "--writebuffer",
                type=int, default=100000,
                help="number of lines to hold in WRITE buffer")
            parser.add_argument(
                "--field-sep", "--fieldsep", nargs='*', default=None,
                choices=['TAB', 'SPACE', 'DBLSPACE', 'COMMA',
                         'MIXED', 'PIPE', 'AT', 'UNDER', 'DBLUNDER'],
                help=("FASTA field separator; assumes "
                      "'>database accession locus' format"))
            parser.add_argument(
                "--contig-field", "--contigfield", type=int,
                help=("When headers are split by --field-sep, "
                      "the 0-based index of the contig id."))
            parser.add_argument(
                "--contig-by-file", "--contigbyfile", action="store_true",
                help=("Contigs are designated by separate files."))
            parser.add_argument(
                "--sample-field", "--samplefield", type=int,
                help=("when headers are split by --field-sep, "
                      "the 0-based index of the sample id"))
            parser.add_argument(
                "--manual-coord", "--manualcoord", nargs='*',
                help=("manually specify reference coordinates "
                      "for each file in the format "
                      "CONTIGID:START..STOP, ..."))
            parser.addarg_overwrite()
            return parser
        parser = generate_argparser()
        args = parser.parse_args(self.arguments[1:])
        fasta2mvf(args)
        return ''

    def ConvertMAFtoMVF(self):
        def generate_argparser():
            parser = MvfArgumentParser("MAF to MVF converter")
            parser.add_argument("--maf", help="input MAF file",
                                required=True, type=os.path.abspath,)
            parser.add_argument("--out", help="output MVF file",
                                type=os.path.abspath, required=True)
            parser.add_argument("--ref-tag", "--reftag",
                                help="old reference tag")
            parser.add_argument(
                "--mvf-ref-label", "--mvfreflabel", default="REF",
                help=("new label for reference sample (default='REF')"))
            parser.add_argument(
                "--sample-tags", "--sampletags", nargs="*", required=True,
                help=("one or more TAG:NEWLABEL or TAG, items, "
                      "if TAG found in sample label, replace with "
                      "NEW (or TAG if NEW not specified) "
                      "NEW and TAG must each be unique."))

            parser.addarg_linebuffer()
            parser.add_overwrite()
            return parser
        parser = generate_argparser()
        args = parser.parse_args(self.arguments[1:])
        maf2mvf(args)
        return ''

    def ConvertMVFtoFASTA(self):
        def generate_argparser():
            parser = MvfArgumentParser(description="MVF to FASTA")
            parser.addarg_mvf()
            parser.add_argument("--out", type=os.path.abspath,
                                help="Output path of FASTA file.",
                                required=True)
            parser.addarg_regions()
            parser.addarg_samples()
            parser.add_argument(
                "--labeltype", choices=('long', 'short'), default='long',
                help=("Long labels with all metadata or short ids"))
            parser.add_argument(
                "--outdata", choices=("dna", "rna", "prot"),
                help="Output dna, rna or prot data.")
            parser.add_argument(
                "--buffer", type=int, default=10,
                help="size (Mbp) of write buffer for each sample")
            parser.add_argument(
                "--tmpdir", default=".",
                help="directory to write temporary fasta files")
            return parser
        parser = generate_argparser()
        args = parser.parse_args(self.arguments[1:])
        mvf2fasta(args)
        return ''

    def MVFtoPhylip(self):
        def generate_argparser():
            parser = MvfArgumentParser()
            parser.addarg_mvf()
            parser.add_argument("--out", type=os.path.abspath,
                                help="Output Phylip file.", required=True)
            parser.addarg_contigs()
            parser.addarg_regions()
            parser.add_argument(
                "--labeltype", choices=('long', 'short'), default='short',
                help="Long labels with all metadata or short ids")
            parser.add_argument(
                "--outdata", choices=("dna", "rna", "prot"),
                help="Output dna, rna or prot data.")
            parser.addarg_samples()
            parser.add_argument(
                "--buffer", type=int, default=100000,
                help="size (bp) of write buffer for each sample")
            parser.add_argument(
                "--tmpdir", default=".",
                help="directory to write temporary fasta files")
            parser.add_argument(
                "--partition", action="store_true",
                help=("Output a CSV partitions file with RAxML"
                      "formatting for use in partitioned "
                      "phylogenetic methods."))
            return parser
        parser = generate_argparser()
        args = parser.parse_args(self.arguments[1:])
        mvf2phy(args)
        return ''

    def CalcSampleCoverage(self):

        def generate_argparser():
            parser = MvfArgumentParser()
            parser.addarg_mvf()
            parser.addarg_out()
            parser.addarg_contigs()
            parser.addarg_samples()
            return parser

        parser = generate_argparser()
        args = parser.parse_args(self.arguments[1:])
        mvfanalysis.calc_sample_coverage(args)
        return ''

    def CalcDstatCombinations(self):

        def generate_argparser():
            parser = MvfArgumentParser()
            parser.addarg_mvf()
            parser.addarg_out()
            parser.addarg_contigs()
            parser.addarg_samples()
            return parser

        parser = generate_argparser()
        args = parser.parse_args(self.arguments[1:])
        mvfanalysis.calc_dstat_combinations(args)
        return ''

    def CalcCountCharacterWindows(self):

        def generate_argparser():
            parser = MvfArgumentParser()
            parser.addarg_mvf()
            parser.addarg_out()
            parser.addarg_contigs()
            parser.addarg_samples()
            parser.addarg_mincoverage()
            parser.addarg_windowsize()
            parser.add_argument(
                "--base-match", "--basematch",
                help="String of bases to match (i.e. numerator).")
            parser.add_argument(
                "--base-total", "--basetotal",
                help="String of bases for total (i.e. denominator).")
            return parser
        parser = generate_argparser()
        args = parser.parse_args(self.arguments[1:])
        mvfanalysis.calc_count_char_window(args)
        return ''

    def CalcPairwiseDistances(self):

        def generate_argparser():
            parser = MvfArgumentParser()
            parser.addarg_mvf()
            parser.addarg_out()
            parser.addarg_windowsize()
            parser.addarg_mincoverage()
            return parser
        parser = generate_argparser()
        args = parser.parse_args(self.arguments[1:])
        mvfanalysis.calc_pairwise_distances(args)
        return ''

    def CalcPatternCount(self):

        def generate_argparser():
            parser = MvfArgumentParser()
            parser.addarg_mvf()
            parser.addarg_out()
            parser.addarg_windowsize()
            parser.addarg_mincoverage()
            return parser
        parser = generate_argparser()
        args = parser.parse_args(self.arguments[1:])
        mvfanalysis.calc_pattern_count(args)
        return ''

    def CheckMVF(self):

        def generate_argparser():
            parser = MvfArgumentParser()
            parser.addarg_mvf()
            return parser

        parser = generate_argparser()
        args = parser.parse_args(self.arguments[1:])
        check_mvf(args)
        return ''

    def InferGroupSpecificAllele(self):

        def generate_argparser():
            parser = MvfArgumentParser()
            parser.addarg_mvf()
            parser.addarg_gff()
            parser.addarg_out()
            parser.addarg_contigs()
            parser.addarg_samples()
            parser.addarg_windowsize()
            parser.addarg_mincoverage()
            parser.add_argument(
                "--allele-groups", "--allelegroups", nargs='*',
                help=("GROUP1:LABEL,LABEL GROUP2:LABEL,LABEL "))
            parser.add_argument(
                "--species-groups", "--speciesgroups", nargs='*')
            parser.add_argument(
                "--chi-test", "--chitest",
                help=("Input two number values for expected "
                      "Nonsynonymous and Synonymous expected values."))
            parser.add_argument(
                "--target", nargs="*",
                help=("Specify the taxa labels that define the "
                      "target lineage-specific branch to be tested."))
            parser.add_argument(
                "--num-target-species", "--targetspec",
                type=int, default=1,
                help=("Specify the minimum number of "
                      "taxa in the target set "
                      "that are required to conduct analysis"))
            parser.add_argument(
                "--output-align", "--outputalign",
                help=("Output alignment to this file path in "
                      "phylip format."))
            parser.add_argument(
                "--outgroup",
                help("Specify sample name with which to root trees."))
            parser.add_argument(
                "--uselabels", "--use_labels", action="store_true",
                help="Use contig labels instead of IDs in output.")
            parser.add_argument(
                "--codemlpath", "--codeml-path", default="codeml",
                type=os.path.abspath,
                help="Full path for PAML codeml executable.")
            parser.add_argument(
                "--raxml-path", "--raxmlpath",
                type=os.path.abspath, default="raxml",
                help="Full path to RAxML program executable.")
            parser.add_argument(
                "--start-contig", "--startcontig", type=int, default=0,
                help="Numerical ID for the starting contig.")
            parser.add_argument(
                "--end-contig", "--endcontig",
                type=int, default=100000000,
                help="Numerical id for the ending contig.")
            parser.add_argument(
                "--pamltmp", "--paml-tmp",
                default="pamltmp", type=os.path.abspath,
                help="path for temporary folder for PAML output files")
            parser.add_argument(
                "--all-sample-trees", "--allsampletrees",
                action="store_true",
                help=("Makes trees from all samples instead of "
                      "only the most complete sequence from "
                      "each species"))
            parser.add_argument(
                "--branchlrt", "--branch-lrt", type=os.path.abspath,
                help=("Specify the output file for and turn on the "
                      "RAxML-PAML format LRT test scan for "
                      "selection on the target branch in addition "
                      "to the basic patterns scan"))
            return parser
        parser = generate_argparser()
        args = parser.parse_args(self.arguments[1:])
        calc_group_unique_allele_window(args)
        return ''

    def InferWindowTree(self):
        def generate_argparser():
            parser = MvfArgumentParser()
            parser.addarg_mvf()
            parser.addarg_out()
            parser.addarg_samples()
            parser.addarg_contigs()
            parser.addarg_windowsize()
            parser.add_argument(
                "--raxml-outgroups", "--raxmloutgroups",
                help=("Comma-separated list of outgroup "
                      "taxon labels to use in RAxML."))
            parser.add_argument(
                "--root-with", "--rootwith",
                help=("Comma-separated list of taxon labels "
                      "to root trees with after RAxML"))
            parser.add_argument(
                "--output-contig-labels", "--outputcontiglabels",
                action="store_true",
                help=("Output will use contig labels "
                      "instead of id numbers."))
            parser.add_argument(
                "--outputempty", action="store_true",
                help=("Include entries of windows "
                      "with no data in output."))
            parser.add_argument(
                "--choose-allele", "--chooseallele", "--hapmode",
                default="none", dest="choose_allele",
                choices=["none", "randomone", "randomboth",
                         "major", "minor", "majorminor"],
                help=("Chooses how heterozygous alleles are "
                      "handled. (none=no splitting (default); "
                      "randomone=pick one allele randomly "
                      "(recommended); randomboth=pick two alleles "
                      "randomly, but keep both; major=pick the "
                      "more common allele; minor=pick the less "
                      "common allele; majorminor= pick the major in "
                      "'a' and minor in 'b'"))
            parser.add_argument(
                "--minsites", type=int, default=100,
                help="minimum number of sites ")
            parser.add_argument(
                "--minseqcoverage", type=float, default=0.1,
                help=("proportion of total alignment a sequence"
                      "must cover to be retianed [0.1]"))
            parser.add_argument(
                "--mindepth", type=int, default=4,
                help=("minimum number of alleles per site"))
            parser.add_argument(
                "--bootstrap", type=int,
                help=("turn on rapid bootstrapping for RAxML and "
                      "perform specified number of replicates"))
            parser.add_argument(
                "--raxml_model", default="GTRGAMMA",
                help=("choose RAxML model"))
            parser.add_argument(
                "--raxml-path", "--raxmlpath", default="raxml",
                help="RAxML path for manual specification.")
            parser.add_argument(
                "--raxmlopts", default="",
                help=("specify additional RAxML arguments as a "
                      "double-quotes encased string"))
            parser.add_argument(
                "--duplicate-seq", "--duplicateseq", default="dontuse",
                choices=["dontuse", "keep", "remove"],
                help=("dontuse=remove duplicate sequences prior to "
                      "RAxML tree inference, then add them to the "
                      "tree manually as zero-branch-length sister "
                      "taxa; keep=keep in for RAxML tree inference "
                      "(may cause errors for RAxML); "
                      "remove=remove entirely from alignment"))
            parser.add_argument(
                "--tempdir", default='./raxmltemp',
                type=os.path.abspath,
                help=("Temporary directory path"))
            parser.add_argument(
                "--tempprefix", default="mvftree",
                help=("Temporary file prefix"))
            return parser
        parser = generate_argparser()
        args = parser.parse_args(self.arguments[1:])
        infer_window_tree(args)
        return ''

    def FilterMVF(self):

        def generate_argparser():
            parser = MvfArgumentParser()
            parser.addarg_mvf()
            parser.addarg_out()
            parser.add_argument(
                "--actions", nargs='*',
                help=("set of actions:args to perform, "
                      "note these are done in order as listed"))
            parser.add_argument(
                "--labels", action="store_true",
                help="use sample labels instead of indices")
            parser.add_argument(
                "--test", help="manually input a line for testing")
            parser.add_argument(
                "--test-nchar", type=int,
                help="total number of samples for test string")
            parser.add_argument(
                "--morehelp", action="store_true",
                help="prints full module list and descriptions")
            parser.addarg_linebuffer()
            parser.add_argument(
                "--verbose", action="store_true",
                help="report every line (for debugging)")
            parser.addarg_overwrite()
            return parser

        parser = generate_argparser()
        args = parser.parse_args(self.arguments[1:])
        filter_mvf(args)
        return ''

    def AnnotateMVF(self):

        def generate_argparser():
            parser = MvfArgumentParser()
            parser.addarg_mvf()
            parser.addarg_out()
            parser.addarg_gff()
            parser.add_argument(
                "--filter-annotation", "--filterannotation",
                help=("Skip entries in the GFF file that "
                      "contain this string in their 'Notes'"))
            parser.add_argument(
                "--nongenic-mode", "--nongenicmode",
                action="store_true",
                help=("Instead of returning annotated genes, "
                      "return the non-genic regions without "
                      "without changing contigs or coordinates"))
            parser.add_argument(
                "--nongenic-margin", "--nongenicmargin",
                type=int, default=0,
                help=("for --unnanotated-mode, only retain "
                      "positions that are this number of bp away "
                      "from an annotated region boundary"))
            parser.addarg_linebuffer()
            return parser

        parser = generate_argparser()
        args = parser.parse_args(self.arguments[1:])
        annotate_mvf(args)
        return ''

    def TranslateMVF(self):

        def generate_argparser():
            parser = MvfArgumentParser()
            parser.addarg_mvf()
            parser.addarg_out()
            parser.add_argument(
                "--gff", type=os.path.abspath,
                help=("Input GFF3 file. If GFF3 not provided, "
                      "alignments are assumed to be "
                      "in-frame coding sequences."))
            parser.add_argument(
                "--outtype", choices=['protein', 'codon'],
                default="codon",
                help=("protein=single data column "
                      "of protein alleles; "
                      "codon=four columns with: "
                      "protein frame1 frame2 frame3"))
            parser.add_argument(
                "--filter-annotation", "--filterannotation",
                help=("skip GFF entries with text "
                      "matching this in their 'Notes' field"))
            parser.addarg_linebuffer()
            parser.addarg_overwrite()
            return parser
        parser = generate_argparser()
        args = parser.parse_args(self.arguments[1:])
        translate_mvf(args)
        return ''

    def JoinMVF(self):

        def generate_argparser():
            parser = MvfArgumentParser()
            parser.add_argument(
                "--mvf", nargs="*", type=os.path.abspath, required=True,
                help="One or more mvf files.")
            parser.addarg_out()
            parser.add_argument(
                "--new-contigs", "--newcontigs", action="store_true",
                help=("By default, contigs are matched between files "
                      "using their text labels in the header. "
                      "Use this option to turn matching off and treat "
                      "each file's contigs as distinct."))
            parser.add_argument(
                "--newsamples", action="store_true",
                help=("By default, samples are matched between files "
                      "using their text labels in the header. "
                      "Use this option to turn matching off and treat "
                      "each file's sample columns as distinct."))
            parser.add_argument(
                "--main_header_file", "--mainheaderfile",
                help=("Output file will use same headers as "
                      "this input file (default=first in list)."))
            parser.addarg_linebuffer()
            parser.addarg_overwrite()
            return parser
        parser = generate_argparser()
        args = parser.parse_args(self.arguments[1:])
        mvf_join(args)
        return ''


if __name__ == "__main__":
    time0 = time()
    MVFcall()
    print("Total time: ", time() - time0)
