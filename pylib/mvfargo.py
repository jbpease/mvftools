# -*- coding: utf-8 -*-
"""
MVFtools: Multisample Variant Format Toolkit
James B. Pease and Ben K. Rosenzweig
http://www.github.org/jbpease/mvftools

If you use this software please cite:
Pease, James B. and Benjamin K. Rosenzweig. 2018.
"Encoding Data Using Biological Principles: the Multisample Variant Format
for Phylogenomics and Population Genomics"
IEEE/ACM Transactions on Computational Biology and Bioinformatics.
15(4) 1231–1238.
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

import os
import argparse


class MvfArgumentParser(argparse.ArgumentParser):
    """Creates common argparse elements for various MVF modules"""

    def __init__(self):
        version = '0.6.2'
        super(MvfArgumentParser, self).__init__()
        setattr(self, 'formatter_class',
                argparse.ArgumentDefaultsHelpFormatter)
        self.add_argument("--version", action="version",
                          version=version,
                          help="Display version information.")
        self.add_argument("--versionx", default=version,
                          help=argparse.SUPPRESS)
        self.add_argument(
            "--quiet", action="store_true",
            help="Suppress screen output.")

    def addarg_contig_ids(self):
        """Contig ids arguments"""
        self.add_argument(
            "--contig-ids", "--contigids", nargs=1,
            help=("Specify comma-separated list of contig short ids. "
                  "Must match exactly. Do not use with --contig-labels."))

    def addarg_contig_labels(self):
        """Contig labels arguments"""
        self.add_argument(
            "--contig-labels", "--contiglabels", nargs=1,
            help=("Specify comma-separated list of contig full labels. "
                  "Must match exactly. Do not use with --contig-ids"))

    def addarg_gff(self):
        """GFF argument"""
        self.add_argument("--gff", type=os.path.abspath,
                          help="Input gff annotation file.")

    def addarg_linebuffer(self):
        """Linebuffer argument"""
        self.add_argument(
            "--line-buffer", "--linebuffer", type=int, default=100000,
            help="Number of entries to store in memory at a time.")

    def addarg_mincoverage(self):
        """Minimum coverage argument"""
        self.add_argument(
            "--mincoverage", type=int,
            help="Mininum sample coverage for sites.")

    def addarg_mvf(self, required=True):
        """MVF file argument"""
        self.add_argument(
            "--mvf", type=os.path.abspath, required=required,
            help="Input MVF file.")

    def addarg_overwrite(self):
        """Overwrite argument"""
        self.add_argument(
            "--overwrite", action="store_true",
            help="USE WITH CAUTION: force overwrite of outputs")

    def addarg_out(self, required=True):
        """Output argument"""
        self.add_argument(
            "--out", help="Output file",
            required=required, type=os.path.abspath)

    def addarg_outgroup_indices(self, nmin=None):
        """Outgroup indices argument"""
        self.add_argument(
            "--outgroup-indices", "--outgroupindices",
            nargs=1,
            help=("Specify comma-separated list of {}outgroup "
                  "sample numerical indices (first column is 0). "
                  "Leave blank for all samples. "
                  "Do not use with --outgroup_labels.".format(
                      str(nmin) + " or more " if nmin is not None
                      else "")))

    def addarg_outgroup_labels(self, nmin=None):
        """Outgroup labels argument"""
        self.add_argument(
            "--outgroup-labels", "--outgrouplabels",
            nargs=1,
            help=("Specify comma-separated list of {}outgroup "
                  "sample labels. "
                  "Labels must be exact (case-sensitive). "
                  "Leave blank for all samples."
                  "Do not use with --outgroup_indicies.".format(
                      str(nmin) + " or more " if nmin is not None
                      else "")))

    def addarg_sample_indices(self, nmin=None):
        """Sample indices argument"""
        self.add_argument(
            "--sample-indices", "--sampleindices",
            nargs=1,
            help=("Specify comma-separated list of {}sample "
                  "numerical indices (first sample is 0). "
                  "Leave blank for all samples. "
                  "Do not use with --sample_labels.".format(
                      str(nmin) + " or more " if nmin is not None
                      else "")))

    def addarg_sample_labels(self, nmin=None):
        """Sample labels argument"""
        self.add_argument(
            "--sample-labels",
            nargs=1,
            help=("Specify comma-separated list of {}sample labels. "
                  "Labels must be exact (case-sensitive). "
                  "Leave blank for all samples."
                  "Do not use with --sample_indicies.".format(
                      str(nmin) + " or more " if nmin is not None
                      else "")))

    def addarg_windowsize(self):
        """Window size argument"""
        self.add_argument(
            "--windowsize", default=100000,
            action=int_range_action(-1, 'Inf'),
            help=("Set integer window size. "
                  "Use 0 for whole file. "
                  "Use -1 for whole contigs. "))

    def addarg_regions(self):
        """Regions argument"""
        self.add_argument(
            "--regions", type=os.path.abspath,
            help=("Path of a plain text file containing one more lines "
                  "with entries 'contigid,stop,start' "
                  "(one per line, inclusive coordinates) "
                  "all data will be returned if left blank."))


def int_range_action(min_value, max_value):
    """Reformats a range of integers"""
    class _IntRangeAction(argparse.Action):

        def __init__(self, option_strings, dest, nargs=None, **kwargs):
            if nargs is not None:
                raise ValueError("nargs not allowed")
            super(_IntRangeAction, self).__init__(
                option_strings, dest, **kwargs)

        def __call__(self, parser, namespace, values, option_string=None):
            values = int(values)
            if min_value != '-Inf':
                if values < min_value:
                    parser.error("Minimum value for {0} is {1}".format(
                        option_string, min_value))
            if max_value != 'Inf':
                if values > max_value:
                    parser.error("Maximum value for {0} is {1}".format(
                        option_string, max_value))
            setattr(namespace, self.dest, values)
    return _IntRangeAction


def mutex_check(args):
    """Check for Mutually Exclusive arguments in an
       argparse parsed ArgumentParser"""
    if "sample_indices" in args and "sample_labels" in args:
        if (args.sample_indices is not None
                and args.sample_labels is not None):
            raise RuntimeError(
                "--sample-indices and --sample-labels should not "
                "be used simulataneously.")
    if "outgroup_indices" in args and "outgroup_labels" in args:
        if (args.outgroup_indices is not None
                and args.outgroup_labels is not None):
            raise RuntimeError(
                "--outgroup-indices and --outgroup-labels should not "
                "be used simulataneously.")
    if "contig_ids" in args and "contig_labels" in args:
        if (args.contig_ids is not None
                and args.contig_labels is not None):
            raise RuntimeError(
                "--contig-ids and --contig-labels should not "
                "be used simulataneously.")
