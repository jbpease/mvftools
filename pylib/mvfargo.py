# -*- coding: utf-8 -*-
"""Custom Argparse Calls"""

import os
import argparse


class MvfArgumentParser(argparse.ArgumentParser):

    def __init__(self):
        super(MvfArgumentParser, self).__init__()
        setattr(self, 'formatter_class',
                argparse.ArgumentDefaultsHelpFormatter)
        self.add_argument("--version", action="version",
                          version="0.8.0",
                          help="Display version information.")
        self.add_argument(
            "--quiet", action="store_true", default=True,
            help="Suppress screen output.")

    def addarg_contigs(self):
        self.add_argument(
            "--contigs",
            help="Specify comma-separated list of contigs.")

    def addarg_gff(self):
        self.add_argument("--gff", type=os.path.abspath,
                          help="Input gff annotation file.")

    def addarg_linebuffer(self):
        self.add_argument(
            "--line-buffer", "--linebuffer", type=int, default=100000,
            help="Number of entries to store in memory at a time.")

    def addarg_mincoverage(self):
        self.add_argument(
                "--mincoverage", type=int,
                help="Mininum sample coverage for sites.")

    def addarg_mvf(self):
        self.add_argument(
            "--mvf", type=os.path.abspath, required=True,
            help="Input MVF file.")

    def addarg_overwrite(self):
        self.add_argument(
            "--overwrite", action="store_true",
            help="USE WITH CAUTION: force overwrite of outputs")

    def addarg_out(self):
        self.add_argument(
            "--out", help="Output file",
            required=True, type=os.path.abspath)

    def addarg_samples(self):
        self.add_argument(
            "--samples",
            help=("Specify comma-separated list of samples, "
                  "Leave blank for all samples."))

    def addarg_windowsize(self):
        self.add_argument(
            "--windowsize", default=100000,
            action=int_range_action(-1, 'Inf'),
            help=("Set integer window size. "
                  "Use 0 for whole file. "
                  "Use -1 for whole contigs. "))

    def addarg_regions(self):
        self.add_argument(
            "--regions", type=os.path.abspath,
            help=("Path of a plain text file containing one more lines "
                  "with entries 'contigid,stop,start' "
                  "(one per line, inclusive coordinates) "
                  "all data will be returned if left blank."))


def int_range_action(min_value, max_value):
    class IntRangeAction(argparse.Action):

        def __init__(self, option_strings, dest, nargs=None, **kwargs):
            if nargs is not None:
                raise ValueError("nargs not allowed")
            super(IntRangeAction, self).__init__(
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
    return IntRangeAction
