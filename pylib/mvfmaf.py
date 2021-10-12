#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This program analyzes a DNA MVF alignment using the modules specified below,
use the --morehelp option for additional module information.

MVFtools: Multisample Variant Format Toolkit
James B. Pease and Ben K. Rosenzweig
http://www.github.org/jbpease/mvftools

If you use this software please cite:
Pease, James B. and Benjamin K. Rosenzweig. 2018.
"Encoding Data Using Biological Principles: the Multisample Variant Format
for Phylogenomics and Population Genomics"
IEEE/ACM Transactions on Computational Biology and Bioinformatics. 15(4) 1231–1238.
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

# TODO: handle multiple contigs in same file -
# default 's SAMPLE.CONTIG START LENGTH etc...', allow user to
# specify alternate formats
# TODO: filter alignment blocks by score (i.e. reject all blocks below score X)
# TODO: discover sample names without user input (either include an option to
# guarantee that all samples are represented in all alignment blocks, or do a
# first-pass over the whole file to find all sample names that will appear.
# This second option should be combined with some sort of indexing scheme
# since it forces us to preprocess the whole file anyway...)

import os
import gzip
import re
from pylib.mvfbase import encode_mvfstring, MultiVariantFile


RE_CONTIG_NAME = re.compile("ID=(.*?),")
RE_CONTIG_LENGTH = re.compile("length=(.*?)>")


class MultiAlignFile(object):
    """Multiple Alignment File (MAF v1) handler

    Attributes:
        alignments: Ordered list of MAF alignment objects
        meta: Dictionary of (key,value) metadata information
        -lengths: Dictionary of total lengths of contigs
        -name_index: Dictionary of names with coordinates of sequence data
        -path: filepath for the MAF
        overwrite=True to overwrite the file
    """

    def __init__(self, args):
        self.path = os.path.abspath(args.maf)
        self.ref = args.ref_tag
        self.metadata = {'sourceformat': 'MAF'}
        self.metadata['labels'] = set()
        self.entrystart = 0
        # Check for Gzip and establish file object
        self.metadata['isgzip'] = (self.path.endswith(".gz") or
                                   args.get('isgzip', False))
        # READ MODE
        filehandler = (self.metadata['isgzip'] and
                       gzip.open(self.path, 'rt') or open(self.path, 'rt'))
        # Process header lines
        line = filehandler.readline()
        if line.startswith("##maf"):
            line = filehandler.readline()
        if line[0] == 'a':
            line = filehandler.readline()
            while line.strip != '' and line[0] != 'a':
                if line[0] == 's':
                    label = line[1:].strip().split()[0]
                    label = label.split('.')[0]
                    self.metadata['labels'].add(label)
                line = filehandler.readline()
        filehandler.close()

    def __iter__(self):
        """Simple entry iterator
           Returns (str(chrom), int(pos), list(allele entries))
        """
        if self.metadata['isgzip']:
            filehandler = gzip.open(self.path, 'rt')
        else:
            filehandler = open(self.path, 'rt')
        filehandler.seek(self.entrystart)
        line = filehandler.readline()
        while line:
            start_pos = -1
            if line[0] == 'a':
                block = {}
                block_length = -1
                line = filehandler.readline()
                while line.strip != '' and line[0] != 'a':
                    if line[0] == 's':
                        (label, start, length,
                         _, _, seq) = (
                             line[1:].strip().split())
                        label = label.split('.')[0]
                        block[label] = seq
                        if label == self.ref:
                            start_pos = int(start)
                            block_length = int(length)
                    line = filehandler.readline()
                self.metadata['labels'].update(block.keys())
                yield (start_pos, block_length, block)
            line = filehandler.readline()
        filehandler.close()

#    def reorient(self, seq, orientation, start, length):
#        if orientation == '+':
#            return seq, start
#        elif orientation == '-':
#            return ''.join(self.WatsonCrick[b] for b in seq[::-1]), start
#        else:
#            raise RuntimeError(
#                "Error: orientation {} not recognized".format(orientation))
#        return "", None


def maf2mvf(args):
    """Main method"""
    # ESTABLISH MAF
    args.qprint("Starting ConvertMAF2MVF")
    maf = MultiAlignFile(args)
    args.qprint("MAF Established")
    # ESTABLISH MVF
    mvf = MultiVariantFile(args.out, 'write', overwrite=args.overwrite)
    args.qprint("MVF output initialized")
    # PROCESS SAMPLE INFO
    contig_translate = {1: 1}
    samplelabels = [s.split(':')[0] for s in args.sample_tags.split(',')]
    args.qprint("Sample tags processed: {}".format(samplelabels))
    if args.ref_tag not in samplelabels:
        raise IndexError("--ref-tag not in the tags listed in --sample-tags")
    samplelabels.remove(args.ref_tag)
    samplelabels.insert(0, args.ref_tag)
    mvf.sample_ids = samplelabels[:]
    mvf.sample_indices = list(range(len(mvf.sample_ids)))
    for i, label in enumerate(samplelabels):
        mvf.sample_data[i] = {'id': label, 'index': i}
    mvf.reset_max_sample()
    mvf.metadata['sourceformat'] = maf.metadata['sourceformat']
    mvf.metadata.notes.append(args.command_string)
    # WRITE MVF HEADER
    mvf.write_data(mvf.get_header())
    args.qprint("MAF Headers Written")
    mvfentries = []
    nentry = 0
    total_entries = 0
    args.qprint("Begin data conversion")
    for pos, length, msa in maf:
        for sname in samplelabels:
            if sname not in msa:
                msa[sname] = '-'*length
        msa['contig'] = 1
        for i in range(length):
            mvf_alleles = encode_mvfstring(
                ''.join(msa[s][i].strip() for s in samplelabels))
            if mvf_alleles:
                mvfentries.append(
                    (contig_translate.get(msa['contig']),
                     pos+i, (mvf_alleles,)))
                nentry += 1
                if nentry == args.line_buffer:
                    total_entries += nentry
                    mvf.write_entries(mvfentries, encoded=True)
                    args.qprint("{} entries written".format(total_entries))
                    mvfentries = []
                    nentry = 0
    if mvfentries:
        total_entries += nentry
        mvf.write_entries(mvfentries)
        args.qprint("{} entries written".format(total_entries))
    args.qprint("Complete.")
    return ''
