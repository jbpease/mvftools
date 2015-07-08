#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
MVFtools: Multisample Variant Format Toolkit
http://www.github.org/jbpease/mvftools

MAF2MVF: MAF to MVF conversion program
@author: James B. Pease
@author: Ben K. Rosenzweig

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
import sys
import argparse
import gzip
import re
from mvfbase import encode_mvfstring, MultiVariantFile

# TODO: handle multiple contigs in same file - default 's SAMPLE.CONTIG START LENGTH etc...', allow user to specify alternate formats
# TODO: filter alignment blocks by score (i.e. reject all blocks below score X)
# TODO: discover sample names without user input (either include an option to guarantee that all samples are represented in all alignment blocks, or do a first-pass over the whole file to find all sample names that will appear.  This second option should be combined with some sort of indexing scheme since it forces us to preprocess the whole file anyway...)

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
        self.ref = args.reftag
        self.metadata = {'sourceformat': 'MAF'}
        self.metadata['labels'] = set()
        # if 'meta' in kwargs:
        #     self.metadata.update(kwargs['meta'])
        self.entrystart = 0
        # Check for Gzip and establish file object
        self.metadata['isgzip'] = (self.path.endswith(".gz") or
                                   args.get('isgzip', False))
        # READ MODE
        # try:
        print self.path
        filehandler = (self.metadata['isgzip']
                       and gzip.open(self.path, 'r')
                       or open(self.path, 'r'))
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
        # except:
        #     raise RuntimeError("file {} cannot be opened".format(self.path))

    def __iter__(self):
        """Simple entry iterator
           Returns (str(chrom), int(pos), list(allele entries))
        """
        if self.metadata['isgzip']:
            filehandler = gzip.open(self.path, 'rb')
        else:
            filehandler = open(self.path, 'r')
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
                        label, start, length, orientation, total_length, seq = (
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

    def reorient(self, seq, orientation, start, length):
        if orientation == '+':
            return seq, start
        elif orientation == '-':
            return ''.join(self.WatsonCrick[b] for b in seq[::-1]), start
        else:
            raise RuntimeError(
                "Error: orientation {} not recognized".format(orientation))
        return "", None


def main(arguments=sys.argv[1:]):
    """Main method for maf2mvf"""
    parser = argparse.ArgumentParser(description="""
    Converts Multiple Alignment Files to MVF file with some filters """)
    parser.add_argument("--maf", help="input MAF file")
    parser.add_argument("--out", help="output MVF file")
    parser.add_argument("--reftag", help="old reference tag")
    parser.add_argument("--mvfreflabel", default="REF",
                        help="new label for reference sample (default='REF')")
    parser.add_argument("--contigids", nargs='*',
                        help=("manually specify one or more contig ids"
                              " as ID:NAME"))
    parser.add_argument("--sampletags", nargs="*",
                        help="""one or more TAG:NEWLABEL or TAG, items,
                                if TAG found in sample label, replace with
                                NEW (or TAG if NEW not specified)
                                NEW and TAG must each be unique""")
    parser.add_argument("--linebuffer", type=int, default=100000,
                        help="number of lines to hold in read/write buffer")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(args=arguments)
    ## ESTABLISH MAF
    maf = MultiAlignFile(args)
    ## ESTABLISH MVF
    mvf = MultiVariantFile(args.out, 'write', overwrite=args.overwrite)
    # PROCESS CONTIG INFO
    # contigs = dict.fromkeys((sorted([x for x in maf.meta['name_index']
    #                             if x.find(args.reftag) > -1])), {})
    # print(contigs)
    # maxcontigid = 0
    # newids = set([])
    # if args.contigids:
    #     for cid, cname in (x.split(':') for x in args.contigids):
    #         for tempid in contigs:
    #             if cname in contigs[tempid]['label']:
    #                 try:
    #                     cid = int(cid)
    #                 except ValueError:
    #                     pass
    #                 mvf.metadata['contigs'][cid] = contigs[tempid].copy()
    #                 del contigs[tempid]
    #                 newids.update([cid])
    #                 break
    #     for cid in newids:
    #         try:
    #             maxcontigid = max([maxcontigid, int(cid) + 1])
    #         except ValueError:
    #             continue
    # tempids = set(contigs.keys()) - newids
    # for tempid, newid in zip(
    #         tempids, xrange(maxcontigid, maxcontigid + len(tempids))):
    #    # mvf.metadata['contigs'][newid] = maf.meta['contigs'][tempid]
    #         pass
    # contig_translate = dict([(mvf.metadata['contigs'][x]['label'], x)
    #                          for x in mvf.metadata['contigs']])
    # PROCESS SAMPLE INFO
    contig_translate = {1: 1}
    samplelabels = [s.split(':')[0] for s in args.sampletags]
    samplelabels.remove(args.reftag)
    samplelabels.insert(0, args.reftag)
    # if args.samplereplace:
    #     newsample = [':' in tuple(x) and x.split(':') or tuple([x,x])
    #                  for x in args.samplereplace]
    mvf.metadata['labels'] = samplelabels[:]
    for i, label in enumerate(samplelabels):
        mvf.metadata['samples'][i] = {'label': label}
    mvf.metadata['ncol'] = len(mvf.metadata['labels'])
    mvf.metadata['sourceformat'] = maf.metadata['sourceformat']
    ## WRITE MVF HEADER
    mvf.write_data(mvf.get_header())
    mvfentries = []
    nentry = 0
    for pos, length, msa in maf:
        for s in samplelabels:
            if s not in msa:
                msa[s] = '-'*length
        msa['contig'] = 1
        for i in range(length):
            mvf_alleles = encode_mvfstring(
                ''.join(msa[s][i].strip() for s in samplelabels))
            if mvf_alleles:
                mvfentries.append(
                    (contig_translate.get(msa['contig']),
                     pos+i, (mvf_alleles,)))
                nentry += 1
                if nentry == args.linebuffer:
                    mvf.write_entries(mvfentries, encoded=True)
                    mvfentries = []
                    nentry = 0
    if mvfentries:
        mvf.write_entries(mvfentries)

    return ''


if __name__ == "__main__":
    main()
