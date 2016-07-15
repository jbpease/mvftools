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

geno2MVF: Simple Genotype Format to MVF conversion program
@author: James B. Pease
@author: Ben K. Rosenzweig

version: 2015-02-01 - First Public Release
version: 2015-09-04 - Style cleanup
@version: 2015-12-31 - Headers and cleanup, changed to clarify this is not a GATK format

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
import gzip
import os
from mvfbase import encode_mvfstring, MultiVariantFile


class GenoFile(object):
    """Processor for a Simple Genotype File"""

    def __init__(self, path, entrystart=None, indexcontigs=True):
        if not path:
            raise IOError(path, " path not found for .geno file")
        self.path = os.path.abspath(path)
        self.metadata = {'contigs': {}, 'samples': []}
        self.metadata['sourceformat'] = '.geno'
        self.entrystart = entrystart or 0
        if path.endswith(".gz"):
            filehandler = gzip.open(self.path, 'rb')
        else:
            filehandler = open(self.path, 'rb')
        headerline = filehandler.readline()
        self.entrystart = filehandler.tell() - 1
        self.metadata['samples'] = headerline.split()[2:]
        if indexcontigs:
            self._index_contigs()

    def _index_contigs(self, fieldsep=" "):
        """Automatically index contig names"""
        tempid = 0
        contigx = {}
        if self.path.endswith('.gz'):
            filehandler = gzip.open(self.path, 'rb')
        else:
            filehandler = open(self.path, 'rb')
        filehandler.seek(self.entrystart)
        for line in filehandler:
            try:
                arr = line.rstrip().split(fieldsep)
                contig = arr[0]
                coord = int(arr[1])
                if contig not in contigx:
                    self.metadata['contigs'][tempid] = {
                        'label': contig, 'length': coord}
                    contigx[contig] = tempid
                    tempid += 1
                else:
                    self.metadata['contigs'][contigx[contig]]['length'] = max(
                        self.metadata['contigs'][contigx[contig]]['length'],
                        coord)
            except:
                continue
        return ''

    def iterentries(self, args):
        """Iterate entries from a .geno file
            Argruments:
                args: passthrough dict of arguments
        """
        if self.path.endswith('.gz'):
            filehandler = gzip.open(self.path, 'rb')
        else:
            filehandler = open(self.path, 'rb')
        nline = 0
        filehandler.seek(self.entrystart)
        linebuffer = []
        for line in filehandler:
            nline += 1
            linebuffer.append(line)
            if nline == args['linebuffer']:
                for line in linebuffer:
                    record = self._parse_entry(line, **args)
                    if record == -1:
                        continue
                    yield record
                linebuffer = []
                nline = 0
        for line in linebuffer:
            record = self._parse_entry(line, **args)
            if record == -1:
                continue
            yield record
            linebuffer = []
            nline = 0
        filehandler.close()

    def _parse_entry(self, entry, **kwargs):

        """Reads and processes a geno multi-sample record
            Arguments:
                entry: single line from geno file
                args: passthrough dict of arguments
        """
        record = {}
        arr = entry.rstrip().split(kwargs.get('fieldsep', None))
        if len(arr) < 3:
            return -1
        record['contig'] = arr[0]
        record['coord'] = int(arr[1])
        record['genotypes'] = [x == 'N' and 'X' or x for x in arr[2:]]
        return record


def main(arguments=sys.argv[1:]):
    """Main method for geno2mvf"""
    parser = argparse.ArgumentParser(description="""
    Converts Simple Genotype Format to MVF file with some filters """)
    parser.add_argument("--geno", help="input .geno file", required=True)
    parser.add_argument("--out", help="output MVF file", required=True)
    parser.add_argument("--contigids", nargs='*',
                        help=("manually specify one or more contig ids"
                              " as ID:NAME"))
    parser.add_argument("--samplereplace", nargs="*",
                        help="""one or more TAG:NEWLABEL or TAG, items,
                                if TAG found in sample label, replace with
                                NEW (or TAG if NEW not specified)
                                NEW and TAG must each be unique""")
    parser.add_argument("--reflabel", default="REF",
                        help="""label of the reference sample
                                (default is first entry)""")
    parser.add_argument("--no_autoindex", action="store_true",
                        help="do not automatically index contigs")
    parser.add_argument("--fieldsep", default="SPACE",
                        choices=['TAB', 'SPACE', 'DBLSPACE', 'COMMA', 'MIXED'],
                        help="""entry field separator (default='SPACE')""")
    parser.add_argument("--linebuffer", type=int, default=100000,
                        help="number of lines to hold in read/write buffer")
    parser.add_argument("--overwrite", action="store_true",
                        help="USE WITH CAUTION: force overwrite of outputs")
    parser.add_argument("--quiet", action="store_true",
                        help="suppress progress meter")
    parser.add_argument("-v", "--version", action="store_true",
                        help="display version information")
    args = parser.parse_args(args=arguments)
    if args.version:
        print("Version 2015-12-31")
        sys.exit()
    sepchars = dict([("TAB", "\t"), ("SPACE", " "), ("DBLSPACE", "  "),
                     ("COMMA", ","), ("MIXED", None)])
    args.fieldsep = sepchars[args.fieldsep]
    # ESTABLISH GENO
    geno = GenoFile(args.geno, indexcontigs=(not args.no_autoindex))
    # ESTABLISH MVF
    mvf = MultiVariantFile(args.out, 'write', overwrite=args.overwrite)
    # PROCESS CONTIG INFO
    contigs = geno.metadata['contigs'].copy()
    maxcontigid = 0
    newids = set([])
    if args.contigids:
        for cid, cname in (x.split(':') for x in args.contigids):
            for tempid in contigs:
                if cname in contigs[tempid]['label']:
                    try:
                        cid = int(cid)
                    except ValueError:
                        pass
                    mvf.metadata['contigs'][cid] = contigs[tempid].copy()
                    del contigs[tempid]
                    newids.update([cid])
                    break
        for cid in newids:
            try:
                maxcontigid = max([maxcontigid, int(cid) + 1])
            except ValueError:
                continue
    tempids = set(contigs.keys()) - newids
    for tempid, newid in sorted(zip(
            tempids, range(maxcontigid, maxcontigid + len(tempids)))):
        mvf.metadata['contigs'][newid] = geno.metadata['contigs'][tempid]
    contig_translate = dict([(mvf.metadata['contigs'][x]['label'], x)
                             for x in mvf.metadata['contigs']])
    # PROCESS SAMPLE INFO
    samplelabels = geno.metadata['samples'][:]
    if args.samplereplace:
        newsample = [':' in tuple(x) and x.split(':') or tuple([x, x])
                     for x in args.samplereplace]
        unmatched = [x for x in enumerate(samplelabels)]
        for old, new in newsample:
            labelmatched = False
            for j, (i, name) in enumerate(unmatched):
                if old in name:
                    samplelabels[i] = new
                    labelmatched = j
                    break
            if labelmatched is not False:
                del unmatched[labelmatched]
    mvf.metadata['labels'] = samplelabels[:]
    for i, label in enumerate(samplelabels):
        mvf.metadata['samples'][i] = {'label': label}
    mvf.metadata['ncol'] = len(mvf.metadata['labels'])
    mvf.metadata['sourceformat'] = geno.metadata['sourceformat']
    # WRITE MVF HEADER
    mvf.write_data(mvf.get_header())
    mvfentries = []
    nentry = 0
    for record in geno.iterentries(vars(args)):
        mvf_alleles = encode_mvfstring(''.join(record['genotypes']))
        if mvf_alleles:
            mvfentries.append(
                (contig_translate.get(record['contig'], record['contig']),
                 record['coord'], mvf_alleles))
            nentry += 1
            if nentry == args.linebuffer:
                mvf.write_entries(mvfentries, encoded=True)
                mvfentries = []
                nentry = 0
    if mvfentries:
        mvf.write_entries(mvfentries)
        mvfentries = []
    return ''


if __name__ == "__main__":
    main()
