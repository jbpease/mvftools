#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
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

import os
import gzip
import re
import sys
from math import log10
from pylib.mvfbase import encode_mvfstring, MultiVariantFile, is_int
from pylib.mvfbiolib import MvfBioLib

MLIB = MvfBioLib()

RE_CONTIG_NAME = re.compile("ID=(.*?),")
RE_CONTIG_LENGTH = re.compile("length=(.*?)>")


class VariantCallFile():
    """Variant Call Format Handler
    Object Structure:
        path = file path (converted to absolute path)
        metadata = dict of associated data including (but not limited to):
            contigs: dict[id] = dict(metadata)
            samples: dict[index] = dict(sample_info)
    """

    def __init__(self, path, indexcontigs=True):
        if not path:
            raise IOError(path, " path not found for VCF file")
        self.path = os.path.abspath(path)
        self.metadata = {'contigs': {}, 'samples': []}
        if path.endswith(".gz"):
            filehandler = gzip.open(self.path, 'rt')
        else:
            filehandler = open(self.path, 'rt')
        header_lines = []
        self.entrystart = filehandler.tell()
        line = filehandler.readline()
        while line:
            if not line.startswith('#'):
                break
            header_lines.append(line.rstrip())
            self.entrystart = filehandler.tell() - 1
            line = filehandler.readline()
        self._process_header(header_lines)
        if not self.metadata['contigs'] and indexcontigs:
            self._index_contigs()

    def _process_header(self, headerlines):
        """Process VCF header information
            Arguments:
                headerlines: list of VCF header strings
        """
        tempid = 0
        for line in headerlines:
            if line.startswith('##fileformat'):
                self.metadata['sourceformat'] = line[line.find('=') + 1:]
            if line.startswith("##contig"):
                contig_name = re.findall(RE_CONTIG_NAME, line)
                contig_name = (contig_name[0] if contig_name else
                               'contig{}'.format(tempid))
                contig_length = re.findall(RE_CONTIG_LENGTH, line)
                contig_length = (int(contig_length[0])
                                 if contig_length
                                 else 0)
                self.metadata['contigs'][tempid] = {
                    'id': str(tempid),
                    'label': contig_name, 'length': contig_length}
                tempid += 1
            elif line.startswith("#CHROM"):
                self.metadata['samples'] = [
                    elem[elem.rfind('/') + 1:] if '/' in elem else elem
                    for elem in line.split()[9:]]
                self.metadata['ncol'] = len(self.metadata['samples'])
        return ''

    def _index_contigs(self, fieldsep="\t"):
        """Go through entries and index contig names
            Arguments:
                fielsep: separator for fields in VCF file
        """
        tempid = 0
        contigndx = {}
        if self.path.endswith('.gz'):
            filehandler = gzip.open(self.path, 'rt')
        else:
            filehandler = open(self.path, 'rt')
        filehandler.seek(self.entrystart)
        for vcfline in filehandler:
            try:
                arr = vcfline.rstrip().split(fieldsep)
                contig = arr[0]
                coord = int(arr[1])
                if contig not in contigndx:
                    self.metadata['contigs'][tempid] = {
                        'label': contig, 'length': coord}
                    contigndx[contig] = tempid
                    tempid += 1
                else:
                    self.metadata['contigs'][contigndx[contig]]['length'] = (
                        max(self.metadata['contigs'][
                            contigndx[contig]]['length'],
                            coord))
            except Exception as exception:
                continue
        return ''

    def iterentries(self, args):
        """Iterate VCF data entries
            Arguments:
                args: dict passthrough object from toplevel argparse options
        """
        if self.path.endswith('.gz'):
            filehandler = gzip.open(self.path, 'rt')
        else:
            filehandler = open(self.path, 'rt')
        nline = 0
        linebuffer = []
        for line in filehandler:
            nline += 1
            linebuffer.append(line)
            if nline == args.line_buffer:
                for xline in linebuffer:
                    if args.verbose is True:
                        print(xline)
                    try:
                        vcfrecord = self._parse_entry(xline, **vars(args))
                    except:
                        print("ERROR ON LINE:", line)
                        print(sys.exc_info()[0])
                        raise Exception
                    # vcfrecord is either (0, DATA) or (1, ERROR_MESSAGE)
                    if vcfrecord[0] == 1:  # 1 indicates error, 0=pass
                        if args.verbose is True:
                            print(vcfrecord)
                        continue
                    yield vcfrecord[1]
                linebuffer = []
                nline = 0
        for line in linebuffer:
            vcfrecord = self._parse_entry(line, **vars(args))
            if vcfrecord[0] == 1:  # 1=error, 0=pass
                continue
            yield vcfrecord[1]
            linebuffer = []
            nline = 0
        filehandler.close()

    def _parse_entry(self, vcfline, **kwargs):
        """Reads and processes a VCF multi-sample record
            Arguments:
                vcfile: string from VCF data entry
                **kwargs: passthrough dict of arguments
        """
        record = {}
        # Check if Depth is Zero
        if 'DP=0' in vcfline:
            return (1, "DP=0")
        # Check for indels
        if "INDEL" in vcfline:
            return (1, "INDEL")
        arr = vcfline.rstrip().split(kwargs.get('fieldsep', '\t'))
        # Check for too few fields
        if len(arr) < 9:
            return (1, "VCF FIELDS < 9")
        # Check for indels by altnerative allele length
        if len(arr[3]) > 1:
            return (1, "ALTERNATE ALLELE LENGTH > 1")
        # Check for no depth where N followed by *
        if arr[3] == 'N' and arr[4] == '*':
            return (1, "N and *, NO DEPTH")
        # Check for a heterozygous deletion call
        if "*" in arr[4]:
            return (1, "HETEROZYGOUS INDEL")
        record['contig'] = arr[0]
        record['coord'] = int(arr[1])
        record['alleles'] = [arr[3]]
        if arr[4] != '.':
            for altbase in arr[4].split(','):
                if len(altbase) > 1:
                    return (1, "AN ALTBASE HAS LEN>1")
                record['alleles'].append(altbase)
        record['tagindex'] = {}
        tags = arr[8].split(':')
        record['samples'] = []
        for elem in arr[9:]:
            record['samples'].append(dict(
                zip(tags, elem.split(':'))))
        if record['alleles'][0] in 'NnXxBbDdHhVv':
            record['genotypes'] = ['X']
            record['qscores'] = ['@']
        else:
            record['genotypes'] = [record['alleles'][0]]
            record['qscores'] = ['h']
        # ADDED CASE FOR WHEN DP IS ENCODED IN INFO FOR SINGLE_COLUMN INSTEAD OF COLUMN 9
        if len(record['samples']) == 1:
            if "DP=" in arr[7] and 'DP' not in arr[8]:
                xdepth = re.findall(r'DP=([0-9]*)', arr[7])
                if xdepth:
                    record['samples'][0]['DP'] = xdepth[0]
        for j in range(len(record['samples'])):
            (allele, quality, _) = (
                self._call_allele(record['samples'][j],
                                  record['alleles'],
                                  **kwargs))
            if kwargs.get("out_flavor") in ("dnaqual", 'dnaqual-indel'):
                record['qscores'].append(chr(min(quality, 40) + 64))
            record['genotypes'].append(allele)
        if kwargs.get("alleles_from"):
            info = dict(field.split('=') for field in arr[7].split(';'))
            record['genotypes'].extend(
                info.get(label, '-') for label in kwargs.get("alleles_from"))
        # Returns zero as non-error code, 0th element is 1 if error
        return (0, record)

    def _call_allele(self, sample, alleles, **kwargs):
        """Determine the allele from a VCF entry
            Arguments:
                sample: list of entry string elements
                alleles: alleles for all samples across lines
                indices: locations in each sample string of information fields
                **kwargs: passthrough arguments
            Returns (allele, quality, depth)
        """
        # phased = False
        if '|' in sample.get('GT', ''):
            # phased = True
            sample['GT'] = sample['GT'].replace('|', '/')
        if list(sample.values())[0][0] == '.':
            return ('-', 0, 0)
        sample_depth = int(sample.get('DP', -1))
        if sample_depth == 0:
            return ('-', 0, 0)
        if -1 < sample_depth < kwargs.get("mask_depth", 1):
            return ('X', -1, sample_depth)
        # Fixed sites
        if kwargs['ploidy'] == 2:
            if sample.get('GT', '') in ('.|.', './.'):
                return ('X', 0, 0)
        elif kwargs['ploidy'] == 4:
            if sample.get('GT', '') in ('.|.|.|.', './././.'):
                return ('X', 0, 0)
        elif kwargs['ploidy'] == 6:
            if sample.get('GT', '') in ('.|.|.|.|.|.', './././././.'):
                return ('X', 0, 0)
        if all(sample.get(x, -1) in (-1, '.')
               for x in ('PL', 'GL', 'GQ', 'GP')):
            quality = -1
            allele = sample.get('GT', 'X')
            if '/' in allele:
                allele = [int(x) for x in allele.split('/')]
                allele = ''.join(list(set(alleles[x] for x in allele)))
                allele = MLIB.joinbases[allele]
            elif '|' in allele:
                allele = [int(x) for x in allele.split('|')]
                allele = ''.join(list(set(alleles[x] for x in allele)))
                allele = MLIB.joinbases[allele]
            else:
                allele = 'X'
        # Invariant sites
        elif len(alleles) == 1:
            allele = alleles[0]
            if sample.get("GQ", -1) == -1:
                quality = -1
            else:
                quality = int(sample['GQ'])
        # Variant site
        elif (sample.get('PL', -1) in (-1, '.')) and (
                sample['GT'] in ('0/0', '1/1')):
            allele = (alleles[0] if sample['GT'] == '0/0' else alleles[1])
            quality = (-1 if sample.get("GQ", -1) == -1 else sample['GQ'])
        elif len(alleles) <= 4:
            if sample.get('PL', -1) == -1 and sample.get('GL', -1) != -1:
                plvalues = [float(x) if x != '.' else -1
                            for x in sample['GL'].split(',')]
            else:
                plvalues = [float(x) if x != '.' else -1
                            for x in sample['PL'].split(',')]
            if all(0 <= x <= 1 for x in plvalues) and sum(plvalues) == 1:
                plvalues = [x == 0 and 1 or x != 1 and
                            int(-10 * log10(x)) or 0 for x in plvalues]
            maxpl = max(plvalues) if 0 not in plvalues else 0
            imaxpl = (-1 if plvalues.count(maxpl) != 1 else
                      plvalues.index(maxpl))
            if imaxpl == -1:
                allele = "X"
            elif kwargs['ploidy'] > 2:
                if kwargs['ploidy'] == 4:
                    alleles = ''.join(list(set(
                        alleles[x] for x in MLIB.vcf_gtcodes_tetra[imaxpl])
                                           - set("-")))
                elif kwargs['ploidy'] == 6:
                    alleles = ''.join(list(set(
                        alleles[x] for x in MLIB.vcf_gtcodes_hex[imaxpl])
                                           - set("-")))
                else:
                    raise RuntimeError("Ploidy is not 2, 4 or 6")
                if alleles == "":
                    allele = "-"
                else:
                    allele = MLIB.joinbasespoly[alleles]
            else:
                alleles = ''.join(list(set(
                    alleles[x] for x in MLIB.vcf_gtcodes[imaxpl])
                                       - set("-")))
                if alleles == "":
                    allele = "-"
                else:
                    allele = MLIB.joinbases[alleles]
            quality = sample['GQ'] if sample.get('GQ', -1) != -1 else -1
        else:
            # Fail-safe check (you should never see a ! in the MVF output)
            allele = '!'
            quality = -1
        quality = int(float(quality)) if quality != '.' else 60
        if -1 < quality < kwargs.get("mask_qual", 10):
            return ('X', quality, sample_depth)
        if (-1 < quality < kwargs.get("low_qual", 20) or
                (-1 < sample_depth < kwargs.get("low_depth", 3))):
            allele = allele.lower()
        if kwargs['ploidy'] > 2:
            if allele in 'NnXx':
                allele = 'X'
        else:
            if allele in 'NnBbDdHhVvXx':
                allele = 'X'
        return (allele, quality, sample_depth)


def vcf2mvf(args=None):
    """Main method for vcf2mvf"""
    sepchars = dict([("TAB", "\t"), ("SPACE", " "), ("DBLSPACE", "  "),
                     ("COMMA", ","), ("MIXED", None)])
    args.fieldsep = sepchars[args.field_sep]
    # ESTABLISH VCF
    args.qprint("Opening input VCF: {}".format(args.vcf))
    vcf = VariantCallFile(args.vcf, indexcontigs=(not args.no_autoindex))
    # ESTABLISH MVF
    args.qprint("Establishing output MVF: {}".format(args.out))
    mvf = MultiVariantFile(args.out, 'write', overwrite=args.overwrite)
    mvf.notes.append(args.command_string)
    mvf.metadata['mvfversion'] = args.versionx
    # PROCESS CONTIG INFO
    args.qprint("Processing VCF headers.")
    vcfcontigs = vcf.metadata['contigs'].copy()
    args.qprint("{} contigs found.".format(len(vcfcontigs)))
    contig_translate = {}
    if args.contig_ids:
        for cid, cvcf, cmvf in (x.split(';') for x in args.contig_ids):
            try:
                cid = int(cid)
            except ValueError:
                pass
            assert cvcf in [vcfcontigs[x]['label'] for x in vcfcontigs]
            for vid in vcfcontigs:
                if vcfcontigs[vid]['label'] == cvcf:
                    contig_translate[cvcf] = [cid, cmvf]
                    if cid in mvf.metadata['contigs']:
                        raise RuntimeError(
                            'Contig id {} is not unique'.format(cid))
                    mvf.metadata['contigs'][cid] = vcfcontigs[vid].copy()
                    if cmvf in mvf.get_contig_labels():
                        raise RuntimeError(
                            'Contig label {} is not unique'.format(cmvf))
                    mvf.metadata['contigs'][cid]['label'] = cmvf[:]
    mvf.reset_max_contig()
    mvf.max_contig_index -= 1
    args.qprint("Processing contigs.")
    static_contig_ids = list(mvf.get_contig_ids())
    for vcid in vcfcontigs:
        vlabel = vcfcontigs[vcid]['label']
        if vlabel not in static_contig_ids:
            newindex = mvf.get_next_contig_index()
            if ((is_int(vlabel) or len(vlabel) < 3) and
                    vlabel not in static_contig_ids):
                newid = vlabel[:]
            else:
                newid = str(newindex)
            mvf.contig_indices.append(newindex)
            mvf.contig_ids.append(newid)
            mvf.contig_data[newindex] = vcfcontigs[vcid].copy()
            static_contig_ids.append(newid)
            contig_translate[vlabel] = [newindex, vlabel]
    mvf.reset_max_contig()
    new_contigs = [(x, mvf.contig_data[x]['label'])
                   for x in mvf.contig_indices]
    if args.skip_contig_label_check is False:
        args.qprint("Checking contigs for label/id overlap errors.")
        xids = [x[0] for x in new_contigs]
        xlabels = [x[1] for x in new_contigs]
        xintersect = set(xids).intersection(xlabels)
        if xintersect:
            for i, (newid, newlabel) in enumerate(new_contigs):
                if i % 100 == 0:
                    args.qprint("{} contigs processed".format(i))
                if newid in xlabels[:i] or newid in xlabels[i+1:]:
                    # if newid in xlabels:
                    # if xlabels.index(newid) != i:
                    raise RuntimeError("Error contig id {} is the same as"
                                       " the label for another contig"
                                       " ({})".format(
                                           newid, xlabels.index(newid)))
                if newlabel in xids[:i] or newlabel in xids[i+1:]:
                    # if newlabel in xids:
                    # if xids.index(newlabel) != i:
                    raise RuntimeError("Error contig label {} is the same"
                                       "as the id for another contig"
                                       "({})".format(
                                           newlabel, xids.index(newlabel)))
    # PROCESS SAMPLE INFO
    args.qprint("Processing samples.")
    samplelabels = [args.ref_label] + vcf.metadata['samples'][:]
    if args.alleles_from:
        args.alleles_from = args.alleles_from.split(':')
        samplelabels += args.alleles_from
    if args.sample_replace:
        newsample = [x.split(':') if ':' in tuple(x) else tuple([x, x])
                     for x in args.sample_replace]
        unmatched = list(enumerate(samplelabels))
        for old, new in newsample:
            labelmatched = False
            for j, (i, name) in enumerate(unmatched):
                if old in name:
                    samplelabels[i] = new
                    labelmatched = j
                    break
            if labelmatched is not False:
                del unmatched[labelmatched]
    mvf.sample_indices = list(range(len(samplelabels)))
    mvf.sample_ids = samplelabels[:]
    for i, label in enumerate(samplelabels):
        mvf.sample_data[i] = {'id': label}
    mvf.metadata['ncol'] = len(mvf.sample_ids)
    mvf.max_sample_index = len(mvf.sample_ids)
    mvf.metadata['sourceformat'] = vcf.metadata['sourceformat']
    # WRITE MVF HEADER
    mvf.write_data(mvf.get_header())
    mvfentries = []
    nentry = 0
    args.qprint("Processing VCF entries.")
    for vcfrecord in vcf.iterentries(args):
        mvfstring = ''.join(vcfrecord['genotypes'])
        if args.filter_nonref_empty is True:
            if all(x in 'Xx-?' for x in mvfstring[1:]):
                continue
        mvf_alleles = encode_mvfstring(mvfstring)
        if args.out_flavor in ('dnaqual',):
            qual_alleles = encode_mvfstring(''.join(vcfrecord['qscores']))
        if mvf_alleles:
            mvfentries.append(
                (contig_translate.get(vcfrecord['contig'])[0],
                 vcfrecord['coord'],
                 ((mvf_alleles, qual_alleles) if
                  args.out_flavor in ('dnaqual',) else
                  (mvf_alleles,))))
            nentry += 1
            if nentry == args.line_buffer:
                mvf.write_entries(mvfentries, encoded=True)
                mvfentries = []
                nentry = 0
    if mvfentries:
        mvf.write_entries(mvfentries)
        mvfentries = []
    return ''
