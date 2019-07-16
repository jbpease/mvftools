#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This program is used to convert a FASTA file into MVF format.
"""

import re
import os
from random import randint
from pylib.mvfbase import encode_mvfstring, is_int
from pylib.mvfbase import MultiVariantFile, fasta_iter


_LICENSE = """
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


def parse_regions_arg(regionfilepath, contigs):
    """Parses the regions into coordinates"""
    fmt_regions = []
    region_max_coord = {}
    if regionfilepath is None:
        fmt_regions = [(x, None, None, None) for x in contigs]
        region_max_coord = dict.fromkeys(contigs, None)
    else:
        with open(regionfilepath) as regfile:
            for line in regfile:
                entry = line.rstrip().split(',')
                if len(entry) > 4 or len(entry) < 1 or len(entry[0]) == 0:
                    print("malformed entry ({}), ignoring...".format(entry))
                    continue
                contig = ''
                if entry[0] in contigs:
                    contig = entry[0][:]
                elif is_int(entry[0]):
                    if int(entry[0]) in contigs:
                        contig = int(entry[0])
                if contig == '':
                    for cid in contigs:
                        if contigs[cid]['label'] == entry[0]:
                            contig = cid
                assert contig in contigs
                if len(entry) == 1:
                    fmt_regions.append((contig, None, None, '+'))
                elif len(entry) == 2:
                    assert int(entry[1]) > 0
                    fmt_regions.append((contig, int(entry[1]), None, '+'))
                else:
                    assert int(entry[1]) > 0
                    assert int(entry[2]) > 0
                    assert int(entry[2]) > int(entry[1])
                    fmt_regions.append((
                        contig, int(entry[1]), int(entry[2]),
                        "+" if int(entry[2]) > int(entry[1]) else "-"))
    fmt_regions.sort()
    for contigid, _, maxcoord, _ in fmt_regions:
        if contigid not in region_max_coord:
            region_max_coord[contigid] = (
                None if maxcoord is None else maxcoord + 0)
        elif maxcoord is None:
            pass
        elif maxcoord > region_max_coord[contigid]:
            region_max_coord[contigid] = maxcoord + 0
    regionlabel = ','.join(["{}{}{}{}{}".format(
        contigs[x[0]]['label'],
        "" if (x[1] == -1 or x[1] == 0 or x[1] is None) else (
            ":{}".format(x[1])),
        "" if (x[2] == -1 or x[2] == 0 or x[2] is None) else '..',
        "" if (x[2] == -1 or x[2] == 0 or x[2] is None) else x[2],
        "" if (x[2] == -1 or x[2] == 0 or x[2] is None) else
        "({})".format(x[3])) for x in fmt_regions])
    return fmt_regions, region_max_coord, regionlabel


def mvf2fasta(args):
    """Main method"""
    mvf = MultiVariantFile(args.mvf, 'read')
    if (mvf.flavor in ("dna", "rna") and args.output_data == "prot") or (
            mvf.flavor == "prot" and args.output_data in ("dna", "rna")):
        raise RuntimeError(
            "--output-data {} incompatiable with '{}' flavor mvf".format(
                args.output_data, mvf.flavor))
    regions, max_region_coord, regionlabel = parse_regions_arg(
        args.regions, mvf.metadata['contigs'])
    sample_labels = mvf.get_sample_labels()
    if args.sample_indices is not None:
        sample_indices = [int(x) for x in
                          args.sample_indices[0].split(",")]
    elif args.sample_labels is not None:
        sample_indices = mvf.get_sample_indices(
            labels=args.sample_labels[0].split(","))
    else:
        sample_indices = mvf.get_sample_indices()
    skipcontig = ''
    tmp_files = dict((fn, open("{}-{}.tmp".format(
        fn, randint(1000000, 9999999)), 'w+', args.buffer))
        for fn in sample_labels)
    labelwritten = dict.fromkeys(sample_labels, False)
    for contig, pos, allelesets in mvf.iterentries(
            contigs=[x for x in max_region_coord], decode=True):
        if contig == skipcontig:
            continue
        if (contig not in max_region_coord) or (
                max_region_coord[contig] is not None and
                pos > max_region_coord[contig]):
            skipcontig = contig[:]
            continue
        inregion = False
        for rcontig, rstart, rstop, _ in regions:
            if contig == rcontig:
                if rstart is None or pos >= rstart:
                    if rstop is None or pos <= rstop:
                        inregion = True
                        break
        if inregion is False:
            continue
        for col, label in zip(sample_indices, sample_labels):
            if not labelwritten[label]:
                if args.label_type == 'long':
                    xlabel = "{} region={}".format(label, regionlabel)
                elif args.label_type == 'short':
                    xlabel = "{}".format(label)
                tmp_files[label].write(">{}\n".format(xlabel))
                labelwritten[label] = True
            if mvf.flavor == 'dna':
                tmp_files[label].write(
                    "N" if allelesets[0][col] == 'X'
                    else allelesets[0][col])
            elif mvf.flavor in ('codon', 'prot') and (
                    args.output_data == 'prot'):
                tmp_files[label].write(allelesets[0][col])
            elif mvf.flavor == 'codon' and args.output_data == 'dna':
                codon = ["N" if allelesets[x][col] == 'X' else
                         allelesets[x][col] for x in (1, 2, 3)]
                tmp_files[label].write(''.join(codon))
    with open(args.out, 'w') as outfile:
        for filehandler in tmp_files.values():
            filehandler.seek(0, 0)
            buff = filehandler.read(args.buffer)
            while len(buff):
                outfile.write(buff)
                buff = filehandler.read(args.buffer)
            outfile.write("\n")
            filehandler.close()
            os.remove(os.path.join(args.temp_dir, filehandler.name))
    return ''


def fasta2mvf(args):
    """Main method"""
    sepchars = dict([("PIPE", "\\|"), ("TAB", "\\t"),
                     ("SPACE", "\\s"), ("DBLSPACE", "\\s\\s"),
                     ("COMMA", "\\,"), ("NONE", None),
                     ("AT", "\\@"), ('UNDER', "\\_"), ("DBLUNDER", "\\_\\_")])
    if args.field_sep is None:
        args.field_sep = ''
    else:
        args.field_sep = re.compile("[{}]".format(
            ''.join([sepchars[x] for x in args.field_sep])))
    if args.manual_coord:
        assert len(args.manual_coord) == len(args.fasta)
        args.manual_coord = [
            (x.split(':')[0], int(x.split(":")[1].split('..')[0]),
             int(x.split(':')[1].split('..')[1]))
            for x in args.manual_coord]
    mvf = MultiVariantFile(args.out, 'write', overwrite=args.overwrite)
    fasta = {}
    current_contig = 0
    fsamples = []
    fcontigs = []
    for ifasta, fastapath in enumerate(args.fasta):
        print("Processing {}".format(fastapath))
        for header, seq in fasta_iter(fastapath):
            if args.field_sep is None:
                header = header[:]
            if args.field_sep != '' and args.field_sep is not None:
                header = [str(x) for x in re.split(args.field_sep, header)]
            if args.contig_by_file is True:
                contig = os.path.basename(fastapath[:])
                if args.sample_field is None:
                    sample = header[:]
                else:
                    sample = header[args.sample_field]
            elif (len(header) < max(args.contig_field if
                                    args.contig_field
                                    is not None else 0,
                                    args.sample_field if
                                    args.sample_field is not None else
                                    0) or
                  args.contig_field is None or args.sample_field is None):
                contig = "UNK{}".format(current_contig)
                sample = header[:]
            elif args.manual_coord:
                contig = args.manual_coord[ifasta][0]
            else:
                contig = header[args.contig_field]
                sample = header[args.sample_field]
            if contig not in fcontigs:
                fcontigs.append(contig)
                fasta[contig] = {}
            if sample not in fsamples:
                fsamples.append(sample)
            fasta[contig][sample] = (len(seq), seq)
    reflabel = None
    if args.ref_label:
        for i, samplename in enumerate(fsamples):
            if args.ref_label in samplename:
                reflabel = i
                break
    if reflabel:
        newref = fsamples.pop(i)
        fsamples = [newref] + fsamples
    for i, contig in enumerate(fcontigs):
        mvf.metadata['contigs'][i] = {
            'label': contig,
            'length': max([fasta[contig][x][0] for x in fasta[contig]])}
    mvf.metadata['labels'] = fsamples[:]
    for i, label in enumerate(fsamples[:]):
        mvf.metadata['samples'][i] = {'label': label}
    mvf.metadata['ncol'] = len(mvf.metadata['labels'])
    mvf.metadata['sourceformat'] = 'fasta'
    mvf.flavor = args.flavor
    # WRITE MVF HEADER
    mvf.write_data(mvf.get_header())
    mvfentries = []
    nentry = 0
    mvf_alleles = {}
    for cind, contig in enumerate(fcontigs):
        for pos in range(mvf.metadata['contigs'][cind]['length']):
            mvf_alleles = encode_mvfstring(
                ''.join(samp not in fasta[contig] and '-' or
                        pos >= fasta[contig][samp][0] and '-' or
                        fasta[contig][samp][1][pos]
                        for samp in fsamples))
            if mvf_alleles:
                if args.flavor == 'dna':
                    mvf_alleles = ''.join(["X" if x in 'NXOBDHVnxobdhv' else x
                                           for x in mvf_alleles])
                mvfentries.append(
                    (cind, pos+1, (mvf_alleles,)))
                nentry += 1
                if nentry == args.write_buffer:
                    mvf.write_entries(mvfentries, encoded=True)
                    mvfentries = []
                    nentry = 0
    if mvfentries:
        mvf.write_entries(mvfentries)
        mvfentries = []
    return ''


def mvf2phy(args):
    """Main method"""
    mvf = MultiVariantFile(args.mvf, 'read')
    if (mvf.flavor in ("dna", "rna") and args.output_data == "prot") or (
            mvf.flavor == "prot" and args.output_data in ("dna", "rna")):
        raise RuntimeError(
            "--outdput-data {} incompatiable with '{}' flavor mvf".format(
                args.output_data, mvf.flavor))
    max_region_coord = dict.fromkeys(mvf.metadata['contigs'], None)
    if args.region is not None:
        _, max_region_coord, _ = parse_regions_arg(
            args.region, mvf.metadata['contigs'])
    sample_labels = mvf.get_sample_labels()
    if args.sample_indices is not None:
        sample_indices = [int(x) for x in
                          args.sample_indices[0].split(",")]
    elif args.sample_labels is not None:
        sample_indices = mvf.get_sample_indices(
            labels=args.sample_labels[0].split(","))
    else:
        sample_indices = mvf.get_sample_indices()
    skipcontig = ''
    tmp_files = dict((fn, open("{}-{}.tmp".format(
        fn, randint(1000000, 9999999)), 'w+', args.buffer))
        for fn in sample_labels)
    labelwritten = dict.fromkeys(sample_labels, False)
    curcontigname = None
    curcontigstart = 1
    curcontigend = 1
    if args.partition is True:
        partprefix = "PROT" if args.output_data == "prot" else "DNA"
        partitionfile = open("{}.part".format(args.out), 'w')
    for contig, _, allelesets in mvf.iterentries(
            contigs=(mvf.metadata['contigs'] if args.region is None else
                     [x for x in max_region_coord]), decode=True):
        if contig == skipcontig:
            continue
        if contig not in max_region_coord:
            skipcontig = contig[:]
            continue
        if curcontigname is None:
            curcontigname = contig[:]
        elif contig != curcontigname:
            if args.partition is True:
                if curcontigend > curcontigstart:
                    partitionfile.write("{}, {} = {}-{}\n".format(
                        partprefix, mvf.get_contig_labels(
                            ids=curcontigname),
                        curcontigstart, curcontigend - 1))
            curcontigname = contig[:]
            # reset start as one position after end of last
            curcontigstart = curcontigend
            curcontigend = curcontigend + 1
        for col, label in zip(sample_indices, sample_labels):
            if not labelwritten[label]:
                if args.label_type == 'long':
                    tmp_files[label].write("{}{}".format(
                        label[:100], " "*(100 - len(label[:100]))))
                elif args.label_type == 'short':
                    tmp_files[label].write("{}{}".format(
                        label[:20], " "*(20 - len(label[:20]))))
                labelwritten[label] = True
            if mvf.flavor == 'dna':
                tmp_files[label].write(
                    allelesets[0][col] == 'X' and
                    'N' or allelesets[0][col])
                if label == sample_labels[0]:
                    curcontigend += 1
            elif ((mvf.flavor == 'codon' and args.output_data == 'prot') or (
                    mvf.flavor == 'prot')):
                tmp_files[label].write(allelesets[0][col])
                if label == sample_labels[0]:
                    curcontigend += 1
            elif mvf.flavor == 'codon':
                codon = ["N" if allelesets[x][col] == 'X' else
                         allelesets[x][col] for x in (1, 2, 3)]
                tmp_files[label].write(''.join(codon))
                if label == sample_labels[0]:
                    curcontigend += 3
    first_file = True
    totalseqlen = 0
    with open(args.out, 'w') as outfile:
        for filehandler in tmp_files.values():
            # read first file to establish sequence length for phylip header
            if first_file is True:
                filehandler.seek(0, 0)
                buff = filehandler.read(args.buffer)
                while buff != '':
                    if " " in buff:
                        totalseqlen += len(buff.strip().split(" ")[-1])
                    else:
                        totalseqlen += len(buff.strip())
                    buff = filehandler.read(args.buffer)
                outfile.write("{} {}\n".format(
                    len(sample_labels), totalseqlen))
                first_file = False
            filehandler.seek(0, 0)
            buff = filehandler.read(args.buffer)
            while buff != '':
                if first_file is True:
                    outfile.write("{} {}\n".format(
                        len(sample_labels), len(buff.split()[1])))
                    first_file = False
                outfile.write(buff)
                buff = filehandler.read(args.buffer)
            outfile.write("\n")
            filehandler.close()
            os.remove(os.path.join(args.temp_dir, filehandler.name))
    if args.partition is True:
        if curcontigend > curcontigstart:
            partitionfile.write("{},{},{},{}\n".format(
                partprefix, mvf.get_contig_labels(ids=curcontigname),
                curcontigstart, curcontigend - 1))
        partitionfile.close()
    return ''
