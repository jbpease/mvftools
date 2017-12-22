# -*- coding: utf-8 -*-
"""
This program analyzes a DNA MVF alignment using the modules specified below,
use the --morehelp option for additional module information.

MVFtools: Multisample Variant Format Toolkit
James B. Pease and Ben K. Rosenzweig
http://www.github.org/jbpease/mvftools

If you use this software please cite:
Pease JB and BK Rosenzweig. 2015.
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

from random import randint
from itertools import combinations
from pylib.mvfbase import MultiVariantFile, OutputFile, zerodiv, same_window
from pylib.mvfbiolib import MvfBioLib

MLIB = MvfBioLib()


def check_mincoverage(val, alleles):
    """Check if site has minimum coverage"""
    if val is not None:
        if sum([int(x not in 'Xx-') for x in alleles]) < val:
            return False
    return True


def pi_diversity(seq):
    """Calculate Pi sequence diversity"""
    seq = ''.join([MLIB.splitbases[x] for x in seq])
    base_count = [seq.count(x) for x in 'ATGC']
    total = float(sum(base_count))
    if not total:
        return 'nodata'
    if any(x == total for x in base_count):
        return 0.0
    else:
        return total / (total - 1) * (1 - sum([(x / total)**2
                                               for x in base_count]))


def calc_sample_coverage(args):
    """Counts the total number of non-gap/ambiguous characters for
      each sample.
      """
    mvf = MultiVariantFile(args.mvf, 'read')
    labels = mvf.get_sample_labels()
    data = {}
    for contig, _, allelesets in mvf.iterentries(
            contigs=args.contigs, subset=args.samples,
            decode=True):
        if contig not in data:
            data[contig] = dict.fromkeys(labels, 0)
            data[contig]['contig'] = contig
        for j, base in enumerate(allelesets[0]):
            data[contig][labels[j]] += int(
                base not in 'Xx-')
        labels = labels
    outfile = OutputFile(path=args.out,
                         headers=(["contig"] + labels))
    for contig in data:
        outfile.write_entry(data[contig])
    return ''


def calc_dstat_combinations(args):
    """Calculate genome-wide D-statstics for
       all possible trio combinations of samples
    """
    mvf = MultiVariantFile(args.mvf, 'read')
    data = {}
    nsamples = len(args.samples) - 1
    samplenames = (args.samples.split(",")
                   if args.samples is not None
                   else mvf.get_sample_labels())
    contigs = (args.contigs.split(",")
               if args.contigs is not None
               else mvf.get_contig_labels())
    for contig, _, allelesets in mvf:
        alleles = mvf.decode(allelesets[0])
        if alleles[-1] in 'X-':
            continue
        for i in range(nsamples - 2):
            if alleles[i] not in 'ATGC':
                continue
            for j in range(i + 1, nsamples - 1):
                if alleles[j] not in 'ATGC':
                    continue
                for k in range(j + 1, nsamples):
                    if alleles[k] not in 'ATGC':
                        continue
                    subset = [alleles[x] for x in [i, j, k, -1]]
                    if subset[-1] not in subset[:3]:
                        continue
                    if len(set(subset)) != 2:
                        continue
                    val = (1 * (subset[j] == subset[-1]) +
                           2 * (subset[k] == subset[-1]))
                    trio = (i, j, k)
                    if trio not in data:
                        data[trio] = {}
                    if contig not in data[trio]:
                        data[trio][contig] = [0, 0, 0]
                    data[trio][contig][val - 1] += 1
    # WRITE OUTPUT
    headers = ['sample0', 'sample1', 'sample2']
    for xcontig in contigs:
        headers.extend(['{}:abba'.format(xcontig),
                        '{}:baba'.format(xcontig),
                        '{}:bbaa'.format(xcontig),
                        '{}:D'.format(xcontig)
                        ])
    outfile = OutputFile(path=args.out, headers=headers)
    for i, j, k in combinations(range(nsamples), 3):
        trio = tuple([i, j, k])
        if trio not in data:
            continue
        entry = dict([('sample{}'.format(i),
                       samplenames[x][0])
                      for i, x in enumerate(trio)])
        for contig in contigs:
            if contig not in data[trio]:
                entry.update(dict().fromkeys([
                    '{}:abba'.format(contig),
                    '{}:baba'.format(contig),
                    '{}:bbaa'.format(contig),
                    '{}:D'.format(contig)],
                                             '0'))
            else:
                [abba, baba, bbaa] = data[trio][contig]
                if abba > baba and abba > bbaa:

                    dstat = zerodiv(baba - bbaa, baba + bbaa)
                elif baba > bbaa and baba > abba:
                    dstat = zerodiv(abba - bbaa, abba + bbaa)
                else:
                    dstat = zerodiv(abba - baba, abba + baba)
                entry.update([('{}:abba'.format(contig), abba),
                              ('{}:baba'.format(contig), baba),
                              ('{}:bbaa'.format(contig), bbaa),
                              ('{}:D'.format(contig), dstat)
                              ])
        outfile.write_entry(entry)
    return ''


def calc_pattern_count(args):
    """Count biallelic patterns spatially along
       chromosomes (e.g,, for use in DFOIL or Dstats
       http://www.github.com/jbpease/dfoil)
    """
    mvf = MultiVariantFile(args.mvf, 'read')
    labels = mvf.get_sample_labels()
    data = {}
    current_contig = None
    current_position = 0
    sitepatterns = {}
    samples = ([labels.index(x) for x in args.samples]
               if args.samples is not None
               else range(len(labels)))
    nsamples = len(samples)
    for contig, pos, allelesets in mvf:
        # Check Minimum Site Coverage
        if check_mincoverage(args.mincoverage, allelesets[0]) is False:
            continue
        # Establish first contig
        if current_contig is None:
            current_contig = contig[:]
            while pos > current_position + args.windowsize - 1:
                current_position += args.windowsize
        # Check if windows are specified.
        if not same_window((current_contig, current_position),
                           (contig, pos), args.windowsize):
            data[(current_contig, current_position)] = dict([
                ('contig', current_contig),
                ('position', current_position)])
            data[(current_contig, current_position)].update(
                sitepatterns)
            sitepatterns = {}
            if contig != current_contig:
                current_position = 0
                current_contig = contig[:]
            else:
                current_position += (0 if args.windowsize == -1
                                     else args.windowsize)
        if len(allelesets[0]) == 1:
            if allelesets[0] in 'ATGC':
                pattern = 'A' * nsamples
            else:
                continue
        elif allelesets[0][1] == '+':
            continue
        else:
            alleles = mvf.decode(allelesets[0])
            alleles = [alleles[x] for x in samples]
            if any(x in alleles for x in 'X-RYKMWS'):
                continue
            if len(set(alleles)) > 2:
                continue
            pattern = ''.join(['A' if x == alleles[-1] else 'B'
                               for x in alleles[:-1]]) + 'A'
        sitepatterns[pattern] = sitepatterns.get(pattern, 0) + 1
    if sitepatterns:
        data[(current_contig, current_position)] = dict([
            ('contig', current_contig),
            ('position', current_position)])
        data[(current_contig, current_position)].update(
            sitepatterns)
    # WRITE OUTPUT
    headers = ['contig', 'position']
    headers.extend(
        [MLIB.abpattern(x, nsamples)
         for x in range(0, 2 ** nsamples, 2)])
    outfile = OutputFile(path=args.out, headers=headers)
    sorted_entries = sorted([(data[k]['contig'],
                              data[k]['position'], k)
                             for k in data])
    for _, _, k in sorted_entries:
        outfile.write_entry(data[k])
    # WRITE LIST OUTPUT
    sorted_entries = sorted([(data[k]['contig'],
                              data[k]['position'], k)
                             for k in data])
    total_counts = {}
    for contig, pos, k in sorted_entries:
        outfilepath = "{}-{}-{}.counts.list".format(
            args.out, contig, pos)
        with open(outfilepath, 'w') as outfile:
            outfile.write("pattern,count\n")
            for pattern, pcount in sorted(data[k].items()):
                if pattern in ['contig', 'position']:
                    continue
                outfile.write("{},{}\n".format(pattern, pcount))
                total_counts[pattern] = (
                    total_counts.get(pattern, 0) + pcount)
    outfilepath = "{}-TOTAL.counts.list".format(args.out)
    with open(outfilepath, 'w') as outfile:
        outfile.write("pattern,count\n")
        for pattern, pcount in sorted(total_counts.items()):
            if pattern in ['contig', 'position']:
                continue
            outfile.write("{},{}\n".format(pattern, pcount))
    return ''


def calc_character_count(args):
    """Count the number of and relative rate of certain bases
       spatially along chromosomes
    """
    mvf = MultiVariantFile(args.mvf, 'read')
    labels = mvf.get_sample_labels()
    data = {}
    current_contig = None
    current_position = 0
    match_counts = dict().fromkeys(labels, 0)
    total_counts = dict().fromkeys(labels, 0)
    all_match = 0
    all_total = 0
    data_in_buffer = 0
    for contig, pos, allelesets in mvf:
        # Check Minimum Site Coverage
        if check_mincoverage(args.mincoverage, allelesets[0]) is False:
            continue
        # Establish first contig
        if current_contig is None:
            current_contig = contig[:]
            while pos > current_position + args.windowsize - 1:
                current_position += args.windowsize
        # Check if windows are specified.
        if not same_window((current_contig, current_position),
                           (contig, pos), args.windowsize):
            data[(current_contig, current_position)] = {
                'contig': current_contig, 'position': current_position}
            for k in match_counts:
                data[(current_contig, current_position)].update([
                    (k + '.match', match_counts[k] + all_match),
                    (k + '.total', total_counts[k] + all_total),
                    (k + '.prop', (
                        (float(match_counts[k] + all_match) /
                         float(total_counts[k] + all_total)) if
                        total_counts[k] + all_total > 0 else 0))])
            if contig != current_contig:
                current_contig = contig[:]
                current_position = 0
            else:
                current_position += (0 if args.windowsize == -1
                                     else args.windowsize)
            match_counts = dict().fromkeys(labels, 0)
            total_counts = dict().fromkeys(labels, 0)
            all_total = 0
            all_match = 0
            data_in_buffer = 0
        else:
            alleles = allelesets[0]
            if len(alleles) == 1:
                if alleles in args.base_match:
                    all_match += 1
                if alleles in args.base_total:
                    all_total += 1
            else:
                alleles = mvf.decode(alleles)
                for i, base in enumerate(alleles):
                    if base in args.base_match:
                        match_counts[labels[i]] += 1
                    if base in args.base_total:
                        total_counts[labels[i]] += 1
            data_in_buffer = 1
    if data_in_buffer:
        data[(current_contig, current_position)] = {
            'contig': current_contig, 'position': current_position}
        for k in match_counts:
            data[(current_contig, current_position)].update([
                (k + '.match', match_counts[k] + all_match),
                (k + '.total', total_counts[k] + all_total),
                (k + '.prop', ((float(match_counts[k] + all_match) /
                                float(total_counts[k] + all_total)) if
                               total_counts[k] + all_total > 0 else 0))])
    # WRITE OUTPUT
    headers = ['contig', 'position']
    for label in labels:
        headers.extend([label + x for x in ('.match', '.total', '.prop')])
    outfile = OutputFile(path=args.out,
                         headers=headers)
    sorted_entries = sorted([(data[k]['contig'],
                              data[k]['position'], k)
                             for k in data])
    for _, _, k in sorted_entries:
        outfile.write_entry(data[k])
    return ''


def calc_pairwise_distances(args):
    """Count the pairwise nucleotide distance between
       combinations of samples in a window
    """
    mvf = MultiVariantFile(args.mvf, 'read')
    data = {}
    labels = mvf.get_sample_labels()
    current_contig = None
    current_position = 0
    data_in_buffer = False
    ncol = mvf.metadata['ncol']
    sample_pairs = [tuple(x) for x in combinations(range(ncol), 2)]
    base_matches = dict([(x, {}) for x in sample_pairs])
    all_match = {}
    for contig, pos, allelesets in mvf:
        # Check Minimum Site Coverage
        if check_mincoverage(args.mincoverage, allelesets[0]) is False:
            continue
        # Establish first contig
        if current_contig is None:
            current_contig = contig[:]
            while pos > current_position + args.windowsize - 1:
                current_position += args.windowsize
        # Check if windows are specified.
        if not same_window((current_contig, current_position),
                           (contig, pos), args.windowsize):
            data[(current_contig, current_position)] = {
                'contig': current_contig, 'position': current_position}
            if mvf.flavor == 'dna':
                all_diff, all_total = pairwise_distance_nuc(all_match)
            elif mvf.flavor == 'prot':
                all_diff, all_total = pairwise_distance_prot(all_match)
            for samplepair in base_matches:
                if mvf.flavor == 'dna':
                    ndiff, ntotal = pairwise_distance_nuc(
                        base_matches[samplepair])
                elif mvf.flavor == 'prot':
                    ndiff, ntotal = pairwise_distance_prot(
                        base_matches[samplepair])
                taxa = "{};{}".format(labels[samplepair[0]],
                                      labels[samplepair[1]])
                data[(current_contig, current_position)].update({
                    '{};ndiff'.format(taxa): ndiff + all_diff,
                    '{};ntotal'.format(taxa): ntotal + all_total,
                    '{};dist'.format(taxa): zerodiv(ndiff + all_diff,
                                                    ntotal + all_total)})
            if contig != current_contig:
                current_contig = contig[:]
                current_position = 0
                while pos > current_position + args.windowsize - 1:
                    current_position += args.windowsize
            else:
                current_position += args.windowsize
            base_matches = dict([(x, {}) for x in sample_pairs])
            all_match = {}
            data_in_buffer = False
        alleles = allelesets[0]
        if len(alleles) == 1:
            all_match["{}{}".format(alleles, alleles)] = (
                all_match.get("{}{}".format(alleles, alleles),
                              0) + 1)
            data_in_buffer = True
            continue
        if alleles[1] == '+':
            if 'X' in alleles or '-' in alleles:
                continue
            samplepair = (0, int(alleles[3:]))
            basepair = "{}{}".format(alleles[0], alleles[2])
            base_matches[samplepair][basepair] = (
                base_matches[samplepair].get(basepair, 0) + 1)
            data_in_buffer = True
            continue
        alleles = mvf.decode(alleles)
        valid_positions = [i for i, x in enumerate(alleles)
                           if x not in 'X-']
        for i, j in combinations(valid_positions, 2):
            samplepair = (i, j)
            basepair = "{}{}".format(alleles[i], alleles[j])
            base_matches[samplepair][basepair] = (
                base_matches[samplepair].get(basepair, 0) + 1)
        data_in_buffer = True
    if data_in_buffer is True:
        # Check whether, windows, contigs, or total
        if args.windowsize == 0:
            current_contig = 'TOTAL'
            current_position = 0
        elif args.windowsize == -1:
            current_position = 0
        data[(current_contig, current_position)] = {
            'contig': current_contig, 'position': current_position}
        if mvf.flavor == 'dna':
            all_diff, all_total = pairwise_distance_nuc(all_match)
        elif mvf.flavor == 'prot':
            all_diff, all_total = pairwise_distance_prot(all_match)
        for samplepair in base_matches:
            if mvf.flavor == 'dna':
                ndiff, ntotal = pairwise_distance_nuc(base_matches[samplepair])
            elif mvf.flavor == 'prot':
                ndiff, ntotal = pairwise_distance_prot(
                    base_matches[samplepair])
            taxa = "{};{}".format(labels[samplepair[0]],
                                  labels[samplepair[1]])
            data[(current_contig, current_position)].update({
                '{};ndiff'.format(taxa): ndiff + all_diff,
                '{};ntotal'.format(taxa): ntotal + all_total,
                '{};dist'.format(taxa): zerodiv(ndiff + all_diff,
                                                ntotal + all_total)})
    headers = ['contig', 'position']
    for samplepair in sample_pairs:
        headers.extend(['{};{};{}'.format(
            labels[samplepair[0]],
            labels[samplepair[1]],
            x) for x in ('ndiff', 'ntotal', 'dist')])
    outfile = OutputFile(path=args.out, headers=headers)
    sorted_entries = sorted([(
        data[k]['contig'], data[k]['position'], k)
                             for k in data])
    for _, _, k in sorted_entries:
        outfile.write_entry(data[k])
    return ''


def pairwise_distance_nuc(basepairs, strict=False):
    """Calculates pairwise distances between two sequences
        strict = only use ATGC if True,
        choose random heterozygous base if False.
    """
    total = 0
    diff = 0
    for pairbases, paircount in basepairs.items():
        base0, base1 = pairbases[0], pairbases[1]
        if base0 not in MLIB.validchars['dna+ambig'] or (
                base1 not in MLIB.validchars['dna+ambig']):
            continue
        if base0 in MLIB.validchars['dnaambig2']:
            if strict is True:
                continue
            base0 = MLIB.splitbases[base0][randint(0, 1)]
        if base1 in MLIB.validchars['dnaambig2']:
            if strict is True:
                continue
            base1 = MLIB.splitbases[base1][randint(0, 1)]
        if base0 != base1:
            diff += paircount
        total += paircount
    return diff, total


def pairwise_distance_prot(basepairs):
    """Calculates pairwise distances between two sequences
        strict = only use ACDEFGHIKLMNPQRSTVXY if True,
        choose random heterozygous base if False.
    """
    total = 0
    diff = 0
    for pairbases, paircount in basepairs.items():
        base0, base1 = pairbases[0], pairbases[1]
        if base0 not in MLIB.validchars['amino'] or (
                base1 not in MLIB.validchars['amino']):
            continue
        if base0 != base1:
            diff += paircount
        total += paircount
    return diff, total
