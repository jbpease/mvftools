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

from itertools import combinations
from pylib.mvfbase import MultiVariantFile, OutputFile, Counter
from pylib.mvfbiolib import MvfBioLib
from pylib.mvfpaml import paml_branchsite
from pylib.mvftranslate import parse_gff_analysis

MLIB = MvfBioLib()


def procarg_speciesgroups(xarg, mvf):
    groups = {}
    for elem in xarg:
        elem = elem.split(':')
        groups[elem[0]] = mvf.get_sample_indices(labels=elem[1].split(','))
    xarg = groups.copy()
    for specgroup in groups:
        ngroup = 0
        for allelegroup in xarg.values():
            if set(allelegroup) & set(groups[specgroup]):
                ngroup += 1
                if ngroup > 1:
                    raise RuntimeError(specgroup, "split across 2+ groups")
    return xarg


def procarg_allelegroups(xarg, mvf):
    groups = {}
    for elem in xarg:
        elem = elem.split(':')
        groups[elem[0]] = mvf.get_sample_indices(labels=elem[1].split(','))
    xarg = groups.copy()
    for grp0, grp1 in combinations(groups, 2):
        if set(groups[grp0]) & set(groups[grp1]):
            raise RuntimeError("Groups contain same element",
                               set(groups[grp0]) & set(groups[grp1]))
    return xarg


def parse_dndsfile(dndsfile):
    entries = {}
    with open(dndsfile, 'r') as infile:
        for line in infile:
            arr = line.rstrip().split()
            entries[arr[0]] = (float(arr[3]), float(arr[4]))
    return entries


def hapgroup(group):
    """Takes codon group and splits it into haplotypes
    """
    new_group = set([])
    for codon in group:
        nhet = sum([int(x in 'RYWKMS') for x in codon])
        if nhet == 0:
            new_group.update((codon,))
        elif nhet > 1:
            continue
        else:
            for frame in (0, 1, 2):
                if codon[frame] in 'RYWKMS':
                    new_group.update((
                        ''.join([
                            ''.join(codon[0:frame]),
                            MLIB.splitbases[codon[frame]][0],
                            ''.join(codon[frame+1:3])]),
                        ''.join([
                            ''.join(codon[0:frame]),
                            MLIB.splitbases[codon[frame]][1],
                            ''.join(codon[frame+1:3])]),))
                    break
    return new_group


def calc_group_unique_allele_window(args):
    """Count the number of and relative rate of uniquely held alleles
       spatially along chromosomes (i.e. Lineage-specific rates)"""
    data = {}
    mvf = MultiVariantFile(args.mvf, 'read')
    annotations = {}
    coordinates = {}
    labels = mvf.get_sample_labels()[:]
    ncol = len(labels)
    current_contig = None
    current_position = 0
    counts = Counter()
    totals = Counter()
    args.start_contig = (
        args.start_contig if args.start_contig is not None else 0)
    args.end_contig = (
        args.end_contig if args.end_contig is not None else 100000000000)

    if args.output_align is True:
        outputalign = []
    if args.gff is not None:
        annotations, coordinates = (parse_gff_analysis(args.gff))
    if args.allele_groups is not None:
        args.allele_groups = procarg_allelegroups(
            args.allele_groups, mvf)
    if args.species_groups is not None:
        args.species_groups = procarg_speciesgroups(
            args.species_groups, mvf)
    fieldtags = [
        'likelihood', 'bgdnds0', 'bgdnds1', 'bgdnds2a', 'bgdnds2b',
        'fgdnds0', 'fgdnds1', 'fgdnds2a', 'fgdnds2b', 'dndstree',
        'errorstate']
    if args.branch_lrt is not None:
        with open(args.branch_lrt, 'w') as branchlrt:
            genealign = []
            branchlrt.write("\t".join(
                ['contig', 'ntaxa', 'alignlength', 'lrtscore'] +
                ["null.{}".format(x) for x in fieldtags] +
                ["test.{}".format(x) for x in fieldtags] +
                ['tree']) + "\n")
    groups = args.allele_groups.values()
    if args.species_groups is not None:
        speciesgroups = args.species_groups.values()
    allsets = set([])
    for group in groups:
        allsets.update(group)
    allsets = list(sorted(allsets))
    speciesnames = args.species_groups.keys()
    speciesrev = {}
    for species in args.species_groups:
        speciesrev.update([
            (x, species) for x in args.species_groups[species]])
    if args.mincoverage is not None:
        if args.mincoverage < len(groups) * 2:
            raise RuntimeError("""
                Error: GroupUniqueAlleleWindow:
                --mincoverage cannot be lower than the twice the number
                of specified groups in --allele-groups
                """)
    for contig, pos, allelesets in mvf:
        if not current_contig:
            current_contig = contig[:]
        if contig != current_contig or (
                args.windowsize > 0 and
                pos > current_position + args.windowsize):
            xkey = (current_contig, current_position,)
            data[xkey] = counts.copy()
            data[xkey].update([
                ('contig', (args.use_labels and
                            mvf.get_contig_label(current_contig))),
                ('position', current_position),
                ('nonsynyonymous_changes',
                 counts.get('nonsynonymous_changes', 0) or 0),
                ('synyonymous_changes',
                 counts.get('synonymous_changes', 0) or 0)
                ])
            data[xkey].update([
                ('ns_ratio', (float(data[xkey].get(
                    'nonsynonymous_changes', 0)) / (
                        data[xkey].get('synonymous_changes', 1.0)))),
                ('annotation',
                 annotations.get(data[xkey]['contig'], '.')),
                ('coordinates',
                 coordinates.get(data[xkey]['contig'], '.'))
                ])
            if genealign:
                if (args.end_contig >= int(current_contig)) and (
                        args.start_contig <= int(current_contig)):
                    (pamlnull, pamltest, tree) = paml_branchsite(
                        genealign, labels[:],
                        species=speciesnames,
                        speciesrev=speciesrev,
                        codemlpath=args.codeml_path,
                        raxmlpath=args.raxml_path,
                        pamltmp=args.paml_tmp,
                        target=args.target,
                        targetspec=args.num_target_species,
                        allsampletrees=args.all_sample_trees,
                        outgroup=args.outgroup)
                    lrtscore = -1
                    if (pamlnull.get('likelihood', -1) != -1 and
                            pamltest.get('likelihood', -1) != -1):
                        lrtscore = 2 * (pamltest['likelihood'] -
                                        pamlnull['likelihood'])
                    with open(args.branch_lrt, 'a') as branchlrt:
                        branchlrt.write("\t".join([str(x) for x in [
                            data[xkey]['contig'],
                            len(genealign),
                            len(genealign[0]) * 3,
                            lrtscore] + [
                                pamlnull.get(y, -1) for y in fieldtags] + [
                                    pamltest.get(y, -1) for y in fieldtags] + [
                                        str(tree).rstrip()]]) + "\n")
            genealign = None
            totals.add('genes_total')
            if counts.get('total_codons', 0) > 0:
                totals.add('genes_tested')
            if counts.get('total_nsyn_codons', 0) > 0:
                totals.add('genes_with_nsyn')
            if contig != current_contig:
                current_contig = contig[:]
                current_position = 0
            elif args.windowsize > 0:
                current_position += args.windowsize
            counts = Counter()
        proteins = allelesets[0]
        codons = allelesets[1:4]
        if len(proteins) == 1 and all(len(x) == 1 for x in codons):
            if proteins == '*' or ''.join(codons) in MLIB.stop_codons:
                continue
            counts.add('total_codons')
            totals.add('total_codons')
            if args.output_align is True:
                if not outputalign:
                    outputalign = [[''.join(codons)]
                                   for x in range(mvf.metadata['ncol'])]
                else:
                    for ialign, xalign in enumerate(outputalign):
                        xalign.append(''.join(codons))
            if args.branch_lrt is not None:
                if not genealign:
                    genealign = [[''.join(codons)]
                                 for x in range(ncol)]
                else:
                    for ialign in range(len(genealign)):
                        genealign[ialign].append(''.join(codons))
            continue
        if len(proteins) > 1:
            if allelesets[0][1] == '+':
                continue
        proteins = mvf.decode(proteins)
        if args.mincoverage is not None:
            if sum([int(x not in 'X-') for x in proteins]) < (
                    args.mincoverage):
                continue
        species_groups = [[proteins[i] for i in x
                           if proteins[i] not in '-X']
                          for x in speciesgroups]
        if any(len(x) == 0 for x in species_groups):
            continue
        xcodons = [mvf.decode(x) for x in codons]
        codons = [''.join(x) for x in zip(*xcodons)]
        if any(codons[x] in MLIB.stop_codons for x in allsets):
            continue
        if any(any(x != species_groups[0][0] for x in y)
               for y in species_groups):
            totals.add('total_nsyn_codons')
            counts.add('total_nsyn_codons')
        totals.add('total_codons')
        totals.add('tested_codons')
        counts.add('total_codons')
        totals.add('variable_codons',
                   val=int(sum([int(len(set(x) - set('X-')) > 1)
                                for x in xcodons]) > 0))
        if args.output_align is not None:
            if not outputalign:
                outputalign = [[x] for x in codons]
            else:
                for ialign in range(len(outputalign)):
                    outputalign[ialign].append(codons[ialign])
        if args.branch_lrt is not None:
            if not genealign:
                genealign = [[x] for x in codons]
            else:
                for ialign in range(len(codons)):
                    genealign[ialign].append(codons[ialign])
        nonsyn_change = False
        synon_change = False
        codon_groups = [
            set([codons[i] for i in x if '-' not in codons[i] and
                 'X' not in codons[i]])
            for x in groups]
        protein_groups = None
        for i in range(len(codon_groups)):
            if any(base in codon for base in 'RYWKMS'
                   for codon in codon_groups[i]):
                codon_groups[i] = hapgroup(codon_groups[i])
        if all(grp1.isdisjoint(grp0) for grp0, grp1 in
               combinations(codon_groups, 2)):
            protein_groups = [set(
                [MLIB.codon_tables['full'][''.join(x)]
                 for x in codon_groups[i]])
                              for i in range(len(codon_groups))]
            if all(grp1.isdisjoint(grp0) for grp0, grp1 in
                   combinations(protein_groups, 2)):
                nonsyn_change = True
            elif all(grp1 == grp0 for grp0, grp1 in combinations(
                    protein_groups, 2)):
                synon_change = True
        if nonsyn_change:
            print('NON', contig, pos, allelesets, codon_groups,
                  protein_groups, groups, mvf.get_contig_label(contig))
            counts.add('nonsynonymous_changes')
            totals.add('nonsynonymous_changes')
        elif synon_change:
            print('SYN', contig, pos, allelesets, codon_groups,
                  protein_groups, groups, mvf.get_contig_label(contig))
            counts.add('synonymous_changes')
            totals.add('synonymous_changes')
    args.totals = totals
    # WRITE OUTPUT
    headers = ["contig", "position", "nonsynonymous_changes",
               "synonymous_changes", "ns_ratio",
               "nonsynonymous_total", "synonymous_total",
               "pvalue",
               "total_codons", "annotation", "coordinates"]
    if args.windowsize == -1:
        headers.remove('position')
    if args.chi_test is None:
        headers.remove('pvalue')
    outfile = OutputFile(path=args.out, headers=headers)
    sorted_entries = sorted([
        (data[k]['ns_ratio'], k)
        for k in data
        if data[k].get('nonsynonymous_changes', 0) > 0],
                            reverse=True)
    for _, k in sorted_entries:
        outfile.write_entry(data[k])
    with open(args.out + '.total', 'w') as totalfile:
        for entry in args.totals.iter_sorted():
            totalfile.write(entry)
    if args.output_align is not None:
        with open(args.output_align, 'w') as alignfile:
            alignfile.write(
                "\n".join([">{}\n{}".format(mvf.metadata['labels'][i],
                                            ''.join(outputalign[i]))
                           for i in range(len(outputalign))]))
    return ''
