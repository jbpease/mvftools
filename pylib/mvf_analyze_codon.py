# -*- coding: utf-8 -*-
"""
MVFtools: Multisample Variant Format Toolkit
James B. Pease and Ben K. Rosenzweig
http://www.github.org/jbpease/mvftools


THIS IS A LEGACY MODULE.
"""

import os
import sys
import argparse
from random import randint
from itertools import combinations
from mvfbase import MultiVariantFile, AnalysisModule, OutputFile, Counter
from mvfbiolib import MvfBioLib  # HAPSPLIT, FULL_CODON_TABLE, AMBIGSTOPS
from mvfpaml import paml_pwcalc_dnds
from mvftranslate import parse_gff_legacy_annotate
MLIB = MvfBioLib()

_LICENSE = """
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



MODULENAMES = ("Coverage", "GroupUniqueAlleleWindow", "PiDiversityWindow",
               "PairwiseNS")


class PairwiseDNDS(AnalysisModule):
    """Count the number of and relative rate of uniquely held alleles
       spatially along chromosomes (i.e. Lineage-specific rates)"""

    def analyze(self, mvf):
        """Analyze Entries for GroupUniqueAlleleWindow Module
        """
        annotations = {}
        coordinates = {}
        if self.params['gff']:
            annotations, coordinates = parse_gff_legacy_annotate(
                self.params['gff'])
        self.params['labels'] = mvf.get_sample_labels()[:]
        ncol = len(self.params['labels'])
        current_contig = None
        current_position = 0
        counts = Counter()
        totals = Counter()
        if self.params['output_align']:
            outputalign = []
        fieldtags = ['likelihood', 'bgdnds0', 'bgdnds1', 'bgdnds2a',
                     'bgdnds2b', 'fgdnds0', 'fgdnds1', 'fgdnds2a', 'fgdnds2b',
                     'dndstree', 'errorstate']
        with open(self.params['branchlrt'], 'w') as branchlrt:
            genealign = []
            branchlrt.write("\t".join(
                ['contig', 'ntaxa', 'alignlength', 'lrtscore'] +
                ["null.{}".format(x) for x in fieldtags] +
                ["test.{}".format(x) for x in fieldtags] +
                ['tree']) + "\n")
        groups = self.params['allele_groups'].values()
        speciesgroups = self.params['speciesgroups'].values()
        allsets = set([])
        for group in groups:
            allsets.update(group)
        allsets = list(sorted(allsets))
        speciesrev = {}
        for species in self.params['speciesgroups']:
            speciesrev.update([(x, species)
                               for x in self.params['speciesgroups'][species]])
        if self.params['mincoverage']:
            if self.params['mincoverage'] < len(groups) * 2:
                raise RuntimeError("""
                    Error: GroupUniqueAlleleWindow:
                    --mincoverage cannot be lower than the twice the number
                    of specified groups in --allele-groups
                    """)
        for contig, pos, allelesets in mvf:
            if not current_contig:
                current_contig = contig[:]
            if contig != current_contig or (
                    self.params['windowsize'] != -1 and
                    pos > current_position + self.params['windowsize']):
                xkey = (current_contig, current_position,)
                self.data[xkey] = counts.copy()
                self.data[xkey].update([
                    ('contig', (self.params['uselabels'] and
                                mvf.get_contig_label(current_contig))),
                    ('position', current_position),
                    ('nonsynyonymous_changes',
                     counts.get('nonsynonymous_changes', 0) or 0),
                    ('synyonymous_changes',
                     counts.get('synonymous_changes', 0) or 0)
                    ])
                self.data[xkey].update([
                    ('ns_ratio', (float(self.data[xkey].get(
                        'nonsynonymous_changes', 0)) / (
                            self.data[xkey].get('synonymous_changes', 1.0)))),
                    ('annotation',
                     annotations.get(self.data[xkey]['contig'], '.')),
                    ('coordinates',
                     coordinates.get(self.data[xkey]['contig'], '.'))
                    ])
                if genealign:
                    if (self.params.get('endcontig', 1000000) >=
                            int(current_contig)) and (
                                self.params.get('startcontig', 0) <=
                                int(current_contig)):
                        # print(current_contig)
                        (dnval, dsval) = paml_pwcalc_dnds(genealign)
                        with open(self.params['branchlrt'],
                                  'a') as branchlrt:
                            branchlrt.write("\t".join([str(x) for x in [
                                self.data[xkey]['contig'],
                                len(genealign),
                                len(genealign[0]) * 3,
                                dnval, dsval]]) + "\n")
                genealign = None
                totals.add('genes_total')
                if counts.get('total_codons', 0) > 0:
                    totals.add('genes_tested')
                if counts.get('total_nsyn_codons', 0) > 0:
                    totals.add('genes_with_nsyn')
                if contig != current_contig:
                    current_contig = contig[:]
                    current_position = 0
                elif self.params['windowsize'] != -1:
                    current_position += self.params['windowsize']
                counts = Counter()
            proteins = allelesets[0]
            codons = allelesets[1:4]
            if len(proteins) == 1 and all(len(x) == 1 for x in codons):
                if proteins == '*' or ''.join(codons) in MLIB.stop_codons:
                    continue
                counts.add('total_codons')
                totals.add('total_codons')
                if self.params['output_align']:
                    if not outputalign:
                        outputalign = [[''.join(codons)]
                                       for x in range(mvf.metadata['ncol'])]
                    else:
                        for ialign, _ in enumerate(outputalign):
                            outputalign[ialign].append(''.join(codons))
                if self.params['branchlrt']:
                    if not genealign:
                        genealign = [[''.join(codons)]
                                     for x in range(ncol)]
                    else:
                        for ialign, _ in enumerate(genealign):
                            genealign[ialign].append(''.join(codons))
                continue
            if len(proteins) > 1:
                if allelesets[0][1] == '+':
                    continue
            proteins = mvf.decode(proteins)
            if self.params['mincoverage']:
                if sum([int(x not in 'X-') for x in proteins]) < (
                        self.params['mincoverage']):
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
            if self.params['output_align']:
                if not outputalign:
                    outputalign = [[x] for x in codons]
                else:
                    for ialign, _ in enumerate(outputalign):
                        outputalign[ialign].append(codons[ialign])
            if self.params['branchlrt']:
                if not genealign:
                    genealign = [[x] for x in codons]
                else:
                    for ialign, _ in enumerate(codons):
                        genealign[ialign].append(codons[ialign])
            nonsyn_change = False
            synon_change = False
            codon_groups = [
                set([codons[i] for i in x if '-' not in codons[i] and
                     'X' not in codons[i]])
                for x in groups]
            protein_groups = None
            for i, icodongrp in enumerate(codon_groups):
                if any(base in codon for base in 'RYWKMS'
                       for codon in icodongrp):
                    codon_groups[i] = hapgroup(icodongrp)
            if all(grp1.isdisjoint(grp0) for grp0, grp1 in
                   combinations(codon_groups, 2)):
                protein_groups = [
                    set([MLIB.codon_tables['full'][''.join(x)]
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
        self.params['totals'] = totals
        self.write()
        if self.params['output_align']:
            with open(self.params['output_align'], 'w') as alignfile:
                alignfile.write(
                    "\n".join([">{}\n{}".format(mvf.metadata['labels'][i],
                                                ''.join(outputalign[i]))
                               for i in range(len(outputalign))]))
        return ''

    def write(self):
        """Write Output"""
        headers = ["contig", "position", "nonsynonymous_changes",
                   "synonymous_changes", "ns_ratio",
                   "nonsynonymous_total", "synonymous_total",
                   "pvalue",
                   "total_codons", "annotation", "coordinates"]
        if self.params['windowsize'] == -1:
            headers.remove('position')
        if not self.params['chi_test']:
            headers.remove('pvalue')
        outfile = OutputFile(path=self.params['out'], headers=headers)
        sorted_entries = sorted([
            (self.data[k]['ns_ratio'], k)
            for k in self.data
            if self.data[k].get('nonsynonymous_changes', 0) > 0],
                                reverse=True)
        for _, k in sorted_entries:
            outfile.write_entry(self.data[k])
        with open(self.params['out'] + '.total', 'w') as totalfile:
            for entry in self.params['totals'].iter_sorted():
                totalfile.write(entry)
        return ''


class PiDiversityWindow(AnalysisModule):
    """Calculate amino acid pi for windows"""

    def pi_diversity(self, alleles):
        """Calculate Pi sequence diversity"""
        base_count = filter(lambda x: x > 0,
                            [alleles.count(y) for y in 'ACDEFGHIKLMNPQRSTVWY'])
        total = float(sum(base_count))
        if not total:
            return 'nodata'
        if any(x == total for x in base_count):
            return 0.0
        else:
            return total / (total - 1) * (1 - sum([(x / total)**2
                                                   for x in base_count]))

    def analyze(self, mvf):
        """Analyze Entries for GroupUniqueAlleleWindow Module"""
        current_contig = None
        current_position = 0
        pidiff = 0.
        nsites = 0
        subpival = 0.
        for contig, pos, allelesets in mvf:
            if not current_contig:
                current_contig = contig[:]
            if contig != current_contig or (
                    self.params['windowsize'] != -1 and
                    pos > current_position + self.params['windowsize']):
                self.data[(current_contig, current_position)] = {
                    'pi': pidiff / float(nsites + 1) if nsites > 0 else 0,
                    'pinumer': pidiff,
                    'nsites': nsites,
                    'contig': (mvf.get_contig_label(current_contig)
                               if self.params['uselabels'] is True
                               else current_contig),
                    'position': current_position}
                if contig != current_contig:
                    current_contig = contig[:]
                    current_position = 0
                elif self.params['windowsize'] != -1:
                    current_position += self.params['windowsize']
                pidiff = 0.
                nsites = 0
                subpival = 0.
            else:
                alleles = allelesets[0]
                if len(alleles) == 1:
                    nsites += 1
                    continue
                if alleles[1] == '+':
                    continue
                if len(alleles) == 2:
                    alleles = mvf.decode(alleles)
                elif alleles[2] == '+':
                    alleles = mvf.decode(alleles)
                if self.params['mincoverage']:
                    if sum([int(x not in 'X-') for x in alleles]) < (
                            self.params['mincoverage']):
                        continue
                if not set(alleles) - set('X-'):
                    continue
                subpival = self.pi_diversity(alleles)
                if subpival != 'nodata':
                    pidiff += subpival
                    nsites += 1
        self.write()
        return ''

    def write(self):
        """Write Output"""
        outfile = OutputFile(path=self.params['out'],
                             headers=('pi', 'pinumer', 'nsites', 'contig',
                                      'position'))
        sorted_entries = sorted([(self.data[k]['pi'], k)
                                 for k in self.data], reverse=True)
        for _, k in sorted_entries:
            outfile.write_entry(self.data[k])
        return ''


class PairwiseProtein(AnalysisModule):
    """Count the number of and relative rate of uniquely held alleles
       spatially along chromosomes"""

    def analyze(self, mvf):
        """Analyze Entries for GroupUniqueAlleleWindow Module"""
        current_contig = None
        current_position = 0
        site_count = 0
        var_count = 0
        total_count = 0
        groups = self.params['allele_groups'].values()
        for contig, pos, allelesets in mvf:
            if not current_contig:
                current_contig = contig[:]
            if contig != current_contig or (
                    self.params['windowsize'] != -1 and
                    pos > current_position + self.params['windowsize']):
                self.data[(current_contig, current_position)] = {
                    'nsites': site_count,
                    'nvar': var_count,
                    'ntotal': total_count,
                    'ratio': (var_count and
                              float(site_count) / float(var_count) or
                              0.),
                    'psites': (total_count and
                               float(site_count)/float(total_count) or
                               0.),
                    'contig': (self.params['uselabels'] and
                               mvf.get_contig_label(current_contig) or
                               current_contig),
                    'position': current_position,
                    }
                if contig != current_contig:
                    current_contig = contig[:]
                    current_position = 0
                elif self.params['windowsize'] != -1:
                    current_position += self.params['windowsize']
                site_count = 0
                var_count = 0
                total_count = 0
            else:
                alleles = allelesets[0]
                if len(alleles) in (1, 2):
                    total_count += 1
                    continue
                if alleles[1] == '+':
                    continue
                var_count += 1
                if alleles[2] == '+':
                    alleles = mvf.decode(alleles)
                if self.params['mincoverage']:
                    if sum([int(x not in 'X-') for x in alleles]) < (
                            self.params['mincoverage']):
                        continue
                allele_groups = [[alleles[i] for i in x
                                  if alleles[i] not in '-X']
                                 for x in groups]
                if not all((len(groups[i]) == 1 and len(x) == 1) or len(x) > 1
                           for i, x in enumerate(allele_groups)):
                    continue
                total_count += 1
                if any(any(x in grp1 for x in grp0) for grp0, grp1 in
                       combinations(allele_groups, 2)):
                    continue
                site_count += 1
        self.write()
        return ''

    def write(self):
        """Write Output"""
        outfile = OutputFile(path=self.params['out'],
                             headers=('ratio', 'psites', 'nsites', 'nvar',
                                      'ntotal', 'contig',
                                      'position'))
        sorted_entries = sorted([(self.data[k]['ratio'], k)
                                 for k in self.data
                                 if self.data[k]['nsites'] > 0],
                                reverse=True)
        for _, k in sorted_entries:
            outfile.write_entry(self.data[k])
        return ''


class PairwiseNS(AnalysisModule):
    """Count the number of synonymous
       and non-synonymous differences pairwise
       """

    def analyze(self, mvf):
        """Analyze Entries for PairwiseNS Module"""
        # invar = 0
        total_count = 0
        dnon = dict.fromkeys([tuple([i, j]) for (i, j) in combinations(
            range(mvf.metadata['ncol']), 2)], 0)
        dsyn = dict.fromkeys([tuple([i, j]) for (i, j) in combinations(
            range(mvf.metadata['ncol']), 2)], 0)
        # invalid = 0
        # ncol = mvf.metadata['ncol'] + 0
        for _, _, allelesets in mvf:
            if any('*' in x or '-' in x for x in allelesets):
                # invalid += 1
                continue
            elif any('X' in x for x in allelesets[1:3]):
                # invalid += 1
                continue
            total_count += 1
            codons = allelesets[1:3]
            if all(len(x) == 1 for x in codons):
                # invar += 1
                continue
            aminoacids = mvf.decode(allelesets[0])
            codons = [mvf.decode(allelesets[x]) for x in (1, 2, 3)]
            for i, j in combinations(range(len(aminoacids)), 2):
                if aminoacids[i] != aminoacids[j]:
                    dnon[(i, j)] += 1
                else:
                    synonchange = False
                    for pos in (0, 1, 2):
                        if (codons[pos][i] in 'ATGC' and
                                codons[pos][j] in 'ATGC'):
                            if codons[pos][i] != codons[pos][j]:
                                synonchange = True
                                break
                        elif (MLIB.splitbases[codons[pos][i]][randint(0, 1)] !=
                              MLIB.splitbases[codons[pos][j]][randint(0, 1)]):
                            synonchange = True
                            break
                    if synonchange:
                        dsyn[(i, j)] += 1
        labels = mvf.get_sample_labels()
        for (i, j) in dnon:
            dkey = labels[i] + ';' + labels[j]
            self.data[dkey] = {'label': dkey, 'Ka': dnon[(i, j)],
                               'Ks': dsyn[(i, j)],
                               'dN': float(dnon[(i, j)]) / float(total_count),
                               'dS': float(dsyn[(i, j)]) / float(total_count)}
        self.write()
        return ''

    def write(self):
        """Write Output"""
        outfile = OutputFile(path=self.params['out'],
                             headers=('label', 'Ka', 'Ks', 'dN', 'dS'))
        for k in self.data:
            outfile.write_entry(self.data[k])
        return ''


def modulehelp(modulenames=MODULENAMES):
    """Prints extra description of modules"""
    for modulename in modulenames:
        print("{}: {}".format(modulename, eval(modulename + ".__doc__")))
    return ''


def generate_argparser():
    """Generate Argparse"""
    parser = argparse.ArgumentParser(
        prog="mvf_analyze_codon.py",
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=_LICENSE)
    parser.add_argument("module", choices=MODULENAMES)
    parser.add_argument("-i", "--mvf", help='Input MVF file.',
                        type=os.path.abspath)
    parser.add_argument("-o", "--out", help="output file",
                        type=os.path.abspath)
    parser.add_argument("-c", "--contigs", nargs='*',
                        help="List of space-separated contig ids.")
    parser.add_argument("-s", "--samples", nargs='*',
                        help="List of space-separated sample names.")
    parser.add_argument("--allele-groups", "--allelegroups", nargs='*',
                        help="""GROUP1:LABEL,LABEL GROUP2:LABEL,LABEL
                                (GroupUniqueAlleleWindow)""")
    parser.add_argument("--species-groups", "--speciesgroups", nargs='*')
    parser.add_argument("-g", "--gff", help="GFF3 file for use in annotation")
    parser.add_argument("-x", "--chi-test", "--chitest",
                        help=("(GroupUniqueAlleleWindow,PairwiseDNDS)"
                              "Input two number values for expected "
                              "Nonsynonymous and Synonymous expected values. "
                              ))
    parser.add_argument("-m", "--mincoverage", type=int,
                        help=("Minimum number of samples with alleles needed "
                              "to use site for analysis."))
    parser.add_argument("--target", nargs="*",
                        help=("(GroupUniqueAlleleWindow) "
                              "Specify the taxa labels that define the "
                              "target lineage-specific branch to be tested."))
    parser.add_argument("--num-target-species", "--targetspec",
                        type=int, default=1,
                        help=("(GroupUniqueAlleleWindow) "
                              "Specify the minimum number of "
                              "taxa in the target set "
                              "that are required to conduct analysis"))
    parser.add_argument("--output-align", "--outputalign",
                        help=("(GroupUniqueAlleleWindow) "
                              "Output alignment to this file path in "
                              "phylip format."))
    parser.add_argument("-O", "--outgroup",
                        help=("(GroupUniqueAlleleWindow) Specify sample "
                              "name with which to root trees."))
    parser.add_argument("-w", "--windowsize", type=int, default=-1,
                        help="Window size in bp, use -1 for whole contig.")
    parser.add_argument("-L", "--uselabels", action="store_true",
                        help="Use contig labels instead of IDs in output.")
    parser.add_argument("-P", "--codemlpath", default="codeml",
                        type=os.path.abspath,
                        help="Full path for PAML codeml executable.")
    parser.add_argument("--raxml-path", "--raxmlpath", type=os.path.abspath,
                        default="raxml",
                        help="Full path to RAxML program executable.")
    parser.add_argument("--start-contig", "--startcontig", type=int, default=0,
                        help="Numerical ID for the starting contig.")
    parser.add_argument("--end-contig", "--endcontig", type=int, 
                        default=100000000,
                        help="Numerical id for the ending contig.")
    parser.add_argument("--pamltmp", default="pamltmp", type=os.path.abspath,
                        help="path for temporary folder for PAML output files")
    parser.add_argument("--all-sample-trees", "--allsampletrees",
                        action="store_true",
                        help=("(GroupUniqueAlleleWindow) Makes trees from "
                              "all samples instead of only the "
                              "most complete sequence from each species"))
    parser.add_argument("--morehelp", action="store_true",
                        help="Get additional information on modules.")
    parser.add_argument("--branchlrt", "--branch-lrt",
                        type=os.path.abspath,
                        help=("(GroupUniqueAlleleWindow) Specify the "
                              "output file for and turn on the "
                              "RAxML-PAML format LRT test scan for "
                              "selection on the target branch in addition "
                              "to the basic patterns scan"))
    parser.add_argument("--version", action="version",
                        version="2017-06-24",
                        help="display version information")
    return parser


def main(arguments=None):
    """Main method"""
    arguments = sys.argv[1:] if arguments is None else arguments
    parser = generate_argparser()
    args = parser.parse_args(args=arguments)
    # HELP MENU
    if args.morehelp:
        modulehelp(MODULENAMES)
        sys.exit()
    # ESTABLISH MVF
    mvf = MultiVariantFile(args.mvf, 'read')
    # Argument Pre-processing
    if args.allele_groups:
        groups = {}
        for elem in args.allele_groups:
            elem = elem.split(':')
            groups[elem[0]] = mvf.get_sample_indices(labels=elem[1].split(','))
        args.allele_groups = groups.copy()
        for grp0, grp1 in combinations(groups, 2):
            if set(groups[grp0]) & set(groups[grp1]):
                raise RuntimeError("Groups contain same element",
                                   set(groups[grp0]) & set(groups[grp1]))
    if args.speciesgroups:
        groups = {}
        for elem in args.speciesgroups:
            elem = elem.split(':')
            groups[elem[0]] = mvf.get_sample_indices(labels=elem[1].split(','))
        args.speciesgroups = groups.copy()
        for specgroup in groups:
            ngroup = 0
            for allelegroup in args.allele_groups.values():
                if set(allelegroup) & set(groups[specgroup]):
                    ngroup += 1
                    if ngroup > 1:
                        raise RuntimeError(specgroup, "split across 2+ groups")
    # MODULES
    if args.module == 'Coverage':
        module = Coverage(params=vars(args))
    elif args.module == 'GroupUniqueAlleleWindow':
        module = GroupUniqueAlleleWindow(params=vars(args))
    elif args.module == 'PiDiversityWindow':
        module = PiDiversityWindow(params=vars(args))
    elif args.module == 'PairwiseNS':
        module = PairwiseNS(params=vars(args))
    # RUN MODULE
    module.analyze(mvf)

    return ''


if __name__ == "__main__":
    if "--morehelp" in sys.argv:
        modulehelp()
        sys.exit()
    main()
