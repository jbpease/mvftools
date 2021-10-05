#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mvfbiolib - Biological sequence object library for use in MVFtools
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

from itertools import combinations, permutations
from random import randint, choice


class MvfBioLib():
    """MVF Biological Information Library Object
    """

    def __init__(self):
        self.codon_tables = {}
        self._populate_codon_tables()
        self.complement_bases = {}
        self._populate_complement()
        self.validchars = {
            'dna': 'ACGT',
            'amino': 'ACDEFGHIKLMNPQRSTVWY',
            'dnaambig2': 'KMRSWY',
            'dnaambig3': 'BDHV',
            'dnaambig23': 'KMRSWYBDHV',
            'dnaambig4': 'NX',
            'dna+ambig2': 'ACGTKMRSWY',
            'dna+ambig3': 'ACGTKMRSWYBDHV',
            'dnaambigall': 'KMRSWYBDHVNX',
            'dna+ambigall': 'ACGTKMRYWSBDHVNX'
            }
        self.stop_codons = (
            'TAA', 'TAG', 'TGA', 'WGA', 'WAG', 'WAA',
            'KGA', 'KAG', 'KAA', 'YGA', 'YAG', 'YAA',
            'TRA', 'TRG', 'TKA', 'TMG', 'TMA', 'TSA',
            'TWG', 'TWA', 'TGW', 'TAS', 'TAW', 'TGR',
            'TAR', 'TGM', 'TAK', 'TAM')

        self.splitbases = {
            'A': 'AA', 'C': 'CC', 'G': 'GG', 'T': 'TT',
            'R': 'AG', 'Y': 'CT', 'M': 'AC', 'K': 'GT',
            'S': 'CG', 'W': 'AT',
            'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG',
            'N': 'XX', '-': '--', 'X': 'XX'}
        # Joinbases for diploid
        self.joinbases = dict(
            [("AC", "M"), ("AG", "R"), ("AT", "W"),
             ("CG", "S"), ("CT", "Y"), ("GT", "K")] +
            [(x, x) for x in 'ATGCRYWSMK'] +
            [(x + x, x) for x in 'ATGCRYWSMK'] +
            [(x[::-1], y) for x, y in [
                ("AC", "M"), ("AG", "R"), ("AT", "W"),
                ("CG", "S"), ("CT", "Y"), ("GT", "K")]] +
            [(x + 'N', 'X') for x in 'ATGCKMRSWYNX'] +
            [('N' + x, 'X') for x in 'ATGCKMRSWYX'] +
            [(x + 'X', 'X') for x in 'ATGCKMRSWYNX'] +
            [('X' + x, 'X') for x in 'ATGCKMRSWYN'] +
            [(x + '*', x) for x in 'ATGCKMRSWY'] +
            [('*' + x, x) for x in 'ATGCKMRSWY'] +
            [('N*', 'X'), ('*N', 'X'), ('**', '-')] +
            [(x + y + z, 'X') for (x, y, z) in permutations('ATGC', 3)] +
            [(w + x + y + z, 'X') for (w, x, y, z) in permutations('ATGC', 4)]
            )
        self.joinbases.update([(k.lower(), v.lower())
                               for k, v in self.joinbases.items()])
        # Joinbases for polyploid
        self.joinbasespoly = self.joinbases.copy()
        self.joinbasespoly.update(
            [("ACG", "V"), ("AGC", "V"), ("GCA", "V"),
             ("GAC", "V"), ("CAG", "V"), ("CGA", "V"),
             ("ACT", "H"), ("ATC", "H"), ("TCA", "H"),
             ("TAC", "H"), ("CAT", "H"), ("CTA", "H"),
             ("GCT", "B"), ("GTC", "B"), ("TCG", "B"),
             ("TGC", "B"), ("CGT", "B"), ("CTG", "B"),
             ("GAT", "D"), ("GTA", "D"), ("TAG", "D"),
             ("TGA", "D"), ("AGT", "D"), ("ATG", "D")] +
            [(w + x + y + z, 'N') for (w, x, y, z) in permutations('ATGC', 4)]
            )
        self.joinbasespoly.update([
            (k.lower(), v.lower())
            for k, v in self.joinbasespoly.items()])
        # VCF genotype index codes for PL order for diploid
        self.vcf_gtcodes = [
            (0, 0), (0, 1), (1, 1), (0, 2), (1, 2),
            (2, 2), (0, 3), (1, 3), (2, 3), (3, 3)]
        # VCF genotype index codes for PL order for tetraploid
        self.vcf_gtcodes_tetra = [tuple(
            [0] * x[0] + [1] * x[1] + [2] * x[2] + [3] * x[3]) for x in [
                # Stanza 0
                (4, 0, 0, 0), (3, 1, 0, 0), (2, 2, 0, 0), (1, 3, 0, 0),
                (0, 4, 0, 0),
                #
                (3, 0, 1, 0), (2, 1, 1, 0), (1, 2, 1, 0), (0, 3, 1, 0),
                (2, 0, 2, 0), (1, 1, 2, 0), (0, 2, 2, 0),
                (1, 0, 3, 0), (0, 1, 3, 0),
                (0, 0, 4, 0),
                # Stanza 1
                (3, 0, 0, 1), (2, 1, 0, 1), (1, 2, 0, 1), (0, 3, 0, 1),
                (2, 0, 1, 1), (1, 1, 1, 1), (0, 2, 1, 1),
                (1, 0, 2, 1), (0, 1, 2, 1),
                (0, 0, 3, 1),
                # Stanza 2
                (2, 0, 0, 2), (1, 1, 0, 2), (0, 2, 0, 2),
                (1, 0, 1, 2), (0, 1, 1, 2),
                (0, 0, 2, 2),
                # Stanza 3
                (1, 0, 0, 3), (0, 1, 0, 3),
                (0, 0, 1, 3),
                # Stanza 4
                (0, 0, 0, 4)]]
        # VCF genotype index codes for PL order for hexaploid
        self.vcf_gtcodes_hex = [tuple(
            [0] * x[0] + [1] * x[1] +
            [2] * x[2] + [3] * x[3]) for x in [
                # Stanza 0
                (6, 0, 0, 0), (5, 1, 0, 0), (4, 2, 0, 0), (3, 3, 0, 0),
                (2, 4, 0, 0), (1, 5, 0, 0), (0, 6, 0, 0),
                #
                (5, 0, 1, 0), (4, 1, 1, 0), (3, 2, 1, 0), (2, 3, 1, 0),
                (1, 4, 1, 0), (0, 5, 1, 0),
                #
                (4, 0, 2, 0), (3, 1, 2, 0), (2, 2, 2, 0), (1, 3, 2, 0),
                (0, 4, 2, 0),
                #
                (3, 0, 3, 0), (2, 1, 3, 0), (1, 2, 3, 0), (0, 3, 3, 0),
                (2, 0, 4, 0), (1, 1, 4, 0), (0, 2, 4, 0),
                (1, 0, 5, 0), (0, 1, 5, 0),
                (0, 0, 6, 0),
                # Stanza 1
                (5, 0, 0, 1), (4, 1, 0, 1), (3, 2, 0, 1), (2, 3, 0, 1),
                (1, 4, 0, 1), (0, 5, 0, 1),
                #
                (4, 0, 1, 1), (3, 1, 1, 1), (2, 2, 1, 1), (1, 3, 1, 1),
                (0, 4, 1, 1),
                #
                (3, 0, 2, 1), (2, 1, 2, 1), (1, 2, 2, 1), (0, 3, 2, 1),
                (2, 0, 3, 1), (1, 1, 3, 1), (0, 2, 3, 1),
                (1, 0, 4, 1), (0, 1, 4, 1),
                (0, 0, 5, 1),
                # Stanza 2
                (4, 0, 0, 2), (3, 1, 0, 2), (2, 2, 0, 2), (1, 3, 0, 2),
                (0, 4, 0, 2),
                #
                (3, 0, 1, 2), (2, 1, 1, 2), (1, 2, 1, 2), (0, 3, 1, 2),
                (2, 0, 2, 2), (1, 1, 2, 2), (0, 2, 2, 2),
                (1, 0, 3, 2), (0, 1, 3, 2),
                (0, 0, 4, 2),
                # Stanza 3
                (3, 0, 0, 3), (2, 1, 0, 3), (1, 2, 0, 3), (0, 3, 0, 3),
                (2, 0, 1, 3), (1, 1, 1, 3), (0, 2, 1, 3),
                (1, 0, 2, 3), (0, 1, 2, 3),
                (0, 0, 3, 3),
                # Stanza 4
                (2, 0, 0, 4), (1, 1, 0, 4), (0, 2, 0, 4),
                (1, 0, 1, 4), (0, 1, 1, 4),
                (0, 0, 2, 4),
                # Stanza 5
                (1, 0, 0, 5), (0, 1, 0, 5),
                (0, 0, 1, 5),
                # Stanza 6
                (0, 0, 0, 6)]]

    def _populate_complement(self):
        compbases = [
            ('A', 'T'), ('C', 'G'), ('K', 'M'), ('R', 'Y'), ('S', 'W'),
            ('B', 'V'), ('D', 'H'), ('N', 'X'), ('-', '-'), ('X', 'X')]

        self.complement_bases = dict(
            compbases +
            [(v, k) for (k, v) in compbases] +
            [(k.lower(), v.lower()) for (k, v) in compbases] +
            [(v.lower(), k.lower()) for (k, v) in compbases]
            )
        return ''

    def _populate_codon_tables(self):
        self.codon_tables['standard'] = {
            "": "", "---": '-',
            "AAA": "K", "AAC": "N", "AAG": "K", "AAT": "N",
            "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
            "AGA": "R", "AGC": "S", "AGG": "R", "AGT": "S",
            "ATA": "I", "ATC": "I", "ATG": "M", "ATT": "I",
            "CAA": "Q", "CAC": "H", "CAG": "Q", "CAT": "H",
            "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
            "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
            "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
            "GAA": "E", "GAC": "D", "GAG": "E", "GAT": "D",
            "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
            "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
            "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
            "TAA": "*", "TAC": "Y", "TAG": "*", "TAT": "Y",
            "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
            "TGA": "*", "TGC": "C", "TGG": "W", "TGT": "C",
            "TTA": "L", "TTC": "F", "TTG": "L", "TTT": "F"
            }
        self.codon_tables['ambig'] = dict(
            [(bases + ambig, aa) for bases, aa in [
                ('AC', 'T'), ('CC', 'P'), ('CG', 'R'), ('CT', 'L'),
                ('GC', 'A'), ('GG', 'G'), ('GT', 'V'), ('TC', 'S')]
             for ambig in 'KMRYWSBDHVX'] +
            [(bases + 'R', aa) for bases, aa in [
                ('AA', 'K'), ('AG', 'R'),
                ('CA', 'Q'), ('GA', 'E'),
                ('TA', '*'), ('TT', 'L')]] +
            [(bases + 'Y', aa) for bases, aa in [
                ('AA', 'N'), ('AG', 'S'),
                ('CA', 'H'), ('GA', 'D'),
                ('TA', 'Y'), ('TG', 'C'), ('TT', 'F')]] +
            [('ATH', 'I'), ('ATW', 'I'), ('ATM', 'I'), ('ATY', 'I'),
             ('MGA', 'R'), ('MGG', 'R'), ('MGR', 'R'),
             ('YTA', 'L'), ('YTG', 'L'), ('YTR', 'L'),
             ('TRA', '*'),
             ])
        self.codon_tables['full'] = self.codon_tables['standard'].copy()
        self.codon_tables['full'].update(self.codon_tables['ambig'])
        return ''

    def complement(self, sequence):
        """Returns sequence complement"""
        return ''.join([self.complement_bases[x] for x in sequence[::-1]])

    def hapsplit(self, alleles, mode):
        """Process Alleles into Haplotypes"""
        if all(x not in 'RYMKWSBHDV' for x in alleles):
            if mode in ['major', 'minor', 'randomone']:
                return alleles
            if mode in ['majorminor', 'randomboth']:
                return ''.join([base*2 for base in alleles])
        if mode in ['major', 'minor', 'majorminor']:
            hapleles = ''.join([self.splitbases[x] for x in alleles])
            counts = sorted([(hapleles.count(x), x) for x in set(hapleles)],
                            reverse=True)
            order = [x[1] for x in counts]
            newalleles = []
            for base in alleles:
                if base in 'RYMKWSBHDV':
                    newalleles.extend(
                        [x for x in order if x in self.splitbases[base]])
                else:
                    newalleles.extend([base, base])
            if mode == 'major':
                alleles = ''.join([x[0] for x in newalleles])
            elif mode == 'minor':
                alleles = ''.join([x[1] for x in newalleles])
            elif mode == 'majorminor':
                alleles = ''.join(newalleles)
        elif mode == 'randomone':
            alleles = ''.join([
                self.splitbases[x][randint(0, 1)] for x in alleles])
        elif mode == 'randomboth':
            randx = randint(0, 1)
            alleles = ''.join([self.splitbases[x][randx] +
                               self.splitbases[x][1 - randx]
                               for x in alleles])
        return alleles

    def make_allele_resolution_table(self):
        """Generate Allele Conflict Table"""
        ambig = [('AG', 'R'), ('CT', 'Y'), ('AC', 'M'), ('GT', 'K'),
                 ('AT', 'W'), ('CG', 'S'), ('AR', 'R'), ('AY', 'X'),
                 ('AM', 'M'), ('AK', 'X'), ('AW', 'W'), ('AS', 'X'),
                 ('CR', 'X'), ('CY', 'Y'), ('CM', 'M'), ('CK', 'X'),
                 ('CW', 'X'), ('CS', 'S'), ('GR', 'R'), ('GY', 'X'),
                 ('GM', 'X'), ('GK', 'K'), ('GW', 'X'), ('GS', 'S'),
                 ('TR', 'X'), ('TY', 'Y'), ('TM', 'X'), ('TK', 'K'),
                 ('TW', 'W'), ('TS', 'X')]
        xtable = dict([(x + x, x) for x in 'AaCcGgKkMmRrSsTtWwXxYy-'] +
                      [(x + y, x) for x in 'AaCcGgKkMmRrSsTtWwYy'
                       for y in 'XxNn-'] +
                      [(x + '-', x) for x in 'XxNn'] +
                      [(x + y, 'X') for x in 'XxNn' for y in 'XxNn'] +
                      [(x + y, 'X')
                       for x, y in combinations(
                           'BbDdHhKkMmNnRrSsVvWwXxYy', 2)] +
                      [(nuc + x, nuc) for x in 'acgkmrstwybdhv'
                       for nuc in 'ACGKMRSTWY'] +
                      ambig +
                      [(x.lower(), y.lower()) for (x, y) in ambig])
        xtable.update([(x[::-1], y) for (x, y) in xtable.items()])
        return xtable

    def make_dna_pwtable(self):
        """Returns pairwise distance class codes
        H=homozygous, h=heterozygous
        0 = H-H match
        1 = H-H mismatch
        2 = H-h partial match
        3 = H-h mismatch
        4 = h-h double match
        5 = h-h partial match
        6 = h-h double mismatch
        """
        pwcode = dict([(x + x, 0) for x in 'ATGC'] +
                      [(x + y, 1) for (x, y) in permutations('ATGC', 2)] +
                      [('A' + x, 2) for x in 'MRW'] +
                      [('T' + x, 2) for x in 'KYW'] +
                      [('G' + x, 2) for x in 'KRS'] +
                      [('C' + x, 2) for x in 'MYS'] +
                      [('A' + x, 3) for x in 'KYS'] +
                      [('T' + x, 3) for x in 'MRS'] +
                      [('G' + x, 3) for x in 'MYW'] +
                      [('C' + x, 3) for x in 'KRW'] +
                      [(x + x, 4) for x in 'KMRSWY'] +
                      [('K' + x, 5) for x in 'RYWS'] +
                      [('M' + x, 5) for x in 'RYWS'] +
                      [('R' + x, 5) for x in 'KMWS'] +
                      [('S' + x, 5) for x in 'KMRY'] +
                      [('W' + x, 5) for x in 'KMRY'] +
                      [('Y' + x, 5) for x in 'KMWS'] +
                      [('KM', 6), ('RY', 6), ('WS', 6)])
        pwcode.update([(k[::-1], v) for (k, v) in pwcode.items()])
        return pwcode

    def merge_bases(self, bases):
        """Merges bases in a list or string"""
        newbase = 'X'
        bases = set(bases) - set('-')
        if not bases:
            newbase = '-'
        elif bases == set('X'):
            newbase = 'X'
        elif len(bases) == 1:
            newbase = list(bases)[0]
        else:
            try:
                newbases = set(''.join([self.splitbases[x]
                                        for x in bases if x not in 'X-']))
                if len(newbases) > 2:
                    newbase = 'X'
                elif len(newbases) == 1:
                    newbase = list(newbases)[0]
                else:
                    newbase = self.joinbases[''.join(sorted(list(newbases)))]
            except:
                print(bases)
        return newbase

    def abpattern(self, num, digits=0):
        """Convert a decimal integer to a binary AB pattern
        """
        return bin(num)[2:].zfill(digits).replace('0', 'A').replace('1', 'B')

    def randomnuc(self, ambigcode, alleles=1):
        """Return random nucleotide for given ambiguity code
        """
        return choice(self.splitbases[ambigcode])



if __name__ == ("__main__"):
    print("""MVF biosequence library, please run one of the
          other MVFtools scripts to access these functions""")
