#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

MVFtools: Multisample Variant Format Toolkit
http://www.github.org/jbpease/mvftools

MVFbiolib: Biological sequence object library for use in MVFtools
@author: James B. Pease
@author: Ben K. Rosenzweig

version: 2015-02-01 - First Public Release
version: 2015-06-11 - Cleanup and Python 3.x compatibility fixes
@version 2015-09-04 - Cleanup and a few fixes

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
from itertools import combinations, permutations

STANDARD_CODON_TABLE = {
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

AMBIG_CODON_TABLE = dict([
    (bases + ambig, aa) for bases, aa in [
        ('AC', 'T'), ('CC', 'P'), ('CG', 'R'), ('CT', 'L'),
        ('GC', 'A'), ('GG', 'G'), ('GT', 'V'), ('TC', 'S')]
    for ambig in 'KMRYWSN'] +
   [(bases + 'R', aa) for bases, aa in [
     ('AA', 'K'), ('AC', 'R'), ('CA', 'Q'),
     ('GA', 'E'), ('TA', '*'), ('TT', 'L')]] +
   [(bases + 'Y', aa) for bases, aa in [
     ('AA', 'N'), ('AG', 'S'), ('AT', 'I'), ('CA', 'H'), ('GA', 'D'),
     ('TA', 'Y'), ('TG', 'C'), ('TT', 'F')]] +
   [('MGA', 'R'), ('MGG', 'R'), ('MGR', 'R'), ('YTA', 'L'),
    ('YTG', 'L'), ('YTR', 'L'), ('TRA', '*')]
   )

FULL_CODON_TABLE = STANDARD_CODON_TABLE.copy()
FULL_CODON_TABLE.update(AMBIG_CODON_TABLE)

AMBIGSTOPS = ('TAA', 'TAG', 'TGA',
              'WGA', 'WAG', 'WAA', 'KGA', 'KAG', 'KAA', 'YGA', 'YAG', 'YAA',
              'TRA', 'TRG', 'TKA', 'TMG', 'TMA', 'TSA', 'TWG', 'TWA',
              'TGW', 'TAS', 'TAW', 'TGR', 'TAR', 'TGM', 'TAK', 'TAM')

# Converts base to dinucleotides
HAPSPLIT = {'A': 'AA', 'C': 'CC', 'G': 'GG', 'T': 'TT',
            'R': 'AG', 'Y': 'CT', 'M': 'AC', 'K': 'GT', 'S': 'CG', 'W': 'AT',
            'N': 'XX', '-': '--', 'X': 'XX'}

# Converts dinucleotide string to ambigous base
HAPJOIN = dict(
    [("AC", "M"), ("AG", "R"), ("AT", "W"),
     ("CG", "S"), ("CT", "Y"), ("GT", "K")] +
    [(x + x, x) for x in 'ATGC'] +
    [(x[::-1], y) for x, y in [
        ("AC", "M"), ("AG", "R"), ("AT", "W"),
        ("CG", "S"), ("CT", "Y"), ("GT", "K")]] +
    [(x + 'N', 'X') for x in 'ATGCKMRSWYNX'] +
    [('N' + x, 'X') for x in 'ATGCKMRSWYX'] +
    [(x + 'X', 'X') for x in 'ATGCKMRSWYNX'] +
    [('X' + x, 'X') for x in 'ATGCKMRSWYN']
)

# VCF genotype index codes for PL order
GTCODES = [(0, 0), (0, 1), (1, 1), (0, 2), (1, 2),
           (2, 2), (0, 3), (1, 3), (2, 3), (3, 3)]


def make_allele_resolution_table():
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
                   for x, y in combinations('BbDdHhKkMmNnRrSsVvWwXxYy', 2)] +
                  [(nuc + x, nuc) for x in 'acgkmrstwybdhv'
                   for nuc in 'ACGKMRSTWY'] +
                  ambig +
                  [(x.lower(), y.lower()) for (x, y) in ambig]
                  )
    xtable.update([(x[::-1], y) for (x, y) in xtable.items()])
    return xtable

COMPBASES = [('A', 'T'), ('C', 'G'), ('K', 'M'), ('R', 'Y'), ('S', 'W'),
             ('B', 'V'), ('D', 'H'), ('N', 'X'), ('-', '-'), ('X', 'X')]

COMPCODE = dict(COMPBASES +
                [(v, k) for (k, v) in COMPBASES] +
                [(k.lower(), v.lower()) for (k, v) in COMPBASES] +
                [(v.lower(), k.lower()) for (k, v) in COMPBASES]
                )


def complement(sequence):
    """Returns sequence complement"""
    return ''.join([COMPCODE[x] for x in sequence[::-1]])


def make_dna_pwtable():
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
    PWCODE = dict([(x + x, 0) for x in 'ATGC'] +
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
    PWCODE.update([(k[::-1], v) for (k, v) in PWCODE.items()])
    return PWCODE

if __name__ == ("__main__"):
    print("""MVF biosequence library v. 2015-06-11, please run one of the
          other MVFtools scripts to access these functions""")