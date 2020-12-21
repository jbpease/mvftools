# -*- coding: utf-8 -*-
"""
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

from random import choices
from mvfbase import encode_mvfstring, decode_mvfstring

NCOL = 10
TEST_STRINGS = [
    'AAAAAAAAAA',
    'X---------',
    'XXXXXXXXXX',
    'ATTTTTTTTT',
    'A---------',
    'ATCCCCCCCC',
    'A-CCCCCCCC',
    'ATGCCCCCCC',
    'AGCGGGGGGG',
    'AT--------',
    'A-T-------',
    'A--T------',
]

NRAND = 100000
RANDOM_STRINGS = [''.join(choices("ATGCX-", k=10)) for _ in range(NRAND)]
for x in TEST_STRINGS:
    print(x)
    y = encode_mvfstring(x)
    print(y)
    z = decode_mvfstring(y, NCOL)
    print(z)
    print(x == z)
    print("==========")

RANDOM_PASS = 0
for x in RANDOM_STRINGS:
    #print(x)
    y = encode_mvfstring(x)
    #print(y)
    z = decode_mvfstring(y, NCOL)
    #print(z)
    #print(x == z)
    #print("==========")
    if x == z:
        RANDOM_PASS += 1
    else:
        print(x, y, z)
print("RANDOM TRIALS PASSED: {}/{}".format(RANDOM_PASS, NRAND))
